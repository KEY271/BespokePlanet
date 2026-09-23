#!/usr/bin/env Rscript

# Produce final-year climate diagnostics for the moist run over Earth terrain.
# This mirrors scripts/analyze_moist_land_sea.R so the analytic-continent run
# and the Earth-terrain run can be compared figure by figure, and adds the
# checks listed in docs/cases/land-sea-earth.md (monsoon reversal, land water
# budget, ocean surface-pressure ripples, spin-up).
# The script uses base R only. Binary layouts and units are defined by the
# case's metadata.json.

args <- commandArgs(trailingOnly = TRUE)
case_dir <- if (length(args) >= 1) args[[1]] else "output/moist_land_sea_earth_t31"
analysis_dir <- if (length(args) >= 2) args[[2]] else file.path(case_dir, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
png_device_type <- if (capabilities("aqua")) "quartz" else "cairo"
# Shared snow and sea-ice map, kept next to this script.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir <- if (length(script_path) == 1) dirname(normalizePath(script_path)) else "scripts"
source(file.path(script_dir, "snow_sea_ice_map.R"))

metadata_path <- file.path(case_dir, "metadata.json")
if (!file.exists(metadata_path)) stop("Missing ", metadata_path)
metadata <- paste(readLines(metadata_path, warn = FALSE), collapse = "\n")

json_number <- function(key) {
  pattern <- sprintf('"%s"\\s*:\\s*([-+0-9.Ee]+)', key)
  hit <- regmatches(metadata, regexec(pattern, metadata, perl = TRUE))[[1]]
  if (length(hit) != 2) stop("Could not read metadata number: ", key)
  as.numeric(hit[[2]])
}

json_array <- function(key) {
  pattern <- sprintf('"%s"\\s*:\\s*\\[([^]]+)\\]', key)
  hit <- regmatches(metadata, regexec(pattern, metadata, perl = TRUE))[[1]]
  if (length(hit) != 2) stop("Could not read metadata array: ", key)
  as.numeric(trimws(strsplit(hit[[2]], ",", fixed = TRUE)[[1]]))
}

mu <- json_array("mu")
nlon <- as.integer(json_array("nlon"))
a_half <- json_array("hybrid_a_half_pa")
b_half <- json_array("hybrid_b_half")
reference_half_pressure <- json_array("reference_half_level_pressure_pa")
gravity <- json_number("gravity_acceleration_m_s-2")
months_per_year <- as.integer(json_number("months_per_year"))
days_per_month <- as.integer(json_number("days_per_month"))

earth_radius <- 6.371e6 # core/src/config/planet_parameters.f90
nlat <- length(mu)
nlev <- length(a_half) - 1L
npoints <- sum(nlon)
latitude <- asin(mu) * 180 / pi
cos_latitude <- sqrt(pmax(0, 1 - mu^2))
reference_full_pressure_hpa <- 0.5 * (
  reference_half_pressure[-length(reference_half_pressure)] +
    reference_half_pressure[-1]
) / 100

read_grid <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  values <- readBin(con, what = "double", n = npoints, size = 8, endian = "little")
  if (length(values) != npoints) stop("Unexpected grid file size: ", path)
  values
}

read_zonal <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  values <- readBin(con, what = "double", n = nlat * nlev, size = 8, endian = "little")
  if (length(values) != nlat * nlev) stop("Unexpected zonal file size: ", path)
  matrix(values, nrow = nlat, ncol = nlev)
}

ring_mean <- function(values, weights = NULL) {
  result <- rep(NA_real_, nlat)
  offset <- 0L
  for (j in seq_len(nlat)) {
    indices <- offset + seq_len(nlon[[j]])
    if (is.null(weights)) {
      result[[j]] <- mean(values[indices])
    } else {
      keep <- is.finite(values[indices]) & weights[indices] > 0
      if (any(keep)) {
        result[[j]] <- weighted.mean(values[indices][keep], weights[indices][keep])
      }
    }
    offset <- offset + nlon[[j]]
  }
  result
}

ring_expand <- function(values) rep(values, times = nlon)

# Coordinates follow the metadata convention: longitude 2*pi*k/nlon[j], with
# the flat files ordered ring-major from south to north.
grid_latitude <- ring_expand(latitude)
grid_longitude <- unlist(lapply(nlon, function(n) 360 * (seq_len(n) - 1) / n), use.names = FALSE)

monthly_files <- list.files(
  case_dir,
  pattern = "^monthly_surface_temperature_m[0-9]{4}\\.bin$",
  full.names = FALSE
)
month_indices <- sort(as.integer(sub(".*_m([0-9]{4})\\.bin$", "\\1", monthly_files)))
if (length(month_indices) < months_per_year) stop("Fewer than one year of monthly outputs")
final_months <- tail(month_indices, months_per_year)
if (!all(diff(final_months) == 1)) stop("Final-year monthly files are not consecutive")

# m0001 is April, so translate each file index to a conventional calendar month.
calendar_month <- ((final_months + 2L) %% 12L) + 1L
calendar_labels <- month.abb[calendar_month]

temperature_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
land_temperature_monthly <- ocean_temperature_monthly <- temperature_monthly
precipitation_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
evaporation_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
cloud_cover_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
surface_pressure_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
surface_water_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
zonal_v_monthly <- array(NA_real_, dim = c(nlat, nlev, months_per_year))
ice_fields <- c("sea_ice_fraction", "sea_ice_volume", "sea_ice_thickness", "sea_ice_temperature")
ice_files <- outer(ice_fields, final_months, function(field, month)
  file.path(case_dir, sprintf("monthly_%s_m%04d.bin", field, month)))
has_sea_ice <- all(file.exists(ice_files))
if (any(file.exists(ice_files)) && !has_sea_ice) stop("Incomplete monthly sea-ice output")
if (has_sea_ice) {
  sea_ice_fraction_monthly <- sea_ice_volume_monthly <-
    sea_ice_thickness_monthly <- sea_ice_temperature_monthly <-
      matrix(NA_real_, nrow = npoints, ncol = months_per_year)
}

for (m in seq_along(final_months)) {
  suffix <- sprintf("m%04d.bin", final_months[[m]])
  temperature_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_surface_temperature_", suffix)))
  land_path <- file.path(case_dir, paste0("monthly_land_temperature_", suffix))
  ocean_path <- file.path(case_dir, paste0("monthly_ocean_temperature_", suffix))
  land_temperature_monthly[, m] <- if (file.exists(land_path)) read_grid(land_path) else temperature_monthly[, m]
  ocean_temperature_monthly[, m] <- if (file.exists(ocean_path)) read_grid(ocean_path) else temperature_monthly[, m]
  precipitation_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_precipitation_", suffix)))
  evaporation_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_evaporation_", suffix)))
  cloud_cover_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_cloud_cover_", suffix)))
  surface_pressure_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_surface_pressure_", suffix)))
  surface_water_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_surface_water_", suffix)))
  zonal_v_monthly[, , m] <- read_zonal(file.path(case_dir, paste0("monthly_zonal_v_", suffix)))
  if (has_sea_ice) {
    sea_ice_fraction_monthly[, m] <- read_grid(ice_files[1, m])
    sea_ice_volume_monthly[, m] <- read_grid(ice_files[2, m])
    sea_ice_thickness_monthly[, m] <- read_grid(ice_files[3, m])
    sea_ice_temperature_monthly[, m] <- read_grid(ice_files[4, m])
  }
}

land_fraction <- read_grid(file.path(case_dir, "land_fraction.bin"))
surface_height <- read_grid(file.path(case_dir, "surface_height.bin"))
land <- land_fraction >= 0.5
temperature_annual_k <- rowMeans(temperature_monthly)
temperature_annual_c <- temperature_annual_k - 273.15
land_temperature_annual_c <- rowMeans(land_temperature_monthly) - 273.15
ocean_temperature_annual_c <- rowMeans(ocean_temperature_monthly) - 273.15
precipitation_annual_mm <- rowSums(precipitation_monthly * days_per_month)
evaporation_annual_mm <- rowSums(evaporation_monthly * days_per_month)
cloud_cover_annual <- rowMeans(cloud_cover_monthly)
surface_pressure_annual_hpa <- rowMeans(surface_pressure_monthly) / 100
surface_water_annual <- rowMeans(surface_water_monthly)

# Gauss-Legendre weights are reconstructed from P_N'(mu) and split evenly over
# the longitudes of each ring, so the point weights sum to one.
legendre_and_previous <- function(n, x) {
  p_nm1 <- 1
  p_n <- x
  if (n == 0) return(c(1, NA_real_))
  if (n == 1) return(c(p_n, p_nm1))
  for (k in 2:n) {
    p_np1 <- ((2 * k - 1) * x * p_n - (k - 1) * p_nm1) / k
    p_nm1 <- p_n
    p_n <- p_np1
  }
  c(p_n, p_nm1)
}
gaussian_weights <- vapply(mu, function(x) {
  polynomials <- legendre_and_previous(nlat, x)
  derivative <- nlat * (x * polynomials[[1]] - polynomials[[2]]) / (x^2 - 1)
  2 / ((1 - x^2) * derivative^2)
}, numeric(1))
point_weight <- ring_expand((gaussian_weights / 2) / nlon)

# Summer is the high-sun half-year of the hemisphere, as in the analytic
# land/sea analysis: April-September in the north, October-March in the south.
northern <- grid_latitude >= 0
summer_columns_north <- which(calendar_month %in% 4:9)
summer_columns_south <- which(calendar_month %in% c(10:12, 1:3))
summer_mask <- matrix(FALSE, nrow = npoints, ncol = months_per_year)
summer_mask[northern, summer_columns_north] <- TRUE
summer_mask[!northern, summer_columns_south] <- TRUE

temperatures_c <- land_temperature_monthly - 273.15
precipitation_mm_month <- precipitation_monthly * days_per_month
annual_temperature_c <- rowMeans(temperatures_c)
warmest_month_c <- apply(temperatures_c, 1, max)
coldest_month_c <- apply(temperatures_c, 1, min)
months_above_10 <- rowSums(temperatures_c > 10)

summer_precipitation_mm <- rowSums(precipitation_mm_month * summer_mask)
summer_fraction <- ifelse(precipitation_annual_mm > 0,
                          summer_precipitation_mm / precipitation_annual_mm, 0.5)
masked_min <- function(values, mask) apply(ifelse(mask, values, Inf), 1, min)
masked_max <- function(values, mask) apply(ifelse(mask, values, -Inf), 1, max)
driest_month_mm <- apply(precipitation_mm_month, 1, min)
driest_summer_mm <- masked_min(precipitation_mm_month, summer_mask)
wettest_summer_mm <- masked_max(precipitation_mm_month, summer_mask)
driest_winter_mm <- masked_min(precipitation_mm_month, !summer_mask)
wettest_winter_mm <- masked_max(precipitation_mm_month, !summer_mask)

# Aridity threshold of Peel et al. (2007), in mm per year.
dryness_adjustment <- ifelse(summer_fraction >= 0.7, 280,
                             ifelse(summer_fraction <= 0.3, 0, 140))
dryness_threshold <- 20 * annual_temperature_c + dryness_adjustment

# First-letter group with the original -3 degC C/D boundary, kept identical to
# the analytic land/sea analysis so the two runs can be compared directly.
koppen_group_of <- function(cold_boundary) {
  group <- rep(NA_character_, npoints)
  group[precipitation_annual_mm < dryness_threshold] <- "B"
  unset <- is.na(group)
  group[unset & coldest_month_c >= 18] <- "A"
  unset <- is.na(group)
  group[unset & warmest_month_c < 10] <- "E"
  unset <- is.na(group)
  group[unset & coldest_month_c > cold_boundary & warmest_month_c >= 10] <- "C"
  group[is.na(group)] <- "D"
  group
}
koppen_group <- koppen_group_of(-3)
koppen_group_zero_boundary <- koppen_group_of(0)
koppen_group[!land] <- NA_character_
koppen_group_zero_boundary[!land] <- NA_character_

# Full Koppen-Geiger type following Peel et al. (2007), which uses the 0 degC
# C/D boundary. Only land points are classified.
koppen_type <- rep(NA_character_, npoints)
for (i in which(land)) {
  group <- koppen_group_zero_boundary[[i]]
  code <- group
  if (group == "B") {
    code <- paste0("B", if (precipitation_annual_mm[[i]] < 0.5 * dryness_threshold[[i]]) "W" else "S",
                   if (annual_temperature_c[[i]] >= 18) "h" else "k")
  } else if (group == "A") {
    code <- if (driest_month_mm[[i]] >= 60) {
      "Af"
    } else if (driest_month_mm[[i]] >= 100 - precipitation_annual_mm[[i]] / 25) {
      "Am"
    } else {
      "Aw"
    }
  } else if (group %in% c("C", "D")) {
    seasonality <- if (driest_summer_mm[[i]] < 40 && driest_summer_mm[[i]] < wettest_winter_mm[[i]] / 3) {
      "s"
    } else if (driest_winter_mm[[i]] < wettest_summer_mm[[i]] / 10) {
      "w"
    } else {
      "f"
    }
    heat <- if (warmest_month_c[[i]] >= 22) {
      "a"
    } else if (months_above_10[[i]] >= 4) {
      "b"
    } else if (group == "D" && coldest_month_c[[i]] < -38) {
      "d"
    } else {
      "c"
    }
    if (group == "D" && heat == "c" && coldest_month_c[[i]] < -38) heat <- "d"
    code <- paste0(group, seasonality, heat)
  } else if (group == "E") {
    code <- if (warmest_month_c[[i]] > 0) "ET" else "EF"
  }
  koppen_type[[i]] <- code
}

climate_table <- data.frame(
  longitude_deg = grid_longitude,
  latitude_deg = grid_latitude,
  land_fraction = land_fraction,
  surface_height_m = surface_height,
  annual_mean_surface_temperature_c = temperature_annual_c,
  annual_mean_land_temperature_c = ifelse(land_fraction > 0, land_temperature_annual_c, NA_real_),
  annual_mean_ocean_temperature_c = ifelse(land_fraction < 1, ocean_temperature_annual_c, NA_real_),
  annual_precipitation_mm = precipitation_annual_mm,
  annual_evaporation_mm = evaporation_annual_mm,
  annual_mean_cloud_cover = cloud_cover_annual,
  annual_mean_surface_water_kg_m2 = surface_water_annual,
  koppen_group = koppen_group,
  koppen_type = koppen_type
)
write.csv(climate_table, file.path(analysis_dir, "final_year_koppen_groups.csv"), row.names = FALSE)

# Real-world reference sites. Each target is snapped to the nearest native cell
# with land fraction >= 0.5; the reference Koppen type is the observed one on
# Earth and is only listed for comparison with the simulated type.
target_points <- data.frame(
  site = c("Manaus (Amazon)", "Kinshasa (Congo)", "Kolkata (India)",
           "Cairo (Sahara)", "Alice Springs (Australia)", "Buenos Aires (Pampas)",
           "Paris (W Europe)", "Chicago (N America)", "Yakutsk (Siberia)",
           "Lhasa (Tibet)", "Greenland interior", "Vostok (Antarctica)"),
  target_lon = c(300.0, 15.3, 88.4, 31.2, 133.9, 301.6, 2.3, 272.4, 129.7, 91.1, 318.0, 106.8),
  target_lat = c(-3.1, -4.3, 22.6, 30.0, -23.7, -34.6, 48.9, 41.9, 62.0, 29.7, 72.0, -78.5),
  earth_koppen = c("Af", "Aw", "Aw", "BWh", "BWh", "Cfa", "Cfb", "Dfa", "Dfd", "Dwb", "EF", "EF"),
  stringsAsFactors = FALSE
)
candidate <- which(land_fraction >= 0.5)
selected <- integer(nrow(target_points))
snap_distance_deg <- numeric(nrow(target_points))
for (s in seq_len(nrow(target_points))) {
  target_lon_rad <- target_points$target_lon[[s]] * pi / 180
  target_lat_rad <- target_points$target_lat[[s]] * pi / 180
  lon_rad <- grid_longitude[candidate] * pi / 180
  lat_rad <- grid_latitude[candidate] * pi / 180
  cosine_distance <- sin(target_lat_rad) * sin(lat_rad) +
    cos(target_lat_rad) * cos(lat_rad) * cos(lon_rad - target_lon_rad)
  best <- which.max(cosine_distance)
  selected[[s]] <- candidate[[best]]
  snap_distance_deg[[s]] <- acos(pmin(1, cosine_distance[[best]])) * 180 / pi
}

site_table <- data.frame(
  site = target_points$site,
  target_longitude_deg = target_points$target_lon,
  target_latitude_deg = target_points$target_lat,
  longitude_deg = grid_longitude[selected],
  latitude_deg = grid_latitude[selected],
  snap_distance_deg = snap_distance_deg,
  land_fraction = land_fraction[selected],
  surface_height_m = surface_height[selected],
  koppen_group = koppen_group[selected],
  koppen_type = koppen_type[selected],
  earth_reference_koppen = target_points$earth_koppen,
  annual_mean_temperature_c = land_temperature_annual_c[selected],
  annual_precipitation_mm = precipitation_annual_mm[selected],
  annual_mean_surface_water_kg_m2 = surface_water_annual[selected]
)
write.csv(site_table, file.path(analysis_dir, "representative_land_sites.csv"), row.names = FALSE)

# Shared map machinery: draw each native octahedral-grid cell with the prime
# meridian at the centre, splitting cells that straddle the -180/180 seam.
latitude_edges <- c(-90, 0.5 * (latitude[-1] + latitude[-nlat]), 90)
map_longitude <- ifelse(grid_longitude >= 180, grid_longitude - 360, grid_longitude)
ocean_colour <- "#b9d9eb"

draw_cells <- function(colour, include = rep(TRUE, npoints)) {
  offset <- 0L
  for (j in seq_len(nlat)) {
    delta_lon <- 360 / nlon[[j]]
    for (i in seq_len(nlon[[j]])) {
      index <- offset + i
      if (!include[[index]] || is.na(colour[[index]])) next
      centre <- map_longitude[[index]]
      left <- centre - delta_lon / 2
      right <- centre + delta_lon / 2
      cell <- colour[[index]]
      if (left < -180) {
        rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = cell, border = NA)
        rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = cell, border = NA)
      } else if (right > 180) {
        rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = cell, border = NA)
        rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = cell, border = NA)
      } else {
        rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = cell, border = NA)
      }
    }
    offset <- offset + nlon[[j]]
  }
}

map_frame <- function(main, background = ocean_colour) {
  plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
       xlab = "Longitude", ylab = "Latitude", main = main, axes = FALSE)
  if (!is.na(background)) rect(-180, -90, 180, 90, col = background, border = NA)
}

map_axes <- function(grid_alpha = 0.65) {
  abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60),
         col = adjustcolor("white", grid_alpha), lwd = 0.7)
  axis(1, at = seq(-180, 180, 60),
       labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
  axis(2, at = seq(-90, 90, 30),
       labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
  box()
}

# Coastline of the model's own land mask, drawn as the edges between land and
# ocean cells so that every map carries the same outline of the continents.
coastline_segments <- local({
  x0 <- c(); y0 <- c(); x1 <- c(); y1 <- c()
  ring_offset <- cumsum(c(0L, nlon))
  cell_index <- function(j, i) ring_offset[[j]] + ((i - 1L) %% nlon[[j]]) + 1L
  for (j in seq_len(nlat)) {
    delta_lon <- 360 / nlon[[j]]
    for (i in seq_len(nlon[[j]])) {
      index <- cell_index(j, i)
      if (!land[[index]]) next
      centre <- map_longitude[[index]]
      left <- centre - delta_lon / 2
      right <- centre + delta_lon / 2
      if (left < -180 || right > 180) next
      if (!land[[cell_index(j, i + 1L)]]) {
        x0 <- c(x0, right); x1 <- c(x1, right)
        y0 <- c(y0, latitude_edges[[j]]); y1 <- c(y1, latitude_edges[[j + 1]])
      }
      if (!land[[cell_index(j, i - 1L)]]) {
        x0 <- c(x0, left); x1 <- c(x1, left)
        y0 <- c(y0, latitude_edges[[j]]); y1 <- c(y1, latitude_edges[[j + 1]])
      }
      # Meridional edges: test the two neighbouring rings at the same longitude.
      for (neighbour in c(j - 1L, j + 1L)) {
        if (neighbour < 1L || neighbour > nlat) next
        delta_neighbour <- 360 / nlon[[neighbour]]
        i_neighbour <- floor(grid_longitude[[index]] / delta_neighbour + 0.5)
        neighbour_index <- cell_index(neighbour, i_neighbour + 1L)
        if (!land[[neighbour_index]]) {
          edge <- if (neighbour < j) latitude_edges[[j]] else latitude_edges[[j + 1]]
          x0 <- c(x0, left); x1 <- c(x1, right)
          y0 <- c(y0, edge); y1 <- c(y1, edge)
        }
      }
    }
  }
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
})

draw_coastline <- function(colour = "#333333", lwd = 0.8) {
  segments(coastline_segments$x0, coastline_segments$y0,
           coastline_segments$x1, coastline_segments$y1, col = colour, lwd = lwd)
}

open_map_png <- function(name, width = 2200, height = 1250, res = 180) {
  png(file.path(analysis_dir, name), width = width, height = height, res = res,
      type = png_device_type)
}

period_label <- sprintf("final year (months %d–%d)", min(final_months), max(final_months))

# --- Sea ice ---------------------------------------------------------------
# A and V are per ocean area. Point weights include the native Gaussian cell
# area; multiplying by ocean fraction converts them to whole-Earth fractions.
if (has_sea_ice) {
  ocean_fraction <- 1 - land_fraction
  ocean_weight <- point_weight * ocean_fraction
  earth_area_m2 <- 4 * pi * earth_radius^2
  ice_area_fraction <- function(indices, month) {
    sum(ocean_weight[indices] * sea_ice_fraction_monthly[indices, month])
  }
  ice_volume_fraction <- function(indices, month) {
    sum(ocean_weight[indices] * sea_ice_volume_monthly[indices, month])
  }
  hemisphere <- list(global = seq_len(npoints), north = which(grid_latitude >= 0),
                     south = which(grid_latitude < 0))
  ice_monthly <- do.call(rbind, lapply(names(hemisphere), function(region) {
    indices <- hemisphere[[region]]
    area <- vapply(seq_len(months_per_year), function(m) ice_area_fraction(indices, m), numeric(1))
    volume <- vapply(seq_len(months_per_year), function(m) ice_volume_fraction(indices, m), numeric(1))
    data.frame(region = region, output_month = final_months,
               calendar_month = calendar_month, month = calendar_labels,
               ice_area_m2 = area * earth_area_m2,
               ice_volume_m3 = volume * earth_area_m2,
               mean_thickness_m = ifelse(area > 0, volume / area, NA_real_))
  }))
  write.csv(ice_monthly, file.path(analysis_dir, "final_year_sea_ice_monthly.csv"), row.names = FALSE)

  ice_annual_fraction <- rowMeans(sea_ice_fraction_monthly)
  ice_annual_volume <- rowMeans(sea_ice_volume_monthly)
  ice_annual_thickness <- ifelse(ice_annual_fraction > 0,
                                 ice_annual_volume / ice_annual_fraction, NA_real_)
  ice_annual_temperature_c <- rowSums(sea_ice_fraction_monthly * sea_ice_temperature_monthly) /
    rowSums(sea_ice_fraction_monthly) - 273.15
  ice_annual_temperature_c[ice_annual_fraction == 0] <- NA_real_
  ice_annual_thickness[ocean_fraction == 0] <- NA_real_
  ice_annual_temperature_c[ocean_fraction == 0] <- NA_real_

  ice_fraction_breaks <- c(-0.001, 0.001, 0.15, 0.3, 0.5, 0.7, 0.85, 1.001)
  ice_fraction_colors <- c("#b9d9eb", hcl.colors(6, "Blues 3", rev = TRUE))
  ice_fraction_labels <- c("0", "0–15%", "15–30%", "30–50%", "50–70%", "70–85%", "85–100%")
  draw_ice_fraction_map <- function(values, title) {
    classes <- cut(pmin(1, pmax(0, values)), breaks = ice_fraction_breaks,
                   labels = FALSE, include.lowest = TRUE)
    map_frame(title, background = "#e6e2d8")
    draw_cells(ice_fraction_colors[classes], ocean_fraction > 0)
    map_axes(grid_alpha = 0.4)
    draw_coastline()
    legend("bottomleft", legend = c(ice_fraction_labels, "Land"),
           fill = c(ice_fraction_colors, "#e6e2d8"), border = "#555555",
           bg = "white", cex = 0.72, ncol = 2, title = "Ocean ice concentration")
  }

  march <- which(calendar_month == 3)
  september <- which(calendar_month == 9)
  png(file.path(analysis_dir, "final_year_sea_ice_seasonal_maps.png"),
      width = 2200, height = 2300, res = 180, type = png_device_type)
  par(mfrow = c(2, 1), mar = c(4.6, 5, 3.4, 1), oma = c(0.4, 0, 2.2, 0), las = 1)
  draw_ice_fraction_map(sea_ice_fraction_monthly[, march], "March monthly-mean sea ice")
  draw_ice_fraction_map(sea_ice_fraction_monthly[, september], "September monthly-mean sea ice")
  mtext(sprintf("Sea-ice concentration, %s", period_label), side = 3,
        outer = TRUE, line = 0.4, cex = 1.12)
  dev.off()

  ice_thickness_breaks <- c(-0.001, 0.25, 0.5, 1, 2, 4, 8, Inf)
  ice_thickness_colors <- hcl.colors(length(ice_thickness_breaks) - 1, "YlOrRd")
  ice_thickness_class <- cut(ice_annual_thickness, breaks = ice_thickness_breaks,
                             labels = FALSE, include.lowest = TRUE)
  open_map_png("final_year_sea_ice_mean_thickness_map.png")
  par(mar = c(6.2, 5, 4.1, 1), las = 1)
  map_frame(sprintf("Ice-area-weighted mean thickness, %s", period_label), background = "#e6e2d8")
  draw_cells(rep(ocean_colour, npoints), ocean_fraction > 0)
  draw_cells(ice_thickness_colors[ice_thickness_class], ice_annual_fraction > 0 & ocean_fraction > 0)
  map_axes(grid_alpha = 0.4)
  draw_coastline()
  legend("bottomleft", legend = c("<0.25", "0.25–0.5", "0.5–1", "1–2", "2–4", "4–8", ">=8", "No ice", "Land"),
         fill = c(ice_thickness_colors, ocean_colour, "#e6e2d8"),
         border = "#555555", bg = "white", cex = 0.72, ncol = 2, title = "Thickness (m)")
  mtext("Mean volume divided by mean ice concentration; only ice-covered ocean cells are coloured",
        side = 1, line = 4.6, cex = 0.72)
  dev.off()

  ice_temperature_breaks <- c(-Inf, -40, -30, -20, -10, -5, 0, Inf)
  ice_temperature_colors <- hcl.colors(length(ice_temperature_breaks) - 1, "Blue-Red 3")
  ice_temperature_class <- cut(ice_annual_temperature_c, breaks = ice_temperature_breaks,
                               labels = FALSE, include.lowest = TRUE)
  open_map_png("final_year_sea_ice_surface_temperature_map.png")
  par(mar = c(6.2, 5, 4.1, 1), las = 1)
  map_frame(sprintf("Ice-area-weighted surface temperature, %s", period_label), background = "#e6e2d8")
  draw_cells(rep(ocean_colour, npoints), ocean_fraction > 0)
  draw_cells(ice_temperature_colors[ice_temperature_class], ice_annual_fraction > 0 & ocean_fraction > 0)
  map_axes(grid_alpha = 0.4)
  draw_coastline()
  legend("bottomleft", legend = c("<-40", "-40–-30", "-30–-20", "-20–-10", "-10–-5", "-5–0", ">=0", "No ice", "Land"),
         fill = c(ice_temperature_colors, ocean_colour, "#e6e2d8"),
         border = "#555555", bg = "white", cex = 0.72, ncol = 2, title = "Ice skin (°C)")
  mtext("Ice-area and time weighted; only ice-covered ocean cells are coloured",
        side = 1, line = 4.6, cex = 0.72)
  dev.off()

  png(file.path(analysis_dir, "final_year_sea_ice_seasonal_cycle.png"),
      width = 1800, height = 1500, res = 180, type = png_device_type)
  par(mfrow = c(2, 1), mar = c(4.4, 5.4, 3, 1.1), oma = c(0, 0, 2.4, 0), las = 1)
  north_ice <- ice_monthly[ice_monthly$region == "north", ][order(ice_monthly$calendar_month[ice_monthly$region == "north"]), ]
  south_ice <- ice_monthly[ice_monthly$region == "south", ][order(ice_monthly$calendar_month[ice_monthly$region == "south"]), ]
  for (quantity in c("ice_area_m2", "ice_volume_m3")) {
    scale <- if (quantity == "ice_area_m2") 1e12 else 1e13
    label <- if (quantity == "ice_area_m2") expression("Ice area (10"^12*" m"^2*")") else
      expression("Ice volume (10"^13*" m"^3*")")
    plot(1:12, north_ice[[quantity]] / scale, type = "o", pch = 16, col = "#2368a2", lwd = 2.5,
         ylim = range(0, north_ice[[quantity]], south_ice[[quantity]]) / scale,
         xaxt = "n", xlab = "Calendar month", ylab = label,
         main = if (quantity == "ice_area_m2") "Sea-ice area" else "Sea-ice volume",
         panel.first = grid(col = "#dddddd"))
    lines(1:12, south_ice[[quantity]] / scale, type = "o", pch = 16, col = "#c26541", lwd = 2.5)
    axis(1, at = 1:12, labels = month.abb)
    legend("topright", c("Northern Hemisphere", "Southern Hemisphere"),
           col = c("#2368a2", "#c26541"), lwd = 2.5, pch = 16, bg = "white", cex = 0.8)
  }
  mtext(sprintf("Sea ice in the %s (native-grid area weighted)", period_label),
        side = 3, outer = TRUE, line = 0.5, cex = 1.12)
  dev.off()
}

# --- Land snow and sea ice ---------------------------------------------------
# Land snow cover and sea ice on one map, for March and September. Runs
# without snow output (older cases) skip the figure.
snow_fraction_monthly <- snow_ice_read_monthly(case_dir, "snow_fraction", final_months, read_grid)
has_snow <- !is.null(snow_fraction_monthly)
if (has_snow) {
  draw_snow_sea_ice_maps(file.path(analysis_dir, "final_year_snow_and_sea_ice_maps.png"),
                         nlon, latitude, grid_longitude, land_fraction, calendar_month,
                         snow_fraction_monthly, if (has_sea_ice) sea_ice_fraction_monthly else NULL,
                         json_number("masking_water_equivalent_kg_m-2"), period_label, png_device_type)
}

# --- Koppen first-letter groups -------------------------------------------
koppen_colors <- c(A = "#2f9e44", B = "#d8b365", C = "#ffd43b", D = "#4dabf7", E = "#f1f3f5")
koppen_labels <- c(A = "A Tropical", B = "B Dry", C = "C Temperate", D = "D Continental", E = "E Polar")

open_map_png("final_year_koppen_climate_zones.png", height = 1400)
par(mar = c(9.6, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Köppen climate groups over land, %s", period_label))
draw_cells(koppen_colors[koppen_group], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(koppen_labels, "Ocean"),
       fill = c(koppen_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.82)
site_map_longitude <- ifelse(site_table$longitude_deg > 180,
                             site_table$longitude_deg - 360,
                             site_table$longitude_deg)
points(site_map_longitude, site_table$latitude_deg, pch = 21, cex = 1.5,
       bg = "white", col = "#111111", lwd = 1.4)
text(site_map_longitude, site_table$latitude_deg, labels = seq_len(nrow(site_table)),
     cex = 0.55, font = 2)
legend_text <- sprintf("%d %s [%s]", seq_len(nrow(site_table)), site_table$site, site_table$koppen_type)
par(xpd = NA)
legend(x = -180, y = -118, legend = legend_text, ncol = 4, bty = "n", cex = 0.68,
       title = "Climograph sites (simulated type)", title.adj = 0)
par(xpd = FALSE)
mtext("Land is land fraction >= 0.5; C/D boundary is -3 degC", side = 1, line = 8.3, cex = 0.72)

dev.off()

# --- Full Koppen-Geiger types ---------------------------------------------
koppen_type_colors <- c(
  Af = "#0000fe", Am = "#0078ff", Aw = "#46aafa",
  BWh = "#ff0000", BWk = "#ff9695", BSh = "#f5a301", BSk = "#ffdb63",
  Csa = "#ffff00", Csb = "#c6c700", Csc = "#969600",
  Cwa = "#96ff96", Cwb = "#63c764", Cwc = "#329633",
  Cfa = "#c6ff4e", Cfb = "#66ff33", Cfc = "#33c701",
  Dsa = "#ff00fe", Dsb = "#c600c7", Dsc = "#963295", Dsd = "#966495",
  Dwa = "#abb1ff", Dwb = "#5a77db", Dwc = "#4c51b5", Dwd = "#320087",
  Dfa = "#00ffff", Dfb = "#37c8ff", Dfc = "#007e7d", Dfd = "#00456e",
  ET = "#b2b2b2", EF = "#686868"
)
present_types <- names(koppen_type_colors)[names(koppen_type_colors) %in% unique(koppen_type[land])]

open_map_png("final_year_koppen_climate_types.png", height = 1400)
par(mar = c(9.6, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Köppen-Geiger types over land, %s", period_label))
draw_cells(koppen_type_colors[koppen_type], land)
map_axes()
draw_coastline()
par(xpd = NA)
legend(x = -180, y = -118, legend = c(present_types, "Ocean"),
       fill = c(koppen_type_colors[present_types], ocean_colour), border = "#555555",
       ncol = 10, bty = "n", cex = 0.7)
par(xpd = FALSE)
mtext("Peel et al. (2007) criteria, i.e. the 0 degC C/D boundary; only the types that occur are listed",
      side = 1, line = 8.3, cex = 0.72)
dev.off()

# --- Terrain ---------------------------------------------------------------
terrain_breaks <- c(-Inf, 250, 500, 1000, 1500, 2000, 2500, Inf)
terrain_colors <- hcl.colors(length(terrain_breaks) - 1, palette = "Terrain 2")
terrain_labels <- c("<250 m", "250-500 m", "500-1000 m", "1000-1500 m",
                    "1500-2000 m", "2000-2500 m", ">=2500 m")
terrain_class <- cut(surface_height, breaks = terrain_breaks, labels = FALSE,
                     include.lowest = TRUE, right = FALSE)

open_map_png("surface_terrain_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame("Surface terrain (smoothed and truncated ETOPO 2022)")
draw_cells(terrain_colors[terrain_class], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(terrain_labels, "Ocean"),
       fill = c(terrain_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.78,
       title = "Surface elevation")
mtext("Land is land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Annual-mean land surface temperature ----------------------------------
temperature_breaks <- c(-Inf, -30, -20, -10, 0, 10, 20, 30, Inf)
temperature_colors <- hcl.colors(length(temperature_breaks) - 1, palette = "Blue-Red 3")
temperature_labels <- c("<-30 °C", "-30 to -20 °C", "-20 to -10 °C", "-10 to 0 °C",
                        "0 to 10 °C", "10 to 20 °C", "20 to 30 °C", ">=30 °C")
temperature_class <- cut(land_temperature_annual_c, breaks = temperature_breaks, labels = FALSE,
                         include.lowest = TRUE, right = FALSE)

open_map_png("final_year_land_surface_temperature_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual-mean land surface temperature, %s", period_label))
draw_cells(temperature_colors[temperature_class], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(temperature_labels, "Ocean"),
       fill = c(temperature_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.74,
       title = "Temperature")
mtext("Land is land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Annual-mean land-bucket water -----------------------------------------
# Fixed bins span the configured 150 kg/m2 capacity and match the analytic
# land/sea analysis, making the two maps directly comparable.
surface_water_breaks <- seq(0, 150, by = 15)
surface_water_colors <- hcl.colors(length(surface_water_breaks) - 1, palette = "YlGnBu")
surface_water_class <- cut(pmin(150, pmax(0, surface_water_annual)),
                           breaks = surface_water_breaks, labels = FALSE, include.lowest = TRUE)
surface_water_labels <- sprintf("%d–%d", head(surface_water_breaks, -1), tail(surface_water_breaks, -1))

open_map_png("final_year_land_surface_water_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual-mean land surface water, %s", period_label))
draw_cells(surface_water_colors[surface_water_class], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(surface_water_labels, "Ocean"),
       fill = c(surface_water_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.70,
       title = expression(W~(kg~m^{-2})))
mtext("Final-year mean; land is land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Annual land precipitation ---------------------------------------------
precipitation_breaks <- c(-Inf, 250, 500, 1000, 1500, 2000, 3000, 4000, Inf)
precipitation_colors <- hcl.colors(length(precipitation_breaks) - 1, palette = "YlGnBu")
precipitation_labels <- c("<250 mm", "250-500 mm", "500-1000 mm", "1000-1500 mm",
                          "1500-2000 mm", "2000-3000 mm", "3000-4000 mm", ">=4000 mm")
precipitation_class <- cut(precipitation_annual_mm, breaks = precipitation_breaks, labels = FALSE,
                           include.lowest = TRUE, right = FALSE)

open_map_png("final_year_land_annual_precipitation_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual land precipitation, %s", period_label))
draw_cells(precipitation_colors[precipitation_class], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(precipitation_labels, "Ocean"),
       fill = c(precipitation_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.74,
       title = "Precipitation (mm/year)")
mtext("Monthly means are rates in mm/day; every model month has 30 days", side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Annual-mean cloud cover -----------------------------------------------
cloud_breaks <- seq(0, 1, by = 0.1)
cloud_colors <- hcl.colors(length(cloud_breaks) - 1, palette = "Blues 3", rev = TRUE)
cloud_class <- cut(pmin(1, pmax(0, cloud_cover_annual)), breaks = cloud_breaks,
                   labels = FALSE, include.lowest = TRUE)
cloud_labels <- sprintf("%d–%d%%", seq(0, 90, 10), seq(10, 100, 10))

open_map_png("final_year_mean_cloud_cover_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual-mean effective cloud cover, %s", period_label), background = NA)
draw_cells(cloud_colors[cloud_class])
map_axes(grid_alpha = 0.55)
draw_coastline(colour = "#a03030", lwd = 0.9)
legend("topleft", legend = cloud_labels, fill = cloud_colors,
       border = "#555555", bg = "white", cex = 0.72, title = "Cloud cover")
mtext("Diagnosed effective column cloud fraction; land and ocean are both shown",
      side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Land water budget P - E ------------------------------------------------
budget_annual_mm <- precipitation_annual_mm - evaporation_annual_mm
budget_breaks <- c(-Inf, -500, -250, -100, 0, 100, 250, 500, 1000, Inf)
budget_colors <- hcl.colors(length(budget_breaks) - 1, palette = "BrBG")
budget_labels <- c("<-500", "-500 to -250", "-250 to -100", "-100 to 0", "0 to 100",
                   "100 to 250", "250 to 500", "500 to 1000", ">=1000")
budget_class <- cut(budget_annual_mm, breaks = budget_breaks, labels = FALSE,
                    include.lowest = TRUE, right = FALSE)

open_map_png("final_year_land_water_budget_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual land P - E, %s", period_label))
draw_cells(budget_colors[budget_class], land)
map_axes()
draw_coastline()
legend("topleft", legend = c(budget_labels, "Ocean"),
       fill = c(budget_colors, ocean_colour), border = "#555555", bg = "white", cex = 0.72,
       title = "P - E (mm/year)")
mtext("Positive values are the runoff that the bucket-free land surface discards",
      side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Monsoon: high-sun minus low-sun season --------------------------------
jja <- which(calendar_month %in% 6:8)
djf <- which(calendar_month %in% c(12, 1, 2))
precipitation_jja <- rowMeans(precipitation_monthly[, jja, drop = FALSE])
precipitation_djf <- rowMeans(precipitation_monthly[, djf, drop = FALSE])
precipitation_season_diff <- precipitation_jja - precipitation_djf
pressure_jja <- rowMeans(surface_pressure_monthly[, jja, drop = FALSE]) / 100
pressure_djf <- rowMeans(surface_pressure_monthly[, djf, drop = FALSE]) / 100
pressure_season_diff <- pressure_jja - pressure_djf

season_precip_breaks <- c(-Inf, -8, -4, -2, -1, -0.25, 0.25, 1, 2, 4, 8, Inf)
season_precip_colors <- hcl.colors(length(season_precip_breaks) - 1, palette = "Blue-Red 3", rev = TRUE)
season_precip_class <- cut(precipitation_season_diff, breaks = season_precip_breaks,
                           labels = FALSE, include.lowest = TRUE, right = FALSE)
season_precip_labels <- c("<-8", "-8 to -4", "-4 to -2", "-2 to -1", "-1 to -0.25",
                          "-0.25 to 0.25", "0.25 to 1", "1 to 2", "2 to 4", "4 to 8", ">=8")
season_pressure_breaks <- c(-Inf, -12, -8, -4, -2, -0.5, 0.5, 2, 4, 8, 12, Inf)
season_pressure_colors <- hcl.colors(length(season_pressure_breaks) - 1, palette = "Purple-Green")
season_pressure_class <- cut(pressure_season_diff, breaks = season_pressure_breaks,
                             labels = FALSE, include.lowest = TRUE, right = FALSE)
season_pressure_labels <- c("<-12", "-12 to -8", "-8 to -4", "-4 to -2", "-2 to -0.5",
                            "-0.5 to 0.5", "0.5 to 2", "2 to 4", "4 to 8", "8 to 12", ">=12")

png(file.path(analysis_dir, "final_year_monsoon_seasonal_contrast.png"),
    width = 2200, height = 2300, res = 180, type = png_device_type)
par(mfrow = c(2, 1), mar = c(4.6, 5.0, 3.4, 1.0), oma = c(0.4, 0, 2.2, 0), las = 1)
map_frame("Precipitation, JJA minus DJF", background = NA)
draw_cells(season_precip_colors[season_precip_class])
map_axes(grid_alpha = 0.5)
draw_coastline(colour = "#222222", lwd = 0.9)
legend("bottomleft", legend = season_precip_labels, fill = season_precip_colors,
       border = "#555555", bg = "white", cex = 0.62, ncol = 2, title = "mm/day")
map_frame("Surface pressure, JJA minus DJF", background = NA)
draw_cells(season_pressure_colors[season_pressure_class])
map_axes(grid_alpha = 0.5)
draw_coastline(colour = "#222222", lwd = 0.9)
legend("bottomleft", legend = season_pressure_labels, fill = season_pressure_colors,
       border = "#555555", bg = "white", cex = 0.62, ncol = 2, title = "hPa")
mtext(sprintf("Seasonal contrast of the %s", period_label), side = 3, outer = TRUE,
      line = 0.4, cex = 1.12)
dev.off()

# --- Ocean surface pressure, ripple check ----------------------------------
ocean <- land_fraction <= 0.01
ocean_pressure_anomaly <- surface_pressure_annual_hpa -
  weighted.mean(surface_pressure_annual_hpa[ocean], point_weight[ocean])

# Ripples from the truncated terrain would appear as a grid-scale zig-zag along
# each longitude ring, so measure the departure from the two ring neighbours.
ring_offset <- cumsum(c(0L, nlon))
zonal_small_scale <- rep(NA_real_, npoints)
for (j in seq_len(nlat)) {
  indices <- ring_offset[[j]] + seq_len(nlon[[j]])
  values <- surface_pressure_annual_hpa[indices]
  left <- c(values[[nlon[[j]]]], values[-nlon[[j]]])
  right <- c(values[-1], values[[1]])
  zonal_small_scale[indices] <- values - 0.5 * (left + right)
}
ripple_breaks <- c(-Inf, -12, -8, -4, -2, -1, 1, 2, 4, 8, 12, Inf)
ripple_colors <- hcl.colors(length(ripple_breaks) - 1, palette = "Blue-Red 3")
ripple_class <- cut(ocean_pressure_anomaly, breaks = ripple_breaks, labels = FALSE,
                    include.lowest = TRUE, right = FALSE)
ripple_labels <- c("<-12", "-12 to -8", "-8 to -4", "-4 to -2", "-2 to -1", "-1 to 1",
                   "1 to 2", "2 to 4", "4 to 8", "8 to 12", ">=12")

open_map_png("final_year_ocean_surface_pressure_map.png")
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
map_frame(sprintf("Annual-mean ocean surface pressure anomaly, %s", period_label),
          background = "#f2f2f2")
draw_cells(ripple_colors[ripple_class], ocean)
map_axes(grid_alpha = 0.35)
draw_coastline(colour = "#333333", lwd = 0.9)
legend("bottomleft", legend = c(ripple_labels, "Land"),
       fill = c(ripple_colors, "#f2f2f2"), border = "#555555", bg = "white", cex = 0.62,
       ncol = 2, title = "hPa from ocean mean")
mtext("Ocean is land fraction <= 0.01; ripples next to high terrain would show up here",
      side = 1, line = 4.6, cex = 0.72)
dev.off()

# --- Zonal means ------------------------------------------------------------
zonal_of <- function(field) rowMeans(vapply(
  seq_len(months_per_year), function(m) ring_mean(field[, m]), numeric(nlat)
))
temperature_zonal_k <- zonal_of(temperature_monthly)
precipitation_zonal_mm_day <- zonal_of(precipitation_monthly)
evaporation_zonal_mm_day <- zonal_of(evaporation_monthly)
cloud_cover_zonal <- zonal_of(cloud_cover_monthly)
zonal_table <- data.frame(
  latitude_deg = latitude,
  annual_mean_surface_temperature_k = temperature_zonal_k,
  annual_mean_precipitation_mm_day = precipitation_zonal_mm_day,
  annual_mean_evaporation_mm_day = evaporation_zonal_mm_day,
  annual_mean_cloud_cover = cloud_cover_zonal
)
write.csv(zonal_table, file.path(analysis_dir, "final_year_zonal_means.csv"), row.names = FALSE)

png(file.path(analysis_dir, "final_year_zonal_temperature_precipitation.png"),
    width = 1800, height = 1400, res = 180, type = png_device_type)
par(mfrow = c(2, 1), mar = c(2.1, 5.2, 3.1, 1.1), oma = c(4.3, 0, 2.4, 0), las = 1)
plot(latitude, temperature_zonal_k - 273.15, type = "l", lwd = 3, col = "#c23b22",
     xlim = c(-90, 90), xaxt = "n", xlab = "", ylab = "Temperature (°C)",
     main = "Annual-mean surface temperature", panel.first = grid(col = "#dddddd"))
abline(h = 0, col = "#3572a5", lty = 2)
plot(latitude, precipitation_zonal_mm_day, type = "l", lwd = 3, col = "#26734d",
     xlim = c(-90, 90), ylim = range(0, precipitation_zonal_mm_day, evaporation_zonal_mm_day),
     xaxt = "n", xlab = "", ylab = expression("Water flux (mm day"^{-1}*")"),
     main = "Annual-mean precipitation and evaporation", panel.first = grid(col = "#dddddd"))
lines(latitude, evaporation_zonal_mm_day, lwd = 2.4, col = "#3572a5", lty = 2)
legend("topright", legend = c("Precipitation", "Evaporation"), col = c("#26734d", "#3572a5"),
       lty = c(1, 2), lwd = c(3, 2.4), bg = "white", cex = 0.8)
axis(1, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
mtext("Latitude", side = 1, outer = TRUE, line = 2.4)
mtext(sprintf("Earth-terrain land-sea simulation, %s", period_label),
      side = 3, outer = TRUE, line = 0.6, cex = 1.12)
dev.off()

# --- Meridional mass streamfunction ----------------------------------------
mass_streamfunction <- matrix(0, nrow = nlat, ncol = nlev)
surface_mass_residual <- matrix(0, nrow = nlat, ncol = months_per_year)
for (m in seq_len(months_per_year)) {
  zonal_surface_pressure <- ring_mean(surface_pressure_monthly[, m])
  delta_pressure <- outer(rep(1, nlat), diff(a_half)) + outer(zonal_surface_pressure, diff(b_half))
  layer_transport <- (2 * pi * earth_radius / gravity) * cos_latitude * zonal_v_monthly[, , m] * delta_pressure
  psi_month <- t(apply(layer_transport, 1, function(x) cumsum(x) - 0.5 * x))
  mass_streamfunction <- mass_streamfunction + psi_month / months_per_year
  surface_mass_residual[, m] <- rowSums(layer_transport)
}

streamfunction_table <- expand.grid(
  latitude_deg = latitude,
  level = seq_len(nlev),
  KEEP.OUT.ATTRS = FALSE
)
streamfunction_table$reference_pressure_hpa <- reference_full_pressure_hpa[streamfunction_table$level]
streamfunction_table$mass_streamfunction_kg_s <- as.vector(mass_streamfunction)
write.csv(streamfunction_table, file.path(analysis_dir, "final_year_mass_streamfunction.csv"), row.names = FALSE)

pressure_order <- rev(seq_len(nlev))
y <- -log10(reference_full_pressure_hpa[pressure_order])
psi_billion <- mass_streamfunction / 1e9
z <- psi_billion[, pressure_order, drop = FALSE]
z_limit <- ceiling(max(abs(z)) / 10) * 10
if (z_limit == 0) z_limit <- 1
contour_levels <- seq(-z_limit, z_limit, length.out = 25)
stream_palette <- hcl.colors(length(contour_levels) - 1, palette = "Blue-Red 3")

png(file.path(analysis_dir, "final_year_mass_streamfunction.png"),
    width = 1900, height = 1200, res = 180, type = png_device_type)
filled.contour(
  x = latitude, y = y, z = z, levels = contour_levels,
  color.palette = function(n) stream_palette,
  xlim = c(-90, 90),
  plot.title = title(
    main = sprintf("Annual-mean meridional mass streamfunction, %s", period_label),
    xlab = "Latitude", ylab = "Reference pressure (hPa)"
  ),
  plot.axes = {
    axis(1, at = seq(-90, 90, 30),
         labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
    pressure_ticks <- c(1000, 850, 700, 500, 300, 200, 100, 50, 10, 2)
    pressure_ticks <- pressure_ticks[
      pressure_ticks >= min(reference_full_pressure_hpa) & pressure_ticks <= max(reference_full_pressure_hpa)
    ]
    axis(2, at = -log10(pressure_ticks), labels = pressure_ticks, las = 1)
    contour(latitude, y, z, levels = contour_levels, add = TRUE, drawlabels = FALSE,
            col = adjustcolor("black", 0.35), lwd = 0.7)
    contour(latitude, y, z, levels = 0, add = TRUE, drawlabels = FALSE,
            col = "black", lwd = 1.6)
  },
  key.title = title(main = expression(10^9 * " kg " * s^{-1}), cex.main = 0.85),
  key.axes = axis(4, las = 1)
)
dev.off()

# --- Climographs ------------------------------------------------------------
# Reorder April-March output into the conventional January-December display.
calendar_order <- order(calendar_month)
selected_temperature_c <- temperatures_c[selected, calendar_order, drop = FALSE]
selected_precipitation_mm <- precipitation_mm_month[selected, calendar_order, drop = FALSE]
temperature_limits <- range(pretty(range(c(selected_temperature_c, 0)), n = 6))
precipitation_limit <- max(pretty(c(0, selected_precipitation_mm), n = 6))

png(file.path(analysis_dir, "final_year_representative_land_climographs.png"),
    width = 2400, height = 2500, res = 180, type = png_device_type)
par(mfrow = c(4, 3), mar = c(3.9, 4.5, 3.6, 4.5), oma = c(1.2, 0.8, 3.2, 0.8), las = 1)
for (s in seq_along(selected)) {
  temperature_c <- selected_temperature_c[s, ]
  precipitation_mm <- selected_precipitation_mm[s, ]

  plot(seq_len(12), temperature_c, type = "n", xlim = c(0.5, 12.5), ylim = temperature_limits,
       xaxt = "n", xlab = "", ylab = "Temperature (°C)",
       main = sprintf("%d. %s — %s", s, site_table$site[[s]], site_table$koppen_type[[s]]),
       panel.first = grid(col = "#e5e5e5"))
  axis(1, at = seq_len(12), labels = substr(month.abb, 1, 1), cex.axis = 0.8)
  abline(h = 0, lty = 3, col = "#777777")

  par(new = TRUE)
  plot(seq_len(12), precipitation_mm, type = "n", xlim = c(0.5, 12.5),
       ylim = c(0, precipitation_limit), axes = FALSE, xlab = "", ylab = "")
  rect(seq_len(12) - 0.31, 0, seq_len(12) + 0.31, pmin(precipitation_mm, precipitation_limit),
       col = adjustcolor("#3572a5", 0.55), border = "#3572a5")
  axis(4, cex.axis = 0.85)
  mtext("Precip. (mm/month)", side = 4, line = 2.7, cex = 0.6, las = 0)

  par(new = TRUE)
  plot(seq_len(12), temperature_c, type = "o", pch = 16, lwd = 2.3, col = "#c23b22",
       xlim = c(0.5, 12.5), ylim = temperature_limits, axes = FALSE, xlab = "", ylab = "")

  mtext(sprintf("%.1f°%s, %.1f°E, %d m  |  %.1f °C, %d mm/yr  |  Earth: %s",
                abs(site_table$latitude_deg[[s]]),
                ifelse(site_table$latitude_deg[[s]] < 0, "S", "N"),
                site_table$longitude_deg[[s]],
                round(site_table$surface_height_m[[s]]),
                site_table$annual_mean_temperature_c[[s]],
                round(site_table$annual_precipitation_mm[[s]]),
                site_table$earth_reference_koppen[[s]]),
        side = 3, line = 0.25, cex = 0.6)
}
mtext(sprintf("Climographs at real-world sites, %s (shown Jan–Dec)", period_label),
      side = 3, outer = TRUE, line = 1.3, cex = 1.15)
mtext("Temperature and precipitation axes are shared across all panels; \"Earth\" is the observed Köppen type at the site",
      side = 1, outer = TRUE, line = 0.1, cex = 0.7)
dev.off()

# --- Spin-up from the daily global diagnostics ------------------------------
daily <- read.csv(file.path(case_dir, "daily_global.csv"), check.names = FALSE)
daily_year <- daily[["simulation_day"]] / 360
toa_imbalance <- daily[["incoming_shortwave_w_m-2"]] - daily[["reflected_shortwave_w_m-2"]] -
  daily[["outgoing_longwave_w_m-2"]]

if (has_sea_ice) {
  required_daily_ice <- c("sea_ice_area_m2", "sea_ice_volume_m3", "mean_sea_ice_thickness_m")
  if (!all(required_daily_ice %in% names(daily))) stop("Missing daily sea-ice diagnostics")
  png(file.path(analysis_dir, "spinup_sea_ice_timeseries.png"),
      width = 1900, height = 1800, res = 180, type = png_device_type)
  par(mfrow = c(3, 1), mar = c(2.3, 5.5, 2.8, 1.2), oma = c(4.1, 0, 2.4, 0), las = 1)
  for (quantity in required_daily_ice) {
    scale <- switch(quantity, sea_ice_area_m2 = 1e12, sea_ice_volume_m3 = 1e13,
                    mean_sea_ice_thickness_m = 1)
    label <- switch(quantity, sea_ice_area_m2 = expression(10^12*" m"^2),
                    sea_ice_volume_m3 = expression(10^13*" m"^3),
                    mean_sea_ice_thickness_m = "m")
    title <- switch(quantity, sea_ice_area_m2 = "Global sea-ice area",
                    sea_ice_volume_m3 = "Global sea-ice volume",
                    mean_sea_ice_thickness_m = "Ice-area-weighted thickness")
    plot(daily_year, daily[[quantity]] / scale, type = "l", lwd = 1.5,
         col = "#2368a2", xlab = "", ylab = label, main = title,
         panel.first = grid(col = "#dddddd"))
  }
  mtext("Simulation year", side = 1, outer = TRUE, line = 2.3)
  mtext("Sea-ice evolution (daily global diagnostics)", side = 3, outer = TRUE,
        line = 0.5, cex = 1.12)
  dev.off()
}

png(file.path(analysis_dir, "spinup_global_timeseries.png"),
    width = 2000, height = 2200, res = 180, type = png_device_type)
par(mfrow = c(4, 1), mar = c(2.2, 5.4, 2.6, 1.2), oma = c(4.0, 0, 2.6, 0), las = 1)
plot(daily_year, daily[["mean_surface_temperature_k"]], type = "l", lwd = 1.6, col = "#c23b22",
     xlab = "", ylab = "Temperature (K)", main = "Global-mean surface temperature",
     ylim = range(daily[["mean_surface_temperature_k"]], daily[["mean_land_surface_temperature_k"]],
                  daily[["mean_ocean_surface_temperature_k"]]),
     panel.first = grid(col = "#dddddd"))
lines(daily_year, daily[["mean_land_surface_temperature_k"]], lwd = 1.4, col = "#8a6d3b")
lines(daily_year, daily[["mean_ocean_surface_temperature_k"]], lwd = 1.4, col = "#3572a5")
legend("bottomright", legend = c("All", "Land", "Ocean"),
       col = c("#c23b22", "#8a6d3b", "#3572a5"), lwd = 1.6, bg = "white", cex = 0.75, ncol = 3)
plot(daily_year, daily[["precipitation_mm_day-1"]], type = "l", lwd = 1.5, col = "#26734d",
     xlab = "", ylab = expression("mm day"^{-1}), main = "Global-mean precipitation and evaporation",
     ylim = range(daily[["precipitation_mm_day-1"]], daily[["evaporation_mm_day-1"]]),
     panel.first = grid(col = "#dddddd"))
lines(daily_year, daily[["evaporation_mm_day-1"]], lwd = 1.4, col = "#3572a5")
legend("bottomright", legend = c("Precipitation", "Evaporation"), col = c("#26734d", "#3572a5"),
       lwd = 1.5, bg = "white", cex = 0.75, ncol = 2)
plot(daily_year, daily[["precipitable_water_kg_m-2"]], type = "l", lwd = 1.6, col = "#3572a5",
     xlab = "", ylab = expression("kg m"^{-2}), main = "Precipitable water",
     panel.first = grid(col = "#dddddd"))
plot(daily_year, toa_imbalance, type = "l", lwd = 1.4, col = "#7048e8",
     xlab = "", ylab = expression("W m"^{-2}), main = "Top-of-atmosphere net radiation",
     panel.first = grid(col = "#dddddd"))
abline(h = 0, lty = 2, col = "#777777")
mtext("Simulation year", side = 1, outer = TRUE, line = 2.2)
mtext("Spin-up of the Earth-terrain land-sea run (daily global means)",
      side = 3, outer = TRUE, line = 0.7, cex = 1.12)
dev.off()

# --- Comparison with the analytic-continent run -----------------------------
# A saved analytic-continent run may predate the sea-ice physics. Mark that
# comparison explicitly so its differences are not attributed to terrain alone.
reference_dir <- if (length(args) >= 3) args[[3]] else "output/moist_land_sea_t31"
reference_zonal_path <- file.path(reference_dir, "analysis", "final_year_zonal_means.csv")
reference_daily_path <- file.path(reference_dir, "daily_global.csv")
if (file.exists(reference_zonal_path) && file.exists(reference_daily_path)) {
  reference_zonal <- read.csv(reference_zonal_path, check.names = FALSE)
  reference_daily <- read.csv(reference_daily_path, check.names = FALSE)
  reference_year <- reference_daily[["simulation_day"]] / 360
  reference_label <- basename(reference_dir)
  reference_metadata_path <- file.path(reference_dir, "metadata.json")
  reference_metadata <- if (file.exists(reference_metadata_path))
    paste(readLines(reference_metadata_path, warn = FALSE), collapse = "\n") else ""
  reference_has_sea_ice <- grepl('"sea_ice"\\s*:\\s*\\{\\s*"enabled"\\s*:\\s*true',
                                  reference_metadata, perl = TRUE)
  comparison_note <- if (has_sea_ice && !reference_has_sea_ice)
    "reference lacks sea ice; physics also differs" else
    "same sea-ice setting"

  png(file.path(analysis_dir, "comparison_with_analytic_continents.png"),
      width = 2000, height = 1800, res = 180, type = png_device_type)
  par(mfrow = c(2, 2), mar = c(4.3, 5.0, 3.2, 1.2), oma = c(0.6, 0, 3.0, 0), las = 1)

  plot(latitude, temperature_zonal_k - 273.15, type = "l", lwd = 2.6, col = "#c23b22",
       xlim = c(-90, 90), xaxt = "n", xlab = "Latitude", ylab = "Temperature (°C)",
       main = "Zonal-mean surface temperature",
       ylim = range(temperature_zonal_k, reference_zonal$annual_mean_surface_temperature_k) - 273.15,
       panel.first = grid(col = "#dddddd"))
  lines(reference_zonal$latitude_deg, reference_zonal$annual_mean_surface_temperature_k - 273.15,
        lwd = 2.2, col = "#555555", lty = 2)
  axis(1, at = seq(-90, 90, 30), labels = c("90S", "60S", "30S", "0", "30N", "60N", "90N"))
  legend("bottom", legend = c("Earth terrain", "Analytic continents"),
         col = c("#c23b22", "#555555"), lty = c(1, 2), lwd = 2.4, bg = "white", cex = 0.75)

  plot(latitude, precipitation_zonal_mm_day, type = "l", lwd = 2.6, col = "#26734d",
       xlim = c(-90, 90), xaxt = "n", xlab = "Latitude",
       ylab = expression("Precipitation (mm day"^{-1}*")"),
       main = "Zonal-mean precipitation",
       ylim = range(0, precipitation_zonal_mm_day, reference_zonal$annual_mean_precipitation_mm_day),
       panel.first = grid(col = "#dddddd"))
  lines(reference_zonal$latitude_deg, reference_zonal$annual_mean_precipitation_mm_day,
        lwd = 2.2, col = "#555555", lty = 2)
  axis(1, at = seq(-90, 90, 30), labels = c("90S", "60S", "30S", "0", "30N", "60N", "90N"))

  plot(daily_year, daily[["mean_surface_temperature_k"]], type = "l", lwd = 1.5, col = "#c23b22",
       xlab = "Simulation year", ylab = "Temperature (K)",
       main = "Global-mean surface temperature",
       ylim = range(daily[["mean_surface_temperature_k"]],
                    reference_daily[["mean_surface_temperature_k"]]),
       panel.first = grid(col = "#dddddd"))
  lines(reference_year, reference_daily[["mean_surface_temperature_k"]], lwd = 1.3, col = "#555555")

  reference_toa <- reference_daily[["incoming_shortwave_w_m-2"]] -
    reference_daily[["reflected_shortwave_w_m-2"]] - reference_daily[["outgoing_longwave_w_m-2"]]
  plot(daily_year, toa_imbalance, type = "l", lwd = 1.3, col = "#7048e8",
       xlab = "Simulation year", ylab = expression("W m"^{-2}),
       main = "Top-of-atmosphere net radiation",
       ylim = range(toa_imbalance, reference_toa), panel.first = grid(col = "#dddddd"))
  lines(reference_year, reference_toa, lwd = 1.2, col = "#555555")
  abline(h = 0, lty = 2, col = "#777777")

  mtext(sprintf("Earth terrain (colour) against %s (grey); %s", reference_label, comparison_note),
        side = 3, outer = TRUE, line = 0.6, cex = 0.9)
  dev.off()
} else {
  reference_zonal <- NULL
  message("Reference run not found, skipping the comparison figure: ", reference_dir)
}


# --- Summary ----------------------------------------------------------------
land_weight <- point_weight * land_fraction
ocean_weight <- point_weight * (1 - land_fraction)
group_fraction <- vapply(names(koppen_colors), function(group) {
  sum(land_weight[koppen_group == group], na.rm = TRUE) / sum(land_weight[land])
}, numeric(1))

area_mean <- function(values, weights) sum(weights * values) / sum(weights)
monsoon_box <- function(lon_min, lon_max, lat_min, lat_max) {
  which(land & grid_latitude >= lat_min & grid_latitude <= lat_max &
          grid_longitude >= lon_min & grid_longitude <= lon_max)
}
south_asia <- monsoon_box(70, 100, 5, 30)
east_asia <- monsoon_box(100, 125, 20, 45)
amazon <- monsoon_box(290, 310, -15, 5)
tibet_upwind <- monsoon_box(80, 95, 20, 28)
tibet_lee <- monsoon_box(80, 95, 36, 45)
andes_west <- monsoon_box(283, 289, -35, -20)
andes_east <- monsoon_box(293, 302, -35, -20)
box_precipitation <- function(indices, columns) {
  area_mean(rowMeans(precipitation_monthly[indices, columns, drop = FALSE]) * 360,
            point_weight[indices])
}

group_change <- sum(point_weight[land & koppen_group != koppen_group_zero_boundary]) /
  sum(point_weight[land])

summary_table <- data.frame(
  metric = c(
    "final_month_start", "final_month_end",
    "global_land_area_fraction",
    "maximum_streamfunction_1e9_kg_s", "minimum_streamfunction_1e9_kg_s",
    "maximum_surface_mass_residual_1e9_kg_s",
    "global_area_mean_cloud_cover", "land_area_mean_cloud_cover", "ocean_area_mean_cloud_cover",
    "global_area_mean_surface_temperature_c",
    "land_area_mean_surface_temperature_c", "ocean_area_mean_surface_temperature_c",
    "land_area_mean_surface_water_kg_m-2",
    "land_precipitation_minus_evaporation_mm_yr",
    "ocean_precipitation_minus_evaporation_mm_yr",
    "minimum_land_annual_temperature_c", "maximum_land_annual_temperature_c",
    "antarctic_land_annual_temperature_c", "arctic_land_annual_temperature_c",
    "ocean_surface_pressure_anomaly_max_abs_hpa",
    "ocean_surface_pressure_anomaly_sd_hpa",
    "ocean_surface_pressure_grid_scale_max_abs_hpa",
    "south_asia_jja_minus_djf_precipitation_mm_yr",
    "east_asia_jja_minus_djf_precipitation_mm_yr",
    "amazon_djf_minus_jja_precipitation_mm_yr",
    "tibet_upwind_annual_precipitation_mm_yr", "tibet_lee_annual_precipitation_mm_yr",
    "andes_west_annual_precipitation_mm_yr", "andes_east_annual_precipitation_mm_yr",
    "land_fraction_changing_group_with_0c_boundary",
    paste0("land_area_fraction_group_", names(group_fraction))
  ),
  value = c(
    min(final_months), max(final_months),
    sum(point_weight * land_fraction),
    max(psi_billion), min(psi_billion),
    max(abs(rowMeans(surface_mass_residual))) / 1e9,
    area_mean(cloud_cover_annual, point_weight),
    area_mean(cloud_cover_annual, land_weight),
    area_mean(cloud_cover_annual, ocean_weight),
    area_mean(temperature_annual_c, point_weight),
    area_mean(land_temperature_annual_c, land_weight),
    area_mean(ocean_temperature_annual_c, ocean_weight),
    area_mean(surface_water_annual, land_weight),
    area_mean(budget_annual_mm, land_weight),
    area_mean(budget_annual_mm, ocean_weight),
    min(land_temperature_annual_c[land]), max(land_temperature_annual_c[land]),
    area_mean(land_temperature_annual_c[land & grid_latitude < -65],
              point_weight[land & grid_latitude < -65]),
    area_mean(land_temperature_annual_c[land & grid_latitude > 65],
              point_weight[land & grid_latitude > 65]),
    max(abs(ocean_pressure_anomaly[ocean])),
    sqrt(area_mean((ocean_pressure_anomaly[ocean] -
                      area_mean(ocean_pressure_anomaly[ocean], point_weight[ocean]))^2,
                   point_weight[ocean])),
    max(abs(zonal_small_scale[ocean])),
    box_precipitation(south_asia, jja) - box_precipitation(south_asia, djf),
    box_precipitation(east_asia, jja) - box_precipitation(east_asia, djf),
    box_precipitation(amazon, djf) - box_precipitation(amazon, jja),
    box_precipitation(tibet_upwind, seq_len(months_per_year)),
    box_precipitation(tibet_lee, seq_len(months_per_year)),
    box_precipitation(andes_west, seq_len(months_per_year)),
    box_precipitation(andes_east, seq_len(months_per_year)),
    group_change,
    group_fraction
  )
)
if (has_sea_ice) {
  ice_summary <- data.frame(
    metric = c("final_year_mean_sea_ice_area_m2", "final_year_mean_sea_ice_volume_m3",
               "final_year_mean_sea_ice_thickness_m", "final_year_max_sea_ice_area_m2"),
    value = c(mean(ice_monthly$ice_area_m2[ice_monthly$region == "global"]),
              mean(ice_monthly$ice_volume_m3[ice_monthly$region == "global"]),
              sum(ocean_weight * ice_annual_volume) / sum(ocean_weight * ice_annual_fraction),
              max(ice_monthly$ice_area_m2[ice_monthly$region == "global"]))
  )
  summary_table <- rbind(summary_table, ice_summary)
}
write.csv(summary_table, file.path(analysis_dir, "analysis_summary.csv"), row.names = FALSE)

notes <- c(
  "Final-year climate analysis of the Earth-terrain land/sea run",
  sprintf("Input months: %d-%d (12 equal 30-day months; April through March).", min(final_months), max(final_months)),
  "Figures mirror scripts/analyze_moist_land_sea.R so the analytic-continent run can be compared directly.",
  "Koppen group map: first letters A-E with the original -3 degC C/D boundary; land is land_fraction >= 0.5.",
  "Koppen type map: full Peel et al. (2007) criteria, which use the 0 degC C/D boundary; the summary lists how much land changes group between the two boundaries.",
  "The B threshold uses annual temperature, annual precipitation, and the hemisphere-specific high-sun-half-year precipitation fraction.",
  "Climograph sites are real-world locations snapped to the nearest cell with land_fraction >= 0.5; the snap distance is in representative_land_sites.csv.",
  "The observed Earth Koppen type of each site is listed for reference only; it is not a model output.",
  "Zonal means include both land and ocean grid cells.",
  "The cloud map shows final-year mean diagnosed effective column cloud fraction over land and ocean in fixed 10-percentage-point bins.",
  "The surface-water map shows final-year mean land-bucket water in fixed 15 kg/m2 bins spanning the 150 kg/m2 capacity; ocean is masked.",
  "P - E over land is the runoff discarded by the land surface; over the ocean it is the net moisture source.",
  "The monsoon figure uses June-August minus December-February monthly means of precipitation and surface pressure.",
  "The ocean surface pressure map shows the departure of the annual mean from the ocean area mean, as the check for terrain-induced ripples.",
  "Mass streamfunction uses monthly zonal-mean v and monthly zonal-mean surface pressure; pressure shown is the reference full-level pressure.",
  "The spin-up figure uses all daily global means in daily_global.csv, not only the final year.",
  "The analytic-continent comparison is not terrain-only when the reference run predates sea ice; the figure states this when detected.",
  "The snow and sea-ice map shows March and September monthly means: land (land_fraction >= 0.5) by snow cover f = S/(S+S_0) with the equivalent snow water S, other cells by sea-ice concentration; cover below 1% counts as none. It is written only when the case has snow output.",
  "All annual means use the final 12 monthly means with equal weights."
)
if (has_sea_ice) notes <- c(notes,
  "Sea-ice seasonal maps show March and September monthly mean concentration A per ocean area, with land masked.",
  "Sea-ice area and volume are integrated over the native Gaussian cell areas multiplied by ocean fraction; the seasonal CSV separates hemispheres.",
  "Mean ice thickness is mean V divided by mean A; ice skin temperature is weighted by monthly ice area and time. Both maps mask ice-free water.",
  "The sea-ice spin-up figure uses daily global area, volume, and ice-area-weighted thickness.")
writeLines(notes, file.path(analysis_dir, "README.txt"))

cat("Analysis written to", analysis_dir, "\n")
print(site_table[, c("site", "longitude_deg", "latitude_deg", "snap_distance_deg",
                     "surface_height_m", "koppen_type", "earth_reference_koppen",
                     "annual_mean_temperature_c", "annual_precipitation_mm",
                     "annual_mean_surface_water_kg_m2")], row.names = FALSE)
print(summary_table, row.names = FALSE)
