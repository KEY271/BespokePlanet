#!/usr/bin/env Rscript

# Produce final-year climate diagnostics for the mixed land/sea moist run.
# The script uses base R only. Binary layouts and units are defined by the
# case's metadata.json.

args <- commandArgs(trailingOnly = TRUE)
case_dir <- if (length(args) >= 1) args[[1]] else "output/moist_land_sea_t31"
analysis_dir <- if (length(args) >= 2) args[[2]] else file.path(case_dir, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
png_device_type <- if (capabilities("aqua")) "quartz" else "cairo"

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
precipitation_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
cloud_cover_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
surface_pressure_monthly <- matrix(NA_real_, nrow = npoints, ncol = months_per_year)
zonal_v_monthly <- array(NA_real_, dim = c(nlat, nlev, months_per_year))

for (m in seq_along(final_months)) {
  suffix <- sprintf("m%04d.bin", final_months[[m]])
  temperature_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_surface_temperature_", suffix)))
  precipitation_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_precipitation_", suffix)))
  cloud_cover_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_cloud_cover_", suffix)))
  surface_pressure_monthly[, m] <- read_grid(file.path(case_dir, paste0("monthly_surface_pressure_", suffix)))
  zonal_v_monthly[, , m] <- read_zonal(file.path(case_dir, paste0("monthly_zonal_v_", suffix)))
}

land_fraction <- read_grid(file.path(case_dir, "land_fraction.bin"))
land <- land_fraction >= 0.5
temperature_annual_k <- rowMeans(temperature_monthly)
temperature_annual_c <- temperature_annual_k - 273.15
precipitation_annual_mm <- rowSums(precipitation_monthly * days_per_month)
precipitation_mean_mm_day <- rowMeans(precipitation_monthly)
cloud_cover_annual <- rowMeans(cloud_cover_monthly)

# Köppen first-letter climate group. We use the original -3 C C/D boundary.
# B is tested first using annual temperature, annual precipitation, and the
# hemisphere-specific high-sun half-year precipitation fraction.
classify_koppen_group <- function(temperatures_k, precipitation_mm_day, latitudes) {
  temperatures_c <- temperatures_k - 273.15
  annual_temperature <- rowMeans(temperatures_c)
  annual_precipitation <- rowSums(precipitation_mm_day * days_per_month)
  warmest <- apply(temperatures_c, 1, max)
  coldest <- apply(temperatures_c, 1, min)

  summer_fraction <- numeric(nrow(temperatures_c))
  for (i in seq_len(nrow(temperatures_c))) {
    summer_calendar <- if (latitudes[[i]] >= 0) 4:9 else c(10:12, 1:3)
    summer_columns <- which(calendar_month %in% summer_calendar)
    total <- annual_precipitation[[i]]
    summer_fraction[[i]] <- if (total > 0) {
      sum(precipitation_mm_day[i, summer_columns] * days_per_month) / total
    } else {
      0.5
    }
  }

  dryness_adjustment <- ifelse(
    summer_fraction >= 0.7, 280,
    ifelse(summer_fraction <= 0.3, 0, 140)
  )
  dryness_threshold <- 20 * annual_temperature + dryness_adjustment

  group <- rep(NA_character_, nrow(temperatures_c))
  group[annual_precipitation < dryness_threshold] <- "B"
  unset <- is.na(group)
  group[unset & apply(temperatures_c, 1, min) >= 18] <- "A"
  unset <- is.na(group)
  group[unset & warmest < 10] <- "E"
  unset <- is.na(group)
  group[unset & coldest > -3 & warmest >= 10] <- "C"
  group[is.na(group)] <- "D"
  group
}

koppen_group <- classify_koppen_group(temperature_monthly, precipitation_monthly, grid_latitude)
koppen_group[!land] <- NA_character_

climate_table <- data.frame(
  longitude_deg = grid_longitude,
  latitude_deg = grid_latitude,
  land_fraction = land_fraction,
  annual_mean_surface_temperature_c = temperature_annual_c,
  annual_precipitation_mm = precipitation_annual_mm,
  annual_mean_cloud_cover = cloud_cover_annual,
  koppen_group = koppen_group
)
write.csv(climate_table, file.path(analysis_dir, "final_year_koppen_groups.csv"), row.names = FALSE)

# Representative land points, chosen near interiors of the three analytic
# continents and the southern polar cap. Snap to the nearest native cell with
# land fraction >= 0.75 using great-circle distance. These are defined before
# plotting so the same exact points appear on both the map and climographs.
target_points <- data.frame(
  site = c("A north interior", "A south interior", "B equatorial interior",
           "B south interior", "C interior", "South polar cap"),
  target_lon = c(45, 45, 280, 280, 150, 180),
  target_lat = c(50, 30, -5, -30, -30, -80)
)
candidate <- which(land_fraction >= 0.75)
selected <- integer(nrow(target_points))
for (s in seq_len(nrow(target_points))) {
  target_lon_rad <- target_points$target_lon[[s]] * pi / 180
  target_lat_rad <- target_points$target_lat[[s]] * pi / 180
  lon_rad <- grid_longitude[candidate] * pi / 180
  lat_rad <- grid_latitude[candidate] * pi / 180
  cosine_distance <- sin(target_lat_rad) * sin(lat_rad) +
    cos(target_lat_rad) * cos(lat_rad) * cos(lon_rad - target_lon_rad)
  selected[[s]] <- candidate[[which.max(cosine_distance)]]
}

site_table <- data.frame(
  site = target_points$site,
  longitude_deg = grid_longitude[selected],
  latitude_deg = grid_latitude[selected],
  land_fraction = land_fraction[selected],
  koppen_group = koppen_group[selected],
  annual_mean_temperature_c = temperature_annual_c[selected],
  annual_precipitation_mm = precipitation_annual_mm[selected]
)
write.csv(site_table, file.path(analysis_dir, "representative_land_sites.csv"), row.names = FALSE)

# Requested map: draw each native octahedral-grid cell with the prime meridian
# at the centre, splitting cells that straddle the -180/180-degree seam.
koppen_colors <- c(A = "#2f9e44", B = "#d8b365", C = "#ffd43b", D = "#4dabf7", E = "#f1f3f5")
koppen_labels <- c(A = "A Tropical", B = "B Dry", C = "C Temperate", D = "D Continental", E = "E Polar")
latitude_edges <- c(-90, 0.5 * (latitude[-1] + latitude[-nlat]), 90)
map_longitude <- ifelse(grid_longitude >= 180, grid_longitude - 360, grid_longitude)

png(file.path(analysis_dir, "final_year_koppen_climate_zones.png"),
    width = 2200, height = 1250, res = 180, type = png_device_type)
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
     xlab = "Longitude", ylab = "Latitude",
     main = sprintf("Köppen climate groups over land, final year (months %d–%d)",
                    min(final_months), max(final_months)), axes = FALSE)
rect(-180, -90, 180, 90, col = "#b9d9eb", border = NA)
offset <- 0L
for (j in seq_len(nlat)) {
  delta_lon <- 360 / nlon[[j]]
  for (i in seq_len(nlon[[j]])) {
    index <- offset + i
    if (!land[[index]]) next
    centre <- map_longitude[[index]]
    left <- centre - delta_lon / 2
    right <- centre + delta_lon / 2
    colour <- koppen_colors[[koppen_group[[index]]]]
    if (left < -180) {
      rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
    } else if (right > 180) {
      rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = colour, border = NA)
    } else {
      rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
    }
  }
  offset <- offset + nlon[[j]]
}
abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60), col = adjustcolor("white", 0.65), lwd = 0.7)
axis(1, at = seq(-180, 180, 60),
     labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
axis(2, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
box()
legend("topleft", legend = c(koppen_labels, "Ocean"),
       fill = c(koppen_colors, "#b9d9eb"), border = "#555555", bg = "white", cex = 0.82)

# Place labels manually over ocean and connect them to their land points.
site_map_longitude <- ifelse(site_table$longitude_deg > 180,
                             site_table$longitude_deg - 360,
                             site_table$longitude_deg)
# Keep the marker for the point exactly on the 180-degree seam fully visible.
site_map_longitude[abs(site_map_longitude - 180) < 1e-10] <- 178.5
label_x <- c(150, 158, -142, -142, 92, 112)
label_y <- c(78, 13, 18, -53, -53, -62)
label_text <- c("1  A north interior", "2  A south interior", "3  B equatorial interior",
                "4  B south interior", "5  C interior", "6  South polar cap")
segments(site_map_longitude, site_table$latitude_deg, label_x, label_y,
         col = "#303030", lwd = 1.1)
points(site_map_longitude, site_table$latitude_deg, pch = 21, cex = 1.45,
       bg = "white", col = "#111111", lwd = 1.4)
text(site_map_longitude, site_table$latitude_deg, labels = seq_len(nrow(site_table)),
     cex = 0.54, font = 2)
for (s in seq_along(label_text)) {
  width <- strwidth(label_text[[s]], cex = 0.67)
  height <- strheight(label_text[[s]], cex = 0.67)
  rect(label_x[[s]] - width / 2 - 1.3, label_y[[s]] - height / 2 - 0.8,
       label_x[[s]] + width / 2 + 1.3, label_y[[s]] + height / 2 + 0.8,
       col = adjustcolor("white", 0.92), border = "#555555")
  text(label_x[[s]], label_y[[s]], label_text[[s]], cex = 0.67)
}
mtext("Land is defined as land fraction >= 0.5; C/D boundary is -3 degC",
      side = 1, line = 4.6, cex = 0.72)
dev.off()

# Companion terrain map without site markers or climate classification.
surface_height <- read_grid(file.path(case_dir, "surface_height.bin"))
terrain_breaks <- c(-Inf, 250, 500, 1000, 1500, 2000, 2500, Inf)
terrain_colors <- hcl.colors(length(terrain_breaks) - 1, palette = "Terrain 2")
terrain_labels <- c("<250 m", "250-500 m", "500-1000 m", "1000-1500 m",
                    "1500-2000 m", "2000-2500 m", ">=2500 m")
terrain_class <- cut(surface_height, breaks = terrain_breaks, labels = FALSE,
                     include.lowest = TRUE, right = FALSE)

png(file.path(analysis_dir, "surface_terrain_map.png"),
    width = 2200, height = 1250, res = 180, type = png_device_type)
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
     xlab = "Longitude", ylab = "Latitude", main = "Surface terrain",
     axes = FALSE)
rect(-180, -90, 180, 90, col = "#b9d9eb", border = NA)
offset <- 0L
for (j in seq_len(nlat)) {
  delta_lon <- 360 / nlon[[j]]
  for (i in seq_len(nlon[[j]])) {
    index <- offset + i
    if (!land[[index]]) next
    centre <- map_longitude[[index]]
    left <- centre - delta_lon / 2
    right <- centre + delta_lon / 2
    colour <- terrain_colors[[terrain_class[[index]]]]
    if (left < -180) {
      rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
    } else if (right > 180) {
      rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = colour, border = NA)
    } else {
      rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
    }
  }
  offset <- offset + nlon[[j]]
}
abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60),
       col = adjustcolor("white", 0.65), lwd = 0.7)
axis(1, at = seq(-180, 180, 60),
     labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
axis(2, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
box()
legend("topleft", legend = c(terrain_labels, "Ocean"),
       fill = c(terrain_colors, "#b9d9eb"), border = "#555555", bg = "white", cex = 0.78,
       title = "Surface elevation")
mtext("Land is defined as land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# Final-year annual-mean land surface temperature map.
temperature_breaks <- c(-Inf, -30, -20, -10, 0, 10, 20, 30, Inf)
temperature_colors <- hcl.colors(length(temperature_breaks) - 1, palette = "Blue-Red 3")
temperature_labels <- c("<-30 °C", "-30 to -20 °C", "-20 to -10 °C", "-10 to 0 °C",
                        "0 to 10 °C", "10 to 20 °C", "20 to 30 °C", ">=30 °C")
temperature_class <- cut(temperature_annual_c, breaks = temperature_breaks, labels = FALSE,
                         include.lowest = TRUE, right = FALSE)

png(file.path(analysis_dir, "final_year_land_surface_temperature_map.png"),
    width = 2200, height = 1250, res = 180, type = png_device_type)
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
     xlab = "Longitude", ylab = "Latitude",
     main = sprintf("Annual-mean land surface temperature, final year (months %d–%d)",
                    min(final_months), max(final_months)),
     axes = FALSE)
rect(-180, -90, 180, 90, col = "#b9d9eb", border = NA)
offset <- 0L
for (j in seq_len(nlat)) {
  delta_lon <- 360 / nlon[[j]]
  for (i in seq_len(nlon[[j]])) {
    index <- offset + i
    if (!land[[index]]) next
    centre <- map_longitude[[index]]
    left <- centre - delta_lon / 2
    right <- centre + delta_lon / 2
    colour <- temperature_colors[[temperature_class[[index]]]]
    if (left < -180) {
      rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
    } else if (right > 180) {
      rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = colour, border = NA)
    } else {
      rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
    }
  }
  offset <- offset + nlon[[j]]
}
abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60),
       col = adjustcolor("white", 0.65), lwd = 0.7)
axis(1, at = seq(-180, 180, 60),
     labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
axis(2, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
box()
legend("topleft", legend = c(temperature_labels, "Ocean"),
       fill = c(temperature_colors, "#b9d9eb"), border = "#555555", bg = "white", cex = 0.74,
       title = "Temperature")
mtext("Land is defined as land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# Final-year annual land precipitation map. Monthly output is a mean rate in
# mm/day, and every model month has 30 days, so their sum gives mm/year.
precipitation_breaks <- c(-Inf, 250, 500, 1000, 1500, 2000, 3000, 4000, Inf)
precipitation_colors <- hcl.colors(length(precipitation_breaks) - 1, palette = "YlGnBu")
precipitation_labels <- c("<250 mm", "250-500 mm", "500-1000 mm", "1000-1500 mm",
                          "1500-2000 mm", "2000-3000 mm", "3000-4000 mm", ">=4000 mm")
precipitation_class <- cut(precipitation_annual_mm, breaks = precipitation_breaks, labels = FALSE,
                           include.lowest = TRUE, right = FALSE)

png(file.path(analysis_dir, "final_year_land_annual_precipitation_map.png"),
    width = 2200, height = 1250, res = 180, type = png_device_type)
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
     xlab = "Longitude", ylab = "Latitude",
     main = sprintf("Annual land precipitation, final year (months %d–%d)",
                    min(final_months), max(final_months)),
     axes = FALSE)
rect(-180, -90, 180, 90, col = "#b9d9eb", border = NA)
offset <- 0L
for (j in seq_len(nlat)) {
  delta_lon <- 360 / nlon[[j]]
  for (i in seq_len(nlon[[j]])) {
    index <- offset + i
    if (!land[[index]]) next
    centre <- map_longitude[[index]]
    left <- centre - delta_lon / 2
    right <- centre + delta_lon / 2
    colour <- precipitation_colors[[precipitation_class[[index]]]]
    if (left < -180) {
      rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
    } else if (right > 180) {
      rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = colour, border = NA)
    } else {
      rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
    }
  }
  offset <- offset + nlon[[j]]
}
abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60),
       col = adjustcolor("white", 0.65), lwd = 0.7)
axis(1, at = seq(-180, 180, 60),
     labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
axis(2, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
box()
legend("topleft", legend = c(precipitation_labels, "Ocean"),
       fill = c(precipitation_colors, "#b9d9eb"), border = "#555555", bg = "white", cex = 0.74,
       title = "Precipitation (mm/year)")
mtext("Land is defined as land fraction >= 0.5", side = 1, line = 4.6, cex = 0.72)
dev.off()

# Final-year annual-mean effective cloud cover over the whole globe. Cloud
# cover is a diagnosed column area fraction in [0, 1], so use fixed 10%-wide
# bins to make this map directly comparable with future runs.
cloud_breaks <- seq(0, 1, by = 0.1)
cloud_colors <- hcl.colors(length(cloud_breaks) - 1, palette = "Blues 3", rev = TRUE)
cloud_class <- cut(pmin(1, pmax(0, cloud_cover_annual)), breaks = cloud_breaks,
                   labels = FALSE, include.lowest = TRUE)
cloud_labels <- sprintf("%d–%d%%", seq(0, 90, 10), seq(10, 100, 10))

png(file.path(analysis_dir, "final_year_mean_cloud_cover_map.png"),
    width = 2200, height = 1250, res = 180, type = png_device_type)
par(mar = c(6.2, 5.0, 4.1, 1.0), las = 1)
plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
     xlab = "Longitude", ylab = "Latitude",
     main = sprintf("Annual-mean effective cloud cover, final year (months %d–%d)",
                    min(final_months), max(final_months)),
     axes = FALSE)
offset <- 0L
for (j in seq_len(nlat)) {
  delta_lon <- 360 / nlon[[j]]
  for (i in seq_len(nlon[[j]])) {
    index <- offset + i
    centre <- map_longitude[[index]]
    left <- centre - delta_lon / 2
    right <- centre + delta_lon / 2
    colour <- cloud_colors[[cloud_class[[index]]]]
    if (left < -180) {
      rect(-180, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(360 + left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
    } else if (right > 180) {
      rect(left, latitude_edges[[j]], 180, latitude_edges[[j + 1]], col = colour, border = NA)
      rect(-180, latitude_edges[[j]], right - 360, latitude_edges[[j + 1]], col = colour, border = NA)
    } else {
      rect(left, latitude_edges[[j]], right, latitude_edges[[j + 1]], col = colour, border = NA)
    }
  }
  offset <- offset + nlon[[j]]
}
abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60),
       col = adjustcolor("white", 0.55), lwd = 0.7)
axis(1, at = seq(-180, 180, 60),
     labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
axis(2, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
box()
legend("topleft", legend = cloud_labels, fill = cloud_colors,
       border = "#555555", bg = "white", cex = 0.72, title = "Cloud cover")
mtext("Diagnosed effective column cloud fraction; land and ocean are both shown",
      side = 1, line = 4.6, cex = 0.72)
dev.off()

# Requested zonal means, including land and ocean in every longitude ring.
temperature_zonal_k <- rowMeans(vapply(
  seq_len(months_per_year), function(m) ring_mean(temperature_monthly[, m]), numeric(nlat)
))
precipitation_zonal_mm_day <- rowMeans(vapply(
  seq_len(months_per_year), function(m) ring_mean(precipitation_monthly[, m]), numeric(nlat)
))
cloud_cover_zonal <- rowMeans(vapply(
  seq_len(months_per_year), function(m) ring_mean(cloud_cover_monthly[, m]), numeric(nlat)
))
zonal_table <- data.frame(
  latitude_deg = latitude,
  annual_mean_surface_temperature_k = temperature_zonal_k,
  annual_mean_precipitation_mm_day = precipitation_zonal_mm_day,
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
     xlim = c(-90, 90), xaxt = "n", xlab = "", ylab = expression("Precipitation (mm day"^{-1}*")"),
     main = "Annual-mean precipitation", panel.first = grid(col = "#dddddd"))
axis(1, at = seq(-90, 90, 30),
     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
mtext("Latitude", side = 1, outer = TRUE, line = 2.4)
mtext(sprintf("Mixed land–sea simulation, final year (months %d–%d)",
              min(final_months), max(final_months)), side = 3, outer = TRUE, line = 0.6, cex = 1.12)
dev.off()

# Requested annual-mean meridional mass streamfunction. Monthly mean v and
# monthly mean surface pressure are the highest-frequency saved inputs.
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
    main = sprintf("Annual-mean meridional mass streamfunction (months %d–%d)",
                   min(final_months), max(final_months)),
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

# Reorder April–March output into the conventional January–December display.
calendar_order <- order(calendar_month)
selected_temperature_c <- temperature_monthly[selected, calendar_order, drop = FALSE] - 273.15
selected_precipitation_mm <- precipitation_monthly[selected, calendar_order, drop = FALSE] * days_per_month
temperature_limits <- range(pretty(range(c(selected_temperature_c, 0)), n = 6))
precipitation_limit <- max(pretty(c(0, selected_precipitation_mm), n = 6))
png(file.path(analysis_dir, "final_year_representative_land_climographs.png"),
    width = 2200, height = 1900, res = 180, type = png_device_type)
par(mfrow = c(3, 2), mar = c(4.1, 4.7, 3.7, 4.7), oma = c(1.0, 0.8, 3.0, 0.8), las = 1)
for (s in seq_along(selected)) {
  index <- selected[[s]]
  temperature_c <- temperature_monthly[index, calendar_order] - 273.15
  precipitation_mm <- precipitation_monthly[index, calendar_order] * days_per_month

  plot(seq_len(12), temperature_c, type = "n", xlim = c(0.5, 12.5), ylim = temperature_limits,
       xaxt = "n", xlab = "Month", ylab = "Temperature (°C)",
       main = sprintf("%s — %s", site_table$site[[s]], site_table$koppen_group[[s]]),
       panel.first = grid(col = "#e5e5e5"))
  axis(1, at = seq_len(12), labels = month.abb, cex.axis = 0.78)
  lines(seq_len(12), temperature_c, type = "o", pch = 16, lwd = 2.3, col = "#c23b22")
  abline(h = 0, lty = 3, col = "#777777")

  par(new = TRUE)
  plot(seq_len(12), precipitation_mm, type = "n", xlim = c(0.5, 12.5), ylim = c(0, precipitation_limit),
       axes = FALSE, xlab = "", ylab = "")
  rect(seq_len(12) - 0.31, 0, seq_len(12) + 0.31, precipitation_mm,
       col = adjustcolor("#3572a5", 0.58), border = "#3572a5")
  axis(4)
  mtext("Precip. (mm/month)", side = 4, line = 2.8, cex = 0.66, las = 0)
  mtext(sprintf("%.1f°%s, %.1f°E  |  Tmean %.1f °C, P %d mm/yr",
                abs(site_table$latitude_deg[[s]]),
                ifelse(site_table$latitude_deg[[s]] < 0, "S", "N"),
                site_table$longitude_deg[[s]],
                site_table$annual_mean_temperature_c[[s]],
                round(site_table$annual_precipitation_mm[[s]])),
        side = 3, line = 0.25, cex = 0.72)
}
mtext(sprintf("Representative land-point climographs, final year (months %d–%d; shown Jan–Dec)",
              min(final_months), max(final_months)), side = 3, outer = TRUE, line = 1.2, cex = 1.15)
mtext("Temperature and precipitation axes are shared across all six panels",
      side = 1, outer = TRUE, line = 0.0, cex = 0.72)
dev.off()

# Compact summary and methodological record alongside the figures.
group_weight <- numeric(npoints)
# Gauss-Legendre weights are reconstructed from P_N'(mu).
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
group_weight <- ring_expand((gaussian_weights / 2) / nlon)
land_weight <- group_weight * land_fraction
group_fraction <- vapply(names(koppen_colors), function(group) {
  sum(land_weight[koppen_group == group], na.rm = TRUE) / sum(land_weight[land])
}, numeric(1))

summary_table <- data.frame(
  metric = c("final_month_start", "final_month_end", "maximum_streamfunction_1e9_kg_s",
             "minimum_streamfunction_1e9_kg_s", "maximum_surface_mass_residual_1e9_kg_s",
             "global_area_mean_cloud_cover", "land_area_mean_cloud_cover",
             "ocean_area_mean_cloud_cover",
             paste0("land_area_fraction_group_", names(group_fraction))),
  value = c(min(final_months), max(final_months), max(psi_billion), min(psi_billion),
            max(abs(rowMeans(surface_mass_residual))) / 1e9,
            sum(group_weight * cloud_cover_annual),
            sum(group_weight * land_fraction * cloud_cover_annual) / sum(group_weight * land_fraction),
            sum(group_weight * (1 - land_fraction) * cloud_cover_annual) /
              sum(group_weight * (1 - land_fraction)),
            group_fraction)
)
write.csv(summary_table, file.path(analysis_dir, "analysis_summary.csv"), row.names = FALSE)

notes <- c(
  "Final-year mixed land/sea climate analysis",
  sprintf("Input months: %d-%d (12 equal 30-day months; April through March).", min(final_months), max(final_months)),
  "Köppen map: first-letter groups A-E from the 12 monthly means; land is land_fraction >= 0.5.",
  "The B threshold uses annual temperature, annual precipitation, and hemisphere-specific high-sun-half-year precipitation.",
  "The C/D boundary is -3 degC. A requires every month >=18 degC; E has warmest month <10 degC.",
  "Zonal means include both land and ocean grid cells.",
  "The cloud map shows final-year mean diagnosed effective column cloud fraction over land and ocean in fixed 10-percentage-point bins.",
  "Mass streamfunction uses monthly zonal-mean v and monthly zonal-mean surface pressure; pressure shown is the reference full-level pressure.",
  "Climographs use snapped native-grid cells with land_fraction >= 0.75, are reordered to January-December, and share both vertical-axis ranges across all panels.",
  "All annual means use the final 12 monthly means with equal weights."
)
writeLines(notes, file.path(analysis_dir, "README.txt"))

cat("Analysis written to", analysis_dir, "\n")
print(site_table, row.names = FALSE)
print(summary_table, row.names = FALSE)
