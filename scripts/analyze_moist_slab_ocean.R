#!/usr/bin/env Rscript

# Diagnose the five-year moist slab-ocean run and plot the final simulation year.
# This script intentionally depends only on base R and ggplot2 (for portability in
# the project environment). Binary layouts and units come from metadata.json.

args <- commandArgs(trailingOnly = TRUE)
case_dir <- if (length(args) >= 1) args[[1]] else "output/moist_slab_ocean"
analysis_dir <- if (length(args) >= 2) args[[2]] else file.path(case_dir, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(case_dir, "metadata.json")
daily_path <- file.path(case_dir, "daily_global.csv")
if (!file.exists(metadata_path) || !file.exists(daily_path)) {
  stop("Expected metadata.json and daily_global.csv under ", case_dir)
}

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
ocean_heat_capacity <- json_number("heat_capacity_j_m-2_k-1")
days_per_year <- as.integer(json_number("days_per_year"))
months_per_year <- as.integer(json_number("months_per_year"))
solar_day_seconds <- json_number("solar_day_seconds")

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

ring_mean <- function(values) {
  result <- numeric(nlat)
  offset <- 0L
  for (j in seq_len(nlat)) {
    indices <- offset + seq_len(nlon[[j]])
    result[[j]] <- mean(values[indices])
    offset <- offset + nlon[[j]]
  }
  result
}

ring_expand <- function(ring_values) rep(ring_values, times = nlon)

# The Gaussian nodes are already provided, but their quadrature weights are not.
# Compute the Gauss-Legendre weights from P_N'(mu) for area-weighted diagnostics.
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
area_ring_weights <- gaussian_weights / 2
area_grid_weights <- ring_expand(area_ring_weights / nlon)
stopifnot(abs(sum(area_ring_weights) - 1) < 1e-12)

monthly_files <- list.files(case_dir, pattern = "^monthly_surface_temperature_m[0-9]{4}\\.bin$", full.names = FALSE)
month_indices <- sort(as.integer(sub(".*_m([0-9]{4})\\.bin$", "\\1", monthly_files)))
if (length(month_indices) < months_per_year) stop("Fewer than one year of monthly outputs")
final_months <- tail(month_indices, months_per_year)
if (!all(diff(final_months) == 1)) stop("Final-year monthly files are not consecutive")

surface_temperature_grid <- numeric(npoints)
precipitation_grid <- numeric(npoints)
surface_temperature_zonal <- numeric(nlat)
precipitation_zonal <- numeric(nlat)
mass_streamfunction <- matrix(0, nrow = nlat, ncol = nlev)
surface_mass_residual <- matrix(0, nrow = nlat, ncol = length(final_months))

for (month_position in seq_along(final_months)) {
  month <- final_months[[month_position]]
  suffix <- sprintf("m%04d.bin", month)
  surface_temperature <- read_grid(file.path(case_dir, paste0("monthly_surface_temperature_", suffix)))
  precipitation <- read_grid(file.path(case_dir, paste0("monthly_precipitation_", suffix)))
  surface_pressure <- read_grid(file.path(case_dir, paste0("monthly_surface_pressure_", suffix)))
  zonal_v <- read_zonal(file.path(case_dir, paste0("monthly_zonal_v_", suffix)))

  surface_temperature_grid <- surface_temperature_grid + surface_temperature / months_per_year
  precipitation_grid <- precipitation_grid + precipitation / months_per_year
  surface_temperature_zonal <- surface_temperature_zonal + ring_mean(surface_temperature) / months_per_year
  precipitation_zonal <- precipitation_zonal + ring_mean(precipitation) / months_per_year

  zonal_surface_pressure <- ring_mean(surface_pressure)
  delta_pressure <- outer(rep(1, nlat), diff(a_half)) + outer(zonal_surface_pressure, diff(b_half))
  layer_transport <- (2 * pi * earth_radius / gravity) * cos_latitude * zonal_v * delta_pressure
  # Psi at each full level is the top-to-bottom cumulative transport evaluated
  # at the arithmetic layer midpoint, consistent with this model's full-level eta.
  psi_month <- t(apply(layer_transport, 1, function(x) cumsum(x) - 0.5 * x))
  mass_streamfunction <- mass_streamfunction + psi_month / months_per_year
  surface_mass_residual[, month_position] <- rowSums(layer_transport)
}

zonal_profiles <- data.frame(
  latitude_deg = latitude,
  surface_temperature_k = surface_temperature_zonal,
  precipitation_mm_day = precipitation_zonal
)
write.csv(zonal_profiles, file.path(analysis_dir, "final_year_zonal_profiles.csv"), row.names = FALSE)

streamfunction_table <- expand.grid(
  latitude_deg = latitude,
  level = seq_len(nlev),
  KEEP.OUT.ATTRS = FALSE
)
streamfunction_table$reference_pressure_hpa <- reference_full_pressure_hpa[streamfunction_table$level]
streamfunction_table$mass_streamfunction_kg_s <- as.vector(mass_streamfunction)
write.csv(streamfunction_table, file.path(analysis_dir, "final_year_mass_streamfunction.csv"), row.names = FALSE)

daily <- read.csv(daily_path, check.names = TRUE, strip.white = TRUE)
if (nrow(daily) < days_per_year) stop("Fewer than one year of daily diagnostics")
final_days <- tail(seq_len(nrow(daily)), days_per_year)
last <- daily[final_days, , drop = FALSE]
daily$toa_net_w_m2 <- daily$incoming_shortwave_w_m.2 - daily$reflected_shortwave_w_m.2 - daily$outgoing_longwave_w_m.2
last$toa_net_w_m2 <- last$incoming_shortwave_w_m.2 - last$reflected_shortwave_w_m.2 - last$outgoing_longwave_w_m.2

mean_last <- function(column) mean(last[[column]])
linear_trend_per_year <- function(column) {
  fit <- lm(last[[column]] ~ last$simulation_day)
  unname(coef(fit)[[2]]) * days_per_year
}

area_below_freezing <- sum(area_grid_weights[surface_temperature_grid < 273.15])
precip_peak_index <- which.max(precipitation_zonal)
psi_billion <- mass_streamfunction / 1e9
nh_tropics <- which(latitude >= 0 & latitude <= 30)
sh_tropics <- which(latitude >= -30 & latitude <= 0)
nh_peak <- which(psi_billion == max(psi_billion[nh_tropics, , drop = FALSE]), arr.ind = TRUE)[1, ]
sh_peak <- which(psi_billion == min(psi_billion[sh_tropics, , drop = FALSE]), arr.ind = TRUE)[1, ]
surface_residual_peak <- max(abs(rowMeans(surface_mass_residual)))
streamfunction_peak <- max(abs(mass_streamfunction))

signed_water_tendency <- (
  last$signed_column_water_kg_m.2[[nrow(last)]] - last$signed_column_water_kg_m.2[[1]]
) / (nrow(last) - 1)
water_flux_difference <- mean(last$evaporation_mm_day.1 - last$precipitation_mm_day.1)
water_budget_residual <- signed_water_tendency - water_flux_difference

annual_blocks <- split(seq_len(nrow(daily)), ceiling(seq_len(nrow(daily)) / days_per_year))
annual_metrics <- do.call(rbind, lapply(seq_along(annual_blocks), function(year) {
  rows <- annual_blocks[[year]]
  data.frame(
    simulation_year = year,
    mean_ocean_temperature_k = mean(daily$mean_ocean_temperature_k[rows]),
    precipitation_mm_day = mean(daily$precipitation_mm_day.1[rows]),
    evaporation_mm_day = mean(daily$evaporation_mm_day.1[rows]),
    precipitable_water_kg_m2 = mean(daily$precipitable_water_kg_m.2[rows]),
    negative_column_water_kg_m2 = mean(daily$negative_column_water_kg_m.2[rows]),
    toa_net_w_m2 = mean(daily$toa_net_w_m2[rows])
  )
}))
write.csv(annual_metrics, file.path(analysis_dir, "annual_global_diagnostics.csv"), row.names = FALSE)

metrics <- data.frame(
  metric = c(
    "final_month_start", "final_month_end", "final_day_start", "final_day_end",
    "surface_temperature_global_k", "surface_temperature_min_zonal_k", "surface_temperature_max_zonal_k",
    "surface_area_below_freezing_fraction", "ocean_temperature_daily_global_k",
    "ocean_temperature_trend_k_per_year", "atmospheric_temperature_daily_global_k",
    "precipitation_global_mm_day", "evaporation_global_mm_day", "precipitation_minus_evaporation_mm_day",
    "precipitation_peak_mm_day", "precipitation_peak_latitude_deg",
    "precipitable_water_kg_m2", "signed_column_water_kg_m2", "negative_column_water_kg_m2",
    "negative_to_positive_water_fraction", "water_budget_residual_mm_day",
    "incoming_shortwave_w_m2", "reflected_shortwave_w_m2", "outgoing_longwave_w_m2", "toa_net_w_m2",
    "toa_net_trend_w_m2_per_year", "maximum_daily_wind_m_s",
    "streamfunction_max_1e9_kg_s", "streamfunction_min_1e9_kg_s",
    "nh_hadley_peak_1e9_kg_s", "nh_hadley_peak_latitude_deg", "nh_hadley_peak_pressure_hpa",
    "sh_hadley_peak_1e9_kg_s", "sh_hadley_peak_latitude_deg", "sh_hadley_peak_pressure_hpa",
    "streamfunction_surface_residual_fraction"
  ),
  value = c(
    min(final_months), max(final_months), min(last$simulation_day), max(last$simulation_day),
    sum(surface_temperature_grid * area_grid_weights), min(surface_temperature_zonal), max(surface_temperature_zonal),
    area_below_freezing, mean_last("mean_ocean_temperature_k"),
    linear_trend_per_year("mean_ocean_temperature_k"), mean_last("mean_atmospheric_temperature_k"),
    mean_last("precipitation_mm_day.1"), mean_last("evaporation_mm_day.1"),
    mean_last("precipitation_mm_day.1") - mean_last("evaporation_mm_day.1"),
    precipitation_zonal[[precip_peak_index]], latitude[[precip_peak_index]],
    mean_last("precipitable_water_kg_m.2"), mean_last("signed_column_water_kg_m.2"),
    mean_last("negative_column_water_kg_m.2"),
    mean_last("negative_column_water_kg_m.2") / mean_last("precipitable_water_kg_m.2"),
    water_budget_residual,
    mean_last("incoming_shortwave_w_m.2"), mean_last("reflected_shortwave_w_m.2"),
    mean_last("outgoing_longwave_w_m.2"), mean(last$toa_net_w_m2), linear_trend_per_year("toa_net_w_m2"),
    max(last$maximum_wind_speed_m_s.1), max(psi_billion), min(psi_billion),
    psi_billion[nh_peak[[1]], nh_peak[[2]]], latitude[nh_peak[[1]]], reference_full_pressure_hpa[nh_peak[[2]]],
    psi_billion[sh_peak[[1]], sh_peak[[2]]], latitude[sh_peak[[1]]], reference_full_pressure_hpa[sh_peak[[2]]],
    surface_residual_peak / streamfunction_peak
  )
)
write.csv(metrics, file.path(analysis_dir, "physical_metrics.csv"), row.names = FALSE)

# Requested figure 1: final-year zonal-mean surface temperature and precipitation.
png(file.path(analysis_dir, "final_year_zonal_surface_temperature_precipitation.png"),
    width = 1800, height = 1400, res = 180, type = "quartz")
par(mfrow = c(2, 1), mar = c(2.0, 5.0, 3.0, 1.2), oma = c(4.0, 0, 2.5, 0),
    las = 1, mgp = c(3.2, 0.8, 0), cex.axis = 0.9)
plot(latitude, surface_temperature_zonal, type = "l", lwd = 3, col = "#c23b22",
     xlim = c(-90, 90), xaxt = "n", xlab = "", ylab = "Surface temperature (K)",
     main = "Surface temperature", panel.first = grid(col = "#dddddd", lty = 1))
abline(h = 273.15, lty = 2, lwd = 1.5, col = "#3572a5")
legend("bottomright", legend = "Water freezing point", lty = 2, lwd = 1.5,
       col = "#3572a5", bty = "n", cex = 0.85)
plot(latitude, precipitation_zonal, type = "l", lwd = 3, col = "#26734d",
     xlim = c(-90, 90), xaxt = "n", xlab = "", ylab = expression("Precipitation (mm day"^{-1}*")"),
     main = "Precipitation", panel.first = grid(col = "#dddddd", lty = 1))
axis(1, at = seq(-90, 90, 30), labels = paste0(abs(seq(-90, 90, 30)), c("°S", "°S", "°S", "°", "°N", "°N", "°N")))
mtext("Latitude", side = 1, outer = TRUE, line = 2.2)
mtext(sprintf("Moist slab-ocean aquaplanet - final simulation year (months %d-%d)",
              min(final_months), max(final_months)), side = 3, outer = TRUE, line = 0.6, cex = 1.15)
dev.off()

# Requested figure 2: annual-mean meridional mass streamfunction.
pressure_order <- rev(seq_len(nlev))
y <- -log10(reference_full_pressure_hpa[pressure_order])
z <- psi_billion[, pressure_order, drop = FALSE]
z_limit <- ceiling(max(abs(z)) / 10) * 10
if (z_limit == 0) z_limit <- 1
levels <- seq(-z_limit, z_limit, length.out = 25)
palette <- hcl.colors(length(levels) - 1, palette = "Blue-Red 3")
png(file.path(analysis_dir, "final_year_mass_streamfunction.png"),
    width = 1900, height = 1200, res = 180, type = "quartz")
filled.contour(
  x = latitude, y = y, z = z, levels = levels, color.palette = function(n) palette,
  xlim = c(-90, 90),
  plot.title = title(
    main = sprintf("Final-year meridional mass streamfunction (months %d-%d)",
                   min(final_months), max(final_months)),
    xlab = "Latitude", ylab = "Reference pressure (hPa)"
  ),
  plot.axes = {
    axis(1, at = seq(-90, 90, 30), labels = paste0(abs(seq(-90, 90, 30)), c("°S", "°S", "°S", "°", "°N", "°N", "°N")))
    pressure_ticks <- c(1000, 850, 700, 500, 300, 200, 100, 50, 10, 2)
    pressure_ticks <- pressure_ticks[pressure_ticks >= min(reference_full_pressure_hpa) &
                                       pressure_ticks <= max(reference_full_pressure_hpa)]
    axis(2, at = -log10(pressure_ticks), labels = pressure_ticks, las = 1)
    contour(latitude, y, z, levels = levels, add = TRUE, drawlabels = FALSE,
            col = adjustcolor("black", 0.35), lwd = 0.7)
    contour(latitude, y, z, levels = 0, add = TRUE, drawlabels = FALSE,
            col = "black", lwd = 1.6)
  },
  key.title = title(main = expression(10^9 * " kg " * s^{-1}), cex.main = 0.85),
  key.axes = axis(4, las = 1)
)
dev.off()

# Supporting diagnostic: annual means reveal whether the run is equilibrated.
png(file.path(analysis_dir, "five_year_global_diagnostics.png"),
    width = 1800, height = 1450, res = 180, type = "quartz")
par(mfrow = c(3, 1), mar = c(2.0, 5.0, 2.8, 1.2), oma = c(4.0, 0, 2.5, 0), las = 1)
plot(annual_metrics$simulation_year, annual_metrics$mean_ocean_temperature_k, type = "o", pch = 16,
     lwd = 2.5, col = "#c23b22", xlab = "", ylab = "Ocean temperature (K)",
     main = "Annual global-mean ocean temperature", panel.first = grid(col = "#dddddd"))
abline(h = 288, lty = 2, col = "#555555")
plot(annual_metrics$simulation_year, annual_metrics$toa_net_w_m2, type = "o", pch = 16,
     lwd = 2.5, col = "#704cb6", xlab = "", ylab = expression("TOA net (W m"^{-2}*")"),
     main = "Top-of-atmosphere net radiation (positive downward)", panel.first = grid(col = "#dddddd"))
abline(h = 0, lty = 2, col = "#555555")
plot(annual_metrics$simulation_year, annual_metrics$precipitation_mm_day, type = "o", pch = 16,
     lwd = 2.5, col = "#26734d", xlab = "", ylab = expression("Water flux (mm day"^{-1}*")"),
     main = "Global precipitation and evaporation", panel.first = grid(col = "#dddddd"))
lines(annual_metrics$simulation_year, annual_metrics$evaporation_mm_day, type = "o", pch = 17,
      lwd = 2.5, col = "#3572a5")
legend("topright", legend = c("Precipitation", "Evaporation"), col = c("#26734d", "#3572a5"),
       pch = c(16, 17), lwd = 2.5, bty = "n")
mtext("Simulation year", side = 1, outer = TRUE, line = 2.2)
mtext("Moist slab-ocean aquaplanet - equilibration diagnostics", side = 3, outer = TRUE, line = 0.6, cex = 1.15)
dev.off()

cat("Analysis written to", analysis_dir, "\n")
print(metrics, row.names = FALSE)
