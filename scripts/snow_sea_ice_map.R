# Shared map of the cryosphere for the land-sea analysis scripts: the land
# snow cover f = S/(S + S_0) (docs/tendency/snow.md) and the sea-ice
# concentration A (docs/tendency/sea-ice.md) on one map per selected month.
#
# Sourced by analyze_moist_land_sea.R and analyze_moist_land_sea_earth.R. The
# grid is the model's native octahedral grid, ring-major from south to north.
# Every name defined here starts with snow_ice_ (or is draw_snow_sea_ice_maps)
# so that the sourcing scripts cannot overwrite it with their own variables.

# Five-step single-hue ramps, each light end at >= 2:1 contrast against the
# colour of its own "no cover" class (checked with the dataviz palette validator).
snow_ice_snow_colours <- c("#9B92DA", "#8175CB", "#6759B2", "#4C3F8F", "#312271")
snow_ice_ice_colours <- c("#659FDA", "#2C86CA", "#006BAC", "#005089", "#00366C")
snow_ice_bare_land_colour <- "#e6e2d8"
snow_ice_open_ocean_colour <- "#dde6ec"
# Cover below 1% counts as none, so that trace amounts do not tint a cell.
snow_ice_snow_breaks <- c(0.01, 0.1, 0.3, 0.5, 0.7, 1.0001)
snow_ice_ice_breaks <- c(0.01, 0.15, 0.3, 0.5, 0.85, 1.0001)

# Reads the monthly field of every selected output month; NULL if any file is missing.
snow_ice_read_monthly <- function(case_dir, field, months, reader) {
  paths <- file.path(case_dir, sprintf("monthly_%s_m%04d.bin", field, months))
  if (!all(file.exists(paths))) {
    if (any(file.exists(paths))) stop("Incomplete monthly output for ", field)
    return(NULL)
  }
  vapply(paths, reader, numeric(length(reader(paths[[1]]))), USE.NAMES = FALSE)
}

# Cell edges of the native grid, with the prime meridian at the centre.
snow_ice_cell_geometry <- function(nlon, latitude, grid_longitude) {
  nlat <- length(nlon)
  list(nlon = nlon,
       latitude_edges = c(-90, 0.5 * (latitude[-1] + latitude[-nlat]), 90),
       centre = ifelse(grid_longitude >= 180, grid_longitude - 360, grid_longitude),
       grid_longitude = grid_longitude)
}

snow_ice_draw_cells <- function(geometry, colour) {
  offset <- 0L
  edges <- geometry$latitude_edges
  for (j in seq_along(geometry$nlon)) {
    delta_lon <- 360 / geometry$nlon[[j]]
    index <- offset + seq_len(geometry$nlon[[j]])
    left <- geometry$centre[index] - delta_lon / 2
    right <- geometry$centre[index] + delta_lon / 2
    cell <- colour[index]
    keep <- !is.na(cell)
    # Cells straddling the -180/180 seam are drawn on both sides.
    rect(pmax(left[keep], -180), edges[[j]], pmin(right[keep], 180), edges[[j + 1]], col = cell[keep], border = NA)
    west <- keep & left < -180
    east <- keep & right > 180
    if (any(west)) rect(360 + left[west], edges[[j]], 180, edges[[j + 1]], col = cell[west], border = NA)
    if (any(east)) rect(-180, edges[[j]], right[east] - 360, edges[[j + 1]], col = cell[east], border = NA)
    offset <- offset + geometry$nlon[[j]]
  }
}

# Edges between land and ocean cells of the land mask.
snow_ice_coastline <- function(geometry, land) {
  nlon <- geometry$nlon
  edges <- geometry$latitude_edges
  ring_offset <- cumsum(c(0L, nlon))
  cell_index <- function(j, i) ring_offset[[j]] + ((i - 1L) %% nlon[[j]]) + 1L
  x0 <- y0 <- x1 <- y1 <- numeric(0)
  for (j in seq_along(nlon)) {
    delta_lon <- 360 / nlon[[j]]
    for (i in seq_len(nlon[[j]])) {
      index <- cell_index(j, i)
      if (!land[[index]]) next
      left <- geometry$centre[[index]] - delta_lon / 2
      right <- geometry$centre[[index]] + delta_lon / 2
      if (left < -180 || right > 180) next
      if (!land[[cell_index(j, i + 1L)]]) {
        x0 <- c(x0, right); x1 <- c(x1, right); y0 <- c(y0, edges[[j]]); y1 <- c(y1, edges[[j + 1]])
      }
      if (!land[[cell_index(j, i - 1L)]]) {
        x0 <- c(x0, left); x1 <- c(x1, left); y0 <- c(y0, edges[[j]]); y1 <- c(y1, edges[[j + 1]])
      }
      for (neighbour in c(j - 1L, j + 1L)) {
        if (neighbour < 1L || neighbour > length(nlon)) next
        i_neighbour <- floor(geometry$grid_longitude[[index]] / (360 / nlon[[neighbour]]) + 0.5)
        if (!land[[cell_index(neighbour, i_neighbour + 1L)]]) {
          edge <- if (neighbour < j) edges[[j]] else edges[[j + 1]]
          x0 <- c(x0, left); x1 <- c(x1, right); y0 <- c(y0, edge); y1 <- c(y1, edge)
        }
      }
    }
  }
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}

# Colour of every cell: land (land_fraction >= 0.5) by its snow cover, every
# other cell by its sea-ice concentration.
snow_ice_cell_colours <- function(land, snow_fraction, ice_fraction) {
  snow_class <- findInterval(snow_fraction, snow_ice_snow_breaks)
  ice_class <- if (is.null(ice_fraction)) rep(0L, length(land)) else findInterval(ice_fraction, snow_ice_ice_breaks)
  snow_colour <- ifelse(snow_class == 0L, snow_ice_bare_land_colour, snow_ice_snow_colours[pmax(snow_class, 1L)])
  ice_colour <- ifelse(ice_class == 0L, snow_ice_open_ocean_colour, snow_ice_ice_colours[pmax(ice_class, 1L)])
  ifelse(land, snow_colour, ice_colour)
}

snow_ice_percent_labels <- function(breaks) {
  lower <- round(100 * head(breaks, -1))
  upper <- pmin(100, round(100 * tail(breaks, -1)))
  sprintf("%d–%d%%", lower, upper)
}

# Snow water equivalent S = S_0 f/(1 - f) at the snow-cover breaks, in kg m^-2 (mm).
snow_ice_snow_water_labels <- function(masking_water_equivalent) {
  f <- snow_ice_snow_breaks
  s <- ifelse(f >= 1, Inf, masking_water_equivalent * f / (1 - f))
  format_s <- function(x) ifelse(is.infinite(x), "", ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", x)))
  lower <- format_s(head(s, -1))
  upper <- format_s(tail(s, -1))
  ifelse(upper == "", sprintf("S ≥ %s", lower), sprintf("S %s–%s", lower, upper))
}

# One map per selected calendar month, stacked vertically, with both legends
# in the right margin so that no polar region is hidden.
draw_snow_sea_ice_maps <- function(path, nlon, latitude, grid_longitude, land_fraction, calendar_month,
                                   snow_fraction_monthly, sea_ice_fraction_monthly, masking_water_equivalent,
                                   period_label, device_type, months = c(3L, 9L)) {
  land <- land_fraction >= 0.5
  geometry <- snow_ice_cell_geometry(nlon, latitude, grid_longitude)
  coastline <- snow_ice_coastline(geometry, land)
  png(path, width = 2400, height = 1150 * length(months) + 150, res = 180, type = device_type)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(length(months), 1), mar = c(4.6, 5, 3.4, 13.5), oma = c(1.4, 0, 2.2, 0), las = 1)
  for (month in months) {
    column <- which(calendar_month == month)
    if (length(column) != 1L) stop("The final year lacks calendar month ", month)
    ice <- if (is.null(sea_ice_fraction_monthly)) NULL else sea_ice_fraction_monthly[, column]
    colours <- snow_ice_cell_colours(land, snow_fraction_monthly[, column], ice)
    plot(NA, xlim = c(-180, 180), ylim = c(-90, 90), xaxs = "i", yaxs = "i",
         xlab = "Longitude", ylab = "Latitude", axes = FALSE,
         main = sprintf("%s monthly mean", month.name[[month]]))
    snow_ice_draw_cells(geometry, colours)
    abline(h = seq(-60, 60, 30), v = seq(-180, 180, 60), col = adjustcolor("white", 0.5), lwd = 0.7)
    segments(coastline$x0, coastline$y0, coastline$x1, coastline$y1, col = "#333333", lwd = 0.8)
    axis(1, at = seq(-180, 180, 60), labels = c("180°W", "120°W", "60°W", "0°", "60°E", "120°E", "180°E"))
    axis(2, at = seq(-90, 90, 30), labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"))
    box()
    legend("topleft", inset = c(1.01, 0), xpd = NA, bty = "n", cex = 0.72,
           title = "Land snow cover f", title.adj = 0,
           legend = c(sprintf("%s  (%s)", rev(snow_ice_percent_labels(snow_ice_snow_breaks)),
                              rev(snow_ice_snow_water_labels(masking_water_equivalent))), "< 1%  (no snow)"),
           fill = c(rev(snow_ice_snow_colours), snow_ice_bare_land_colour), border = "#555555")
    legend("bottomleft", inset = c(1.01, 0), xpd = NA, bty = "n", cex = 0.72,
           title = if (is.null(ice)) "Sea ice (not in output)" else "Sea-ice concentration A", title.adj = 0,
           legend = c(rev(snow_ice_percent_labels(snow_ice_ice_breaks)), "< 1%  (open water)"),
           fill = c(rev(snow_ice_ice_colours), snow_ice_open_ocean_colour), border = "#555555")
  }
  mtext(sprintf("Land snow cover and sea ice, %s", period_label), side = 3, outer = TRUE, line = 0.4, cex = 1.12)
  mtext(sprintf(paste("Land is land fraction >= 0.5 and is coloured by f = S/(S + S_0), S_0 = %g kg m^-2;",
                      "other cells by the ice concentration per ocean area; S in kg m^-2 (mm water equivalent)"),
                masking_water_equivalent), side = 1, outer = TRUE, line = 0.2, cex = 0.7)
  invisible(path)
}
