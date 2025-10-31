#' Plot simulated compositional rasters
#'
#' @description This function generates a plot of simulated compositional data
#' stored in a `SpatRaster` object. It extracts the relevant layers from the
#' raster object and visualizes them as a series of plots.
#'
#' @param sim A list containing the simulated data. The `sim$data` element must
#'   be a `SpatRaster` or `PackedSpatRaster` object. The `sim$d` element
#'   specifies the number of compositional components to plot.
#'
#' @return Invisibly returns `sim`. Produces a multi-panel plot as side effect.
#'
#' @examples
#' \donttest{
#' r <- terra::rast(system.file("ex/elev.tif", package = "terra"))
#' sim <- simulate_compositional_data(r, d = 3, as_raster = TRUE)
#' plot_compositional_data(sim)
#' }
#'
#' @seealso [simulate_compositional_data()]
#' @export
plot_compositional_data <- function(sim) {
  if (inherits(sim$data, "PackedSpatRaster")) {
    sim$data <- terra::unwrap(sim$data)
  }
  stopifnot(inherits(sim$data, "SpatRaster"))
  sim$data |>
    terra::subset(subset = stringr::str_c(".p_sim", seq_len(sim$d))) |>
    terra::plot(range = c(0, 1), nr = 1)
  invisible(sim)
}
