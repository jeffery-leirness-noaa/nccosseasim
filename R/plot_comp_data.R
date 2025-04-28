#' Plot simulated compositional data
#'
#' This function generates a plot of simulated compositional data stored in a `SpatRaster` object. 
#' It extracts the relevant layers from the raster object and visualizes them as a series of plots.
#'
#' @param sim A list containing the simulated data. The `sim$data` element must be a `SpatRaster` 
#'   or `PackedSpatRaster` object. The `sim$d` element specifies the number of compositional 
#'   components to plot.
#'
#' @return A plot of the compositional data layers. The function does not return a value; 
#'   it generates a plot as a side effect.
#' 
#' @details If `sim$data` is a `PackedSpatRaster` object, it is first unpacked using 
#'   `terra::unwrap()`. The function ensures that `sim$data` is a `SpatRaster` object 
#'   before proceeding. The layers to be plotted are identified by their names, which 
#'   are expected to follow the pattern `"p_sim_1"`, `"p_sim_2"`, ..., `"p_sim_d"`, 
#'   where `d` is the number of components.
#'
#' @export
plot_comp_data <- function(sim) {

  # if sim$data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(sim$data, "PackedSpatRaster")) sim$data <- terra::unwrap(sim$data)

  # error if sim$data is not a SpatRaster object
  stopifnot(inherits(sim$data, "SpatRaster"))

  # create plot
  sim$data |>
    terra::subset(subset = stringr::str_c("p_sim_", 1:sim$d)) |>
    terra::plot(range = c(0, 1), nr = 1)
}
