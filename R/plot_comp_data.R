#' Plot simulated compositional data
#'
#' @param sim description
#'
#' @return description
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
