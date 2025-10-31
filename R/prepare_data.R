#' Orthogonal polynomial raster features
#'
#' @description Expand each raster layer into an orthogonal polynomial basis
#' (degree `poly_degree`) using [stats::poly()] (centered & scaled, non-raw).
#'
#' @param x A `SpatRaster` or `PackedSpatRaster`.
#' @param poly_degree Integer >= 1; number of polynomial terms per original layer.
#'
#' @details
#' Steps:
#' 1. Unwrap (if packed).
#' 2. Extract values matrix.
#' 3. For each layer compute orthogonal polynomial basis (NAs preserved).
#' 4. Assemble into new raster and wrap.
#'
#' Output layer naming: `<original>_poly1`, ..., `_poly<poly_degree>`.
#'
#' @return A wrapped `SpatRaster` containing expanded polynomial layers.
#'
#' @seealso [simulate_compositional_data()], [stats::poly()]
#'
#' @examples
#' \donttest{
#' r <- terra::rast(system.file("ex/elev.tif", package = "terra"))
#' poly_r <- prepare_data(r, poly_degree = 2)
#' names(poly_r)
#' }
#'
#' @export
prepare_data <- function(x, poly_degree = 1) {
  if (inherits(x, "PackedSpatRaster")) {
    x <- terra::unwrap(x)
  }
  if (poly_degree < 1) {
    stop("poly_degree must be >= 1")
  }
  vals <- terra::values(x, mat = TRUE)
  layer_names <- names(x)
  poly_list <- vector("list", ncol(vals))
  for (i in seq_along(poly_list)) {
    idx_i <- !is.na(vals[, i])
    p <- matrix(NA_real_, nrow = nrow(vals), ncol = poly_degree)
    p[idx_i, ] <- stats::poly(
      vals[idx_i, i],
      degree = poly_degree,
      raw = FALSE,
      simple = TRUE
    )
    colnames(p) <- paste0(layer_names[i], "_poly", seq_len(ncol(p)))
    poly_list[[i]] <- p
  }
  vals_poly <- do.call(cbind, poly_list)
  out <- terra::rast(x, nlyr = ncol(vals_poly))
  terra::values(out) <- vals_poly
  names(out) <- colnames(vals_poly)
  terra::wrap(out)
}
