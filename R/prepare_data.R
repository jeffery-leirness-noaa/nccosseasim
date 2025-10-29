#' Generate polynomial feature rasters from a (Packed) SpatRaster
#'
#' @description
#' Expands each layer of a `SpatRaster` (or `PackedSpatRaster`) into
#' `poly_degree` polynomial features using
#' `stats::poly(raw = FALSE, simple = TRUE)`.
#' Returns a `PackedSpatRaster` whose layers are named
#' `<original_layer>_poly1`, ..., `<original_layer>_poly<poly_degree>`.
#'
#' @details
#' The function:
#' 1. Unwraps a `PackedSpatRaster` if supplied.
#' 2. Extracts values as a dense matrix (`terra::values(mat = TRUE)`).
#' 3. For each original layer, computes an orthogonal polynomial basis
#'    (centered & scaled) of the specified degree. NAs are preserved.
#' 4. Reassembles all polynomial columns into a new `SpatRaster` and wraps.
#'
#' If raw powers (x, x^2, ...) are preferred (for direct coefficient
#' interpretation), consider an alternative implementation with `raw = TRUE`.
#'
#' Memory: Allocates one matrix of size `ncell(x) * (nlyr(x) * poly_degree)`
#' plus temporary per-layer matrices; for very large rasters consider chunked
#' processing.
#'
#' @param x A `SpatRaster` or `PackedSpatRaster` whose layers will be expanded.
#' @param poly_degree Integer >= 1; number of polynomial terms per layer.
#'
#' @return A `PackedSpatRaster` containing the polynomial feature layers.
#'
#' @seealso [terra::rast()], [stats::poly()]
#'
#' @examples
#' \donttest{
#' library(terra)
#' f <- system.file("ex/elev.tif", package = "terra")
#' r <- rast(f)
#' poly_r <- prepare_data(r, poly_degree = 3)
#' names(poly_r)
#' }
#'
#' @importFrom terra unwrap values rast wrap
#' @importFrom stats poly
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
    p <- matrix(NA, nrow = nrow(vals), ncol = poly_degree)
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
