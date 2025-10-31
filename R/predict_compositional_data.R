#' Predict compositional probabilities
#'
#' @description Generate expected composition probabilities (.p_hat*) from a
#' fitted `DirichletRegModel` over new covariate data (data frame or raster).
#'
#' @param object A fitted `DirichletRegModel`.
#' @param new_data A data frame or (Packed) `SpatRaster` containing predictor
#' variables used in model fitting.
#'
#' @return If `new_data` is a raster: wrapped `SpatRaster` with `.p_hat1`,
#' `.p_hat2`, ... layers. Else a tibble with the same columns.
#'
#' @details Missing rows (in predictors) yield NA for all components.
#'
#' @seealso [fit_compositional_data()], [metrics_comp()]
#'
#' @examples
#' \donttest{
#' # object <- fit_compositional_data(df, y ~ x1 + x2)
#' # preds <- predict_compositional_data(object, df)
#' }
#'
#' @export
predict_compositional_data <- function(object, new_data) {
  if (inherits(new_data, "PackedSpatRaster")) {
    new_data <- terra::unwrap(new_data)
  }
  d <- object$dims
  if (inherits(new_data, "SpatRaster")) {
    df <- terra::as.data.frame(new_data, cells = TRUE, na.rm = FALSE)
  } else {
    df <- new_data
  }
  valid <- stats::complete.cases(df)
  p_hat <- matrix(NA_real_, nrow = nrow(df), ncol = d)
  p_hat[valid, ] <- DirichletReg::predict.DirichletRegModel(
    object,
    newdata = df[valid, ],
    mu = TRUE
  )
  if (inherits(new_data, "SpatRaster")) {
    out <- terra::rast(terra::subset(new_data, 1), nlyr = d)
    terra::values(out) <- p_hat
    names(out) <- paste0(".p_hat", seq_len(d))
    terra::wrap(out)
  } else {
    tibble::as_tibble(
      p_hat,
      .name_repair = ~ stringr::str_c(".p_hat", seq_len(d))
    )
  }
}
