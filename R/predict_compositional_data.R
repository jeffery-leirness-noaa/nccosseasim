#' Predict substrate composition probabilities
#'
#' This function predicts substrate composition probabilities based on a fitted
#' Dirichlet regression model.
#'
#' @param object A fitted Dirichlet regression model object of class `DirichletRegModel` (from `DirichletReg`).
#' @param new_data A data frame or a `PackedSpatRaster`/`SpatRaster` object containing the input data.
#'
#' @return A modified version of the input `data` object, including predicted
#'   substrate composition probabilities and model-estimated alpha values.
#'   If the input is a `SpatRaster`, the output is a wrapped raster object.
#' @export
predict_compositional_data <- function(object, new_data) {
  # if new_data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(new_data, "PackedSpatRaster")) {
    new_data <- terra::unwrap(new_data)
  }

  # predict substrate composition probabilities for entire study area
  # calculate expected values from model-estimated alpha values
  d <- object$dims
  if (inherits(new_data, "SpatRaster")) {
    df <- terra::as.data.frame(new_data, cells = TRUE, na.rm = FALSE)
  } else {
    df <- new_data
  }
  valid <- complete.cases(df)
  p_hat <- matrix(NA_real_, nrow = nrow(df), ncol = d)
  p_hat[valid, ] <- predict(object, newdata = df[valid, ], mu = TRUE)

  if (inherits(new_data, "SpatRaster")) {
    out <- terra::rast(terra::subset(new_data, 1), nlyr = d)
    terra::values(out) <- p_hat
    names(out) <- paste0(".p_hat", seq_len(d))
    terra::wrap(out)
  } else {
    tibble::as_tibble(
      p_hat,
      .name_repair = ~ stringr::str_c(".p_hat", 1:d)
    )
  }
}
