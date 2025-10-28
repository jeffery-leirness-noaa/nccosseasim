#' Fit a Dirichlet regression model
#'
#' This function fits a Dirichlet regression model using the `DirichletReg::DirichReg()`
#' function.
#'
#' @param data A data frame containing the input data. Must include columns with names
#'   starting with "alpha_sim" for the Dirichlet parameters.
#' @param formula A formula specifying the regression model to be fitted.
#'
#' @return A fitted Dirichlet regression model object.
#' @export
fit_comp_dirichlet <- function(data, formula) {
  # simulate observed substrate compositions
  alpha_sub <- data |>
    dplyr::select(dplyr::starts_with("alpha_sim")) |>
    as.matrix()
  y_sub <- DirichletReg::rdirichlet(nrow(data), alpha = alpha_sub)
  y_sub <- DirichletReg::DR_data(y_sub)

  # fit the model
  data$y <- y_sub
  DirichletReg::DirichReg(substitute_formula, data = data) |>
    substitute(env = list("substitute_formula" = formula)) |>
    eval()
}


#' Sample data and fit a Dirichlet regression model
#'
#' This function samples the input data and fits a Dirichlet regression model using
#' the `fit_comp_dirichlet()` function.
#'
#' @param data A data frame or a `PackedSpatRaster` object containing the input data.
#' @param n Integer. The number of samples to draw.
#' @param method Character. The sampling method to use (e.g., "random", "stratified").
#' @param bias_var Character. The name of the variable to use for biased sampling (optional).
#' @param bias_thresh Numeric. A threshold for biased sampling (optional).
#' @param clh_var Character. The name of the variable for conditional Latin hypercube sampling (optional).
#' @param clh_iter Integer. The number of iterations for conditional Latin hypercube sampling (optional).
#' @param strata_var Character. The name of the variable to use for stratified sampling (optional).
#' @param formula A formula specifying the regression model to be fitted.
#'
#' @return A modified version of the input `data` object, including model coefficients
#'   and predicted substrate composition probabilities.
#' @export
sample_fit_comp_dirichlet <- function(
  data,
  n,
  method,
  bias_var = NULL,
  bias_thresh = NULL,
  clh_var = NULL,
  clh_iter = NULL,
  strata_var = NULL,
  formula
) {
  # if data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(data, "PackedSpatRaster")) {
    data <- terra::unwrap(data)
  }

  # sample the data
  df_sub <- data$data |>
    rastersample::spatial_sample(
      n = n,
      method = method,
      bias_var = bias_var,
      bias_thresh = bias_thresh,
      clh_var = clh_var,
      clh_iter = clh_iter,
      strata_var = strata_var
    )

  # fit dirichlet regression model
  mod <- fit_comp_dirichlet(df_sub, formula = formula)

  # model-estimated fixed effects
  # data$coef_hat <- purrr::map(mod$summary_fixed, .f = \(x) x$mean) |>
  #   unlist()

  # predict substrate composition probabilities for entire study area
  data$data <- predict_comp_dirichlet(mod, new_data = data$data)

  # return object `data`
  data
}

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
predict_comp_dirichlet <- function(object, new_data) {
  # if new_data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(new_data, "PackedSpatRaster")) {
    new_data <- terra::unwrap(new_data)
  }

  # predict substrate composition probabilities for entire study area
  # calculate expected values from model-estimated alpha values
  if (inherits(new_data, "SpatRaster")) {
    df <- new_data |>
      terra::as.data.frame(cells = TRUE, na.rm = TRUE) |>
      tibble::as_tibble()
  } else {
    df <- new_data
  }
  p_hat <- predict(object, newdata = df, mu = TRUE)

  if (inherits(new_data, "SpatRaster")) {
    r_temp <- terra::subset(new_data, stringr::str_c("p_sim_", 1:n_cat))
    terra::values(r_temp) <- NA
    r_p_hat <- r_temp
    for (i in 1:n_cat) {
      r_p_hat[[i]][df$cell] <- p_hat[, i]
    }
    names(r_p_hat) <- stringr::str_c("p_hat_", 1:n_cat)
    c(new_data, r_p_hat) |>
      terra::wrap()
  } else {
    new_data |>
      dplyr::bind_cols(tibble::as_tibble(
        p_hat,
        .name_repair = ~ stringr::str_c("p_hat_", 1:n_cat)
      ))
  }
}
