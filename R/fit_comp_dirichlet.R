#' Fit dirichlet regression model (via dirinla::dirinlareg() function)
#'
#' @param data description
#' @param formula description
#' @param use_dirinla description
#' @param tol0 description
#' @param verbose description
#'
#' @return description
#' @export
fit_comp_dirichlet <- function(data, formula, use_dirinla = FALSE, tol0 = NULL,
                               verbose = FALSE) {
  # simulate observed substrate compositions
  alpha_sub <- data |>
    dplyr::select(dplyr::starts_with("alpha_sim")) |>
    as.matrix()
  y_sub <- DirichletReg::rdirichlet(nrow(data), alpha = alpha_sub)
  # may need to use this data transformation if dirinla::dirinlareg() fails
  # note that the models ability to estimate very low substrate composition
  # probabilities may be affected by this!
  y_sub <- DirichletReg::DR_data(y_sub)
  # colnames(y_sub) <- paste0("y", 1:ncol(y_sub))

  # fit the model
  if (use_dirinla) {
    if (is.null(tol0)) {
      dirinla::dirinlareg(formula, y = y_sub, data.cov = data, k0 = 100,
                          prec = 0.0001, verbose = verbose)
    } else {
      dirinla::dirinlareg(formula, y = y_sub, data.cov = data, tol0 = tol0,
                          k0 = 100, prec = 0.0001, verbose = verbose)
    }
  } else {
    data$y <- y_sub
    DirichletReg::DirichReg(substitute_formula, data = data) |>
      substitute(env = list("substitute_formula" = formula)) |>
      eval()
  }
}


#' Sample the data and fit dirichlet regression model (via dirinla:: function)
#'
#' @param data description
#' @param n description
#' @param method description
#' @param bias_var description
#' @param bias_thresh description
#' @param clh_var description
#' @param clh_iter description
#' @param strata_var description
#' @param formula description
#'
#' @return description
#' @export
sample_fit_comp_dirichlet <- function(data, n, method, bias_var = NULL,
                                      bias_thresh = NULL, clh_var = NULL,
                                      clh_iter = NULL, strata_var = NULL,
                                      formula) {

  # if data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(data, "PackedSpatRaster")) data <- terra::unwrap(data)

  # sample the data
  df_sub <- data$data |>
    rastersample::spatial_sample(n = n, method = method, bias_var = bias_var,
                                 bias_thresh = bias_thresh, clh_var = clh_var,
                                 clh_iter = clh_iter, strata_var = strata_var)

  # fit dirichlet regression model
  mod <- fit_comp_dirichlet(df_sub, formula = formula)

  # model-estimated fixed effects
  data$coef_hat <- purrr::map(mod$summary_fixed, .f = \(x) x$mean) |>
    unlist()

  # predict substrate composition probabilities for entire study area
  data$data <- predict_comp_dirichlet(data$data, model = mod)

  # return object `data`
  data

}

#' Predict substrate composition values based on dirichlet regression model
#'
#' @param data description
#' @param model description
#'
#' @return description
#' @export
predict_comp_dirichlet <- function(data, model) {

  # if data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(data, "PackedSpatRaster")) data <- terra::unwrap(data)

  # predict substrate composition probabilities for entire study area
  # calculate expected values from model-estimated alpha values
  if (inherits(data, "SpatRaster")) {
    df <- data |>
      terra::as.data.frame(cells = TRUE, na.rm = TRUE) |>
      tibble::as_tibble()
  } else {
    df <- data
  }
  if (inherits(model, "dirinlaregmodel")) {
    n_cat <- model$ncat
    coef_hat <- purrr::map(model$summary_fixed, .f = \(x) x$mean) |>
      unlist()
  } else if (inherits(model, "DirichletRegModel")) {
    n_cat <- length(model$n.vars)
    coef_hat <- model$coefficients
  }
  ds_hat <- dirinla::data_stack_dirich(
    y = rep.int(NA, times = nrow(df) * n_cat),
    covariates = dirinla::formula_list(model$formula),
    data = df,
    d = n_cat,
    n = nrow(df)
  )
  eta_hat <- ds_hat %*% coef_hat
  alpha_hat <- eta_hat |>
    exp() |>
    matrix(ncol = n_cat, byrow = TRUE)
  p_hat <- alpha_hat / rowSums(alpha_hat)  # expected values of Dirichlet distribution with parameters equal to alpha_hat

  if (inherits(data, "SpatRaster")) {
    r_temp <- terra::subset(data, stringr::str_c("p_sim_", 1:n_cat))
    terra::values(r_temp) <- NA
    r_alpha_hat <- r_p_hat <- r_temp
    for (i in 1:n_cat) {
      r_alpha_hat[[i]][df$cell] <- alpha_hat[, i]
      r_p_hat[[i]][df$cell] <- p_hat[, i]
    }
    names(r_alpha_hat) <- stringr::str_c("alpha_hat_", 1:n_cat)
    names(r_p_hat) <- stringr::str_c("p_hat_", 1:n_cat)
    c(data, r_alpha_hat, r_p_hat) |>
      terra::wrap()
  } else {
    data |>
      dplyr::bind_cols(tibble::as_tibble(
        alpha_hat, .name_repair = ~ stringr::str_c("alpha_hat_", 1:n_cat))
      ) |>
      dplyr::bind_cols(tibble::as_tibble(
        p_hat, .name_repair = ~ stringr::str_c("p_hat_", 1:n_cat))
      )
  }

}
