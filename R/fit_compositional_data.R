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
fit_compositional_data <- function(data, formula) {
  # simulate observed substrate compositions
  alpha_sub <- data |>
    dplyr::select(dplyr::starts_with(".alpha_sim")) |>
    as.matrix()
  y_sub <- DirichletReg::rdirichlet(nrow(data), alpha = alpha_sub)
  y_sub <- DirichletReg::DR_data(y_sub)

  # fit the model
  data$y <- y_sub
  DirichletReg::DirichReg(substitute_formula, data = data) |>
    substitute(env = list("substitute_formula" = formula)) |>
    eval()
}
