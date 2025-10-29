#' Fit a Dirichlet regression model
#'
#' @description Generates a compositional response from simulated Dirichlet
#' parameters (.alpha_sim*) and fits a Dirichlet regression using
#' [DirichletReg::DirichReg()].
#'
#' @details
#' Steps:
#' 1. Extract columns starting with `.alpha_sim` and form an alpha matrix.
#' 2. Draw one Dirichlet composition per row.
#' 3. Attach the composition as column `y` (a `DirichletRegData` object).
#' 4. Fit the supplied formula with `DirichletReg::DirichReg()`.
#'
#' The formula is passed unmodified. Provide a multi-part (|) formula for
#' component-specific predictors if desired.
#'
#' @param formula A model formula acceptable to [DirichletReg::DirichReg()].
#' @param data A data frame containing columns `.alpha_sim1`, `.alpha_sim2`, ... .
#'
#' @return A `DirichletRegModel` object.
#'
#' @seealso [simulate_compositional_data()], [predict_compositional_data()]
#'
#' @examples
#' \donttest{
#' # Assume sim is result of simulate_compositional_data(as_raster = FALSE)
#' # with columns .alpha_sim1:.alpha_sim3
#' # formula with common predictors for all components:
#' # y ~ x1 + x2
#' # (y will be created internally from .alpha_sim*)
#' # fit_compositional_data(sim, y ~ 1)
#' }
#'
#' @export
fit_compositional_data <- function(formula, data) {
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
