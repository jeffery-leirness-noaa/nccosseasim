#' Performance metrics for compositional data
#'
#' @param truth description
#' @param estimate description
#' @param metric description
#' @param summarize description
#'
#' @return description
#' @examples
metrics_comp <- function(truth, estimate, metric = "all", summarize = FALSE) {
  if (metric == "rmse") {
    comp_rmse(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "mae") {
    comp_mae(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "rho") {
    comp_rho(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "mlr") {
    comp_mlr(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "all") {
    dplyr::bind_rows(
      comp_rmse(truth, estimate = estimate, summarize = summarize),
      comp_mae(truth, estimate = estimate, summarize = summarize),
      comp_rho(truth, estimate = estimate, summarize = summarize),
      comp_mlr(truth, estimate = estimate, summarize = summarize)
    )
  }
}
# root mean square error
rmse_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (truth - estimate) ^ 2 |>
    apply(MARGIN = 2, FUN = mean) |>
    sqrt()
  if (summarize) {
    tibble::tibble(.metric = "rmse", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "rmse", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = id, values_from = .estimate,
                         names_prefix = ".estimate_")
  }
}
# mean absolute error
mae_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (truth - estimate) |>
    abs() |>
    apply(MARGIN = 2, FUN = mean)
  if (summarize) {
    tibble::tibble(.metric = "mae", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "mae", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = id, values_from = .estimate,
                         names_prefix = ".estimate_")
  }
}
# spearman rank correlation
rho_comp <- function(truth, estimate, summarize = FALSE) {
  m <- stats::cor(truth, estimate, use = "complete.obs", method = "spearman") |>
    diag()
  if (summarize) {
    tibble::tibble(.metric = "rho", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "rho", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = id, values_from = .estimate,
                         names_prefix = ".estimate_")
  }
}
# mean log ratio
# caution: this is experimental, has not been tested/verified, nor has any real theory behind its use
mlr_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (estimate / truth) |>
    log() |>
    apply(MARGIN = 2, FUN = mean)
  if (summarize) {
    tibble::tibble(.metric = "mlr", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "mlr", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = id, values_from = .estimate,
                         names_prefix = ".estimate_")
  }
}


#' Spearman rank correlation (implemented as a custom performance metric in the [yardstick](https://yardstick.tidymodels.org/) package)
#'
#' Calculate the Spearman rank correlation via stats::cor(truth, estimate, use = "complete.obs", method = "spearman").
#'
#' @param data A `data.frame` containing the columns specified by the truth and estimate arguments.
#' @param ... Not currently used.
#' @param truth The column identifier for the true results (that is `numeric`). This should be an unquoted column name although this argument is passed by expression and supports quasiquotation (you can unquote column names). For `_vec()`⁠ functions, a `numeric` vector. This value is passed the `x` argument of `stats::cor()`.
#' @param estimate The column identifier for the predicted results (that is also `numeric`). As with `truth` this can be specified different ways but the primary method is to use an unquoted variable name. For `_vec()`⁠ functions, a numeric `vector`. This value is passed the `y` argument of `stats::cor()`.
#' @param na_rm A `logical` value indicating whether `NA` values should be stripped before the computation proceeds.
#' @param case_weights The optional column identifier for case weights. This should be an unquoted column name that evaluates to a numeric column in `data`. For `_vec()` functions, a `numeric` vector.
#'
#' @return
#' @examples
rho <- function(data, ...) {
  UseMethod("rho")
}
rho <- yardstick::new_numeric_metric(rho, direction = "maximize")

rho.data.frame <- function(data, truth, estimate, na_rm = TRUE,
                           case_weights = NULL, ...) {
  yardstick::numeric_metric_summarizer(
    name = "rho",
    fn = rho_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    na_rm = na_rm,
    case_weights = !!rlang::enquo(case_weights)
  )
}

rho_vec <- function(truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  yardstick::check_numeric_metric(truth, estimate, case_weights)

  if (na_rm) {
    result <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)

    truth <- result$truth
    estimate <- result$estimate
    case_weights <- result$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }

  rho_impl(truth, estimate, case_weights = case_weights)
}

rho_impl <- function(truth, estimate, case_weights = NULL) {
  stats::cor(truth, estimate, use = "complete.obs", method = "spearman")
}

#' Median true result when estimated value is greater than or equal to some threshold value (implemented as a custom performance metric in the [yardstick](https://yardstick.tidymodels.org/) package)
#'
#' @param data A `data.frame` containing the columns specified by the truth and estimate arguments.
#' @param ... Not currently used.
#' @param truth The column identifier for the true results (that is `numeric`). This should be an unquoted column name although this argument is passed by expression and supports quasiquotation (you can unquote column names). For `_vec()`⁠ functions, a `numeric` vector.
#' @param estimate The column identifier for the predicted results (that is also `numeric`). As with `truth` this can be specified different ways but the primary method is to use an unquoted variable name. For `_vec()`⁠ functions, a numeric `vector`.
#' @param thresh Threshold value
#' @param na_rm A `logical` value indicating whether `NA` values should be stripped before the computation proceeds.
#' @param case_weights The optional column identifier for case weights. This should be an unquoted column name that evaluates to a numeric column in `data`. For `_vec()` functions, a `numeric` vector.
#'
#' @return
#' @examples
mtr <- function(data, ...) {
  UseMethod("mtr")
}
mtr <- yardstick::new_numeric_metric(mtr, direction = "minimize")

mtr.data.frame <- function(data, truth, estimate, thresh = 0, na_rm = TRUE,
                           case_weights = NULL, ...) {
  yardstick::numeric_metric_summarizer(
    name = "mtr",
    fn = mtr_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    na_rm = na_rm,
    case_weights = !!rlang::enquo(case_weights),
    fn_options = list(thresh = thresh)
  )
}

mtr_vec <- function(truth, estimate, thresh = 0, na_rm = TRUE,
                    case_weights = NULL, ...) {
  yardstick::check_numeric_metric(truth, estimate, case_weights)

  if (na_rm) {
    result <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)

    truth <- result$truth
    estimate <- result$estimate
    case_weights <- result$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }

  mtr_impl(truth, estimate, thresh, case_weights = case_weights)
}

mtr_impl <- function(truth, estimate, thresh = 0, case_weights = NULL) {
  median(truth[estimate >= thresh])
}
