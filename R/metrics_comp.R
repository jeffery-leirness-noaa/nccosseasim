#' Metrics for compositional predictions
#'
#' @description Compute performance metrics between true and estimated
#' compositions (rows sum to 1): RMSE, MAE, Spearman correlation (per component),
#' Kullback–Leibler divergence (component-wise), and Mean Log Ratio.
#'
#' @details
#' Input matrices/data frames must have identical dimensions and component order.
#' Each metric is computed component-wise; when `summarize = TRUE`, the mean
#' across components is returned. No re-scaling to enforce unity is performed.
#'
#' NA handling: rows containing NA in either `truth` or `estimate` for a given
#' component contribute NA to that component's metric before aggregation.
#'
#' Experimental: KLD and MLR are marked experimental (see notes).
#'
#' @param truth A matrix or data frame of true values for the compositional data.
#' @param estimate A matrix or data frame of estimated values for the compositional data.
#' @param metric A character string specifying the metric to compute. Options are:
#'   - `"rmse"`: Root mean square error.
#'   - `"mae"`: Mean absolute error.
#'   - `"rho"`: Spearman rank correlation.
#'   - `"kld"`: Kullback-Leibler divergence.
#'   - `"mlr"`: Mean log ratio.
#'   - `"all"`: Compute all metrics.
#' @param summarize Logical. If `TRUE`, returns a single summarized value for each metric.
#'   If `FALSE`, returns detailed results for each component.
#'
#' @return A tibble containing the computed metric(s). If `summarize = TRUE`, the tibble
#'   contains one row per metric. If `summarize = FALSE`, the tibble contains detailed
#'   results for each component.
#'
#' @seealso [predict_compositional_data()]
#'
#' @examples
#' truth <- matrix(c(0.2, 0.3, 0.5, 0.25, 0.25, 0.5), nrow = 2, byrow = TRUE)
#' est <- truth + matrix(c(0, 0.05, -0.05, 0.02, -0.02, 0), nrow = 2, byrow = TRUE)
#' metrics_comp(truth, est, metric = "all", summarize = FALSE)
#'
#' @export
metrics_comp <- function(truth, estimate, metric = "all", summarize = FALSE) {
  if (metric == "rmse") {
    rmse_comp(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "mae") {
    mae_comp(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "rho") {
    rho_comp(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "kld") {
    kld_comp(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "mlr") {
    mlr_comp(truth, estimate = estimate, summarize = summarize)
  } else if (metric == "all") {
    dplyr::bind_rows(
      rmse_comp(truth, estimate = estimate, summarize = summarize),
      mae_comp(truth, estimate = estimate, summarize = summarize),
      rho_comp(truth, estimate = estimate, summarize = summarize),
      kld_comp(truth, estimate = estimate, summarize = summarize),
      mlr_comp(truth, estimate = estimate, summarize = summarize)
    )
  } else {
    stop("metric must be one of: rmse, mae, rho, kld, mlr, all")
  }
}

rmse_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (truth - estimate)^2 |>
    apply(2, mean, na.rm = TRUE) |>
    sqrt()
  if (summarize) {
    tibble::tibble(.metric = "rmse", .estimate = mean(m, na.rm = TRUE))
  } else {
    tibble::tibble(.metric = "rmse", .estimate = m) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::pivot_wider(
        names_from = id,
        values_from = .estimate,
        names_prefix = ".estimate_"
      )
  }
}

#' Mean Absolute Error (MAE) for compositional data
#'
#' This function computes the MAE between true and estimated values for compositional data.
#'
#' @param truth A matrix or data frame of true values.
#' @param estimate A matrix or data frame of estimated values.
#' @param summarize Logical. If `TRUE`, returns a single summarized MAE value.
#'   If `FALSE`, returns MAE values for each component.
#'
#' @return A tibble containing the MAE values. If `summarize = TRUE`, the tibble contains
#'   one row with the mean MAE. If `summarize = FALSE`, the tibble contains detailed
#'   MAE values for each component.
#' @importFrom rlang .data
mae_comp <- function(truth, estimate, summarize = FALSE) {
  m <- abs(truth - estimate) |>
    apply(2, mean, na.rm = TRUE)
  if (summarize) {
    tibble::tibble(.metric = "mae", .estimate = mean(m, na.rm = TRUE))
  } else {
    tibble::tibble(.metric = "mae", .estimate = m) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::pivot_wider(
        names_from = id,
        values_from = .estimate,
        names_prefix = ".estimate_"
      )
  }
}

#' Spearman Rank Correlation (rho) for compositional data
#'
#' This function computes the Spearman rank correlation between true and estimated values
#' for compositional data.
#'
#' @param truth A matrix or data frame of true values.
#' @param estimate A matrix or data frame of estimated values.
#' @param summarize Logical. If `TRUE`, returns a single summarized correlation value.
#'   If `FALSE`, returns correlation values for each component.
#'
#' @return A tibble containing the correlation values. If `summarize = TRUE`, the tibble
#'   contains one row with the mean correlation. If `summarize = FALSE`, the tibble contains
#'   detailed correlation values for each component.
#' @importFrom rlang .data
rho_comp <- function(truth, estimate, summarize = FALSE) {
  m <- stats::cor(
    truth,
    estimate,
    use = "complete.obs",
    method = "spearman"
  ) |>
    diag()
  if (summarize) {
    tibble::tibble(.metric = "rho", .estimate = mean(m, na.rm = TRUE))
  } else {
    tibble::tibble(.metric = "rho", .estimate = m) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::pivot_wider(
        names_from = id,
        values_from = .estimate,
        names_prefix = ".estimate_"
      )
  }
}

#' Kullback-Leibler Divergence (KLD) for compositional data
#'
#' This function computes the KLD between true and estimated values for compositional data.
#' **Note**: This is experimental and has not been thoroughly tested or verified.
#'
#' @param truth A matrix or data frame of true values.
#' @param estimate A matrix or data frame of estimated values.
#' @param summarize Logical. If `TRUE`, returns a single summarized KLD value.
#'   If `FALSE`, returns KLD values for each component.
#'
#' @return A tibble containing the KLD values. If `summarize = TRUE`, the tibble contains
#'   one row with the mean KLD. If `summarize = FALSE`, the tibble contains detailed
#'   KLD values for each component.
#' @importFrom rlang .data
kld_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (truth * log(truth / estimate)) |>
    apply(2, sum, na.rm = TRUE)
  if (summarize) {
    tibble::tibble(.metric = "kld", .estimate = mean(m, na.rm = TRUE))
  } else {
    tibble::tibble(.metric = "kld", .estimate = m) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::pivot_wider(
        names_from = id,
        values_from = .estimate,
        names_prefix = ".estimate_"
      )
  }
}

#' Mean Log Ratio (MLR) for compositional data
#'
#' This function computes the MLR between true and estimated values for compositional data.
#' **Note**: This is experimental, has not been thoroughly tested or verified, and lacks
#' strong theoretical justification.
#'
#' @param truth A matrix or data frame of true values.
#' @param estimate A matrix or data frame of estimated values.
#' @param summarize Logical. If `TRUE`, returns a single summarized MLR value.
#'   If `FALSE`, returns MLR values for each component.
#'
#' @return A tibble containing the MLR values. If `summarize = TRUE`, the tibble contains
#'   one row with the mean MLR. If `summarize = FALSE`, the tibble contains detailed
#'   MLR values for each component.
#' @importFrom rlang .data
mlr_comp <- function(truth, estimate, summarize = FALSE) {
  m <- log(estimate / truth) |>
    apply(2, mean, na.rm = TRUE)
  if (summarize) {
    tibble::tibble(.metric = "mlr", .estimate = mean(m, na.rm = TRUE))
  } else {
    tibble::tibble(.metric = "mlr", .estimate = m) |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidyr::pivot_wider(
        names_from = id,
        values_from = .estimate,
        names_prefix = ".estimate_"
      )
  }
}
