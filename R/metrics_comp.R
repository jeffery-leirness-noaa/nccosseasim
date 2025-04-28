#' Performance metrics for compositional data
#'
#' This function computes various performance metrics for compositional data, such as 
#' root mean square error (RMSE), mean absolute error (MAE), Spearman rank correlation (rho), 
#' Kullback-Leibler divergence (KLD), and mean log ratio (MLR). It can compute a single metric 
#' or all metrics at once.
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
    stop("Invalid metric argument. Must be one of: 'rmse', 'mae', 'rho', 'kld', 'mlr', or 'all'")
  }
}

#' Root Mean Square Error (RMSE) for compositional data
#'
#' This function computes the RMSE between true and estimated values for compositional data.
#'
#' @param truth A matrix or data frame of true values.
#' @param estimate A matrix or data frame of estimated values.
#' @param summarize Logical. If `TRUE`, returns a single summarized RMSE value. 
#'   If `FALSE`, returns RMSE values for each component.
#'
#' @return A tibble containing the RMSE values. If `summarize = TRUE`, the tibble contains 
#'   one row with the mean RMSE. If `summarize = FALSE`, the tibble contains detailed 
#'   RMSE values for each component.
#' @importFrom rlang .data
rmse_comp <- function(truth, estimate, summarize = FALSE) {
  m <- (truth - estimate) ^ 2 |>
    apply(MARGIN = 2, FUN = mean) |>
    sqrt()
  if (summarize) {
    tibble::tibble(.metric = "rmse", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "rmse", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = .data$id, values_from = .data$.estimate,
                         names_prefix = ".estimate_")
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
  m <- (truth - estimate) |>
    abs() |>
    apply(MARGIN = 2, FUN = mean)
  if (summarize) {
    tibble::tibble(.metric = "mae", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "mae", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = .data$id, values_from = .data$.estimate,
                         names_prefix = ".estimate_")
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
  m <- stats::cor(truth, estimate, use = "complete.obs", method = "spearman") |>
    diag()
  if (summarize) {
    tibble::tibble(.metric = "rho", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "rho", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = .data$id, values_from = .data$.estimate,
                         names_prefix = ".estimate_")
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
    apply(MARGIN = 2, FUN = sum, na.rm = TRUE)
  if (summarize) {
    tibble::tibble(.metric = "kld", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "kld", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = .data$id, values_from = .data$.estimate,
                         names_prefix = ".estimate_")
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
  m <- (estimate / truth) |>
    log() |>
    apply(MARGIN = 2, FUN = mean)
  if (summarize) {
    tibble::tibble(.metric = "mlr", .estimate = mean(m))
  } else {
    tibble::tibble(.metric = "mlr", .estimate = m) |>
      dplyr::mutate(id = 1:dplyr::n()) |>
      tidyr::pivot_wider(names_from = .data$id, values_from = .data$.estimate,
                         names_prefix = ".estimate_")
  }
}
