#' Performance metrics for compositional data
#'
#' @param truth description
#' @param estimate description
#' @param metric description
#' @param summarize description
#'
#' @return description
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
# root mean square error
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
# mean absolute error
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
# spearman rank correlation
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
# kullback-leibler divergence
# caution: this is experimental and has not been tested/verified
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
# mean log ratio
# caution: this is experimental, has not been tested/verified, nor has any real theory behind its use
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
