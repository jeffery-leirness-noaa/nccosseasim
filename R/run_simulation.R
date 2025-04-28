#' Run simulation
#'
#' This function runs a simulation to evaluate the performance of compositional data models. 
#' It uses spatial sampling, fits a Dirichlet regression model, predicts compositional probabilities, 
#' and computes performance metrics such as RMSE, MAE, Spearman rank correlation, and KLD.
#'
#' @param sim_data A `SpatRaster` or `PackedSpatRaster` object containing the simulated data.
#'   If it is a `PackedSpatRaster`, it will be unpacked using `terra::unwrap()`.
#' @param sites An optional `sf` object containing spatial site locations. If provided, 
#'   it will be transformed to match the coordinate reference system of `sim_data`.
#' @param formula A formula specifying the regression model to be fitted.
#' @param n Integer. The number of samples to draw during spatial sampling.
#' @param method Character. The sampling method to use (e.g., "random", "stratified").
#' @param strata_var Character. The name of the variable to use for stratified sampling (optional).
#' @param num_sim Integer. The number of simulation iterations to run.
#' @param use_dirinla Logical. If `TRUE`, uses the `dirinla::dirinlareg()` function for model fitting; 
#'   otherwise, uses `DirichletReg::DirichReg()`.
#' @param tol0 Numeric. A tolerance parameter for `dirinla::dirinlareg()`. Ignored if `use_dirinla` is `FALSE`.
#' @param verbose Logical. If `TRUE`, prints additional information during model fitting and simulation.
#' @param parallel Logical. If `TRUE`, runs simulations in parallel.
#' @param n_cores Integer. The number of cores to use for parallel processing. Defaults to 1.
#'
#' @return A simulation object containing the results of the simulation, including performance metrics 
#'   (e.g., RMSE, MAE, Spearman rank correlation, KLD) for each iteration. The results are stored in 
#'   a structured format for further analysis.
#'
#' @details The function performs the following steps:
#'   1. Unpacks `PackedSpatRaster` objects and ensures coordinate reference system consistency.
#'   2. Initializes a simulation object using the `SimEngine` package.
#'   3. Configures simulation levels and parameters, including the number of iterations and sampling method.
#'   4. Defines a simulation script that performs spatial sampling, model fitting, prediction, and metric computation.
#'   5. Runs the simulation and returns the results.
#'
#' @importFrom rlang .data
#' @export
run_simulation <- function(sim_data, sites = NULL, formula, n, method,
                           strata_var = NULL, num_sim, use_dirinla = TRUE,
                           tol0 = NULL, verbose = FALSE, parallel = FALSE,
                           n_cores = 1L) {

  # if sim_data is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(sim_data, "PackedSpatRaster")) sim_data <- terra::unwrap(sim_data)

  # ensure sites and sim_data are using the same coordinate reference system
  if (!is.null(sites)) {
    sites <- sf::st_transform(sites, crs = terra::crs(sim_data))
  }

  # initialize simulation object
  sim <- SimEngine::new_sim()

  # set simulation levels
  sim <- sim |> SimEngine::set_levels(
    n = n,
    method = method
  )

  # set simulation configuration
  sim <- sim |> SimEngine::set_config(
    num_sim = num_sim,
    parallel = parallel,
    n_cores
  )

  # specify simulation script
  # sim_script <- function() {
  #   samp_str <- sim_data |>
  #     rastersample::spatial_sample(n = L$n, method = L$method,
  #                                  strata_var = L$strata_var, drop_na = TRUE)
  #   m_samp <- samp_str |>
  #     fit_comp_dirichlet(formula = formula, use_dirinla = use_dirinla,
  #                        tol0 = tol0, verbose = verbose)
  #   pred_samp <- predict_comp_dirichlet(sim_data, model = m_samp) |>
  #     terra::unwrap()
  #   res_samp <- terra::as.data.frame(pred_samp) |>
  #     tibble::as_tibble() |>
  #     dplyr::select(dplyr::starts_with("p_")) |>
  #     tidyr::drop_na() |>
  #     pipebind::bind(._, metrics_comp(truth = dplyr::select(._, dplyr::starts_with("p_sim")),
  #                                     estimate = dplyr::select(._, dplyr::starts_with("p_hat")),
  #                                     summarize = FALSE)) |>
  #     dplyr::rowwise() |>
  #     dplyr::mutate(.estimate_mean = mean(dplyr::c_across(dplyr::starts_with(".estimate"))),
  #                   .estimate_min = min(dplyr::c_across(dplyr::starts_with(".estimate"))),
  #                   .estimate_max = max(dplyr::c_across(dplyr::starts_with(".estimate"))))
  #   list(
  #     "rmse_mean" = res_samp |>
  #       dplyr::filter(.data$.metric == "rmse") |>
  #       dplyr::pull(.data$.estimate_mean),
  #     "rmse_min" = res_samp |>
  #       dplyr::filter(.data$.metric == "rmse") |>
  #       dplyr::pull(.data$.estimate_min),
  #     "rmse_max" = res_samp |>
  #       dplyr::filter(.data$.metric == "rmse") |>
  #       dplyr::pull(.data$.estimate_max),
  #     "mae_mean" = res_samp |>
  #       dplyr::filter(.data$.metric == "mae") |>
  #       dplyr::pull(.data$.estimate_mean),
  #     "mae_min" = res_samp |>
  #       dplyr::filter(.data$.metric == "mae") |>
  #       dplyr::pull(.data$.estimate_min),
  #     "mae_max" = res_samp |>
  #       dplyr::filter(.data$.metric == "mae") |>
  #       dplyr::pull(.data$.estimate_max),
  #     "rho_mean" = res_samp |>
  #       dplyr::filter(.data$.metric == "rho") |>
  #       dplyr::pull(.data$.estimate_mean),
  #     "rho_min" = res_samp |>
  #       dplyr::filter(.data$.metric == "rho") |>
  #       dplyr::pull(.data$.estimate_min),
  #     "rho_max" = res_samp |>
  #       dplyr::filter(.data$.metric == "rho") |>
  #       dplyr::pull(.data$.estimate_max),
  #     ".complex" = list("performance_metrics" = res_samp)
  #   )
  # }

  # set simulation script
  sim <- sim |> SimEngine::set_script(
    function() {
      samp_str <- sim_data |>
        rastersample::spatial_sample(n = L$n, method = L$method,
                                     strata_var = L$strata_var, drop_na = TRUE)
      m_samp <- samp_str |>
        nccosseasim::fit_comp_dirichlet(formula = formula,
                                        use_dirinla = use_dirinla, tol0 = tol0,
                                        verbose = verbose)
      pred_samp <- nccosseasim::predict_comp_dirichlet(sim_data, model = m_samp) |>
        terra::unwrap()
      res_samp <- terra::as.data.frame(pred_samp) |>
        tibble::as_tibble() |>
        dplyr::select(dplyr::starts_with("p_")) |>
        tidyr::drop_na() |>
        pipebind::bind(._, nccosseasim::metrics_comp(
          truth = dplyr::select(._, dplyr::starts_with("p_sim")),
          estimate = dplyr::select(._, dplyr::starts_with("p_hat")),
          summarize = FALSE)
        ) |>
        dplyr::rowwise() |>
        dplyr::mutate(.estimate_mean = mean(dplyr::c_across(dplyr::starts_with(".estimate"))),
                      .estimate_min = min(dplyr::c_across(dplyr::starts_with(".estimate"))),
                      .estimate_max = max(dplyr::c_across(dplyr::starts_with(".estimate"))))
      list(
        "rmse_mean" = res_samp |>
          dplyr::filter(.data$.metric == "rmse") |>
          dplyr::pull(.data$.estimate_mean),
        "rmse_min" = res_samp |>
          dplyr::filter(.data$.metric == "rmse") |>
          dplyr::pull(.data$.estimate_min),
        "rmse_max" = res_samp |>
          dplyr::filter(.data$.metric == "rmse") |>
          dplyr::pull(.data$.estimate_max),
        "mae_mean" = res_samp |>
          dplyr::filter(.data$.metric == "mae") |>
          dplyr::pull(.data$.estimate_mean),
        "mae_min" = res_samp |>
          dplyr::filter(.data$.metric == "mae") |>
          dplyr::pull(.data$.estimate_min),
        "mae_max" = res_samp |>
          dplyr::filter(.data$.metric == "mae") |>
          dplyr::pull(.data$.estimate_max),
        "rho_mean" = res_samp |>
          dplyr::filter(.data$.metric == "rho") |>
          dplyr::pull(.data$.estimate_mean),
        "rho_min" = res_samp |>
          dplyr::filter(.data$.metric == "rho") |>
          dplyr::pull(.data$.estimate_min),
        "rho_max" = res_samp |>
          dplyr::filter(.data$.metric == "rho") |>
          dplyr::pull(.data$.estimate_max),
        "kld_mean" = res_samp |>
          dplyr::filter(.data$.metric == "kld") |>
          dplyr::pull(.data$.estimate_mean),
        "kld_min" = res_samp |>
          dplyr::filter(.data$.metric == "kld") |>
          dplyr::pull(.data$.estimate_min),
        "kld_max" = res_samp |>
          dplyr::filter(.data$.metric == "kld") |>
          dplyr::pull(.data$.estimate_max),
        ".complex" = list("performance_metrics" = res_samp)
      )
    }
  )

  # run simulations
  sim |> SimEngine::run()

}
