#' Run compositional model simulation
#'
#' @description Execute repeated spatial sampling, model fitting, prediction,
#' and metric computation using [SimEngine]. Returns a `SimEngine` simulation
#' object with summarized metrics per iteration.
#'
#' @param data (Packed) `SpatRaster` containing simulated `.p_sim*` layers
#' and covariates used for fitting.
#' @param sites Optional `sf` point object for site locations (transformed to
#' `data` CRS).
#' @param formula Model formula passed to [DirichletReg::DirichReg()].
#' @param n Sample size per iteration.
#' @param method Sampling method passed to [rastersample::spatial_sample()].
#' @param strata_var Optional stratification variable name (character).
#' @param replications Number of simulation replications.
#' @param verbose Logical; if TRUE prints fitting progress (may slow runs).
#' @param parallel Logical; if TRUE enables parallel execution via SimEngine.
#' @param n_cores Integer; number of cores (ignored if `parallel = FALSE`).
#'
#' @return A `SimEngine` simulation result (after `SimEngine::run()`), containing
#' metric summaries (`rmse_*`, `mae_*`, `rho_*`, `kld_*`) and a `.complex`
#' list with full per-component metric tibble.
#'
#' @details
#' Each iteration:
#' 1. Sample spatial locations.
#' 2. Fit Dirichlet regression (using internal fit/predict wrappers).
#' 3. Predict compositions for full raster.
#' 4. Compute metrics via [metrics_comp()] (component-wise + summary stats).
#'
#' @seealso [simulate_compositional_data()], [fit_compositional_data()],
#' [predict_compositional_data()], [metrics_comp()]
#'
#' @examples
#' \donttest{
#' # Requires SimEngine and rastersample packages.
#' # sim_obj <- run_simulation(sim_data = sim$data, formula = y ~ 1,
#' #                           n = 50, method = "random",
#' #                           replications = 10, parallel = FALSE)
#' }
#'
#' @export
run_simulation_compositional_data <- function(
  data,
  sites = NULL,
  formula,
  n,
  method,
  strata_var = NULL,
  replications,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 1L
) {
  if (inherits(data, "SpatRaster")) {
    data <- terra::wrap(data)
  }
  if (!is.null(sites)) {
    crs_ref <- terra::crs(terra::unwrap(data))
    sites <- sf::st_transform(sites, crs = crs_ref)
  }
  design <- if (is.null(strata_var)) {
    SimDesign::createDesign(n = n, method = method)
  } else {
    SimDesign::createDesign(n = n, method = method, strata_var = strata_var)
  }
  generate_data <- function(condition, fixed_objects) {
    n <- condition$n
    method <- condition$method
    data <- fixed_objects$data
    if (inherits(data, "PackedSpatRaster")) {
      data <- terra::unwrap(data)
    }
    rastersample::spatial_sample(
      data,
      n = n,
      method = method,
      strata_var = strata_var,
      drop_na = TRUE
    )
  }
  analyse_data <- function(condition, dat, fixed_objects) {
    new_data <- fixed_objects$data
    if (inherits(new_data, "PackedSpatRaster")) {
      new_data <- terra::unwrap(new_data)
    }
    if (inherits(new_data, "SpatRaster")) {
      new_data <- terra::as.data.frame(new_data, na.rm = TRUE)
    }
    sample_fit <- fit_compositional_data(fixed_objects$formula, data = dat)
    pred <- predict_compositional_data(sample_fit, new_data = new_data)
    metrics <- metrics_comp(
      truth = dplyr::select(new_data, dplyr::starts_with(".p_sim")),
      estimate = pred,
      summarize = FALSE
    ) |>
      tidyr::pivot_longer(
        cols = dplyr::starts_with(".estimate_"),
        names_to = ".class",
        names_prefix = ".estimate_",
        values_to = ".estimate"
      )
    adist <- robCompositions::aDist(
      dplyr::select(new_data, tidyselect::starts_with(".p_sim")),
      y = pred
    )
    ret <- c(metrics$.estimate, adist)
    names(ret) <- c(paste0(metrics$.metric, "_", metrics$.class), "adist")
    ret
  }
  summarise_data <- function(condition, results, fixed_objects) {
    c(mean = SimDesign::bias(results, parameter = 0))
  }
  if (parallel) {
    parallel <- "future"
    future::plan(future.mirai::mirai_multisession, workers = n_cores)
  }
  SimDesign::runSimulation(
    design,
    replications = replications,
    generate = generate_data,
    analyse = analyse_data,
    summarise = summarise_data,
    fixed_objects = list(data = data, formula = formula),
    save = FALSE,
    parallel = parallel
  )
  # SimEngine::set_script(
  #   function() {
  #     samp_str <- sim_data |>
  #       rastersample::spatial_sample(
  #         n = L$n,
  #         method = L$method,
  #         strata_var = L$strata_var,
  #         drop_na = TRUE
  #       )
  #     m_samp <- fit_compositional_data(
  #       samp_str,
  #       formula = formula
  #     )
  #     pred_samp <- predict_compositional_data(
  #       m_samp,
  #       sim_data
  #     ) |>
  #       terra::unwrap()
  #     res_samp <- terra::as.data.frame(pred_samp) |>
  #       tibble::as_tibble() |>
  #       dplyr::select(dplyr::starts_with("p_")) |>
  #       tidyr::drop_na() |>
  #       pipebind::bind(
  #         ._,
  #         metrics_comp(
  #           truth = dplyr::select(._, dplyr::starts_with("p_sim")),
  #           estimate = dplyr::select(._, dplyr::starts_with("p_hat")),
  #           summarize = FALSE
  #         )
  #       ) |>
  #       dplyr::rowwise() |>
  #       dplyr::mutate(
  #         .estimate_mean = mean(dplyr::c_across(dplyr::starts_with(
  #           ".estimate"
  #         ))),
  #         .estimate_min = min(dplyr::c_across(dplyr::starts_with(
  #           ".estimate"
  #         ))),
  #         .estimate_max = max(dplyr::c_across(dplyr::starts_with(
  #           ".estimate"
  #         )))
  #       )
  #     list(
  #       rmse_mean = res_samp |>
  #         dplyr::filter(.metric == "rmse") |>
  #         dplyr::pull(.estimate_mean),
  #       rmse_min = res_samp |>
  #         dplyr::filter(.metric == "rmse") |>
  #         dplyr::pull(.estimate_min),
  #       rmse_max = res_samp |>
  #         dplyr::filter(.metric == "rmse") |>
  #         dplyr::pull(.estimate_max),
  #       mae_mean = res_samp |>
  #         dplyr::filter(.metric == "mae") |>
  #         dplyr::pull(.estimate_mean),
  #       mae_min = res_samp |>
  #         dplyr::filter(.metric == "mae") |>
  #         dplyr::pull(.estimate_min),
  #       mae_max = res_samp |>
  #         dplyr::filter(.metric == "mae") |>
  #         dplyr::pull(.estimate_max),
  #       rho_mean = res_samp |>
  #         dplyr::filter(.metric == "rho") |>
  #         dplyr::pull(.estimate_mean),
  #       rho_min = res_samp |>
  #         dplyr::filter(.metric == "rho") |>
  #         dplyr::pull(.estimate_min),
  #       rho_max = res_samp |>
  #         dplyr::filter(.metric == "rho") |>
  #         dplyr::pull(.estimate_max),
  #       kld_mean = res_samp |>
  #         dplyr::filter(.metric == "kld") |>
  #         dplyr::pull(.estimate_mean),
  #       kld_min = res_samp |>
  #         dplyr::filter(.metric == "kld") |>
  #         dplyr::pull(.estimate_min),
  #       kld_max = res_samp |>
  #         dplyr::filter(.metric == "kld") |>
  #         dplyr::pull(.estimate_max),
  #       .complex = list(performance_metrics = res_samp)
  #     )
  #   }
  # )
  # SimEngine::run(sim)
}
