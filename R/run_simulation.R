#' Run compositional model simulation
#'
#' @description Execute repeated spatial sampling, model fitting, prediction,
#' and metric computation using [SimEngine]. Returns a `SimEngine` simulation
#' object with summarized metrics per iteration.
#'
#' @param sim_data (Packed) `SpatRaster` containing simulated `.p_sim*` layers
#' and covariates used for fitting.
#' @param sites Optional `sf` point object for site locations (transformed to
#' `sim_data` CRS).
#' @param formula Model formula passed to [DirichletReg::DirichReg()].
#' @param n Sample size per iteration.
#' @param method Sampling method passed to [rastersample::spatial_sample()].
#' @param strata_var Optional stratification variable name (character).
#' @param num_sim Number of simulation iterations.
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
#' #                           num_sim = 10, parallel = FALSE)
#' }
#'
#' @export
run_simulation <- function(
  sim_data,
  sites = NULL,
  formula,
  n,
  method,
  strata_var = NULL,
  num_sim,
  verbose = FALSE,
  parallel = FALSE,
  n_cores = 1L
) {
  if (inherits(sim_data, "PackedSpatRaster")) {
    sim_data <- terra::unwrap(sim_data)
  }
  if (!is.null(sites)) {
    sites <- sf::st_transform(sites, crs = terra::crs(sim_data))
  }
  sim <- SimEngine::new_sim()
  sim <- sim |>
    SimEngine::set_levels(
      n = n,
      method = method,
      strata_var = strata_var
    ) |>
    SimEngine::set_config(
      num_sim = num_sim,
      parallel = parallel,
      n_cores = n_cores
    ) |>
    SimEngine::set_script(
      function() {
        samp_str <- sim_data |>
          rastersample::spatial_sample(
            n = L$n,
            method = L$method,
            strata_var = L$strata_var,
            drop_na = TRUE
          )
        m_samp <- nccosseasim::fit_compositional_data(
          samp_str,
          formula = formula
        )
        pred_samp <- nccosseasim::predict_compositional_data(
          m_samp,
          sim_data
        ) |>
          terra::unwrap()
        res_samp <- terra::as.data.frame(pred_samp) |>
          tibble::as_tibble() |>
          dplyr::select(dplyr::starts_with("p_")) |>
          tidyr::drop_na() |>
          pipebind::bind(
            ._,
            nccosseasim::metrics_comp(
              truth = dplyr::select(._, dplyr::starts_with("p_sim")),
              estimate = dplyr::select(._, dplyr::starts_with("p_hat")),
              summarize = FALSE
            )
          ) |>
          dplyr::rowwise() |>
          dplyr::mutate(
            .estimate_mean = mean(dplyr::c_across(dplyr::starts_with(
              ".estimate"
            ))),
            .estimate_min = min(dplyr::c_across(dplyr::starts_with(
              ".estimate"
            ))),
            .estimate_max = max(dplyr::c_across(dplyr::starts_with(
              ".estimate"
            )))
          )
        list(
          rmse_mean = res_samp |>
            dplyr::filter(.metric == "rmse") |>
            dplyr::pull(.estimate_mean),
          rmse_min = res_samp |>
            dplyr::filter(.metric == "rmse") |>
            dplyr::pull(.estimate_min),
          rmse_max = res_samp |>
            dplyr::filter(.metric == "rmse") |>
            dplyr::pull(.estimate_max),
          mae_mean = res_samp |>
            dplyr::filter(.metric == "mae") |>
            dplyr::pull(.estimate_mean),
          mae_min = res_samp |>
            dplyr::filter(.metric == "mae") |>
            dplyr::pull(.estimate_min),
          mae_max = res_samp |>
            dplyr::filter(.metric == "mae") |>
            dplyr::pull(.estimate_max),
          rho_mean = res_samp |>
            dplyr::filter(.metric == "rho") |>
            dplyr::pull(.estimate_mean),
          rho_min = res_samp |>
            dplyr::filter(.metric == "rho") |>
            dplyr::pull(.estimate_min),
          rho_max = res_samp |>
            dplyr::filter(.metric == "rho") |>
            dplyr::pull(.estimate_max),
          kld_mean = res_samp |>
            dplyr::filter(.metric == "kld") |>
            dplyr::pull(.estimate_mean),
          kld_min = res_samp |>
            dplyr::filter(.metric == "kld") |>
            dplyr::pull(.estimate_min),
          kld_max = res_samp |>
            dplyr::filter(.metric == "kld") |>
            dplyr::pull(.estimate_max),
          .complex = list(performance_metrics = res_samp)
        )
      }
    )
  SimEngine::run(sim)
}
