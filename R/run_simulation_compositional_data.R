#' Run repeated spatial sampling and compositional model simulation
#'
#' @description
#' Performs a Monte Carlo simulation using [SimDesign::runSimulation()] that,
#' for each replication:
#' 1. Samples spatial locations from a (Packed) `SpatRaster` via
#'    [rastersample::spatial_sample()].
#' 2. Generates one compositional response per sampled row from `.alpha_sim*`
#'    layers (see [fit_compositional_data()]) and fits a Dirichlet regression
#'    with [DirichletReg::DirichReg()].
#' 3. Predicts composition probabilities `.p_hat*` over the full raster (all
#'    cells converted to a data frame for prediction).
#' 4. Computes component-wise performance metrics against ground truth
#'    probabilities `.p_sim*` using [metrics_comp()] (RMSE, MAE, Spearman rho,
#'    KLD, MLR) plus the Aitchison distance (robCompositions::aDist).
#'
#' @details
#' Expected raster layers:
#' - `.alpha_sim1`, `.alpha_sim2`, ... : Dirichlet parameters used to simulate
#'   the sampled response internally (not used directly in prediction).
#' - `.p_sim1`, `.p_sim2`, ... : Ground-truth composition probabilities used
#'   for metric calculation.
#' - Predictor (covariate) layers referenced in the model `formula`.
#'
#' The returned object is the list produced by [SimDesign::runSimulation()],
#' containing per-variable summaries (currently only the mean/bias relative to
#' 0 via `SimDesign::bias()`) for each metric-component combination:
#' `rmse_1`, `rmse_2`, ..., `mae_1`, ..., `rho_1`, ..., `kld_1`, ..., `mlr_1`,
#' ... and `adist` (Aitchison distance). Each variable's replication-level
#' values are the raw metric outputs; summarization reduces them to a single
#' `mean` column. (No additional distribution summaries are computed.)
#'
#' Parallel execution (if `parallel = TRUE`) uses the `"future"` backend and
#' sets a multisession plan from `future.mirai`.
#'
#' Note: The `sites` argument is currently only reprojected to match `data`
#' but otherwise unused in the simulation (placeholder for future stratified or
#' fixed-location designs).
#'
#' @param data A (Packed) `SpatRaster` containing `.alpha_sim*`, `.p_sim*` and
#' predictor layers required by `formula`.
#' @param sites Optional `sf` POINT layer; reprojected to `data` CRS (currently
#' unused in sampling logic).
#' @param formula A Dirichlet regression formula passed to
#' [DirichletReg::DirichReg()]. Response (`y`) is constructed internally.
#' @param n Integer; sample size per replication.
#' @param method Sampling method for [rastersample::spatial_sample()].
#' @param bias_var Optional character; name of a layer used for biased sampling
#' (passed through to `rastersample::spatial_sample()`).
#' @param bias_thresh Optional numeric; threshold for biased sampling
#' (passed through to `rastersample::spatial_sample()`).
#' @param clh_var Optional character; layer(s) used for conditioned Latin
#' hypercube sampling (passed through to `rastersample::spatial_sample()`).
#' @param clh_iter Optional integer; number of iterations for conditioned Latin
#' hypercube sampling (passed through to `rastersample::spatial_sample()`).
#' @param strata_var Optional character; name of a layer used for stratified
#' sampling (passed through to `rastersample::spatial_sample()`).
#' @param replications Integer; number of Monte Carlo replications.
#' @param verbose Logical; reserved (currently not producing messages).
#' @param parallel Logical; if `TRUE`, enables future-based parallelism.
#' @param n_cores Integer; number of workers (ignored if `parallel = FALSE`).
#'
#' @return A list (class `"simulation"`) from [SimDesign::runSimulation()] with:
#' - `$results`: replication-level metric-component values.
#' - `$summaries`: a data frame of per-variable mean (bias vs 0).
#' - Other bookkeeping elements from SimDesign.
#'
#' @seealso
#' [simulate_compositional_data()],
#' [fit_compositional_data()],
#' [predict_compositional_data()],
#' [metrics_comp()],
#' [rastersample::spatial_sample()],
#' [SimDesign::runSimulation()]
#'
#' @examples
#' \dontrun{
#' f <- system.file("ex/elev.tif", package = "terra")
#' r <- terra::rast(f)
#' r <- c(r, MultiscaleDTM::BPI(r, w = c(2, 4)))
#' r_prep <- prepare_data(r, poly_degree = 3)
#' sim_data <- simulate_compositional_data(
#'   x = r_prep,
#'   d = 3,
#'   n_cov_sim = 3,
#'   as_raster = TRUE,
#'   seed = 123
#' )
#' dat <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data))
#' res <- run_simulation_compositional_data(
#'   dat,
#'   formula = y ~ 1 + elevation_poly1 + elevation_poly2,
#'   n = c(20, 50, 100, 500, 1000, 5000, 10000),
#'   method = "random",
#'   replications = 10,
#'   parallel = FALSE
#' )
#' ggplot2::ggplot(res, mapping = ggplot2::aes(x = n, y = mean.adist)) +
#'   ggplot2::geom_point() +
#'   ggplot2::geom_line()
#' }
#'
#' @export
run_simulation_compositional_data <- function(
  data,
  sites = NULL,
  formula,
  n,
  method,
  bias_var = NULL,
  bias_thresh = NULL,
  clh_var = NULL,
  clh_iter = NULL,
  strata_var = NULL,
  replications = 1L,
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
      bias_var = bias_var,
      bias_thresh = bias_thresh,
      clh_var = clh_var,
      clh_iter = clh_iter,
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
        names_to = ".comp",
        names_prefix = ".estimate_",
        values_to = ".estimate"
      )
    adist <- robCompositions::aDist(
      dplyr::select(new_data, tidyselect::starts_with(".p_sim")),
      y = pred
    )
    # adist <- adist / nrow(new_data) # average Aitchison distance
    ret <- c(metrics$.estimate, adist)
    names(ret) <- c(paste0(metrics$.metric, "_", metrics$.comp), "adist")
    ret
  }
  summarise_data <- function(condition, results, fixed_objects) {
    c(mean = SimDesign::bias(results, parameter = 0))
  }
  if (parallel) {
    parallel <- "future"
    if (n_cores == 1L) {
      n_cores <- future::availableCores()
    }
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
