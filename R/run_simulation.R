#' Run simulation
#'
#' @param sim_data description
#' @param sites description
#' @param formula description
#' @param n description
#' @param method description
#' @param strata_var description
#' @param num_sim description
#' @param use_dirinla description
#' @param tol0 description
#' @param verbose description
#'
#' @return description
#' @importFrom rlang .data
#' @export
run_simulation <- function(sim_data, sites = NULL, formula, n, method,
                           strata_var = NULL, num_sim, use_dirinla = TRUE,
                           tol0 = NULL, verbose = FALSE) {

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
    num_sim = num_sim
  )

  # specify simulation script
  sim_script <- function() {
    samp_str <- sim_data |>
      rastersample::spatial_sample(n = L$n, method = L$method,
                                   strata_var = L$strata_var, drop_na = TRUE)
    m_samp <- samp_str |>
      fit_comp_dirichlet(formula = formula, use_dirinla = use_dirinla,
                         tol0 = tol0, verbose = verbose)
    pred_samp <- predict_comp_dirichlet(sim_data, model = m_samp) |>
      terra::unwrap()
    res_samp <- terra::as.data.frame(pred_samp) |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::starts_with("p_")) |>
      tidyr::drop_na() |>
      pipebind::bind(._, metrics_comp(truth = dplyr::select(._, dplyr::starts_with("p_sim")),
                                      estimate = dplyr::select(._, dplyr::starts_with("p_hat")),
                                      summarize = FALSE)) |>
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
      ".complex" = list("performance_metrics" = res_samp)
    )
  }

  # set simulation script
  sim <- sim |> SimEngine::set_script(sim_script())

  # run simulations
  sim |> SimEngine::run()

}
