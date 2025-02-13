#' Simulate compositional data from a Dirichlet likelihood
#'
#' @param x description
#' @param d description
#' @param layers description
#' @param poly_degree description
#' @param n_cov_sim description
#' @param as_raster description
#' @param seed description
#'
#' @return description
#' @export
simulate_comp_data <- function(x, d, layers = NULL, poly_degree = NULL,
                               n_cov_sim = NULL, as_raster = FALSE,
                               seed = NULL) {

  # set random number generator seed
  if (is.numeric(seed)) set.seed(seed)

  # if x is a PackedSpatRaster object, use terra::unwrap() to unpack it
  if (inherits(x, "PackedSpatRaster")) x <- terra::unwrap(x)

  # create covariate polynomial rasters
  if (is.null(layers) & !is.null(poly_degree)) {
    if (poly_degree > 1) {
      x_add <- purrr::map(2:poly_degree, \(deg) x^deg |>
                            stats::setNames(nm = stringr::str_c(names(x), "^", deg))) |>
        terra::rast()
      x <- c(x, x_add)
    }
    if (!is.null(n_cov_sim)) {
      layers <- sample.int(terra::nlyr(x), size = n_cov_sim) |>
        sort()
    } else {
      layers <- 1:terra::nlyr(x)
    }
  }

  # specify model formula structure that generates the simulated data
  if (is.character(layers)) {
    layer_names <- layers
  } else if (is.numeric(layers)) {
    layer_names <- names(x)[layers]
  }
  form_sim <- layer_names |>
    stringr::str_flatten(collapse = " + ") |>
    pipebind::bind(._, stringr::str_c("1 + ", ._)) |>
    rep.int(times = d) |>
    stringr::str_flatten(collapse = " | ") |>
    pipebind::bind(._, stringr::str_c("y ~ ", ._)) |>
    stats::as.formula()


  # simulate compositional data ---------------------------------------------

  # # method 1: create SpatRaster object containing simulated substrate composition
  # # probabilities
  # # linear effects formulation
  # x_sim <- stats::rnorm(d * 2, mean = 0, sd = 2) |>
  #   matrix(nrow = d) |>
  #   tibble::as_tibble(.name_repair = ~ paste0("beta_", 0:1))
  #
  # # create simulated substrate composition probabilities as SpatRaster objects
  # r_sim <- terra::subset(x, subset = layers) * x_sim$beta_1 + x_sim$beta_0
  # alpha_sim <- exp(r_sim)
  # names(alpha_sim) <- paste0("alpha_sim", 1:terra::nlyr(alpha_sim))
  # r_sim <- alpha_sim / sum(alpha_sim)
  # names(r_sim) <- paste0("p_sim", 1:terra::nlyr(r_sim))
  #
  # # create SpatRaster object containing substrate composition probabilities and covariates
  # r <- list(r_sim, alpha_sim, terra::subset(x, subset = layers)) |>
  #   terra::rast()
  #
  # # create data frame (tibble) containing substrate composition probabilities and covariates
  # df <- rastersample_prep(r, cells = TRUE)

  # method 2: create dgCMatrix object via data_stack_dirich() function
  # this makes it easier to automate the matrix algebra based on the structure of
  # form_sim
  df <- x |>
    terra::as.data.frame(cells = TRUE, na.rm = TRUE) |>
    tibble::as_tibble()
  ds_sim <- dirinla::data_stack_dirich(
    y = rep.int(NA, times = nrow(df) * d),
    covariates = dirinla::formula_list(form_sim),
    data = df,
    d = d,
    n = nrow(df)
  )
  x_sim <- ncol(ds_sim) |>
    stats::rnorm(mean = 0, sd = 2)
  eta_sim <- ds_sim %*% x_sim
  alpha_sim <- eta_sim |>
    exp() |>
    matrix(ncol = d, byrow = TRUE)
  p_sim <- alpha_sim / rowSums(alpha_sim)  # expected values of Dirichlet distribution with parameters equal to alpha

  if (inherits(x, "SpatRaster") & as_raster) {
    r_temp <- terra::subset(x, subset = 1)
    terra::values(r_temp) <- NA
    r_alpha_sim <- r_p_sim <- rep(r_temp, d)
    for (i in 1:d) {
      r_alpha_sim[[i]][df$cell] <- alpha_sim[, i]
      r_p_sim[[i]][df$cell] <- p_sim[, i]
    }
    names(r_alpha_sim) <- stringr::str_c("alpha_sim_", 1:d)
    x <- c(x, r_alpha_sim)
    names(r_p_sim) <- stringr::str_c("p_sim_", 1:d)
    x <- c(x, r_p_sim)
    tibble::lst(d = d, form_sim = form_sim, coef_sim = x_sim,
                data = terra::wrap(x))
  } else {
    df <- df |>
      dplyr::bind_cols(
        tibble::as_tibble(alpha_sim, .name_repair = ~ paste0("alpha_sim_", 1:d))
      ) |>
      dplyr::bind_cols(
        tibble::as_tibble(p_sim, .name_repair = ~ paste0("p_sim_", 1:d))
      )
    tibble::lst(d = d, form_sim = form_sim, coef_sim = x_sim, data = df)
  }
}
