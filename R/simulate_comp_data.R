#' Simulate compositional data from a Dirichlet likelihood
#'
#' This function simulates compositional data based on a Dirichlet likelihood. It allows for 
#' the generation of covariate rasters, polynomial transformations, and compositional probabilities 
#' for a specified number of components. The output can be returned as a raster object or a data frame.
#'
#' @param x A `SpatRaster` or `PackedSpatRaster` object containing covariate data. If it is a 
#'   `PackedSpatRaster`, it will be unpacked using `terra::unwrap()`.
#' @param d Integer. The number of compositional components to simulate.
#' @param layers A character vector or numeric vector specifying the layers of `x` to use as covariates. 
#'   If `NULL`, all layers are used or a subset is selected based on `n_cov_sim`.
#' @param poly_degree Integer. The degree of polynomial transformations to apply to the covariates. 
#'   If greater than 1, additional layers are generated as powers of the original covariates.
#' @param n_cov_sim Integer. The number of covariate layers to randomly select for simulation. 
#'   Ignored if `layers` is specified.
#' @param as_raster Logical. If `TRUE`, the output includes simulated data as a `SpatRaster` object. 
#'   If `FALSE`, the output is a data frame.
#' @param seed Integer. A random seed for reproducibility. If `NULL`, no seed is set.
#'
#' @return A list containing the following elements:
#'   - `d`: The number of compositional components.
#'   - `layers`: The layers used for simulation.
#'   - `poly_degree`: The degree of polynomial transformations applied.
#'   - `n_cov_sim`: The number of covariate layers used.
#'   - `as_raster`: Whether the output includes raster data.
#'   - `seed`: The random seed used for simulation.
#'   - `form_sim`: The formula used to generate the simulated data.
#'   - `coef_sim`: The coefficients used in the simulation.
#'   - `data`: A `SpatRaster` object (if `as_raster = TRUE`) or a data frame containing the simulated 
#'     compositional probabilities (`p_sim_*`) and Dirichlet parameters (`alpha_sim_*`).
#'
#' @details The function performs the following steps:
#'   1. Unpacks `PackedSpatRaster` objects if necessary.
#'   2. Generates polynomial transformations of covariates if `poly_degree > 1`.
#'   3. Constructs a formula for simulating compositional data.
#'   4. Simulates Dirichlet parameters (`alpha_sim_*`) and compositional probabilities (`p_sim_*`).
#'   5. Returns the simulated data as a raster or data frame, depending on the `as_raster` parameter.
#'
#' @examples
#' # Example usage
#' library(terra)
#' covariates <- rast(matrix(runif(100), 10, 10))
#' names(covariates) <- "cov1"
#' sim_data <- simulate_comp_data(x = covariates, d = 3, poly_degree = 2, as_raster = TRUE)
#'
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
    tibble::lst(d = d, layers = layers, poly_degree = poly_degree,
                n_cov_sim = n_cov_sim, as_raster = as_raster, seed = seed,
                form_sim = form_sim, coef_sim = x_sim,
                data = terra::wrap(x))
  } else {
    df <- df |>
      dplyr::bind_cols(
        tibble::as_tibble(alpha_sim, .name_repair = ~ paste0("alpha_sim_", 1:d))
      ) |>
      dplyr::bind_cols(
        tibble::as_tibble(p_sim, .name_repair = ~ paste0("p_sim_", 1:d))
      )
    tibble::lst(d = d, layers = layers, poly_degree = poly_degree,
                n_cov_sim = n_cov_sim, as_raster = as_raster, seed = seed,
                form_sim = form_sim, coef_sim = x_sim, data = df)
  }
}
