#' Simulate compositional data from a Dirichlet likelihood
#'
#' @description
#' Simulates compositional probabilities (`p_sim*`) for `d` components from
#' selected covariate raster layers. For each component a coefficient vector
#' (intercept + covariates) is drawn i.i.d. `N(0, 2^2)`, linear predictors are
#' exponentiated (`.alpha_sim*`), then normalized to the simplex (`.p_sim*`).
#'
#' @details
#' Steps:
#' 1. Unwraps a `PackedSpatRaster` if supplied.
#' 2. Selects covariate layers (all, or random subset when `n_cov_sim` given).
#' 3. Extracts covariate values (`terra::values(mat = TRUE)`).
#' 4. Builds feature matrix: intercept + covariates.
#' 5. Draws a coefficient matrix (`n_features x d`).
#' 6. Computes linear predictors for rows without missing covariates.
#' 7. Applies `exp` to get `.alpha_sim*` and row-normalizes to obtain `.p_sim*`.
#' 8. Returns both `.alpha_sim*` and `.p_sim*` as a wrapped raster or tibble.
#'
#' Rows with any `NA` in covariates yield `NA` in all component probabilities.
#'
#' Memory: Allocates feature matrix (`ncell * n_features`) and probability
#' matrix (`ncell * d`).
#'
#' @param x A `SpatRaster` or `PackedSpatRaster` containing covariate layers.
#' @param d Integer (>= 2). Number of compositional components.
#' @param layers Optional character or numeric vector of layer names/indices to
#'   use. If `NULL`, all layers are used or random subset is chosen via
#'   `n_cov_sim`.
#' @param n_cov_sim Optional integer: number of layers randomly sampled when
#'   `layers` is `NULL`.
#' @param as_raster Logical; if `TRUE` returns a wrapped `SpatRaster` with
#'   probability layers. Otherwise a tibble.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' * `d` Number of components.
#' * `layers` Selected layers.
#' * `n_cov_sim` Number sampled (if applicable).
#' * `as_raster` Input flag.
#' * `seed` Seed used (or `NULL`).
#' * `form_sim` Formula describing the component-specific linear predictors.
#' * `coef_sim` Numeric matrix (`n_features x d`).
#' * `data` Raster (wrapped) or tibble with columns/layers `.alpha_sim*` and
#'     `p_sim*`.
#'
#' @seealso [prepare_data()], [terra::values()], [stats::rnorm()]
#'
#' @examples
#' \donttest{
#' library(terra)
#' r <- rast(system.file("ex/elev.tif", package = "terra"))
#' sim <- simulate_comp_data(x = r, d = 3, as_raster = TRUE, seed = 123)
#' sim$coef_sim
#' names(sim$data)
#' }
#'
#' @importFrom terra unwrap subset values rast wrap nlyr
#' @importFrom stats rnorm as.formula
#' @importFrom tibble tibble lst as_tibble
#' @importFrom dplyr bind_cols
#' @export
simulate_compositional_data <- function(
  x,
  d,
  layers = NULL,
  n_cov_sim = NULL,
  as_raster = FALSE,
  seed = NULL
) {
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  if (inherits(x, "PackedSpatRaster")) {
    x <- terra::unwrap(x)
  }
  if (!is.integer(d)) {
    d <- as.integer(d)
  }
  if (d < 2L) {
    stop("d must be >= 2")
  }

  n_all <- terra::nlyr(x)
  if (is.null(layers)) {
    if (!is.null(n_cov_sim)) {
      if (n_cov_sim < 1 || n_cov_sim > n_all) {
        stop("n_cov_sim out of range")
      }
      layers <- sort(sample.int(n_all, size = n_cov_sim))
    } else {
      layers <- seq_len(n_all)
    }
  }
  layer_names <- if (is.numeric(layers)) names(x)[layers] else layers

  vals <- terra::values(terra::subset(x, layers), mat = TRUE)
  n_cells <- nrow(vals)

  features <- cbind(Intercept = rep(1, n_cells), vals)
  valid <- rowSums(is.na(features[, -1L, drop = FALSE])) == 0
  n_feat <- ncol(features)

  betas <- matrix(
    stats::rnorm(
      n_feat * d,
      mean = stats::rnorm(1, mean = 0, sd = 10),
      sd = 5
    ),
    nrow = n_feat,
    ncol = d
  )
  rownames(betas) <- c("Intercept", layer_names)
  colnames(betas) <- paste0("comp", seq_len(d))

  eta_sim <- matrix(NA_real_, nrow = n_cells, ncol = d)
  eta_sim[valid, ] <- features[valid, , drop = FALSE] %*% betas

  alpha_sim <- exp(eta_sim)
  row_sums <- rowSums(alpha_sim, na.rm = TRUE)
  p_sim <- alpha_sim / row_sums

  base_terms <- paste(c("1", layer_names), collapse = " + ")
  form_sim <- stats::as.formula(
    paste("y ~", paste(rep(base_terms, d), collapse = " | "))
  )

  if (as_raster) {
    out <- terra::rast(terra::subset(x, layers), nlyr = 2 * d)
    terra::values(out) <- cbind(alpha_sim, p_sim)
    names(out) <- c(
      paste0(".alpha_sim", seq_len(d)),
      paste0(".p_sim", seq_len(d))
    )
    data_out <- terra::wrap(out)
  } else {
    df <- tibble::tibble(cell = seq_len(n_cells)) |>
      dplyr::bind_cols(
        tibble::as_tibble(
          alpha_sim,
          .name_repair = ~ paste0(".alpha_sim", seq_len(d))
        ),
        tibble::as_tibble(p_sim, .name_repair = ~ paste0(".p_sim", seq_len(d)))
      )
    data_out <- df
  }

  tibble::lst(
    d = d,
    layers = layers,
    n_cov_sim = n_cov_sim,
    as_raster = as_raster,
    seed = seed,
    form_sim = form_sim,
    coef_sim = betas,
    data = data_out
  )
}
