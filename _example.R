devtools::load_all()

# Example usage
library(terra)
f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)

r <- c(r, MultiscaleDTM::BPI(r, w = c(2, 4)))

r_prep <- prepare_data(r, poly_degree = 3)
terra::plot(terra::unwrap(r_prep))
sim_data <- simulate_compositional_data(
  x = r_prep,
  d = 3,
  n_cov_sim = 3,
  as_raster = TRUE,
  seed = 123
)
plot_compositional_data(sim_data)
df <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data)) |>
  rastersample::spatial_sample(n = 50, method = "random")
fit <- fit_compositional_data(
  y ~ 1 + elevation_poly1 + elevation_poly2,
  data = df
)
temp <- predict_compositional_data(fit, new_data = r_prep)
plot_compositional_data(sim_data)
terra::plot(terra::unwrap(temp), range = c(0, 1), nr = 1)

pred <- terra::unwrap(temp) |>
  terra::values(mat = TRUE, na.rm = TRUE)
layers <- stringr::str_replace_all(
  colnames(pred),
  pattern = "_hat",
  replacement = "_sim"
)
sim <- terra::unwrap(sim_data$data) |>
  terra::subset(layers) |>
  terra::values(mat = TRUE, na.rm = TRUE)
robCompositions::aDist(sim, y = pred)


# test simulation framework (using SimDesign) with compositional data ---------------
res <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data)) |>
  run_simulation_compositional_data(
    formula = y ~ 1 + elevation_poly1 + elevation_poly2,
    n = c(20, 50, 100, 500, 1000, 5000, 10000),
    method = "random",
    replications = 10,
    parallel = FALSE
  )
res <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data)) |>
  run_simulation_compositional_data(
    formula = y ~ 1 + elevation_poly1 + elevation_poly2,
    n = c(20, 50, 100, 500, 1000, 5000, 10000),
    method = "random",
    replications = 10,
    parallel = TRUE,
    n_cores = 4
  )
ggplot2::ggplot(res, mapping = ggplot2::aes(x = n, y = mean.adist)) +
  ggplot2::geom_point() +
  ggplot2::geom_line()


df <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data)) |>
  terra::as.data.frame(na.rm = TRUE) |>
  tibble::as_tibble()
fixed_objects <- list(
  df = df,
  formula = y ~ 1 + elevation_poly1 + elevation_poly2
)

design <- SimDesign::createDesign(
  n = c(20, 50, 100, 500, 1000, 5000),
  method = "random"
)

generate_data <- function(condition, fixed_objects) {
  n <- condition$n
  method <- condition$method
  fixed_objects$df |>
    rastersample::spatial_sample(n = n, method = method, drop_na = TRUE)
}

analyse_data <- function(condition, dat, fixed_objects) {
  sample_fit <- nccosseasim::fit_compositional_data(
    fixed_objects$formula,
    data = dat
  )
  pred <- nccosseasim::predict_compositional_data(
    sample_fit,
    new_data = fixed_objects$df
  )
  adist <- robCompositions::aDist(
    dplyr::select(fixed_objects$df, tidyselect::starts_with(".p_sim")),
    y = dplyr::select(pred, tidyselect::starts_with(".p_hat"))
  )
  ret <- c(adist = adist)
  ret
}

summarise_data <- function(condition, results, fixed_objects) {
  ret <- c(mean = SimDesign::bias(results, parameter = 0))
  ret
}

res <- SimDesign::runSimulation(
  design,
  replications = 5,
  generate = generate_data,
  analyse = analyse_data,
  summarise = summarise_data,
  fixed_objects = fixed_objects,
  save = FALSE
)

future::plan(future.mirai::mirai_multisession, workers = 4)
res <- SimDesign::runSimulation(
  design,
  replications = 5,
  generate = generate_data,
  analyse = analyse_data,
  summarise = summarise_data,
  fixed_objects = fixed_objects,
  save = FALSE,
  parallel = "future"
)
