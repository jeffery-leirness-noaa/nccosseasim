devtools::load_all()

# Example usage
library(terra)
f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)

r_prep <- prepare_data(r, poly_degree = 5)
terra::plot(terra::unwrap(r_prep))
sim_data <- simulate_compositional_data(
  x = r_prep,
  d = 3,
  as_raster = TRUE,
  seed = 123
)
plot_compositional_data(sim_data)
df <- c(terra::unwrap(r_prep), terra::unwrap(sim_data$data)) |>
  rastersample::spatial_sample(n = 100, method = "random")
fit <- fit_compositional_data(
  df,
  formula = y ~ elevation_poly1 + elevation_poly2
)
ck <- predict_compositional_data(fit, new_data = r_prep)
plot_compositional_data(sim_data)
terra::plot(terra::unwrap(ck), range = c(0, 1), nr = 1)
