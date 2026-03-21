test_that("sim_factorial generates correct designs", {
  d1 = sim_factorial(n = 50, design = "mixed", seed = 1)
  expect_equal(ncol(d1), 4)
  expect_true("catalyst" %in% names(d1))
  expect_true(is.factor(d1$catalyst))
  expect_true(is.numeric(d1$temp))
  expect_equal(nrow(d1), 50)

  d2 = sim_factorial(n = 50, design = "continuous", seed = 1)
  expect_true(all(c("temp", "pressure", "speed") %in% names(d2)))

  d3 = sim_factorial(n = 50, design = "categorical", seed = 1)
  expect_true(all(c("catalyst", "operator", "shift") %in% names(d3)))
})

test_that("bin_continuous produces correct number of bins", {
  x = rnorm(100)
  b3 = bin_continuous(x, n_bins = 3, method = "quantile")
  expect_true(is.factor(b3))
  expect_lte(nlevels(b3), 3)

  b5 = bin_continuous(x, n_bins = 5, method = "equal")
  expect_true(is.factor(b5))
  expect_lte(nlevels(b5), 5)
})

test_that("detect_factor_type correctly classifies", {
  dat = sim_factorial(n = 50, design = "mixed", seed = 1)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  expect_equal(detect_factor_type(fit, "temp"), "continuous")
  expect_equal(detect_factor_type(fit, "pressure"), "continuous")
  expect_equal(detect_factor_type(fit, "catalyst"), "categorical")
})

test_that("make_prediction_grid builds correct dimensions with categorical facet_by", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  result = make_prediction_grid(fit, x = "temp", y = "pressure",
                                facet_by = "catalyst", n = 10)

  # 10 temp x 10 pressure x 3 catalyst = 300
  expect_equal(nrow(result$grid), 10 * 10 * 3)
  expect_equal(length(result$binned_by), 0)  # no binning needed
})

test_that("make_prediction_grid bins continuous facet_by", {
  dat = sim_factorial(n = 200, design = "continuous", seed = 42)
  fit = lm(y ~ temp * pressure * speed, data = dat)

  result = make_prediction_grid(fit, x = "temp", y = "pressure",
                                facet_by = "speed", n = 10, n_bins = 3)

  # speed is binned to 3 levels: 10 * 10 * 3 = 300
  expect_equal(nrow(result$grid), 10 * 10 * 3)
  expect_true("speed" %in% names(result$binned_by))
  expect_equal(length(result$binned_by$speed$midpoints), 3)
})

test_that("interaction_surface returns plotly with correct metadata", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure",
                          facet_by = "catalyst", n = 10)

  expect_s3_class(p, "plotly")

  meta = attr(p, "ixsurface_meta")
  expect_equal(meta$x, "temp")
  expect_equal(meta$y, "pressure")
  expect_equal(meta$n_surfaces, 3)
})

test_that("interaction_surface works with show_points and color coding", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure",
                          facet_by = "catalyst", n = 10,
                          show_points = TRUE)

  expect_s3_class(p, "plotly")
})

test_that("interaction_surface works with show_contour", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure",
                          facet_by = "catalyst", n = 10,
                          show_crossings = TRUE, show_contour = TRUE)

  expect_s3_class(p, "plotly")
})

test_that("interaction_surface works with GLM", {
  dat = sim_factorial(n = 200, design = "mixed", seed = 42)
  dat$success = rbinom(nrow(dat), 1, plogis(scale(dat$y)))
  gfit = glm(success ~ temp * pressure * catalyst,
             data = dat, family = binomial)

  p = interaction_surface(gfit, x = "temp", y = "pressure",
                          facet_by = "catalyst", n = 10)

  expect_s3_class(p, "plotly")

  # predictions should be on probability scale
  meta = attr(p, "ixsurface_meta")
  z_vals = unlist(lapply(meta$z_matrices, c))
  expect_true(all(z_vals >= 0 & z_vals <= 1, na.rm = TRUE))
})

test_that("interaction_surface works with continuous facet_by (binned)", {
  dat = sim_factorial(n = 200, design = "continuous", seed = 42)
  fit = lm(y ~ temp * pressure * speed, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure",
                          facet_by = "speed", n = 10, n_bins = 4)

  expect_s3_class(p, "plotly")
  meta = attr(p, "ixsurface_meta")
  expect_true("speed" %in% names(meta$binned_by))
})

test_that("interaction_surface accepts labs argument", {
  dat = sim_factorial(n = 50, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure",
                          facet_by = "catalyst", n = 10,
                          labs = list(x = "Temperature (C)",
                                     y = "Pressure (psi)",
                                     z = "Yield (%)"))

  expect_s3_class(p, "plotly")
})

test_that("find_crossings returns a data.frame", {
  dat = sim_factorial(n = 200, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  cx = find_crossings(fit, "temp", "pressure", "catalyst", n = 20)
  expect_true(is.data.frame(cx))
  expect_true(all(c("cx", "cy", "cz", "pair_label") %in% names(cx)))
})

test_that("interaction_surface_grid returns named list", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  plots = interaction_surface_grid(fit, n = 10)
  expect_true(is.list(plots))
  expect_equal(length(plots), 3)  # C(3,2) = 3 pairs
  expect_true("temp__pressure" %in% names(plots))
})

test_that("interaction_surface with no facet_by works", {
  dat = sim_factorial(n = 100, design = "mixed", seed = 42)
  fit = lm(y ~ temp * pressure * catalyst, data = dat)

  p = interaction_surface(fit, x = "temp", y = "pressure", n = 10)
  expect_s3_class(p, "plotly")
  meta = attr(p, "ixsurface_meta")
  expect_equal(meta$n_surfaces, 1)
})

test_that("categorical x categorical surface works", {
  dat = sim_factorial(n = 200, design = "categorical", seed = 42)
  fit = lm(y ~ catalyst * operator * shift, data = dat)

  p = interaction_surface(fit, x = "catalyst", y = "operator",
                          facet_by = "shift", n = 10)
  expect_s3_class(p, "plotly")
})
