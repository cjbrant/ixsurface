#!/usr/bin/env Rscript
# ixsurface demo — interaction surface visualization
# source("demo.R") in RStudio, or Rscript demo.R

library(ixsurface)

# =============================================================================
# 1. Mixed design: 2 continuous + 1 categorical (the bread and butter)
# =============================================================================
cat("=== Mixed Design ===\n")
dat_mixed = sim_factorial(n = 300, design = "mixed", seed = 42)
fit_mixed = lm(y ~ temp * pressure * catalyst, data = dat_mixed)

# surfaces for each catalyst level, with color-coded observed data
p1 = interaction_surface(
  fit_mixed,
  x = "temp", y = "pressure", facet_by = "catalyst",
  alpha = 0.5, show_points = TRUE, show_contour = TRUE,
  labs = list(x = "Temperature (C)", y = "Pressure (psi)", z = "Yield"),
  title = "Mixed: temp x pressure | catalyst"
)

# programmatic crossing detection
cx1 = find_crossings(fit_mixed, "temp", "pressure", "catalyst")
cat(sprintf("Crossings: %d points, %d surface pairs\n",
            nrow(cx1), length(unique(cx1$pair_label))))

# =============================================================================
# 2. Fully continuous — automatic binning of facet_by variable
# =============================================================================
cat("\n=== Continuous Design (binned facet_by) ===\n")
dat_cont = sim_factorial(n = 300, design = "continuous", seed = 42)
fit_cont = lm(y ~ temp * pressure * speed, data = dat_cont)

# speed is continuous -> gets binned into 4 quantile groups
p2 = interaction_surface(
  fit_cont,
  x = "temp", y = "pressure", facet_by = "speed",
  n = 30, n_bins = 4, bin_method = "quantile",
  show_crossings = TRUE, show_contour = TRUE,
  title = "Continuous: temp x pressure | speed (4 quantile bins)"
)

# =============================================================================
# 3. Fully categorical design
# =============================================================================
cat("\n=== Categorical Design ===\n")
dat_cat = sim_factorial(n = 300, design = "categorical", seed = 42)
fit_cat = lm(y ~ catalyst * operator * shift, data = dat_cat)

p3 = interaction_surface(
  fit_cat,
  x = "catalyst", y = "operator", facet_by = "shift",
  show_points = TRUE,
  title = "Categorical: catalyst x operator | shift"
)

# =============================================================================
# 4. GLM example — logistic regression
# =============================================================================
cat("\n=== GLM (Logistic) Demo ===\n")
dat_glm = sim_factorial(n = 500, design = "mixed", seed = 42)
dat_glm$success = rbinom(nrow(dat_glm), 1, plogis((dat_glm$y - 50) / 5))
gfit = glm(success ~ temp * pressure * catalyst,
           data = dat_glm, family = binomial)

p4 = interaction_surface(
  gfit,
  x = "temp", y = "pressure", facet_by = "catalyst",
  n = 30, alpha = 0.5,
  labs = list(z = "P(success)"),
  title = "GLM: P(success) ~ temp x pressure | catalyst"
)

# =============================================================================
# 5. Grid of all pairwise plots
# =============================================================================
cat("\n=== Grid (all pairwise) ===\n")
all_plots = interaction_surface_grid(fit_mixed, n = 20, show_points = TRUE)
cat(sprintf("Generated %d plots: %s\n",
            length(all_plots), paste(names(all_plots), collapse = ", ")))

cat("\nDone. Objects: p1, p2, p3, p4, all_plots\n")
cat("Type any plot name in RStudio to view interactively.\n")
