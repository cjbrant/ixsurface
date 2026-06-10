# ixsurface 0.2.0

## New features

* `model_surface()` for comparing prediction surfaces from multiple
  fitted models on the same predictor pair. Three modes: `"surface"`
  (smooth surfaces, one per model), `"fitted"` (predictions at observed
  points), `"residual"` (vertical residual segments for one model).

* `model_dominance()` for visualizing regional model performance as a
  2D heatmap colored by the locally best-performing model. Uses
  cross-validated residuals when available (via `caret::train()` with
  `savePredictions = "final"`), with explicit warnings when falling
  back to training residuals.

* `extract_residuals()` utility for honest residual extraction.
  Dispatches on model class: returns CV residuals from caret `train`
  objects when available, training residuals with a warning otherwise.

## Notes

* New functions complement `interaction_surface()`: that function
  compares factor levels within a single model; the new functions
  compare multiple models on shared predictors.

# ixsurface 0.1.0

* Initial release.
* `interaction_surface()` for interactive 3D surface plots from fitted models.
* `plot_crossings()` for standalone crossing-point visualization.
* `find_crossings()` for programmatic crossing detection.
* `interaction_surface_grid()` for all-pairwise exploration.
* `sim_factorial()` for generating test data with known interactions.
* Support for `lm`, `glm`, `gam`, and any model with a `predict()` method.
* Automatic binning of continuous conditioning variables.
* Categorical axis support with tick label mapping.
