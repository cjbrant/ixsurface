#' Generate All Pairwise Interaction Surface Plots
#'
#' For a model with multiple factors, generates interaction surface plots for
#' every pair of focal variables, using remaining variables as conditioning
#' factors (\code{facet_by}). Useful for exploratory analysis.
#'
#' @param model A fitted model object.
#' @param factors Character vector or NULL. Variables to consider. If NULL,
#'   uses all predictors.
#' @param facet_max Integer. Maximum number of \code{facet_by} variables per
#'   plot (default 2). Extra variables are held at central values.
#' @param n Integer. Grid resolution (default 30, lower for speed).
#' @param n_bins Integer. Bins for continuous facet_by variables.
#' @param alpha Numeric. Surface opacity.
#' @param ... Passed to \code{interaction_surface}.
#'
#' @return A named list of plotly objects. Names use pattern \code{"x__y"}.
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(design = "mixed", seed = 42)
#' fit = lm(y ~ temp * pressure * catalyst, data = dat)
#' plots = interaction_surface_grid(fit)
#' plots$temp__pressure
#' }
#'
#' @export
interaction_surface_grid = function(model, factors = NULL, facet_max = 2,
                                    n = 30, n_bins = 3,
                                    alpha = 0.6, ...) {

  mf = model.frame(model)
  all_predictors = names(mf)[-1]

  if (is.null(factors)) {
    factors = all_predictors
  } else {
    bad = setdiff(factors, all_predictors)
    if (length(bad) > 0) {
      stop(sprintf("Variables not found in model: %s", paste(bad, collapse = ", ")))
    }
  }

  if (length(factors) < 2) {
    stop("Need at least 2 factors to generate interaction surface plots.")
  }

  pairs = utils::combn(factors, 2, simplify = FALSE)
  plots = list()

  for (pair in pairs) {
    x_var = pair[1]
    y_var = pair[2]
    remaining = setdiff(factors, pair)

    facet_by = if (length(remaining) > 0) {
      head(remaining, facet_max)
    } else {
      NULL
    }

    plot_name = paste0(x_var, "__", y_var)
    plots[[plot_name]] = interaction_surface(
      model, x = x_var, y = y_var, facet_by = facet_by,
      n = n, n_bins = n_bins, alpha = alpha, ...
    )
  }

  plots
}
