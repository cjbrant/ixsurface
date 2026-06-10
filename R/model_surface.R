#' Visualize and compare prediction surfaces from multiple models
#'
#' Renders interactive 3D plotly surfaces showing how multiple fitted
#' models predict across a 2D predictor grid. Where surfaces diverge,
#' models disagree about the response — indicating regions where model
#' class affects predictions. Companion to [interaction_surface()], which
#' compares factor levels within a single model; this function compares
#' models on the same predictors.
#'
#' Three rendering modes are available:
#'
#' * `"surface"` shows smooth prediction surfaces, one per model.
#'   Answers: what does each model think the function looks like?
#' * `"fitted"` shows each model's predictions at observed `(x, y)`
#'   locations as colored points. Answers: where does each model land
#'   relative to the data?
#' * `"residual"` shows vertical segments from each observation to one
#'   model's fitted value. Answers: how big are this model's errors,
#'   and where? Only one model at a time in this mode.
#'
#' @param models Named list of fitted model objects. Each must have a
#'   `predict` method recognized by [safe_predict()].
#' @param x Character. Variable name mapped to the x-axis. Must be a
#'   predictor in all models.
#' @param y Character. Variable name mapped to the y-axis. Must be a
#'   predictor in all models.
#' @param mode Character. Visualization mode: `"surface"` (default),
#'   `"fitted"`, or `"residual"`.
#' @param n Integer. Grid resolution per continuous axis (default 50).
#' @param alpha Numeric in `[0, 1]`. Surface opacity (default 0.5).
#' @param show_observed Logical. Overlay observed data points (default `TRUE`).
#' @param labs Named list for axis labels, e.g.
#'   `list(x = "Weight", y = "Year", z = "MPG")`.
#' @param title Character or `NULL`. Plot title.
#' @param theme Character. Color theme: `"default"`, `"viridis"`, or
#'   `"grey"`.
#' @param ... Reserved for future use.
#'
#' @return A `plotly` htmlwidget object.
#'
#' @details
#' Non-focal predictors are held at their median (continuous) or mode
#' (categorical), using [make_prediction_grid()]. The same convention
#' as [interaction_surface()].
#'
#' For categorical predictors, integer-coded axes are used internally
#' with level names as tick labels, matching [interaction_surface()]
#' behavior.
#'
#' All models in `models` must share the same response variable (extracted
#' from the first model via `model.frame()`). Mixing models with different
#' responses produces nonsense comparisons and is not checked.
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(n = 300, design = "continuous", seed = 42)
#' fit_lm = lm(y ~ temp * pressure, data = dat)
#' fit_gam = mgcv::gam(y ~ s(temp) + s(pressure), data = dat)
#'
#' model_surface(
#'   models = list(linear = fit_lm, gam = fit_gam),
#'   x = "temp", y = "pressure"
#' )
#' }
#'
#' @seealso [interaction_surface()] for within-model factor comparison,
#'   [model_dominance()] for regional performance comparison.
#'
#' @export
model_surface = function(models, x, y,
                         mode = c("surface", "fitted", "residual"),
                         n = 50,
                         alpha = 0.5,
                         show_observed = TRUE,
                         labs = NULL,
                         title = NULL,
                         theme = "default",
                         ...) {

  mode = match.arg(mode)

  if (!is.list(models) || is.null(names(models)) || length(models) == 0) {
    stop("`models` must be a non-empty named list of fitted model objects.")
  }

  # use first model as reference for data, response, predictor types
  ref_model = models[[1]]
  mf = model.frame(ref_model)
  predictor_names = names(mf)[-1]
  response_name = names(mf)[1]

  for (v in c(x, y)) {
    if (!v %in% predictor_names) {
      stop(sprintf(
        "Variable '%s' not found among model predictors: %s",
        v, paste(predictor_names, collapse = ", ")
      ))
    }
  }

  # axis labels (mirroring interaction_surface)
  x_lab = if (!is.null(labs$x)) labs$x else x
  y_lab = if (!is.null(labs$y)) labs$y else y
  z_lab = if (!is.null(labs$z)) {
    labs$z
  } else if (mode == "residual") {
    paste0(response_name, " / fitted")
  } else {
    paste0("predicted ", response_name)
  }

  # build prediction grid using existing package utility
  grid_result = make_prediction_grid(ref_model, x, y, facet_by = NULL,
                                     n = n)
  grid = grid_result$grid

  # detect predictor types for axis handling
  x_type = detect_factor_type(ref_model, x)
  y_type = detect_factor_type(ref_model, y)

  x_tickvals = NULL; x_ticktext = NULL
  y_tickvals = NULL; y_ticktext = NULL

  if (x_type %in% c("categorical", "quasi_categorical")) {
    x_levels = if (is.factor(mf[[x]])) levels(mf[[x]]) else sort(unique(mf[[x]]))
    grid$.x_num = match(as.character(grid[[x]]), as.character(x_levels))
    x_tickvals = seq_along(x_levels)
    x_ticktext = as.character(x_levels)
  } else {
    grid$.x_num = grid[[x]]
  }

  if (y_type %in% c("categorical", "quasi_categorical")) {
    y_levels = if (is.factor(mf[[y]])) levels(mf[[y]]) else sort(unique(mf[[y]]))
    grid$.y_num = match(as.character(grid[[y]]), as.character(y_levels))
    y_tickvals = seq_along(y_levels)
    y_ticktext = as.character(y_levels)
  } else {
    grid$.y_num = grid[[y]]
  }

  x_vals = sort(unique(grid$.x_num))
  y_vals = sort(unique(grid$.y_num))

  # one color per model
  model_colors = surface_palette(length(models))

  p = plotly::plot_ly()

  # dispatch by mode
  if (mode == "surface") {
    p = .add_model_surfaces(p, models, grid, x, y, x_vals, y_vals,
                            model_colors, alpha, x_lab, y_lab, z_lab)

  } else if (mode == "fitted") {
    p = .add_model_fitted(p, models, mf, x, y, x_type, y_type,
                          x_levels = if (exists("x_levels")) x_levels else NULL,
                          y_levels = if (exists("y_levels")) y_levels else NULL,
                          response_name, model_colors,
                          x_lab, y_lab, z_lab)

  } else if (mode == "residual") {
    if (length(models) > 1) {
      message("Residual mode uses only the first model: '",
              names(models)[1], "'")
    }
    p = .add_model_residuals(p, models[1], mf, x, y, x_type, y_type,
                             x_levels = if (exists("x_levels")) x_levels else NULL,
                             y_levels = if (exists("y_levels")) y_levels else NULL,
                             response_name, model_colors,
                             x_lab, y_lab, z_lab)
  }

  # add observed data overlay (always except in residual mode where
  # it's handled inside .add_model_residuals)
  if (show_observed && mode != "residual") {
    p = .add_observed_overlay(p, mf, x, y, x_type, y_type,
                              x_levels = if (exists("x_levels")) x_levels else NULL,
                              y_levels = if (exists("y_levels")) y_levels else NULL,
                              response_name, x_lab, y_lab)
  }

  # title
  auto_title = if (is.null(title)) {
    titles_by_mode = c(
      surface  = paste0("model surfaces: ", z_lab),
      fitted   = paste0("fitted values: ", z_lab),
      residual = paste0("residuals: ", names(models)[1])
    )
    titles_by_mode[mode]
  } else {
    title
  }

  # scene with optional categorical tick handling
  scene = list(
    xaxis = list(title = x_lab),
    yaxis = list(title = y_lab),
    zaxis = list(title = z_lab),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))
  )
  if (!is.null(x_tickvals)) {
    scene$xaxis$tickvals = x_tickvals
    scene$xaxis$ticktext = x_ticktext
  }
  if (!is.null(y_tickvals)) {
    scene$yaxis$tickvals = y_tickvals
    scene$yaxis$ticktext = y_ticktext
  }

  p = p %>% plotly::layout(
    title = list(text = auto_title, font = list(size = 14)),
    scene = scene,
    legend = list(orientation = "h", yanchor = "bottom", y = -0.15)
  )

  attr(p, "ixsurface_meta") = list(
    type = "model_surface",
    mode = mode,
    x = x, y = y,
    n_models = length(models),
    model_names = names(models),
    x_vals = x_vals, y_vals = y_vals
  )

  p
}


# helper: render smooth surface per model
# @keywords internal
.add_model_surfaces = function(p, models, grid, x, y, x_vals, y_vals,
                               colors, alpha, x_lab, y_lab, z_lab) {

  for (i in seq_along(models)) {
    grid$.pred = safe_predict(models[[i]], grid)

    # pivot to matrix
    z_mat = matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
    for (r in seq_len(nrow(grid))) {
      xi = match(grid$.x_num[r], x_vals)
      yi = match(grid$.y_num[r], y_vals)
      if (!is.na(xi) && !is.na(yi)) z_mat[yi, xi] = grid$.pred[r]
    }

    base_col = colors[i]
    surf_colorscale = list(
      list(0, lighten_color(base_col, 0.4)),
      list(1, base_col)
    )

    p = p %>% plotly::add_surface(
      x = x_vals, y = y_vals, z = z_mat,
      opacity = alpha,
      colorscale = surf_colorscale,
      showscale = FALSE,
      name = names(models)[i],
      legendgroup = names(models)[i],
      showlegend = TRUE,
      hovertemplate = paste0(
        x_lab, ": %{x}<br>",
        y_lab, ": %{y}<br>",
        z_lab, ": %{z:.3f}<br>",
        names(models)[i],
        "<extra></extra>"
      )
    )
  }
  p
}


# helper: render fitted-value markers per model
# @keywords internal
.add_model_fitted = function(p, models, mf, x, y, x_type, y_type,
                             x_levels, y_levels, response_name,
                             colors, x_lab, y_lab, z_lab) {

  obs = mf

  if (x_type %in% c("categorical", "quasi_categorical")) {
    obs$.xn = match(as.character(obs[[x]]), as.character(x_levels))
  } else {
    obs$.xn = obs[[x]]
  }
  if (y_type %in% c("categorical", "quasi_categorical")) {
    obs$.yn = match(as.character(obs[[y]]), as.character(y_levels))
  } else {
    obs$.yn = obs[[y]]
  }

  for (i in seq_along(models)) {
    fitted = safe_predict(models[[i]], obs)
    obs$.fitted = fitted

    p = p %>% plotly::add_markers(
      data = obs,
      x = ~.xn, y = ~.yn, z = ~.fitted,
      marker = list(
        size = 4,
        color = colors[i],
        opacity = 0.75,
        line = list(width = 0.5, color = "black")
      ),
      name = names(models)[i],
      legendgroup = names(models)[i],
      hovertemplate = paste0(
        x_lab, ": %{x}<br>",
        y_lab, ": %{y}<br>",
        z_lab, ": %{z:.3f}<br>",
        names(models)[i],
        "<extra></extra>"
      )
    )
  }
  p
}


# helper: render residual segments for one model
# @keywords internal
.add_model_residuals = function(p, model_list, mf, x, y, x_type, y_type,
                                x_levels, y_levels, response_name,
                                colors, x_lab, y_lab, z_lab) {

  model = model_list[[1]]
  model_name = names(model_list)[1]
  obs = mf
  obs$.z_obs = obs[[response_name]]
  obs$.z_fit = safe_predict(model, obs)

  if (x_type %in% c("categorical", "quasi_categorical")) {
    obs$.xn = match(as.character(obs[[x]]), as.character(x_levels))
  } else {
    obs$.xn = obs[[x]]
  }
  if (y_type %in% c("categorical", "quasi_categorical")) {
    obs$.yn = match(as.character(obs[[y]]), as.character(y_levels))
  } else {
    obs$.yn = obs[[y]]
  }

  n_obs = nrow(obs)
  seg_x = as.vector(rbind(obs$.xn, obs$.xn, rep(NA, n_obs)))
  seg_y = as.vector(rbind(obs$.yn, obs$.yn, rep(NA, n_obs)))
  seg_z = as.vector(rbind(obs$.z_obs, obs$.z_fit, rep(NA, n_obs)))

  p = p %>% plotly::add_trace(
    x = seg_x, y = seg_y, z = seg_z,
    type = "scatter3d", mode = "lines",
    line = list(color = colors[1], width = 2),
    opacity = 0.5,
    name = paste0(model_name, " residuals"),
    showlegend = TRUE
  )

  # observed and fitted point overlays for context
  p = p %>% plotly::add_markers(
    data = obs,
    x = ~.xn, y = ~.yn, z = ~.z_obs,
    marker = list(size = 3, color = "black", opacity = 0.5),
    name = "observed",
    hovertemplate = paste0(
      "observed<br>",
      x_lab, ": %{x}<br>",
      y_lab, ": %{y}<br>",
      response_name, ": %{z:.3f}",
      "<extra></extra>"
    )
  )

  p = p %>% plotly::add_markers(
    data = obs,
    x = ~.xn, y = ~.yn, z = ~.z_fit,
    marker = list(size = 3, color = colors[1], opacity = 0.7),
    name = paste0(model_name, " fitted"),
    hovertemplate = paste0(
      "fitted (", model_name, ")<br>",
      x_lab, ": %{x}<br>",
      y_lab, ": %{y}<br>",
      response_name, ": %{z:.3f}",
      "<extra></extra>"
    )
  )

  p
}


# helper: observed data overlay (used by surface and fitted modes)
# @keywords internal
.add_observed_overlay = function(p, mf, x, y, x_type, y_type,
                                 x_levels, y_levels, response_name,
                                 x_lab, y_lab) {

  obs = mf
  obs$.z = obs[[response_name]]

  if (x_type %in% c("categorical", "quasi_categorical")) {
    obs$.xn = match(as.character(obs[[x]]), as.character(x_levels))
  } else {
    obs$.xn = obs[[x]]
  }
  if (y_type %in% c("categorical", "quasi_categorical")) {
    obs$.yn = match(as.character(obs[[y]]), as.character(y_levels))
  } else {
    obs$.yn = obs[[y]]
  }

  p %>% plotly::add_markers(
    data = obs,
    x = ~.xn, y = ~.yn, z = ~.z,
    marker = list(size = 3, color = "black", opacity = 0.55),
    name = "observed",
    hovertemplate = paste0(
      "observed<br>",
      x_lab, ": %{x}<br>",
      y_lab, ": %{y}<br>",
      response_name, ": %{z:.3f}",
      "<extra></extra>"
    )
  )
}
