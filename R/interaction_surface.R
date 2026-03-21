#' Interactive 3D Surface Plot for Multi-Factor Interactions
#'
#' Generates an interactive plotly surface plot from a fitted model. Two focal
#' variables are mapped to the x and y aesthetics, with the predicted response
#' on z. Additional conditioning factors (\code{facet_by}) generate separate
#' surfaces — where surfaces cross indicates interaction effects.
#'
#' @param model A fitted model object with a \code{predict} method (e.g., from
#'   \code{lm}, \code{aov}, \code{glm}, \code{gam}).
#' @param x Character. Variable name mapped to the x-axis.
#' @param y Character. Variable name mapped to the y-axis.
#' @param facet_by Character vector or NULL. Variable(s) whose levels generate
#'   separate surfaces. Continuous variables are automatically binned.
#' @param n Integer. Grid resolution per continuous axis (default 50).
#'   Analogous to \code{n} in \code{geom_contour}.
#' @param n_bins Integer. Number of bins for continuous \code{facet_by}
#'   variables (default 3).
#' @param bin_method Character. Binning method for continuous \code{facet_by}:
#'   "quantile" (default), "equal", or "pretty".
#' @param alpha Numeric in \code{[0, 1]}. Surface opacity (default 0.6).
#' @param show_points Logical. If TRUE, overlays the observed data points,
#'   color-coded by \code{facet_by} level. Analogous to adding
#'   \code{geom_point}.
#' @param show_crossings Logical. If TRUE (default), marks regions where
#'   surfaces cross with red markers.
#' @param show_contour Logical. If FALSE (default), no contour projection.
#'   If TRUE, projects crossing curves onto the x-y floor of the plot as a
#'   2D summary of interaction regions.
#' @param contour_z Numeric or NULL. The z-value at which to draw the contour
#'   projection. If NULL, uses the minimum z in the plot.
#' @param labs Named list for axis labels, e.g.,
#'   \code{list(x = "Temperature", y = "Pressure", z = "Yield")}.
#'   Analogous to \code{ggplot2::labs()}.
#' @param title Character or NULL. Plot title.
#' @param theme Character. Color theme: "default", "viridis", or "grey".
#' @param ... Additional arguments (reserved for future use).
#'
#' @return A \code{plotly} htmlwidget object.
#'
#' @details
#' \strong{Geometric interpretation:}
#' \itemize{
#'   \item Parallel surfaces \eqn{\to} no interaction between \code{facet_by}
#'     and the focal variables
#'   \item Crossing surfaces \eqn{\to} interaction present
#'   \item Twisted/warped surfaces \eqn{\to} higher-order or nonlinear interaction
#' }
#'
#' For categorical focal variables, the surface is constructed over an integer
#' grid with axis tick labels showing level names.
#'
#' Non-focal, non-\code{facet_by} variables are held at their median
#' (continuous) or mode (categorical).
#'
#' For GLM family models, predictions are returned on the response scale via
#' \code{predict(..., type = "response")}.
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(design = "mixed", seed = 42)
#' fit = lm(y ~ temp * pressure * catalyst, data = dat)
#' interaction_surface(fit, x = "temp", y = "pressure", facet_by = "catalyst")
#'
#' # with GLM
#' dat$success = rbinom(nrow(dat), 1, plogis(scale(dat$y)))
#' gfit = glm(success ~ temp * pressure * catalyst, data = dat, family = binomial)
#' interaction_surface(gfit, x = "temp", y = "pressure", facet_by = "catalyst")
#' }
#'
#' @export
interaction_surface = function(model, x, y, facet_by = NULL,
                               n = 50, n_bins = 3,
                               bin_method = "quantile",
                               alpha = 0.6,
                               show_points = FALSE,
                               show_crossings = TRUE,
                               show_contour = FALSE,
                               contour_z = NULL,
                               labs = NULL,
                               title = NULL,
                               theme = "default", ...) {

  # --- validate ---
  mf = model.frame(model)
  predictor_names = names(mf)[-1]
  response_name = names(mf)[1]

  for (v in c(x, y, facet_by)) {
    if (!v %in% predictor_names) {
      stop(sprintf("Variable '%s' not found among model predictors: %s",
                   v, paste(predictor_names, collapse = ", ")))
    }
  }

  # --- axis labels ---
  x_lab = if (!is.null(labs$x)) labs$x else x
  y_lab = if (!is.null(labs$y)) labs$y else y
  z_lab = if (!is.null(labs$z)) labs$z else paste0("predicted ", response_name)

  # --- build prediction grid ---
  grid_result = make_prediction_grid(model, x, y,
                                     facet_by = facet_by, n = n,
                                     n_bins = n_bins, bin_method = bin_method)
  grid = grid_result$grid
  binned_by = grid_result$binned_by
  grid$.pred = safe_predict(model, grid)

  # --- axis type handling ---
  x_type = detect_factor_type(model, x)
  y_type = detect_factor_type(model, y)

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

  # --- split by facet_by ---
  if (!is.null(facet_by) && length(facet_by) > 0) {
    by_combos = unique(grid[, facet_by, drop = FALSE])
    n_surfaces = nrow(by_combos)
  } else {
    by_combos = data.frame(.dummy = 1)
    n_surfaces = 1
  }

  colors = surface_palette(n_surfaces)
  x_vals = sort(unique(grid$.x_num))
  y_vals = sort(unique(grid$.y_num))

  # --- build plotly ---
  p = plotly::plot_ly()
  z_matrices = list()
  surface_labels = character(n_surfaces)

  for (i in seq_len(n_surfaces)) {
    if (!is.null(facet_by) && length(facet_by) > 0) {
      combo = by_combos[i, , drop = FALSE]
      label = make_by_label(combo, facet_by, binned_by)

      mask = rep(TRUE, nrow(grid))
      for (v in facet_by) {
        mask = mask & (as.character(grid[[v]]) == as.character(combo[[v]]))
      }
      sub = grid[mask, ]
    } else {
      label = "predicted"
      sub = grid
    }

    surface_labels[i] = label

    # pivot to z matrix: rows = y, cols = x
    z_mat = matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
    for (r in seq_len(nrow(sub))) {
      xi = match(sub$.x_num[r], x_vals)
      yi = match(sub$.y_num[r], y_vals)
      if (!is.na(xi) && !is.na(yi)) {
        z_mat[yi, xi] = sub$.pred[r]
      }
    }
    z_matrices[[i]] = z_mat

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
      name = label,
      legendgroup = label,
      hovertemplate = paste0(
        x_lab, ": %{x}<br>",
        y_lab, ": %{y}<br>",
        z_lab, ": %{z:.3f}<br>",
        label,
        "<extra></extra>"
      )
    )
  }

  # --- crossings ---
  if (show_crossings && n_surfaces >= 2) {
    crossing_data = compute_crossings(z_matrices, x_vals, y_vals,
                                      by_combos, facet_by, binned_by,
                                      tolerance = NULL)
    if (nrow(crossing_data) > 0) {
      p = p %>% plotly::add_markers(
        data = crossing_data,
        x = ~cx, y = ~cy, z = ~cz,
        marker = list(size = 3, color = "red", symbol = "cross",
                      line = list(width = 1, color = "darkred")),
        name = "crossings",
        showlegend = TRUE,
        hovertemplate = paste0(
          "CROSSING<br>",
          x_lab, ": %{x}<br>",
          y_lab, ": %{y}<br>",
          z_lab, ": %{z:.3f}<br>",
          "%{text}",
          "<extra></extra>"
        ),
        text = ~pair_label
      )

      # --- contour projection of crossings onto floor ---
      if (show_contour) {
        all_z = unlist(lapply(z_matrices, c))
        floor_z = if (!is.null(contour_z)) contour_z else min(all_z, na.rm = TRUE)

        p = p %>% plotly::add_markers(
          data = crossing_data,
          x = ~cx, y = ~cy,
          z = rep(floor_z, nrow(crossing_data)),
          marker = list(size = 2, color = "red", opacity = 0.4,
                        symbol = "circle"),
          name = "crossing contour",
          showlegend = TRUE,
          hovertemplate = paste0(
            "CROSSING (projected)<br>",
            x_lab, ": %{x}<br>",
            y_lab, ": %{y}<br>",
            "%{text}",
            "<extra></extra>"
          ),
          text = ~pair_label
        )
      }
    }
  }

  # --- observed data points, color-coded by facet_by ---
  if (show_points) {
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

    if (!is.null(facet_by) && length(facet_by) > 0) {
      # assign each observation to the nearest surface (for coloring)
      obs$.group = assign_obs_to_surface(obs, facet_by, by_combos, binned_by)

      for (i in seq_len(n_surfaces)) {
        label = surface_labels[i]
        sub_obs = obs[obs$.group == i, ]
        if (nrow(sub_obs) == 0) next

        p = p %>% plotly::add_markers(
          data = sub_obs,
          x = ~.xn, y = ~.yn, z = ~.z,
          marker = list(size = 4, color = colors[i], opacity = 0.7,
                        line = list(width = 0.5, color = "black")),
          name = paste0("obs: ", label),
          legendgroup = label,
          showlegend = TRUE,
          hovertemplate = paste0(
            "OBSERVED<br>",
            x_lab, ": %{x}<br>",
            y_lab, ": %{y}<br>",
            response_name, ": %{z:.3f}<br>",
            label,
            "<extra></extra>"
          )
        )
      }
    } else {
      p = p %>% plotly::add_markers(
        data = obs,
        x = ~.xn, y = ~.yn, z = ~.z,
        marker = list(size = 4, color = "black", opacity = 0.7),
        name = "observed",
        hovertemplate = paste0(
          "OBSERVED<br>",
          x_lab, ": %{x}<br>",
          y_lab, ": %{y}<br>",
          response_name, ": %{z:.3f}",
          "<extra></extra>"
        )
      )
    }
  }

  # --- layout ---
  auto_title = if (is.null(title)) {
    if (!is.null(facet_by)) {
      paste0(z_lab, " ~ ", x, " * ", y, " | ",
             paste(facet_by, collapse = " + "))
    } else {
      paste0(z_lab, " ~ ", x, " * ", y)
    }
  } else {
    title
  }

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

  # attach metadata for downstream use
  attr(p, "ixsurface_meta") = list(
    x = x, y = y, facet_by = facet_by,
    n_surfaces = n_surfaces,
    surface_labels = surface_labels,
    z_matrices = z_matrices,
    x_vals = x_vals, y_vals = y_vals,
    binned_by = binned_by
  )

  p
}


#' Assign observations to the nearest surface group
#'
#' For color-coding observed data by facet_by level. Handles both categorical
#' and binned continuous facet_by variables.
#'
#' @param obs Data frame of observations.
#' @param by_vars Character vector of facet_by variable names.
#' @param by_combos Data frame of surface-level combinations.
#' @param binned_by Named list of bin info.
#' @return Integer vector of surface indices (1-indexed).
#' @keywords internal
assign_obs_to_surface = function(obs, by_vars, by_combos, binned_by) {
  n_obs = nrow(obs)
  n_surf = nrow(by_combos)
  assignments = integer(n_obs)

  for (k in seq_len(n_obs)) {
    best_i = 1
    best_dist = Inf

    for (i in seq_len(n_surf)) {
      dist = 0
      match_ok = TRUE

      for (v in by_vars) {
        obs_val = obs[[v]][k]
        combo_val = by_combos[[v]][i]

        if (v %in% names(binned_by)) {
          # continuous: distance from midpoint
          dist = dist + abs(as.numeric(obs_val) - as.numeric(combo_val))
        } else {
          # categorical: exact match or penalty
          if (as.character(obs_val) != as.character(combo_val)) {
            match_ok = FALSE
            break
          }
        }
      }

      if (match_ok && dist < best_dist) {
        best_dist = dist
        best_i = i
      }
    }
    assignments[k] = best_i
  }
  assignments
}
