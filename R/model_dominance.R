#' Map which model performs best across the predictor space
#'
#' Renders a 2D heatmap of the predictor space, colored by which model
#' has the smallest mean absolute residual in each region. Uses honest
#' (cross-validated) residuals when available via [extract_residuals()],
#' falling back to training residuals with a warning otherwise.
#'
#' For each grid cell, the function finds nearby observed points and
#' identifies the model with the smallest mean absolute residual at
#' those points. The heatmap is clipped to the convex hull of the
#' observed data by default, since regions far from observations don't
#' have meaningful "winners".
#'
#' Companion to [model_surface()] (which visualizes model behavior) —
#' this function visualizes model performance.
#'
#' @param models Named list of fitted model objects. Each must have a
#'   `predict` method. For honest comparison, use `caret::train()`
#'   objects fit with `savePredictions = "final"` in `trainControl()`.
#' @param x Character. Variable name mapped to the x-axis. Must be a
#'   predictor in all models and a numeric column in `data`.
#' @param y Character. Variable name mapped to the y-axis. Must be a
#'   predictor in all models and a numeric column in `data`.
#' @param response Character. Name of the response variable.
#' @param data Data frame containing observed predictor values and
#'   response. If `NULL` (default), `model.frame(models[[1]])` is used.
#' @param n Integer. Grid resolution per axis (default 60).
#' @param k_nearest Integer. Number of nearest observed points to average
#'   over for each grid cell (default 5). Higher = smoother regions.
#'   Acts as a bandwidth parameter.
#' @param show_points Logical. Overlay observed data points on the
#'   heatmap (default `TRUE`).
#' @param show_ties Logical. Mark grid cells where two or more models
#'   have nearly equal residuals as `"tied"` instead of picking a
#'   winner (default `TRUE`).
#' @param tie_threshold Numeric. Two models are considered tied if their
#'   mean residuals differ by less than this fraction of the response
#'   standard deviation (default 0.05).
#' @param clip_to_hull Logical. Restrict heatmap to the convex hull of
#'   observed data (default `TRUE`).
#' @param palette Character. Color palette from RColorBrewer (default
#'   `"Set2"`).
#' @param labs Named list for axis labels.
#' @param title Character or `NULL`. Plot title.
#' @param verbose Logical. Emit warnings about residual sources
#'   (default `TRUE`).
#'
#' @return A `ggplot2` plot object. The grid data with winners is
#'   accessible via `plot$data`.
#'
#' @details
#' This function does *not* compute predictions on a grid. Instead, it
#' uses each model's residuals at the observed data points (CV residuals
#' when available) to estimate regional performance via nearest-neighbor
#' averaging. This is more honest than grid-based prediction because
#' it doesn't extrapolate model behavior to regions where no data exists.
#'
#' The `k_nearest` parameter acts as a smoothing bandwidth:
#'
#' * Small `k_nearest` (1-5): high resolution, more noise-sensitive
#' * Medium `k_nearest` (5-25): typical use, balances local detail
#'   with smoothing
#' * Large `k_nearest` (\eqn{\to n}): converges toward global
#'   comparison, regional structure is lost
#'
#' @examples
#' \donttest{
#' tr = caret::trainControl(method = "cv", number = 10,
#'                          savePredictions = "final")
#'
#' fit_lm = caret::train(mpg ~ weight + year, data = Auto,
#'                       method = "lm", trControl = tr)
#' fit_knn = caret::train(mpg ~ weight + year, data = Auto,
#'                        method = "knn", trControl = tr,
#'                        tuneGrid = data.frame(k = 20))
#'
#' model_dominance(
#'   models = list(OLS = fit_lm, kNN = fit_knn),
#'   x = "weight", y = "year",
#'   response = "mpg", data = Auto
#' )
#' }
#'
#' @seealso [model_surface()] for visualizing model behavior,
#'   [extract_residuals()] for the underlying residual extraction.
#'
#' @export
model_dominance = function(models, x, y, response,
                           data = NULL,
                           n = 60,
                           k_nearest = 5,
                           show_points = TRUE,
                           show_ties = TRUE,
                           tie_threshold = 0.05,
                           clip_to_hull = TRUE,
                           palette = "Set2",
                           labs = NULL,
                           title = NULL,
                           verbose = TRUE) {

  if (!is.list(models) || is.null(names(models)) || length(models) < 2) {
    stop("`models` must be a named list of at least two fitted models.")
  }

  # default data: use the model frame from the first model
  if (is.null(data)) {
    data = model.frame(models[[1]])
  }

  for (v in c(x, y, response)) {
    if (!v %in% names(data)) {
      stop(sprintf("Variable '%s' not found in data.", v))
    }
  }

  # extract honest residuals for each model
  resid_info = lapply(models, function(m) {
    extract_residuals(m, data, response, verbose = verbose)
  })

  # warn loudly if residual sources are mixed
  sources = vapply(resid_info, function(r) r$source, character(1))
  if (length(unique(sources)) > 1) {
    warning(
      "Models use mixed residual sources (",
      paste(names(models), "=", sources, collapse = ", "),
      "). Comparison may be misleading.",
      call. = FALSE
    )
  }

  # absolute residuals as a matrix: n_obs x n_models
  abs_resids = vapply(resid_info, function(r) abs(r$residuals),
                      numeric(nrow(data)))

  # axis labels
  x_lab = if (!is.null(labs$x)) labs$x else x
  y_lab = if (!is.null(labs$y)) labs$y else y

  # build grid
  x_vals = seq(min(data[[x]]), max(data[[x]]), length.out = n)
  y_vals = seq(min(data[[y]]), max(data[[y]]), length.out = n)
  grid = expand.grid(setNames(list(x_vals, y_vals), c(x, y)))

  # standardize for distance calculation
  obs_x_z = (data[[x]] - mean(data[[x]])) / stats::sd(data[[x]])
  obs_y_z = (data[[y]] - mean(data[[y]])) / stats::sd(data[[y]])

  resid_tol = tie_threshold * stats::sd(data[[response]], na.rm = TRUE)

  # for each grid cell: nearest-neighbor mean residuals per model
  winners = character(nrow(grid))
  win_margin = numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    gx_z = (grid[i, 1] - mean(data[[x]])) / stats::sd(data[[x]])
    gy_z = (grid[i, 2] - mean(data[[y]])) / stats::sd(data[[y]])

    dists = sqrt((obs_x_z - gx_z)^2 + (obs_y_z - gy_z)^2)
    nn_idx = order(dists)[seq_len(min(k_nearest, nrow(data)))]

    mean_resids = colMeans(abs_resids[nn_idx, , drop = FALSE], na.rm = TRUE)
    sorted = sort(mean_resids)
    best = which.min(mean_resids)

    if (show_ties && length(mean_resids) >= 2 &&
        (sorted[2] - sorted[1]) < resid_tol) {
      winners[i] = "tied"
    } else {
      winners[i] = names(models)[best]
    }

    win_margin[i] = if (length(mean_resids) >= 2) sorted[2] - sorted[1] else NA
  }

  grid$winner = factor(winners,
                       levels = c(names(models),
                                  if (show_ties) "tied" else NULL))
  grid$margin = win_margin

  # clip to convex hull
  if (clip_to_hull) {
    hull_idx = grDevices::chull(data[[x]], data[[y]])
    hull = data[hull_idx, c(x, y)]
    inside = .point_in_polygon(grid[[x]], grid[[y]],
                               hull[[x]], hull[[y]])
    grid = grid[inside, , drop = FALSE]
  }

  # subtitle reflects residual source
  source_label = if (all(sources == "cv")) {
    "cross-validated residuals"
  } else if (all(sources == "training")) {
    "training residuals (biased — refit with savePredictions='final' for honest comparison)"
  } else {
    paste0("mixed: ",
           paste(names(models), sources, sep = "=", collapse = ", "))
  }

  auto_title = if (is.null(title)) "regional model dominance" else title
  subtitle = paste0(
    "color = model with smallest mean |residual| over ",
    k_nearest, " nearest observations; ", source_label
  )

  # build heatmap
  p = ggplot2::ggplot(
    grid,
    ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[["winner"]])
  ) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::scale_fill_brewer(palette = palette, name = "best model",
                               na.translate = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = auto_title, subtitle = subtitle,
                  x = x_lab, y = y_lab)

  if (show_points) {
    p = p + ggplot2::geom_point(
      data = data,
      ggplot2::aes(x = .data[[x]], y = .data[[y]]),
      inherit.aes = FALSE,
      color = "black", size = 0.6, alpha = 0.5
    )
  }

  attr(p, "ixsurface_meta") = list(
    type = "model_dominance",
    x = x, y = y, response = response,
    n_models = length(models),
    model_names = names(models),
    residual_sources = sources,
    k_nearest = k_nearest,
    tie_threshold = tie_threshold
  )

  p
}


#' Point-in-polygon test using ray casting
#'
#' @keywords internal
.point_in_polygon = function(px, py, vx, vy) {
  n_pts = length(px)
  n_vert = length(vx)
  inside = logical(n_pts)

  for (i in seq_len(n_pts)) {
    j = n_vert
    in_poly = FALSE
    for (k in seq_len(n_vert)) {
      if (((vy[k] > py[i]) != (vy[j] > py[i])) &&
          (px[i] < (vx[j] - vx[k]) * (py[i] - vy[k]) /
             (vy[j] - vy[k]) + vx[k])) {
        in_poly = !in_poly
      }
      j = k
    }
    inside[i] = in_poly
  }
  inside
}
