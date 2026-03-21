#' Compute approximate surface crossing points
#'
#' For each pair of surfaces, finds grid cells where the z-difference changes
#' sign, indicating the surfaces cross in that region.
#'
#' @param z_matrices List of z matrices (one per surface).
#' @param x_vals Numeric vector of x grid values.
#' @param y_vals Numeric vector of y grid values.
#' @param by_combos Data frame of facet_by combinations.
#' @param by_vars Character vector of facet_by variable names.
#' @param binned_by Named list of bin info for continuous by-variables.
#' @param tolerance Numeric or NULL. Adaptive if NULL (2 percent of z-range).
#' @return A data.frame with columns cx, cy, cz, pair_label.
#' @keywords internal
compute_crossings = function(z_matrices, x_vals, y_vals,
                             by_combos, by_vars, binned_by = list(),
                             tolerance = NULL) {

  n_surfaces = length(z_matrices)
  empty = data.frame(cx = numeric(0), cy = numeric(0),
                     cz = numeric(0), pair_label = character(0),
                     stringsAsFactors = FALSE)
  if (n_surfaces < 2) return(empty)

  all_z = unlist(lapply(z_matrices, c))
  z_range = diff(range(all_z, na.rm = TRUE))
  if (is.null(tolerance)) {
    tolerance = z_range * 0.02
  }

  results = list()
  pair_idx = 0

  pairs = utils::combn(n_surfaces, 2)
  for (p in seq_len(ncol(pairs))) {
    i = pairs[1, p]
    j = pairs[2, p]

    diff_mat = z_matrices[[i]] - z_matrices[[j]]
    nr = nrow(diff_mat)
    nc = ncol(diff_mat)

    if (!is.null(by_vars) && length(by_vars) > 0) {
      lab_i = make_by_label(by_combos[i, , drop = FALSE], by_vars, binned_by)
      lab_j = make_by_label(by_combos[j, , drop = FALSE], by_vars, binned_by)
    } else {
      lab_i = paste0("surface_", i)
      lab_j = paste0("surface_", j)
    }
    pair_lab = paste0(lab_i, " vs ", lab_j)

    crossing_cells = matrix(FALSE, nrow = nr, ncol = nc)

    # scan x-direction
    for (row in seq_len(nr)) {
      for (col in seq_len(nc - 1)) {
        v1 = diff_mat[row, col]
        v2 = diff_mat[row, col + 1]
        if (!is.na(v1) && !is.na(v2)) {
          if (v1 * v2 < 0 || abs(v1) < tolerance) {
            crossing_cells[row, col] = TRUE
          }
        }
      }
    }

    # scan y-direction
    for (row in seq_len(nr - 1)) {
      for (col in seq_len(nc)) {
        v1 = diff_mat[row, col]
        v2 = diff_mat[row + 1, col]
        if (!is.na(v1) && !is.na(v2)) {
          if (v1 * v2 < 0 || abs(v1) < tolerance) {
            crossing_cells[row, col] = TRUE
          }
        }
      }
    }

    for (row in seq_len(nr)) {
      for (col in seq_len(nc)) {
        if (crossing_cells[row, col]) {
          z_avg = mean(c(z_matrices[[i]][row, col],
                         z_matrices[[j]][row, col]), na.rm = TRUE)
          pair_idx = pair_idx + 1
          results[[pair_idx]] = data.frame(
            cx = x_vals[col],
            cy = y_vals[row],
            cz = z_avg,
            pair_label = pair_lab,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(results) == 0) return(empty)
  do.call(rbind, results)
}


#' Find Crossing Regions Between Interaction Surfaces
#'
#' Identifies where predicted response surfaces cross for different levels of
#' conditioning factors. Returns a data frame of approximate crossing locations.
#'
#' @param model A fitted model object.
#' @param x Character. First focal variable.
#' @param y Character. Second focal variable.
#' @param facet_by Character vector. Conditioning variable(s).
#' @param n Integer. Grid resolution (default 50).
#' @param n_bins Integer. Bins for continuous facet_by (default 3).
#' @param bin_method Character. Binning method.
#' @param tolerance Numeric or NULL. Crossing detection tolerance.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{cx}{x-coordinate of crossing}
#'   \item{cy}{y-coordinate of crossing}
#'   \item{cz}{predicted response at crossing (average of both surfaces)}
#'   \item{pair_label}{which surface pair crosses}
#' }
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(design = "mixed", seed = 42)
#' fit = lm(y ~ temp * pressure * catalyst, data = dat)
#' crossings = find_crossings(fit, "temp", "pressure", "catalyst")
#' head(crossings)
#' }
#'
#' @export
find_crossings = function(model, x, y, facet_by,
                          n = 50, n_bins = 3,
                          bin_method = "quantile",
                          tolerance = NULL) {

  grid_result = make_prediction_grid(model, x, y, facet_by = facet_by,
                                     n = n, n_bins = n_bins,
                                     bin_method = bin_method)
  grid = grid_result$grid
  binned_by = grid_result$binned_by
  grid$.pred = safe_predict(model, grid)

  mf = model.frame(model)
  x_type = detect_factor_type(model, x)
  y_type = detect_factor_type(model, y)

  if (x_type %in% c("categorical", "quasi_categorical")) {
    x_levels = if (is.factor(mf[[x]])) levels(mf[[x]]) else sort(unique(mf[[x]]))
    grid$.x_num = match(as.character(grid[[x]]), as.character(x_levels))
  } else {
    grid$.x_num = grid[[x]]
  }

  if (y_type %in% c("categorical", "quasi_categorical")) {
    y_levels = if (is.factor(mf[[y]])) levels(mf[[y]]) else sort(unique(mf[[y]]))
    grid$.y_num = match(as.character(grid[[y]]), as.character(y_levels))
  } else {
    grid$.y_num = grid[[y]]
  }

  by_combos = unique(grid[, facet_by, drop = FALSE])
  n_surfaces = nrow(by_combos)
  x_vals = sort(unique(grid$.x_num))
  y_vals = sort(unique(grid$.y_num))

  z_matrices = list()
  for (i in seq_len(n_surfaces)) {
    combo = by_combos[i, , drop = FALSE]
    mask = rep(TRUE, nrow(grid))
    for (v in facet_by) {
      mask = mask & (as.character(grid[[v]]) == as.character(combo[[v]]))
    }
    sub = grid[mask, ]

    z_mat = matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
    for (r in seq_len(nrow(sub))) {
      xi = match(sub$.x_num[r], x_vals)
      yi = match(sub$.y_num[r], y_vals)
      if (!is.na(xi) && !is.na(yi)) {
        z_mat[yi, xi] = sub$.pred[r]
      }
    }
    z_matrices[[i]] = z_mat
  }

  compute_crossings(z_matrices, x_vals, y_vals, by_combos, facet_by,
                    binned_by, tolerance)
}


#' Plot Crossing Regions as a Standalone 3D Scatter
#'
#' Visualizes only the crossing points between interaction surfaces as an
#' interactive 3D scatter plot, color-coded by surface pair. This isolates
#' where interaction effects are strongest, without the surfaces themselves.
#'
#' @inheritParams find_crossings
#' @param labs Named list for axis labels, e.g.,
#'   \code{list(x = "Temperature", y = "Pressure", z = "Yield")}.
#' @param title Character or NULL. Plot title.
#' @param marker_size Numeric. Marker size (default 3).
#' @param marker_opacity Numeric in \code{[0, 1]}. Marker opacity (default 0.7).
#'
#' @return A \code{plotly} htmlwidget object. If no crossings are found,
#'   returns an empty plot with a "No crossings detected" annotation.
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(design = "mixed", seed = 42)
#' fit = lm(y ~ temp * pressure * catalyst, data = dat)
#' plot_crossings(fit, "temp", "pressure", "catalyst")
#'
#' # with custom labels
#' plot_crossings(fit, "temp", "pressure", "catalyst",
#'                labs = list(x = "Temp (C)", y = "Press (psi)", z = "Yield"))
#' }
#'
#' @export
plot_crossings = function(model, x, y, facet_by,
                          n = 50, n_bins = 3,
                          bin_method = "quantile",
                          tolerance = NULL,
                          labs = NULL, title = NULL,
                          marker_size = 3,
                          marker_opacity = 0.7) {

  cx = find_crossings(model, x, y, facet_by,
                      n = n, n_bins = n_bins,
                      bin_method = bin_method,
                      tolerance = tolerance)

  x_lab = if (!is.null(labs$x)) labs$x else x
  y_lab = if (!is.null(labs$y)) labs$y else y
  z_lab = if (!is.null(labs$z)) labs$z else paste0("predicted ", names(model.frame(model))[1])

  auto_title = if (is.null(title)) {
    paste0("Crossings: ", x, " x ", y, " | ",
           paste(facet_by, collapse = " + "))
  } else {
    title
  }

  scene = list(
    xaxis = list(title = x_lab),
    yaxis = list(title = y_lab),
    zaxis = list(title = z_lab),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))
  )

  if (nrow(cx) == 0) {
    p = plot_ly() %>%
      layout(
        title = list(text = auto_title, font = list(size = 14)),
        scene = scene,
        annotations = list(list(
          text = "No crossings detected",
          showarrow = FALSE, xref = "paper", yref = "paper",
          x = 0.5, y = 0.5, font = list(size = 16, color = "grey")
        ))
      )
    return(p)
  }

  pairs = unique(cx$pair_label)
  colors = surface_palette(length(pairs))
  color_map = setNames(colors, pairs)

  p = plot_ly()

  for (i in seq_along(pairs)) {
    sub = cx[cx$pair_label == pairs[i], ]
    p = p %>% add_markers(
      data = sub,
      x = ~cx, y = ~cy, z = ~cz,
      marker = list(size = marker_size, color = colors[i],
                    opacity = marker_opacity,
                    line = list(width = 0.5, color = "black")),
      name = pairs[i],
      hovertemplate = paste0(
        x_lab, ": %{x}<br>",
        y_lab, ": %{y}<br>",
        z_lab, ": %{z:.3f}<br>",
        pairs[i],
        "<extra></extra>"
      )
    )
  }

  p = p %>% layout(
    title = list(text = auto_title, font = list(size = 14)),
    scene = scene,
    legend = list(orientation = "h", yanchor = "bottom", y = -0.15)
  )

  p
}
