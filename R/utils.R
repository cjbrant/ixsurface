#' Detect whether a variable in a model is continuous or categorical
#'
#' @param model A fitted model object.
#' @param varname Character string naming the variable.
#' @return Character: "continuous" or "categorical".
#' @keywords internal
detect_factor_type = function(model, varname) {
  mf = model.frame(model)
  if (!varname %in% names(mf)) {
    stop(sprintf("Variable '%s' not found in model frame.", varname))
  }
  v = mf[[varname]]
  if (is.factor(v) || is.character(v) || is.logical(v)) {
    return("categorical")
  }
  if (is.numeric(v) && length(unique(v)) <= 6) {
    return("quasi_categorical")
  }
  return("continuous")
}


#' Bin a continuous variable into discrete levels
#'
#' Used when a continuous variable appears in \code{facet_by} to produce a
#' manageable number of surfaces.
#'
#' @param x Numeric vector to bin.
#' @param n_bins Integer. Number of bins (default 3).
#' @param method Character. One of \code{"quantile"} (equal-count),
#'   \code{"equal"} (equal-width), or \code{"pretty"} (round breakpoints).
#' @return A factor with descriptive level labels.
#'
#' @examples
#' x = rnorm(100, mean = 50, sd = 15)
#' bin_continuous(x, n_bins = 3, method = "quantile")
#' bin_continuous(x, n_bins = 4, method = "equal")
#'
#' @export
bin_continuous = function(x, n_bins = 3, method = c("quantile", "equal", "pretty")) {
  method = match.arg(method)

  breaks = switch(method,
    quantile = {
      probs = seq(0, 1, length.out = n_bins + 1)
      unique(quantile(x, probs, na.rm = TRUE))
    },
    equal = {
      seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n_bins + 1)
    },
    pretty = {
      pretty(x, n = n_bins)
    }
  )

  if (length(breaks) < 2) {
    return(factor(rep("all", length(x))))
  }

  cut(x, breaks = breaks, include.lowest = TRUE, dig.lab = 3)
}


#' Safe prediction wrapper with GLM/GAM support
#'
#' Dispatches \code{predict()} with \code{type = "response"} for GLMs and GAMs
#' to return predictions on the response scale (e.g., probabilities for logistic
#' regression). Falls back to plain \code{predict()} for lm and others.
#'
#' @param model A fitted model object.
#' @param newdata A data.frame for prediction.
#' @return Numeric vector of predictions.
#' @keywords internal
safe_predict = function(model, newdata) {
  if (inherits(model, "glm") || inherits(model, "gam")) {
    predict(model, newdata = newdata, type = "response")
  } else {
    predict(model, newdata = newdata)
  }
}


#' Build a prediction grid over two focal variables
#'
#' For continuous variables, generates a regular sequence across the observed
#' range. For categorical/factor variables, uses all observed levels. Continuous
#' \code{facet_by} variables are binned, and predictions use bin midpoints.
#'
#' @param model A fitted model object.
#' @param x Character. First focal variable (mapped to x-axis).
#' @param y Character. Second focal variable (mapped to y-axis).
#' @param facet_by Character vector or NULL. Conditioning variable(s) whose
#'   levels generate separate surfaces.
#' @param n Integer. Grid density for continuous axes.
#' @param n_bins Integer. Number of bins for continuous \code{facet_by} variables.
#' @param bin_method Character. Binning method: "quantile", "equal", or "pretty".
#' @return A list with components:
#'   \describe{
#'     \item{grid}{data.frame suitable for \code{predict()}}
#'     \item{binned_by}{named list of bin info for continuous facet_by variables}
#'   }
#'
#' @examples
#' \donttest{
#' dat = sim_factorial(design = "mixed", seed = 42)
#' fit = lm(y ~ temp * pressure * catalyst, data = dat)
#' result = make_prediction_grid(fit, x = "temp", y = "pressure",
#'                               facet_by = "catalyst", n = 10)
#' str(result$grid)
#' }
#'
#' @export
make_prediction_grid = function(model, x, y, facet_by = NULL, n = 50,
                                n_bins = 3, bin_method = "quantile") {
  mf = model.frame(model)
  predictor_vars = names(mf)[-1]
  by_vars = if (!is.null(facet_by)) facet_by else character(0)

  val_list = list()
  binned_by = list()

  for (v in predictor_vars) {
    col = mf[[v]]

    if (v %in% c(x, y)) {
      if (is.factor(col) || is.character(col)) {
        val_list[[v]] = if (is.factor(col)) levels(col) else sort(unique(col))
      } else {
        val_list[[v]] = seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE),
                            length.out = n)
      }

    } else if (v %in% by_vars) {
      if (is.factor(col) || is.character(col)) {
        val_list[[v]] = if (is.factor(col)) levels(col) else sort(unique(col))
      } else {
        # continuous facet_by: bin and use midpoints
        binned = bin_continuous(col, n_bins = n_bins, method = bin_method)
        bin_levels = levels(binned)
        midpoints = vapply(bin_levels, function(lev) {
          mean(col[binned == lev], na.rm = TRUE)
        }, numeric(1))

        binned_by[[v]] = list(
          levels = bin_levels,
          midpoints = midpoints,
          factor_col = binned
        )
        val_list[[v]] = midpoints
      }

    } else {
      # background: hold at median/mode
      if (is.factor(col)) {
        tab = table(col)
        val_list[[v]] = factor(names(tab)[which.max(tab)], levels = levels(col))
      } else if (is.character(col)) {
        tab = table(col)
        val_list[[v]] = names(tab)[which.max(tab)]
      } else {
        val_list[[v]] = median(col, na.rm = TRUE)
      }
    }
  }

  grid = expand.grid(val_list, stringsAsFactors = TRUE)
  list(grid = grid, binned_by = binned_by)
}


#' Generate a label for a facet_by combination
#'
#' @param row A single-row data.frame.
#' @param by_vars Character vector of by-variable names.
#' @param binned_by Named list of bin info for continuous by-variables.
#' @return A character label.
#' @keywords internal
make_by_label = function(row, by_vars, binned_by = list()) {
  parts = vapply(by_vars, function(v) {
    val = row[[v]]
    if (v %in% names(binned_by)) {
      info = binned_by[[v]]
      idx = which.min(abs(info$midpoints - as.numeric(val)))
      paste0(v, ": ", info$levels[idx])
    } else {
      paste0(v, "=", val)
    }
  }, character(1))
  paste(parts, collapse = ", ")
}


#' Get a qualitative color palette for surfaces
#'
#' @param n Number of colors needed.
#' @return Character vector of hex colors.
#' @keywords internal
surface_palette = function(n) {
  base_colors = c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
  )
  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  }
  grDevices::colorRampPalette(base_colors)(n)
}


#' Lighten a hex color
#'
#' @param hex_color Character hex color.
#' @param amount Numeric in \code{[0, 1]}. How much to lighten.
#' @return Character hex color.
#' @keywords internal
lighten_color = function(hex_color, amount = 0.3) {
  rgb_vals = grDevices::col2rgb(hex_color)[, 1] / 255
  lightened = rgb_vals + (1 - rgb_vals) * amount
  grDevices::rgb(lightened[1], lightened[2], lightened[3])
}
