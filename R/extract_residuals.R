#' Extract honest residuals from a fitted model
#'
#' Extracts residuals for use in performance-based visualizations. The
#' source of residuals depends on the model object type:
#'
#' * For `caret` `train` objects with `savePredictions = "final"`,
#'   returns cross-validated residuals.
#' * For `caret` `train` objects without saved predictions, returns
#'   training residuals with a warning.
#' * For other model objects (base `lm`, `glm`, `rlm`, `gam`, etc.),
#'   returns training residuals with a warning.
#'
#' The returned residuals are signed (observed - predicted), aligned with
#' the original data row order.
#'
#' @param model A fitted model object.
#' @param data The original data frame the model was fit on.
#' @param response Character. Name of the response variable in `data`.
#' @param verbose Logical. If `TRUE` (default), emit warnings about the
#'   residual source.
#'
#' @return A list with:
#'   \describe{
#'     \item{`residuals`}{Numeric vector of signed residuals, length
#'       `nrow(data)`.}
#'     \item{`source`}{Character: `"cv"` or `"training"`.}
#'     \item{`description`}{Human-readable description of the
#'       residual source.}
#'   }
#'
#' @details
#' For honest regional model comparison, cross-validated residuals are
#' required. Training residuals on flexible models (like kNN) are
#' systematically smaller than true generalization error, which would
#' bias performance comparisons toward flexible models.
#'
#' To get CV residuals from caret, fit with `savePredictions = "final"`
#' in `trainControl()`.
#'
#' @examples
#' \donttest{
#' tr = caret::trainControl(method = "cv", number = 10,
#'                          savePredictions = "final")
#' fit = caret::train(mpg ~ weight + year, data = Auto,
#'                    method = "lm", trControl = tr)
#' result = extract_residuals(fit, Auto, "mpg")
#' result$source  # "cv"
#' }
#'
#' @export
extract_residuals = function(model, data, response, verbose = TRUE) {

  # caret train object: check for stored CV predictions
  if (inherits(model, "train")) {
    if (!is.null(model$pred) && nrow(model$pred) > 0) {
      return(.extract_caret_cv(model, data, response))
    }
    if (verbose) {
      warning(
        "caret train object has no saved CV predictions. ",
        "Using training residuals (biased toward flexible models). ",
        "For honest comparison, refit with savePredictions = 'final' ",
        "in trainControl().",
        call. = FALSE
      )
    }
    preds = safe_predict(model, data)
    return(list(
      residuals = data[[response]] - preds,
      source = "training",
      description = "training residuals (caret object without saved CV predictions)"
    ))
  }

  # any other model with a predict method
  if (verbose) {
    warning(
      "Model is not a caret train object. ",
      "Using training residuals, which understate true generalization error. ",
      "Consider fitting via caret::train() with savePredictions = 'final'.",
      call. = FALSE
    )
  }

  preds = tryCatch(
    safe_predict(model, data),
    error = function(e) {
      stop("Could not extract predictions from model: ", e$message,
           call. = FALSE)
    }
  )

  list(
    residuals = data[[response]] - preds,
    source = "training",
    description = paste("training residuals from", class(model)[1], "object")
  )
}


#' Extract CV residuals from a caret train object
#'
#' @keywords internal
.extract_caret_cv = function(model, data, response) {
  cv_preds = model$pred

  # filter to best-tune row (if hyperparameter search was performed)
  best_tune = model$bestTune
  if (!is.null(best_tune) && nrow(best_tune) == 1) {
    for (param in names(best_tune)) {
      cv_preds = cv_preds[cv_preds[[param]] == best_tune[[param]], ]
    }
  }

  cv_preds = cv_preds[order(cv_preds$rowIndex), ]

  # handle repeated CV: average predictions per observation
  if (any(duplicated(cv_preds$rowIndex))) {
    cv_preds = stats::aggregate(
      pred ~ rowIndex + obs,
      data = cv_preds,
      FUN = mean
    )
    cv_preds = cv_preds[order(cv_preds$rowIndex), ]
  }

  residuals = rep(NA_real_, nrow(data))
  residuals[cv_preds$rowIndex] = cv_preds$obs - cv_preds$pred

  n_missing = sum(is.na(residuals))
  if (n_missing > 0) {
    warning(
      n_missing, " observation(s) not covered by any CV fold; ",
      "their residuals will be NA.",
      call. = FALSE
    )
  }

  cv_method = if (!is.null(model$control$method)) {
    model$control$method
  } else "cv"

  list(
    residuals = residuals,
    source = "cv",
    description = paste0("cross-validated residuals (", cv_method, ")")
  )
}
