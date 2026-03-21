#' Simulate a Multi-Factor Experimental Dataset
#'
#' Generates synthetic data from a factorial design with known interaction
#' structure. Useful for demonstrating and testing \code{interaction_surface}.
#'
#' @param n Integer. Total number of observations (default 200).
#' @param design Character. One of:
#' \describe{
#'   \item{"mixed"}{Two continuous (temp, pressure) + one categorical (catalyst).}
#'   \item{"continuous"}{Three continuous factors (temp, pressure, speed).}
#'   \item{"categorical"}{Three categorical factors (catalyst, operator, shift).}
#' }
#' @param noise Numeric. Standard deviation of Gaussian noise (default 0.5).
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A data.frame with factor columns and a response column \code{y}.
#'
#' @details
#' The data generating process includes main effects for all factors, two-way
#' interactions between the first two factors, and a three-way interaction
#' (weaker) involving all factors. This makes it straightforward to verify that
#' \code{interaction_surface} correctly detects the embedded structure.
#'
#' @examples
#' dat = sim_factorial(design = "mixed", seed = 42)
#' head(dat)
#'
#' @export
sim_factorial = function(n = 200, design = c("mixed", "continuous", "categorical"),
                         noise = 0.5, seed = NULL) {

  design = match.arg(design)
  if (!is.null(seed)) set.seed(seed)

  if (design == "continuous") {
    temp = runif(n, 150, 250)
    pressure = runif(n, 10, 50)
    speed = runif(n, 100, 500)

    t_s = (temp - 200) / 50
    p_s = (pressure - 30) / 20
    s_s = (speed - 300) / 200

    y = 50 +
      5 * t_s +
      3 * p_s +
      2 * s_s +
      4 * t_s * p_s +
      1.5 * t_s * s_s +
      1.0 * t_s * p_s * s_s +
      rnorm(n, 0, noise)

    data.frame(temp = temp, pressure = pressure, speed = speed, y = y)

  } else if (design == "categorical") {
    catalyst = factor(sample(c("A", "B", "C"), n, replace = TRUE))
    operator = factor(sample(c("Op1", "Op2"), n, replace = TRUE))
    shift = factor(sample(c("Day", "Night"), n, replace = TRUE))

    cat_eff = ifelse(catalyst == "A", -1, ifelse(catalyst == "B", 0.5, 0.5))
    cat_eff2 = ifelse(catalyst == "A", 0, ifelse(catalyst == "B", -1, 1))
    op_eff = ifelse(operator == "Op1", -1, 1)
    sh_eff = ifelse(shift == "Day", -1, 1)

    y = 50 +
      3 * cat_eff +
      2 * cat_eff2 +
      2 * op_eff +
      1 * sh_eff +
      3 * cat_eff * op_eff +
      1.5 * op_eff * sh_eff +
      0.8 * cat_eff * op_eff * sh_eff +
      rnorm(n, 0, noise)

    data.frame(catalyst = catalyst, operator = operator, shift = shift, y = y)

  } else {
    # mixed: 2 continuous + 1 categorical
    temp = runif(n, 150, 250)
    pressure = runif(n, 10, 50)
    catalyst = factor(sample(c("A", "B", "C"), n, replace = TRUE))

    t_s = (temp - 200) / 50
    p_s = (pressure - 30) / 20
    cat_eff = ifelse(catalyst == "A", -2, ifelse(catalyst == "B", 0.5, 1.5))

    y = 50 +
      5 * t_s +
      3 * p_s +
      cat_eff +
      4 * t_s * p_s +
      2.5 * t_s * cat_eff +
      1.0 * t_s * p_s * cat_eff +
      rnorm(n, 0, noise)

    data.frame(temp = temp, pressure = pressure, catalyst = catalyst, y = y)
  }
}
