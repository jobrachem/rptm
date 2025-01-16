

#' Create equidistant knots for an onion spline of order 3.
#'
#'
#' @param a_inner Inner left boundary
#' @param b_inner Inner right boundary
#' @param p Desired number of parameters.
#' @param extend_range_factor Slightly extends the range over `min(x)` and `max(x)`,
#'      by setting `xlow <- min(x) - (max(x) - min(x)) * extend_range_factor`
#'
#' @return List of attributes.
#' @export
onion_knots <- function(a_inner,
                        b_inner,
                        p,
                        extend_range_factor = 0.0) {
    l <- 3

    self <- list()
    self$a_inner <- a_inner
    self$b_inner <- b_inner
    self$nparam <- p
    self$order <- l
    self$nfixed_left <- 3
    self$nfixed_right <- 3
    self$n_inner_knots <- (p + 3 + 3 + 1) - l + 1

    knots <- equidistant_knots(
        x = c(a_inner, b_inner),
        p = p + 1,
        l = l,
        extend_range_factor = extend_range_factor
    )

    self$step <- mean(diff(knots))

    leftmost_knot <- knots[1]
    rightmost_knot <- knots[length(knots)]

    left_additional_knots <- seq(
        from = leftmost_knot - (self$nfixed_left * self$step),
        to = leftmost_knot - self$step,
        length.out = self$nfixed_left
    )

    right_additional_knots <- seq(
        from = rightmost_knot + self$step,
        to = rightmost_knot + (self$nfixed_right * self$step),
        length.out = self$nfixed_right
    )

    self$knots <- c(left_additional_knots, knots, right_additional_knots)

    self$intercept_knot <- self$knots[3]
    self$nparam_full_domain <- self$nparam + self$nfixed_left + self$nfixed_right + 1
    self$left <- self$knots[self$order + self$nfixed_left - 2]
    self$right <- self$knots[length(self$knots) - (self$order + self$nfixed_right - 1)]

    class(self) <- "OnionKnots"

    self
}

#' Computes the constant shift for onion spline coefficients.
#'
#' Subtracting this shift from unconstrained log increments ensures an average
#' slope of one for the onion spline.
#'
#' @param log_increments Unconstrained onion spline log increments.
#' @param knots A list defining the knots for the onion spline, most easily
#'      created with [onion_knots()].
#'
#' @returns A numeric scalar, the onion shift.
#' @export
#'
#' @examples
#' knots <- onion_knots(-4, 4, p = 10)
#' log_inc <- rnorm(10)
#' shift <- onion_shift(log_inc, knots)
onion_shift <- function(log_increments, knots) {
    l <- knots$knots[5]
    r <- rev(knots$knots)[5]
    shift <- log(sum(exp(log_increments))) - log((r - l - 2 * knots$step))
    shift
}


#' Computes the coefficient increments of an onion spline.
#'
#' Note that these are the increments of the spline coefficients, not the
#' full coefficients. To evaluate the onion spline, you have two options:
#' 1)
#'
#' @param log_increments Unconstrained onion spline log increments.
#' @param knots A list defining the knots for the onion spline, most easily
#'      created with [onion_knots()].
#' @param shift Numeric scalar for normalizing the onion spline's slope. Most
#'      easily created with [onion_shift()]. If `NULL` (default), the shift is
#'      computed automatically to normalize the spline's average slope to one.
#'
#' @returns Numeric vector of coefficient increments.
#' @export
#'
#' @examples
#' knots <- onion_knots(-4, 4, p = 10)
#' log_inc <- rnorm(10)
#' coef <- onion_coef(log_inc, knots)
#'
#' x <- runif(100, -4, 4)
#' B <- splines::splineDesign(knots$knots, x)
#'
#' # Option 1: Compute cumulative sum of coefficients.
#' fx <- B %*% cumsum(coef)
#'
#' # Option 2: Absorb cumulative sum into design matrix.
#'
onion_coef <- function(log_increments, knots, shift = NULL) {
    if (is.null(shift)) {
        shift <- onion_shift(log_increments = log_increments, knots = knots)
    }

    log_step <- rep(log(knots$step), times = 3)

    constrained_log_increments <- c(
        log_step,
        log_increments - shift,
        log_step
    )

    coef <- c(knots$intercept_knot, exp(constrained_log_increments))

    coef
}


#' Generates samples from a random walk of order 1.
#'
#' Assumes the starting condition to be zero.
#'
#' @param n Positive integer, the number of samples.
#' @param p Positive integer, the number of parameters.
#' @param sd The standard deviation of the random walk.
#'
#' @returns A numeric matrix of dimension `(n, p)`.
#' @export
#'
#' @examples
#' samples <- rrw1(1, p = 10)
rrw1 <- function(n, p, sd = 1) {
    samples <- rnorm(n * p) |> matrix(nrow = p)

    L <- matrix(0, p, p)
    L[lower.tri(L, diag = TRUE)] <- 1

    matrix(sd, nrow = n, ncol = p) * t(L %*% samples)
}
