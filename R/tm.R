#' The Transformation Distribution.
#'
#' Density and cumulative distribution function of the transformation
#' distribution.
#'
#' @param z,z_deriv Numeric arrays of transformed values and transformation
#'   derivative evaluations.
#' @param reference_density,reference_cdf Density and cumulative distribution
#'   function of the reference distribution.
#' @param log Logical; if `TRUE`, log densities are returned.
#' @param ... Additional arguments passed on to `reference_density`
#'
#' @returns Numeric array of density evaluations.
#' @export
#'
#' @examples
#' knots <- equidistant_knots(c(-4, 4), p = 10)
#' approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)
#' spline <- get_extrap_bspline_and_deriv(approx_basis_and_deriv)
#'
#' # evaluate spline
#' x <- rnorm(n = 100)
#' coef <- rnorm(n = 10)
#' result <- spline(x, coef)
#'
#' # evaluate transformation density
#' pdf_evaluations <- dtramo(result$fx, result$fx_deriv)
#'
#' # evaluate transformation CDF
#' cdf_evaluations <- ptramo(result$fx)
dtramo <- function(z, z_deriv, reference_density = stats::dnorm, log = FALSE, ...) {

    log_prob <- reference_density(z, log = TRUE, ...) + log(z_deriv)
    if (!log) return(exp(log_prob))

    log_prob
}

#' @rdname dtramo
#' @export
ptramo <- function(z, reference_cdf = stats::pnorm, ...) {
    reference_cdf(z, ...)
}
