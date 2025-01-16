#' Create equidistant knots.
#'
#' @param x Input vector
#' @param p Desired number of parameters.
#' @param l Order of the spline. Defaults to `l=3` for cubic splines.
#' @param extend_range_factor Slightly extends the range over `min(x)` and `max(x)`,
#'      by setting `xlow <- min(x) - (max(x) - min(x)) * extend_range_factor`
#'
#' @return Vector of knots. Includes the outer knots. Will be of length `p + l + 1`.
#' @export
equidistant_knots <- function(x, p = 10, l = 3, extend_range_factor = 0.01) {
    # m: Defines the order of the spline for which the knots are created.
    #    m = 2 means cubic splines.
    # k: Number of parameters.
    # extend_range_factor: Factor to stretch the range of x a little for
    #                      compatibility with splines of order 0

    m <- l - 1

    xrange <- max(x) - min(x)
    xlow <- min(x) - xrange * extend_range_factor
    xup <- max(x) + xrange * extend_range_factor

    dx <- (xup - xlow) / (p - m - 1)
    n_knots <- p + m + 2
    knots <- seq(xlow - dx * (m + 1), xup + dx * (m + 1), length.out = n_knots)

    knots
}


#' Generates functions to quickly approximate B-spline bases and their
#' derivatives.
#'
#' The functions work on a grid approximation, returning a linear interpolation
#' between the basis evaluations at the two closest grid points.
#'
#' The approximation functions will return `NA` for `x` outside the range
#' supported by the given knots, that is for those values `x` for which
#' `x < knots[order + 1]` or `rev(knots)[order + 1] < x` holds.
#'
#'
#' @param knots A numeric vector of knots.
#' @param ngrid A positive integer, indicating how many elements the
#'   approximation grid should have. Higher numbers improve the accuracy of the
#'   approximation.
#' @param postmultiply_by A numeric matrix with `p` rows, where `p` is the
#'   number of basis function evaluations, i.e. the number of columns of the
#'   basis matrix.
#' @param order A non-negative integer, giving the order of the spline. For a
#'   cubic spline, we have `order=3`.
#'
#' @returns A function that takes a single argument `x` (a numeric vector) and
#'   returns a list with four named elements:
#'
#' * `basis`: Approximated numeric matrix of basis function evaluations
#' * `basis_deriv`: Approximated numeric matrix of basis function derivative evaluations.
#' * `min_knot`: The minimum of the grid used to create the bases.
#' * `max_knot`: The maximum of the grid used to create the bases.
#'
#'
#' @export
#'
#' @examples
#' knots <- equidistant_knots(c(-4, 4), p = 10)
#' approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)
#'
#' # evaluate approximation
#' x <- rnorm(n = 10)
#' result <- approx_basis_and_deriv(x)
get_approx_basis_and_deriv_fun <- function(knots, ngrid, postmultiply_by = NULL, order = 3) {
    min_knot <- knots[order + 1]
    max_knot <- rev(knots)[order + 1]

    grid <- seq(min_knot, max_knot, length.out = ngrid)
    grid_step <- grid |>
        diff() |>
        mean()

    grid_basis <- splines::splineDesign(knots, grid, ord = order + 1)

    if (is.null(postmultiply_by)) postmultiply_by <- diag(ncol(grid_basis))

    grid_basis <- grid_basis %*% postmultiply_by

    grid_basis_deriv <- splines::splineDesign(knots, grid, ord = order + 1, derivs = 1)
    grid_basis_deriv <- grid_basis_deriv %*% postmultiply_by

    approx_basis_and_deriv <- function(x) {
        i <- findInterval(x, grid, all.inside = TRUE)
        weight <- (x - grid[i]) / grid_step

        basis <- (1 - weight) * grid_basis[i, , drop = FALSE] + weight * grid_basis[i + 1, , drop = FALSE]
        basis_deriv <- (1 - weight) * grid_basis_deriv[i, , drop = FALSE] + weight * grid_basis_deriv[i + 1, , drop = FALSE]

        # Set basis function evaluations for values outside the support
        # explicitly to zero.
        # Maybe this step is not necessary here because it will get taken care
        # of when evaluating the transformation function.
        basis[weight < 0 | weight > 1, ] <- NA_real_
        basis_deriv[weight < 0 | weight > 1, ] <- NA_real_

        list(basis = basis, basis_deriv = basis_deriv, min_knot = min_knot, max_knot = max_knot)
    }

    approx_basis_and_deriv
}


#' Generates a function to evaluate a B-spline and its derivative with
#' extrapolation.
#'
#' @param basis_and_deriv_fn A function that follows the signature of the
#'      function returned by [get_approx_basis_and_deriv_fun()].
#' @param target_left,target_right Strings, indicating how to extrapolate for
#'      `x` outside the range supported by the spline's knots. Supported values:
#'      * `"identity"`: Extrapolates with the identity function.
#'
#' @returns A function that takes to arguments: `x`, a numeric vector of values
#'      at which to evaluate the spline, and `coef`, a numeric vector of spline
#'      coefficients. `length(coef)` must equal the number of columns of the
#'      basis matrix returned by `basis_and_deriv_fn`. The function returns a
#'      list of two numeric vectors:
#'      * `fx`: Numeric vector of spline evaluations.
#'      * `fx_deriv`: Numeric vector of spline derivative evaluations.
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
get_extrap_bspline_and_deriv <- function(
        basis_and_deriv_fn,
        target_left = c("identity"),
        target_right = c("identity")
    ) {
    target_left <- match.arg(target_left)
    target_right <- match.arg(target_right)

    bdot_and_deriv_fun <- function(x, coef) {
        basis_and_deriv <- basis_and_deriv_fn(x)
        min_knot <- basis_and_deriv$min_knot
        max_knot <- basis_and_deriv$max_knot

        outside_left <- x <= min_knot
        outside_right <- max_knot <= x
        in_range <- !outside_left & !outside_right

        fx <- basis_and_deriv$basis %*% coef
        fx_deriv <- basis_and_deriv$basis_deriv %*% coef

        fx[outside_left] <- x[outside_left]
        fx[outside_right] <- x[outside_right]

        fx_deriv[outside_left] <- 1
        fx_deriv[outside_right] <- 1

        list(fx = fx, fx_deriv = fx_deriv)
    }

    bdot_and_deriv_fun
}


#' Computes the average slope of a B-spline based on its coefficients.
#'
#' Assumes equidistant knots.
#'
#' @param knots Vector of spline knots.
#' @param coef Vector of spline coefficients.
#' @param l Order of the spline. The default is `l = 3` for a cubic spline.
#'
#' @return Average slope of the spline, a scalar.
#' @export
avg_slope_of_bspline <- function(knots, coef, l = 3) {
    n_unique_dists <- knots |>
        diff() |>
        round(5) |>
        unique() |>
        length()
    stopifnot("Knots must be equidistant" = n_unique_dists == 1)

    dk <- knots |>
        diff() |>
        mean()

    p <- length(coef)
    dcoef <- diff(coef)

    numerator <- sum(dcoef[-c(1, 2, p - 2, p - 1)]) +
        sum(dcoef[c(1, p - 1)] / 6) +
        sum(5 * dcoef[c(2, p - 2)] / 6)

    numerator / (dk * (p - l))
}
