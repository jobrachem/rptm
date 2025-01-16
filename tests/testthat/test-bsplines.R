test_that("equidistant_knots create equidistant knots", {
    x <- rnorm(100)
    knots <- equidistant_knots(x)
    expect_equal(stats::sd(diff(knots)), 0)
})

test_that("equidistant_knots leads to expected number of knots", {
    x <- rnorm(100)

    for (p in 6:15) {
        for (l in 0:5) {
            knots <- equidistant_knots(x, p = p, l = l)
            expect_equal(length(knots), p + l + 1)
        }
    }
})

test_that("equidistant_knots leads to expected number of parameters", {
    x <- rnorm(100)
    for (p in 6:15) {
        for (l in 0:5) {
            knots <- equidistant_knots(x, p = p, l = l)
            B <- splines::splineDesign(knots, x, ord = l + 1)
            expect_equal(ncol(B), p)
        }
    }
})


test_that("get_approx_basis_and_deriv_fun approximate bases successfully", {
    knots <- equidistant_knots(c(-4, 4), p = 10)

    approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)

    x <- seq(-4, 4, length.out = 10)
    basis_and_deriv <- approx_basis_and_deriv(x)
    basis_computed <- splines::splineDesign(knots, x)
    basis_deriv_computed <- splines::splineDesign(knots, x, derivs = 1)

    expect_equal(basis_and_deriv$basis, basis_computed, tolerance=1e-4)
    expect_equal(basis_and_deriv$basis_deriv, basis_deriv_computed, tolerance=1e-4)

})

test_that("get_approx_basis_and_deriv_fun approximate bases successfully with x outside the supported range.", {
    knots <- equidistant_knots(c(-4, 4), p = 10)
    min_knot <- knots[4]
    max_knot <- rev(knots)[4]

    approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)

    x <- seq(-10, 10, length.out = 100)
    basis_and_deriv <- approx_basis_and_deriv(x)

    outside_left <- x <= min_knot
    outside_right <- max_knot <= x
    in_range <- !outside_left & !outside_right

    expect_true(all(is.na(basis_and_deriv$basis[outside_left | outside_right])))
    expect_true(all(is.na(basis_and_deriv$basis_deriv[outside_left | outside_right])))
    expect_true(!any(is.na(basis_and_deriv$basis[in_range])))
    expect_true(!any(is.na(basis_and_deriv$basis_deriv[in_range])))

})

test_that("get_extrap_bspline_and_deriv runs as expected with extrapolation on ordered x", {
    knots <- equidistant_knots(c(-4, 4), p = 10)
    min_knot <- knots[4]
    max_knot <- rev(knots)[4]

    approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)
    bdot_and_deriv_fun <- get_extrap_bspline_and_deriv(approx_basis_and_deriv)

    x <- seq(-10, 10, length.out = 100)
    coef <- rnorm(10)

    res <- bdot_and_deriv_fun(x, coef)

    outside_left <- x <= min_knot
    outside_right <- max_knot <= x
    in_range <- !outside_left & !outside_right

    expect_equal(res$fx[outside_left], x[outside_left])
    expect_equal(res$fx[outside_right], x[outside_right])

    expect_equal(res$fx_deriv[outside_left], rep(1, times = sum(outside_left)))
    expect_equal(res$fx_deriv[outside_right], rep(1, times = sum(outside_right)))
})

test_that("get_extrap_bspline_and_deriv runs as expected with extrapolation on unordered x", {
    knots <- equidistant_knots(c(-4, 4), p = 10)
    min_knot <- knots[4]
    max_knot <- rev(knots)[4]

    approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)
    bdot_and_deriv_fun <- get_extrap_bspline_and_deriv(approx_basis_and_deriv)

    x <- seq(-10, 10, length.out = 100)
    x <- sample(x)
    x <- -10
    coef <- rnorm(10)

    res <- bdot_and_deriv_fun(x, coef)

    outside_left <- x <= min_knot
    outside_right <- max_knot <= x
    in_range <- !outside_left & !outside_right

    expect_equal(res$fx[outside_left], x[outside_left])
    expect_equal(res$fx[outside_right], x[outside_right])

    expect_equal(res$fx_deriv[outside_left], rep(1, times = sum(outside_left)))
    expect_equal(res$fx_deriv[outside_right], rep(1, times = sum(outside_right)))
})

test_that("get_extrap_bspline_and_deriv run without error without extrapolation", {
    knots <- equidistant_knots(c(-4, 4), p = 10)

    approx_basis_and_deriv <- get_approx_basis_and_deriv_fun(knots, ngrid = 1000)
    bdot_and_deriv_fun <- get_extrap_bspline_and_deriv(approx_basis_and_deriv)

    x <- seq(-4, 4, length.out = 100)
    coef <- rnorm(10)

    expect_no_error({bdot_and_deriv_fun(x, coef)})
})
