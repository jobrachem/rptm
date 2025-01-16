test_that("create_onion_knots runs", {
    oknots <- onion_knots(-4, 4, p = 10)
    expect_true(!is.null(oknots))
})

test_that("onion_knots has p+7 parameters", {
    oknots <- onion_knots(-4, 4, p = 10)
    expect_equal(oknots$nparam_full_domain, 17)
})

test_that("onion_coef spline has average slope one", {
    oknots <- onion_knots(-4, 4, p = 10)
    for (i in 1:10) {
        li <- rnorm(10)
        coef <- onion_coef(li, oknots) |> cumsum()

        expect_equal(avg_slope_of_bspline(oknots$knots, coef), 1)
    }

})
