if (require("testthat") &&
    require("ggeffects") &&
    require("emmeans") &&
    require("gee")) {
  data(warpbreaks)
  m1 <- suppressMessages(gee(breaks ~ tension, id = wool, data = warpbreaks, silent = TRUE))

  test_that("ggpredict", {
    p <- ggpredict(m1, "tension")
    expect_equal(p$predicted[1], 36.38889, tolerance = 1e-3)
  })

  test_that("ggemmeans", {
    p <- ggemmeans(m1, "tension")
    expect_equal(p$predicted[1], 36.38889, tolerance = 1e-3)
  })
}
