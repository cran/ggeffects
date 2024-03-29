skip_on_os(c("mac", "solaris"))
skip_if_not_installed("AER")
unloadNamespace("VGAM")

data("Affairs", package = "AER")
m1 <- AER::tobit(affairs ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)

test_that("ggpredict, tobit", {
  pr <- ggpredict(m1, "yearsmarried")
  expect_equal(pr$predicted[1], -10.15089, tolerance = 1e-4)
})

test_that("ggeffect, tobit", {
  skip_if_not_installed("effects")
  expect_null(ggeffect(m1, "yearsmarried"))
})

test_that("ggemmeans, tobit", {
  skip_if_not_installed("emmeans")
  pr <- ggemmeans(m1, "yearsmarried")
  expect_equal(pr$predicted[1], -10.15089, tolerance = 1e-4)
})
