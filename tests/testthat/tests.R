context("mvp.matrix")

test_that("mvp works correctly", {
  expect_equal(mvp.matrix(3), rbind(c(1, 0, 0, 1, 1, 0),
                                    c(0, 1, 0, 1, 0, 1),
                                    c(0, 0, 1, 0, 1, 1)))
})



