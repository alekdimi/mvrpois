context("mvp.matrix")

test_that("mvp.matrix", {
  expect_equal(mvp.matrix(3), rbind(c(1, 0, 0, 1, 1, 0),
                                    c(0, 1, 0, 1, 0, 1),
                                    c(0, 0, 1, 0, 1, 1)))
})

test_that("mvp.prob", {
  expect_equal(mvp.prob(c(5, 6, 7), c(1, 2, 3, 4, 5, 6), logarithm = FALSE, method = 'analytical'),
               mvp.prob(c(5, 6, 7), c(1, 2, 3, 4, 5, 6), logarithm = FALSE, method = 'recursive'))

  expect_equal(mvp.prob(c(5, 6, 7), c(1, 2, 3, 4, 5, 6), logarithm = TRUE, method = 'analytical'),
               mvp.prob(c(5, 6, 7), c(1, 2, 3, 4, 5, 6), logarithm = TRUE, method = 'recursive'))

})
