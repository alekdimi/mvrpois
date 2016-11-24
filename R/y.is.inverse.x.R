y.is.inverse.x <- function(x, y) identical(t(t(x)), mvp.matrix(length(x)) %*% y)
