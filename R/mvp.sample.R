#' mvp.sample
#'
#' @param n           number of samples to be drawn
#' @param theta.means underlying Poisson rates for the mean
#' @param theta.covs  underlying Poisson rates for the covariance
#'
#' @return            n-by-m matrix of n m-variate Poisson samples, where m is of size y1
#' @export
#'
#' @examples
mvp.sample <- function(n, theta.means, theta.covs) {
  m <- length(theta.means)
  k <- m + length(theta.covs)
  A <- mvp.matrix(m)
  return(t(A %*% matrix(rpois(n * k, c(theta.means, theta.covs)), nrow = k)))
}
