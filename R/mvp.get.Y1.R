#' mvp.get.Y1
#'
#' @param X n-by-m matrix of n m-variate Poisson samples
#' @param Y2 n-by-m(m-1)/2 matrix of latent Poisson covariance samples
#' @param A m-by-m(m+1)/2 matrix from X = AY
#'
#' @return n-by-m matrix of latent Poisson mean samples
#' @export
#'
#' @examples
mvp.get.Y1 <- function(X, Y2, A) {
 return(X - Y2 %*% t(A[, (ncol(X) + 1):ncol(A)]))
}
