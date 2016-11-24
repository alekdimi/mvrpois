#' mvp.matrix
#'
#' @param m dimension of the m-variate Poisson
#'
#' @return  the MVP matrix A from X = AY (Karlis 2005), where X is MVP-distributed, and Y is a vector of independent Poissons.
#' @export
#'
#' @examples
mvp.matrix <- function(m) {
  if (m < 2) {
    stop("Multivariate Poisson is only defined for m > 2. Are you looking for pois, the univariate Poisson distribution?")
  } else {
    combs <- cbind(rbind(1:m, 1:m), combn(1:m, 2))
    k <- ncol(combs)
    A <- matrix(0, m, k)

    for (i in 1:m)
      for (j in 1:k)
        if (combs[1, j] == i | combs[2, j] == i)
          A[i, j] = 1

    return(A)
  }
}
