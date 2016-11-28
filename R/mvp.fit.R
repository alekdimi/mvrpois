#' mvp.fit
#'
#' @param X n-by-m matrix of n m-variate Poisson samples
#' @param Z n-by-m matrix of covariates for the mean latent rates
#' @param method "HMC" or "Gibbs"
#'
#' @return results of the sampling method
#' @export
#'
#' @examples
mvp.fit <- function(method, X, Z = NULL, verbose = FALSE, ...) {
  switch(method,
         HMC = {
           model <- ifelse(is.null(Z), '../Stan/karlis.stan', '../Stan/karlis.covariates.stan')
           A <- mvp.matrix(ncol(X))
           data <- list(N = nrow(X), M = nrow(A), K = ncol(A),
                        X_full = X, A2 = A[, (nrow(A) + 1):ncol(A)], Z = Z)
           fit <- stan(file = model, data = data, ...)

           if (verbose) {
             print(fit)
             plot(fit)
             plot(fit, plotfun = "trace")
             pairs(fit, pars = c("theta[1]", "theta[2]", "theta[3]"))
             pairs(fit, pars = c("theta[4]", "theta[5]", "theta[6]"))
           }
           return(fit)
         },
         Gibbs = {
           samples <- gibbs(X, Z, ...)
           return(samples)
         },
         stop('Unknown method. Please use "HMC" or "Gibbs".')
  )
}

gibbs <- function(X, Z, N1, N2) {
  A <- mvp.matrix(ncol(X))
  Y2 <- matrix(0, nrow(X), ncol(A) - nrow(A)) # 2 lines deleted here, if error check back

  if (is.null(Z)) {

    samples <- matrix(NA, N1, ncol(A))
    theta <- rep(1, ncol(A)) # c(3, 4, 5, 2, 8, 5)

    cov


    for (n in 1:N1) {
      cat(paste0(n," "))
      for (i in 1:nrow(X)) {
        ss <- met.hast.disc(function(y) sum(dpois(c(mvp.get.Y1(t(X[i, ]), t(y), A), y), theta, log = T)), Y2[i, ], N2)
        Y2[i, ] <- ss[nrow(ss), ]
      }
      Y1 <- mvp.get.Y1(X, Y2, A)
      Y <- cbind(Y1, Y2)
      theta <- rgamma(length(theta), 1 + colSums(Y), 1 + rep(nrow(Y), ncol(Y)))
      samples[n, ] <- theta
    }
  } else {

    beta <- matrix(0, ncol(Z), ncol(X))
    bdims <- c(ncol(Z), ncol(X))
    theta1 <- exp(Z %*% beta)

    #theta2 <- rep(1, ncol(A) - ncol(X))
    #thetas <- cbind(theta1, matrix(rep(theta2, nrow(thetas)), nrow(theta1), byrow = T))
    theta2 <- mvp.get.Y1(X, theta1, mvp.matrix(ncol(X)))
    thetas <- cbind(theta1, theta2) # TEST THIS TOMORROW!

    samples <- matrix(NA, N, length(beta) + length(theta2))
    for (n in 1:N1) {
      cat(paste0(n," "))

      for (i in 1:nrow(X)) {
        ss <- met.hast.disc(function(y) sum(dpois(c(mvp.get.Y1(t(X[i, ]), t(y), A), y), thetas[i, ], log = T)), Y2[i, ], N2)
        Y2[i, ] <- ss[nrow(ss), ]
      }
      Y1 <- mvp.get.Y1(X, Y2, A)
      Y <- cbind(Y1, Y2)

      for (j in 1:ncol(X)) {
        ss <- mcmc::metrop(function(b) sum(Y1[, j] * (Z %*% b) - exp(Z %*% b)), beta[, j], N2 * 1)$batch
        beta[, j] <- ss[nrow(ss), ]
      }

      theta1 <- exp(Z %*% beta)
      theta2 <- rgamma(length(theta2), 1 + colSums(Y2), 1 + rep(nrow(Y2), ncol(Y2)))
      thetas <- cbind(theta1, matrix(rep(theta2, nrow(thetas)), nrow(theta1), byrow = T))
      samples[n, ] <- c(as.vector(beta), theta2)
    }
  }
  return(samples)
}

############################## PRIVATE HELPER FUNCTIONS ##############################

#' @export
met.hast.disc <- function(p, x, N = 10000) {
  D <- length(x)
  samples <- matrix(NA, N, D)
  for (n in 1:N) {
    xn <- x + ifelse(runif(D) > 0.5, 1, -1) * sample(c(1, 2, 3, 4), 1, prob = c(0.5, 0.25, 0.15, 0.1))
    if (runif(1) < min(1, exp(p(xn) - p(x)))) x <- xn
    samples[n, ] <- x
  }
  return(samples)
}

met.hast.cont <- function(p, N = 10000) {
  samples <- rep(0, N)
  xi <- 1
  for (n in 1:N) {
    xj <- rnorm(1, xi, 1)
    if (runif(1) < min(1, exp(p(xj) - p(xi)))) xi <- xj
    samples[n] <- xi
  }
  return(samples)
}
