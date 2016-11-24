.onLoad <- function(libname = find.package("mvrpois"), pkgname = 'mvrpois') {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}
