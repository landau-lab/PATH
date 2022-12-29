pend_len2 <- function(lambda, mu, rho) {
  mu2 <- mu - lambda*(1-rho)
  lambda2 <- rho*lambda
  
  (mu2 + (lambda2 - mu2)*log(1- mu2/lambda2))/mu2^2
}
