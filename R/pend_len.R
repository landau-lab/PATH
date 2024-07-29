pend_len <- function(gamma, xi) {
  # gamma := (birth - death), or net growth rate
  # xi := birth * sampling
  
  (xi - gamma + gamma * log(gamma/xi)) / (gamma - xi)^2
}
