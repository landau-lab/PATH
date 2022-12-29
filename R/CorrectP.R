#' CorrectP
#'
#' This function normalizes rows of the transition probability matrix P to sum to 1.
#'
CorrectP <- function(P) {
  Pc <- t(apply(P, 1, function(x) ifelse(x<0, 0, x) ))
  Pc <- t(apply(Pc, 1, function(x) x/sum(x)))
  return(Pc)
}
