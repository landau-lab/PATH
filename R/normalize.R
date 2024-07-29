#' Normalize vector to sum to 1.
#' 
#' @param x input vector with all positive elements. 
#' @return normalized vector. 
normalize <- function(x) {
  if( length(which(x < 0)) > 0)
  {
    stop("\ninput contains negative elements.")
  }
  x/sum(x)
}
