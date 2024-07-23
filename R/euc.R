#' Euclidean distance
#' 
#' Computed Euclidean distance between x and y. 
#' 
#' @export
euc <- function(x,y) {
  sqrt(sum((x-y)^2))
}

#' Constrain cell state-specific proliferation rates using total proliferation
#' 
#' @param random_values Guess.
#' @param total_prolif Total proliferation rate. 
#' @param state_freqs Vector cell state frequencies. 
constrain_prolif <- function(random_values, total_prolif, state_freqs) {
  normalize(random_values * state_freqs) * total_prolif / state_freqs
}

#' Make Q matrix from vector of parameters
#' 
#' @param parameters Vector of parameters. 
#' @param dimensions Dimension of square output matrix.
makeQ <- function(parameters, dimensions) {
  parameters <- abs(parameters)
  mat <- matrix(NA, dimensions, dimensions)
  seq1 <- seq(1, (dimensions-1) * dimensions, dimensions - 1)
  seq2 <- seq(dimensions - 1, (dimensions - 1) * dimensions, dimensions - 1)
  for(i in 1:dimensions) {
    mat[i, -i] <- parameters[seq1[i]:seq2[i]]
    mat[i, i] <- -sum(parameters[seq1[i]:seq2[i]])
  }
  return(mat)
}

#' Get weight matrices for d depth
#'
#' @import data.table
getWt <- function(W, d, N) {
  #df <- W[edge == 2*d]
  df <- with(W, W[edge == 2 * d])
  wd <- Matrix::sparseMatrix(df$Var1, df$Var2, dims=c(N, N))
  td <- mean(df$pat, na.rm = T)
  return(list("w"=wd, "t"=td))
}
