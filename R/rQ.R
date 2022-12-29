#' Generate random transition matrix
#'
#' \code{rQ} generates a random transition rate matrix $Q$.
#' \code{rP} generates a random transition probability matrix $P$.
#'
#' @param n Number of distinct states (dimension of transition probability matrix)
#' @param max.r Parameter used for randomly generating off-diagonal entries of $Q$. Entries chosen from uniform distribution between 0 and max.r.
#' @export
rQ <- function(n, max.r=1) {
    mat <- matrix(NA, n,n) 
    for(i in 1:n) {
        y <- runif(n-1, 0, max.r)
        mat[i,-i] <- y
        mat[i,i] <- -sum(y)
    }
    return(mat)
}
