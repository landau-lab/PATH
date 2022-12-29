#' Generate random transition matrix
#'
#' \code{rQ} generates a random transition rate matrix $Q$.
#' \code{rP} generates a random transition probability matrix $P$.
#'
#' @param n Number of distinct states (dimension of transition probability matrix)
#' @param max.r Parameter used for randomly generating off-diagonal entries of $Q$. Entries chosen from uniform distribution between 0 and max.r.
#' @export
rP <- function(n, max.r=0.1) {
	expm::expm(rQ(n, max.r))
	}
