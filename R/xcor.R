#' Compute cross-correlations.
#'
#' This function computes both auto- and cross-correlations for
#' arguments: `data`, an \eqn{N \times n} dimensional matrix, and
#' `weight.matrix`, an \eqn{N \times N} dimensional matrix.
#' In the context of computing phylogenetic correlations,
#' \eqn{N} is the number of cells (terminal phylogenetic nodes or leaves),
#' and \eqn{n} is the number of observations (e.g., RNA counts for \eqn{n} genes).
#'
#' @param data \eqn{N \times n} matrix of trait/state data, where \eqn{N} is number of cells (phylogeny leaves) and \eqn{n} is number of traits.
#' @param weight.matrix \eqn{N \times N} phylogenetic weight matrix. 
#' @examples
#' phy <- ran_tree_with_states(matrix(c(-0.1, 0.1, 0.1, -0.1), 2, 2))
#' X <- catMat(phy$states)
#' W <- one_node_tree_dist(phy)
#' xcor(X, W)
#' @return A list containing:\tabular{ll}{
#'  \code{phy_cor} \tab Phylogenetic cross-correlations. Diagonal elements are phylogenetic auto-correlations. \cr
#'  \code{Z.score} \tab A matrix of analytical z scores for cross-correlations. \cr
#'  \code{one.sided.pvalue} \tab One-sided p-values derived from z scores. \cr
#' }
#' @export
#' @importFrom methods as
#' @importFrom Matrix t
#' @importFrom qlcMatrix corSparse
#' @importFrom stats pnorm
xcor <- function(data, weight.matrix) {
  # data format: one row per cell and one column per measurement (e.g. one column per gene expression score)
  # The ij'th element of the weight.matrix should reflect how closely related cells i and j are, 
  # with larger numbers corresponding to closer relationships. The diagonal elements should be 0.

  t <- Matrix::t

  W <- methods::as(weight.matrix, "sparseMatrix")
  n <- nrow(data)
  
  
  d0 <- apply(as.matrix(data), 2, function(xx) xx - mean(xx))
  
  d1 <- (t(d0)%*%d0)/n
  d1.2 <- diag(d1)%*%t(diag(d1))
  d2 <- (t(d0^2)%*%(d0^2))/n
  
  x <- rep(1, n)
  w <- as.numeric(t(x)%*%W%*%x)
  S3 <- as.numeric(t(x)%*%(W * t(W))%*%x)
  S4 <- as.numeric(t(x)%*%(W * W)%*%x)
  S5 <- as.numeric(t(x)%*%W%*%W%*%x)
  S6 <- as.numeric(t(W%*%x)%*%(W%*%x) + t(t(W)%*%x)%*%(t(W)%*%x))
  S1 <- S4 + S3
  S2 <- 2*S5 + S6
  
  I <- as.matrix((t(d0)%*%W%*%d0)/(w*sqrt(d1.2)))
  
  E.I2 <- (-qlcMatrix::corSparse(data)/(n-1))
  
  V00 <- as.matrix((  ((d1^2 * n)/ (d1.2) ) * (2*(w^2 - S2 + S1) + (2*S3 - 2*S5)*(n-3) +
					       S3*(n-2)*(n-3) )  +  (( -d2 )/( d1.2 ))*(6*(w^2 - S2 + S1) +
				       (4*S1-2*S2)*(n-3) + S1*(n-2)*(n-3)  ) 
                      + n*((w^2 - S2 + S1) + (2*S4 -S6)*(n-3) + S4*(n-2)*(n-3) )  ) /
  ((n-1)*(n-2)*(n-3)*(w^2) ) - 
  (  ( ( d1^2 )/( d1.2 ) ) * ( 1/((n-1)^2)  )  ))
  
  Z0 <- (I-E.I2)/sqrt(V00)
  pval <- stats::pnorm(-abs(Z0))*2
  
  return(list("phy_cor" = I, "Z.score" = Z0, "one.sided.pvalue" = pval))
  
}
