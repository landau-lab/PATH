#' Infer cell state transition probabilities from a phylogeny.
#'
#' This function infers cell state transitions from the distribution
#' of states on a phylogeny.
#' @param tree A phylogeny.
#' @param cell_states Location of cell state vector.
#' @param nstates Number of possible cell states.
#' @param impute_branches Logical. Default FALSE. Use and impute phylogenetic branch lengths.
#' @param sample_rate_est Default NULL. Used for branch length imputation if impute_branches=TRUE.
#' @param birth Default 1. Used if impute_branches=TRUE.
#' @param deah Default 0. Used if impute_branches=TRUE.
#' @return A list containing:\tabular{ll}{
#' \code{P} \tab PATH inference of transition probability matrix P at time t=1. \cr
#' \code{Pt} \tab PATH inference of P for time proportional to the mean branch length distance at a node depth 1. \cr
#' }
#' @export
PATH.inference <- function(tree, cell_states = "states",
			   nstates=NULL, impute_branches=FALSE, 
			   sample_rate_est=NULL, birth=1, death=0) {
  t <- SparseM::t

  Ntips <- length(tree$tip.label)
  
  z <- state2mat.sparse(tree[[cell_states]], n=nstates)
  
  w <- fast_dist(tree)
  W <- w$W
  bl <- w$mean.pat

  if(impute_branches==TRUE & is.null(sample_rate_est)==FALSE) {
	  md <- pend_len2(lambda=birth, mu=death, rho=sample_rate_est)*2
  } else {
	  md <- bl
  }

  Freq <- as.matrix(t(z)%*%(W/sum(W))%*%z)
  p0 <- rowSums(Freq)
  
  Pt <- diag(1/p0)%*%Freq
  P <- tryCatch( CorrectP(expm(logm(Pt, method = "Eigen") / (md))),
                 error=function(e) {return(matrix(NA,ncol(z), ncol(z)))} )
  
  
  return(list("P"=P, "Pt"=Pt))
  
}
