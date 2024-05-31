#' Infer transition probability matrix using optimization. 
#'
#' @param n number of cell states
#' @param f State pair frequencies from phylogeny. 
#' @param p Mean patristic distances between cell pairs. 
op_Pinf <- function(n, f, p) {
  loss_func <- function(pars, dims, freq, pat_dist) {
    Qguess <- makeQ(pars, dims)
    D <- diag(rowSums(freq))
    Pguess <- expm::expm(Qguess*(pat_dist/2))
    X <- t(Pguess)%*%D%*%Pguess
    out <- euc(freq, X)
    return(list("out"=out, "pars"=pars))
  }
  op_search <- nlminb(runif(n^2, 0, 1), 
                      function(guess_pars) loss_func(pars=guess_pars, 
                                                     dims=n,
                                                     freq=f, pat_dist=p)$out, 
                      control = list(abs.tol=1e-2, rel.tol=1e-2, x.tol=1e-2))
  inf <- expm::expm(makeQ(op_search$par, n))
  return(inf)
}

#' PATH inference for given patristic distances, cell state, and weight matrices. 
#'
#' @param Z Cell state matrix. 
#' @param W Phylogenetic weight matrix for specific depth. 
#' @param mspd Mean patristic distances between cells. 
P_inf.zw <- function(Z, W, mspd) {
 
  t <- Matrix::t
  W <- rowNorm(W)
  W <- W/sum(W)
  Freq <- as.matrix(t(Z)%*%W%*%Z)
  u <- rowSums(Freq)
  
  Pt <- diag(1/u)%*%Freq
  rownames(Pt) <- colnames(Pt)
  
  Pinf <- tryCatch(CorrectP(expm::expm(expm::logm(Pt, method="Eigen")/mspd)),
                   error=function(e){op_Pinf(ncol(Freq), Freq, mspd)})
  
  return(list("Pt"=Pt, "P"=Pinf, "Q"=expm::logm(Pinf)))
} 

#' PATH inference given cell state vector, phylogenetic weight matrix, and mean patristic distances
#'
#' @param state_vector Vector of cell states. 
#' @param weight_matrix Phylogenetic weight matrix. 
#' @param mspd Mean patristic distances between cells represented in weight_matrix. 
P_inf <- function(state_vector, weight_matrix, mspd) {
  
  Z <- catMat(state_vector)
  W <- rowNorm(weight_matrix)
  W <- W/sum(W)
  
  P_inf.zw(Z, W, mspd)
} 

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
#' \code{Q} \tab PATH inference of transition rate matrix Q = logm(P). 
#' }
#' @export

PATH.inference <- function(tree, cell_states = "states", 
                      nstates = NULL, impute_branches = FALSE, sample_rate_est = NULL, 
                      birth = 1, death = 0) {
  w <- one_node_depth(tree)
  W <- w$W
  bl <- w$mean.pat
  
  if (impute_branches == TRUE & is.null(sample_rate_est) == FALSE) {
    mspd <- pend_len2(lambda = birth, mu = death, rho = sample_rate_est) * 2
  }
  else {
    mspd <- bl
  }
  
  P_inf(state_vector = tree[[cell_states]], weight_matrix = W, mspd = mspd)
  
}
