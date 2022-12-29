#' Using a birth--death model, simulate a sampled somatic evolutionary process
#'
#' This function first uses the castor function \code{generate_tree_hbd_reverse}
#' to simulate a birth--death phylogeny in reverse. Next, cell state transitions
#' are simulated as a Markov chain.
#' @param Q0 Transition rate matrix.
#' @param N Number of cells in the phylogeny.
#' @param rho Sampling rate (number of cells in the phylogeny / number of extant cells in the population).
#' @param lambda Cell birth rate. 
#' @param mu Cell death rate.
#' @export
ran_tree_with_states <- function(Q0, N=1000, rho=1, lambda=1, mu=0) {
  tr1 <- ran_tree(N=N, rho=rho, lambda=lambda, mu=mu)
  trs <- castor::simulate_mk_model(tr1, Q = Q0, include_nodes = F)
  tr1$states <- trs$tip_states
  return(tr1)
}

