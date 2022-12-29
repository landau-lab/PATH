#' Simulate an idealized phylogeny
#'
#' First, this function uses \code{pbtree} from phytools to generate an idealized phylogeny. 
#' Here, an idealized phylogeny, is a discrete-time, completely sampled, bifurcating, and perfectly balanced phylogeny.
#' Second, cell state transitions are simulated as a Markov chain with transition rate matrix Q.
#' @param g Number of generations. Phylogeny will have 2^g cells. 
#' @param Q Transition rate matrix. 
#' @export
idealized_phylogeny <- function(g, Q) {
  phy <- phytools::pbtree(b = 1, d = 0, 
                n = 2^g, type = "discrete")
  mk <- castor::simulate_mk_model(tree = phy, Q = Q, 
                include_nodes = FALSE)
  phy$states <- mk$tip_states
  return(phy)
}
