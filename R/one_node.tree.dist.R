#' Phylogenetic weight matrix for a node-depth of one. 
#'
#' This function computes the sum normalized phylogenetic weight matrix
#' using the node-depth of one only weighting function. This function
#' returns an N by N dimensional matrix with a value in the ijth position
#' if cell i and j are separated by one node on the phylogeny and a 0 otherwise.
#' @param tree A phylogeny.
#' @param norm Logical. Default TRUE. If true, the phylogenetic weight matrix is first row normalized before it is sum normalized. 
#' @export
one_node.tree.dist <- function(tree, norm=TRUE) {
  w <- one_node_depth(tree)$W
  if(norm==TRUE) {
    w <- rowNorm(w)
  }
  w <- w/sum(w)
  return(w)
}
