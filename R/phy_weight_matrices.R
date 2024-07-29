#' Phylogenetic weight matrix for a node-depth of one. 
#'
#' This function computes the sum normalized phylogenetic weight matrix
#' using the node-depth of one only weighting function. This function
#' returns an N by N dimensional matrix with a value in the ijth position
#' if cell i and j are separated by one node on the phylogeny and a 0 otherwise.
#' @param tree A phylogeny.
#' @param norm Logical. Default TRUE. If true, the phylogenetic weight matrix is first row normalized before it is sum normalized. 
#' @export
one_node_tree_dist <- function(tree, norm = TRUE) {
  w <- one_node_depth(tree)$W
  if(norm == TRUE) {
    w <- rowNorm(w)
  }
  w <- w/sum(w)
  return(w)
}

#' Get the inverse tree distances between all cells
#'
#' This function gets the inverse tree or phylogenetic distances, in terms of node distances or branch lengths, between all cell pairs. Inverse tree distances are normalized to sum to one.  
#' @param tree A phylogeny.
#' @param node Logical. Default TRUE. Tree distances are measured in terms of nodes. If set to FALSE distances in terms of branch lengths.
#' @param norm Logical. Default TRUE. If TRUE, iverse tree distances are row normalized before sum normalizing.
#' @return A matrix of inverse phylogenetic distances between cells. 
#' @export
#' @importFrom castor get_all_pairwise_distances
inv_tree_dist <- function(tree, node = TRUE, norm = TRUE) {
    w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
    if(node == TRUE) {
        w <- w - 1
    }
    w <- 1/w
    diag(w) <- 0
    if(norm == TRUE) {
        w <- rowNorm(w)
    }
    return(w/sum(w))
}


#' Get the (inverse) exponential tree distances between all cells
#'
#' This function gets the (inverse) exponential tree or phylogenetic distances, in terms of node distances or branch lengths, between all cell pairs. Exponential tree distances are normalized to sum to one.  
#' @param tree A phylogeny.
#' @param node Logical. Default TRUE. Tree distances are measured in terms of nodes. If set to FALSE distances in terms of branch lengths.
#' @param norm Logical. Default TRUE. If TRUE, iverse tree distances are row normalized before sum normalizing.
#' @return A matrix of inverse phylogenetic distances between cells. 
#' @export exp_tree_dist
#' @usage exp_tree_dist(tree, node=TRUE, norm=TRUE)
#' @importFrom castor get_all_pairwise_distances
exp_tree_dist <- function(tree, node = TRUE, norm = TRUE) {
  w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
  w <- exp(-w)
  diag(w) <- 0
  if(norm == TRUE) {
    w <- rowNorm(w)
  }
  return(w/sum(w))
}

