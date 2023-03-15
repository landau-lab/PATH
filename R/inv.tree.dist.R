#' Get the inverse tree distances between all cells
#'
#' This function gets the inverse tree or phylogenetic distances, in terms of node distances or branch lengths, between all cell pairs. Inverse tree distances are normalized to sum to one.  
#' @param tree A phylogeny.
#' @param node Logical. Default TRUE. Tree distances are measured in terms of nodes. If set to FALSE distances in terms of branch lengths.
#' @param norm. Logical. Default TRUE. If TRUE, iverse tree distances are row normalized before sum normalizing.
#' @return A matrix of inverse phylogenetic distances between cells. 
#' @export
inv.tree.dist <- function(tree, node=TRUE, norm=TRUE) {
    w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
    if(node == TRUE) {
        w <- w - 1
    }
    w <- 1/w
    diag(w) <- 0
    if(norm==TRUE) {
        w <- rowNorm(w)
    }
    return(w/sum(w))
}
