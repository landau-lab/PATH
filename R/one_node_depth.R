#' Get weight matrix of cells that share a most recent common ancestor at a node depth of 1. 
#'
#' @param tree A phylogeny.
#' @return A list containing:\tabular{ll}{
#' \code{W} \tab weight matrix. \cr
#' \code{mean.pat} \tab mean patristic distance separating cell pairs. \cr
#' }
#' @export
#' @importFrom magrittr %>% %$%
one_node_depth <- function(tree) {
  
  bifur <- castor::is_bifurcating(tree)
  Ntips <- length(tree$tip.label)
  tipnodes <- which(tree$edge[,2] %in% 1:Ntips)
  edges <- tree$edge[tipnodes,]
  df <- data.frame("node"=edges[,1], "tip"=edges[,2], 
                   "pat"=tree$edge.length[tipnodes])
  df2 <- df %>% dplyr::group_by(node) %>% 
    dplyr::mutate("clades"=list(tip)) %>% 
    dplyr::rowwise() %>% dplyr::filter(length(clades) >= 2)
  indices <- df2 %>% dplyr::select(-tip, -pat) %>% dplyr::distinct() %>%
    dplyr::pull(clades) %>% 
    {if(!bifur)
      lapply(., function(x) t(utils::combn(x, 2))) %>% do.call(rbind, .)
      else do.call(rbind, .) }
  mat <- Matrix::sparseMatrix(rbind(indices[,1], indices[,2]), 
                      rbind(indices[,2], indices[,1]), 
                      dims = c(Ntips, Ntips))
  mp <- df2 %$% pat %>% mean * 2
  return(list("W"=mat, "mean.pat"=mp, "pats"=df2$pat))
}
