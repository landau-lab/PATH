#' Gene Set Enrichment Analysis (GSEA) of genes ranked by phylogenetic auto-correlation.
#' 
#' @param gene.data Matrix of single-cell gene expression data. Each row represents an individual cell and each column a gene. 
#' @param weight.matrix A normalized phylogenetic weight matrix. 
#' @param pathways Pathways to test for enrichment.
#' @param maxS Max pathway size.
#' @param minS Min pathway size.
#' @param nperm Number of permutations. 
#' @param ncores Number of cores to use. 
#' @return GSEA results using \code{fgsea::fgsea}.
#' @export
#' @importFrom dplyr as_tibble group_by select arrange desc
#' @import magrittr
  xcor_gsea <- function(gene.data, weight.matrix, pathways, maxS = Inf, minS = 5, nperm = 10000, ncores = 2) {
  z <- parM(gene.data, weight.matrix, core_num = ncores)
  z$gene <- rownames(z)
  
  z <- z %>% as_tibble() %>% group_by(gene) %>% select(Z) %>% arrange(desc(Z)) %$% set_names(Z, gene)
  
  out <- fgsea::fgsea(pathways, z, maxSize = maxS, minSize = minS, scoreType = "std", nPermSimple = nperm)
  return(out)
}
