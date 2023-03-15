#' Gene Set Enrichment Analysis (GSEA) of genes ranked by phylogenetic auto-correlation.
#' 
#' @param gene.data Matrix of single-cell gene expression data. Each row represents an individual cell and each column a gene. 
#' @param weight.matrix A normalized phylogenetic weight matrix. 
#' @param pathways Pathways to test for enrichment.
#' @export
  xcor_gsea <- function(gene.data, weight.matrix, pathways, maxS=Inf, minS=5, nperm=10000) {
  z <- parM(gene.data, weight.matrix)
  z$gene <- rownames(z)
  
  z <- z %>% as_tibble() %>% group_by(gene) %>% select(Z) %>% arrange(desc(Z)) %$% set_names(Z, gene)
  
  out <- fgsea(pathways, z, maxSize=maxS, minSize=minS, scoreType="std", nPermSimple=nperm)
  return(out)
}
