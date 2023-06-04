fast_dist.polytomy <- function(tree) {
    
    combn3 <- function(v) {
        if(length(v) > 1) {
            out <- combinat::combn2(v)
        } else {
            out <- t(rep(v, 2))
        }
        return(out)
    }
    
    Ntips <- length(tree$tip.label)
    g <- which(tree$edge[,2] %in% c(1:Ntips))
    g2 <- tree$edge[g,]
    d <- NULL #matrix(NA, length(unique(g2[,1])), 2)
    count <- 1
    for(i in unique(g2[,1])) {
        d <- rbind(as.matrix(combn3((g2[which(g2[,1]==i),2,drop=FALSE]))), d)
    }
    #d <- t(apply(d, 1, sort))
    m <- Matrix::sparseMatrix(rbind(d[,1], d[,2]), rbind(d[,2], d[,1]), dims=c(Ntips, Ntips))
    diag(m) <- FALSE
    md <- mean(tree$edge.length[g])*2
    return(list("W"=m, "mean.edge"=md))
}
