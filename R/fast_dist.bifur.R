fast_dist.bifur <- function(tree) {
    Ntips <- length(tree$tip.label)
    g <- which(tree$edge[,2] %in% c(1:Ntips))
    g2 <- tree$edge[g,]
    el <-  tree$edge.length[g]
    d <- matrix(NA, length(unique(g2[,1])), 2)
    dl <- matrix(NA, length(unique(g2[,1])), 2)
    count <- 1
    for(i in unique(g2[,1])) {
        d[count,] <- t(g2[which(g2[,1]==i),2,drop=FALSE])
        el0 <- el[which(g2[,1]==i)]
        if(length(el0)==2) {
            dl[count,] <- el[which(g2[,1]==i)]
        }
        count <- count + 1
    }
    d <- t(apply(d, 1, sort))
    m <- Matrix::sparseMatrix(d[,1],d[,2], symmetric=TRUE )
    #m <- Matrix::forceSymmetric(Matrix::sparseMatrix(d[,1], d[,2], dims = c(Ntips, Ntips)))
    diag(m) <- FALSE
    #md <- mean(tree$edge.length[g])*2
    md <- mean(dl, na.rm=TRUE) * 2
    return(list("W"=m, "mean.pat"=md))
}
