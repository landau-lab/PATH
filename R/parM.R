#' Calculate heritability in parallel
#'
#' @param d Gene data matrix.
#' @param w Phylogenetic weight matrix.
#' @param break.length Number of genes to measure at a time.
#' @param core_num Number of cores for parallel processing. 
#' @importFrom parallel mcmapply detectCores
parM <- function(d, w, break.length = 100, core_num = 2) {
    y <- c(seq(1, ncol(d), break.length), (ncol(d) + 1))
    y2 <- y[-1] - 1
    y1 <- y[-length(y)]
    mfunc <- function(b, e, d1 = d, w1 = w) {
        M <- xcor(d1[,b:e ], w1)
        x1 <- diag(M$Z)
        x2 <- diag(M$phy_cor)
        out <- data.frame("Z" = x1, "I" = x2)
        return(out)
    }
    out <- do.call(rbind,
                   parallel::mcmapply(function(b1,e1) mfunc(b = b1, e = e1, d, w), b1 = y1, e1 = y2,
                            mc.cores = core_num, SIMPLIFY = F))
    return(out)
}
