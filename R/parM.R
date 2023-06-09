parM <- function(d,w, break.length=100) {
    
  
    y <- c(seq(1,ncol(d), break.length), (ncol(d)+1) )
    y2 <- y[-1] - 1
    y1 <- y[-length(y)]
    
    mfunc <- function(b, e, d1=d, w1=w) {
        M <- xcor(d1[,b:e ], w1)
        x1 <- diag(M$Z)
        x2 <- diag(M$Morans.I)
        out <- data.frame("Z"=x1, "I"=x2)
        return(out)
    }
    
    out <- do.call(rbind,
                   mcmapply(function(b1,e1) mfunc(b=b1,e=e1, d, w), b1=y1, e1=y2,
                            mc.cores = detectCores(), SIMPLIFY = F))
    return(out)
}
