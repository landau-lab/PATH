rowNorm <- function(W) {
  s <- Matrix::rowSums(W)
  W/ifelse(s==0, 1, s)
}
