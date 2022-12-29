fast_dist <- function(tree) {
  if(castor::is_bifurcating(tree)==TRUE) {
    out <- fast_dist.bifur(tree)
  } else {
    out <- fast_dist.polytomy(tree)
  }
  return(out)
}
