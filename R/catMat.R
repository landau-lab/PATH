#' Convert categorical cell state vector into categorical matrix.
#'
#' @param state_vector A vector of cell states from a phylogeny.
#' @param label_order (optional) A vector of cell state names in order.
#' @param n (optional) number of possible cell states. 
#' @param return_sparse (default=TRUE) returns sparse state matrix.
catMat0 <- function(state_vector, label_order=NULL, n=NULL,
                    return_sparse=TRUE) {
  
  if(is.null(n) & is.numeric(state_vector)==FALSE) {
    n <- max(length(unique(state_vector)), length(unique(label_order)))
  } else if(is.numeric(state_vector)==TRUE){
    n <- max(state_vector, n)
  }
  nI <- n
  if(is.null(label_order)) {
    if(is.numeric(state_vector)==FALSE) {
      x <- factor(state_vector, ordered = TRUE)
      nI <- length(levels(x))
    } else {
      if(min(state_vector)==0) {
        min <- 0
        nI <- nI + 1
        
      } else {
        min <- 1
      }
      x <- factor(state_vector, ordered = T, levels = seq(min,n, 1))
    }
  } else {
    x <- factor(state_vector, ordered = TRUE, levels = label_order)
    nI <- length(levels(x))
  }
  z <- Matrix::.sparseDiagonal(nI)[as.numeric(x),]
  tryCatch({colnames(z)[1:length(levels(x))] <- levels(x)},
           error=function(e) {warning("Column names not assigned", call. = F); return(z)})
  if(is.null(names(state_vector))==FALSE) {
    rownames(z) <- make.names(names(x), unique = TRUE)
  }
  if(return_sparse==TRUE) {
    return(z)
  } else {
    return(as.matrix(z))
  }
}


#' Convert categorical cell state vector into categorical matrix.
#'
#' @param cell_states A vector of cell states from a phylogeny.
#' @param num_states (optional) number of possible cell states. 
#' @param state_order (optional) A vector of cell state names in order.
#' @param sparse (default=TRUE) returns sparse state matrix.
#' @export
catMat <- function(cell_states, 
                   num_states=NULL, state_order=NULL,
                   sparse=TRUE, unformatted=FALSE) {
  # Converts a vector of categorical states into a categorical matrix.
  if(unformatted==FALSE) {
    catMat0(state_vector = cell_states, label_order = state_order, n=num_states, return_sparse = sparse)
  } else {
    if(is.null(num_states)) {
      stop("if unformatted=TRUE, num_states cannot be NULL", call. = F)
    }
    diag(num_states)[factor(cell_states),]
  }
}
