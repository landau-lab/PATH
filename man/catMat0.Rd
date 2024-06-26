% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catMat.R
\name{catMat0}
\alias{catMat0}
\title{Convert categorical cell state vector into categorical matrix.}
\usage{
catMat0(state_vector, label_order = NULL, n = NULL, return_sparse = TRUE)
}
\arguments{
\item{state_vector}{A vector of cell states from a phylogeny.}

\item{label_order}{(optional) A vector of cell state names in order.}

\item{n}{(optional) number of possible cell states.}

\item{return_sparse}{(default=TRUE) returns sparse state matrix.}
}
\description{
Convert categorical cell state vector into categorical matrix.
}
