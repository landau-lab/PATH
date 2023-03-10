% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ran_tree_with_states.R
\name{ran_tree_with_states}
\alias{ran_tree_with_states}
\title{Using a birth--death model, simulate a sampled somatic evolutionary process}
\usage{
ran_tree_with_states(Q0, N = 1000, rho = 1, lambda = 1, mu = 0)
}
\arguments{
\item{Q0}{Transition rate matrix.}

\item{N}{Number of cells in the phylogeny.}

\item{rho}{Sampling rate (number of cells in the phylogeny / number of extant cells in the population).}

\item{lambda}{Cell birth rate.}

\item{mu}{Cell death rate.}
}
\description{
This function first uses the castor function \code{generate_tree_hbd_reverse}
to simulate a birth--death phylogeny in reverse. Next, cell state transitions
are simulated as a Markov chain.
}
