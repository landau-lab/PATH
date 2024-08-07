% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PATHpro.R
\name{PATHpro}
\alias{PATHpro}
\title{PATHpro}
\usage{
PATHpro(
  tree,
  cell_states,
  depth = 1,
  prolif = NULL,
  total_prolif = NULL,
  guess_list = NULL,
  alpha = NULL,
  nstates = NULL,
  cell_state_order = NULL
)
}
\arguments{
\item{tree}{Phylogeny.}

\item{cell_states}{A vector of cell states.}

\item{depth}{Max node depth to use for computing phylogenetic correlations.}

\item{prolif}{Vector of cell state-specific proliferation rates, if known (optional).}

\item{total_prolif}{Total proliferation rate, if known (optional).}

\item{guess_list}{List of initial parameter guess vectors (optional).}

\item{alpha}{Regulizer penalty weight.}

\item{nstates}{Number of possible cell states (optional).}

\item{cell_state_order}{Cell state ordering (optional).}
}
\value{
A list containing:\tabular{ll}{
\code{out} \tab Parameter fit; the sum of Euclidean distances between input and inference-based state pair frequencies, plus regulizer penalty.  \cr
\code{Q} \tab Inferred transition rate matrix, \eqn{\mathbf{Q}}. \cr
\code{P} \tab Inferred transition probability matrix, \eqn{\mathbf{P}(t=1) = e^{\mathbf{Q}t}}. \eqn{P_{ij}} is the probability a cell in state \eqn{i} transitions to state \eqn{j} after \eqn{t=1} time, where \eqn{i} is the row, \eqn{j} is the column number, and the units of \eqn{t} are determined by the branch lengths of \code{tree}. \cr
\code{gamma} \tab Inferred cell state-specific proliferation rates (if not provided to \code{prolif}), \eqn{\mathbf{\gamma}}. If \code{total_prolif} is not provided, \code{gamma} is not scaled and represents the differences in state-specific proliferation rates, with the lowest rate set to 0.
}
}
\description{
PATHpro (PATH proliferation) infers cell state transition and proliferation rates from phylogenetic correlations.
}
