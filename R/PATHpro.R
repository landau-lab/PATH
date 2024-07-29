#' PATHpro loss function for optimization
#'
#' @param guess Vector of parameter value guesses.
#' @param pair_freqs List of cell state pair frequencies from phylogeny.
#' @param times Vector of times/patristic distances corresponding to pair_freqs.
#' @param state_freqs Vector of cell state frequencies.
#' @param total_prolif Total proliferation rate.
#' @param prolif_rates Vector of cell state proliferation rates. 
#' @param alpha Regulizer penalty weight to use when comparing pair_freqs with those using guessed parameters. 
#' @return out Sum of Euclidean distances between guessed and input state pair frequencies plus regulizer penalty. 
#' @return P Transition probability matrix from guessed parameters. 
#' @return Q Transition rate matrix from guessed parameters. 
#' @return gamma Proliferation rates from guessed parameters or prolif_rates if not null.
#' @return A list containing:\tabular{ll}{
#' \code{out} \tab Sum of Euclidean distances between guessed and input state pair frequencies plus regulizer penalty.  \cr
#' \code{Q} \tab Transition rate matrix from guessed parameters. \cr
#' \code{P} \tab Transition probability matrix from guessed parameters. \cr
#' \code{gamma} \tab Proliferation rates from guessed parameters or prolif_rates if not null.
#' }
#' @importFrom expm expm

PATHpro_loss <- function(guess, pair_freqs, times, state_freqs, total_prolif, prolif_rates, alpha = 0) {
  n <- length(state_freqs)
  Q <- makeQ(parameters = guess, dimensions = n)
  seq1 <- seq(1, (n-1)*n, n-1)
  seq2 <- seq(n-1, (n-1)*n, n-1)
  if(is.null(prolif_rates)) {
    r <- constrain_prolif(abs(guess[(rev(seq2)[1]+1):(rev(seq2)[1]+n)]), 
                          total_prolif, state_freqs)
    if(is.null(total_prolif)) {
      r <- abs(guess[(rev(seq2)[1]+1):(rev(seq2)[1]+n)])
      r <- r - min(r)
    }
  } else {
    r <- prolif_rates
  }
  Gamma <- diag(r)
  lH <- Gamma + Q
  out <- 0
  for(i in 1:length(pair_freqs)) {
    ti <- times[i]/2
    nu <- c(expm::expm(lH*ti)%*%matrix(1, n))
    A <- diag(normalize(state_freqs*nu))
    IP <- rowNorm(expm::expm(lH*ti))
    out <- out + euc(pair_freqs[[i]], t(IP)%*%A%*%IP)
  }
  out <- out + (alpha * sum((Q[upper.tri(Q)|lower.tri(Q)])^2))
  P <- expm::expm(Q)
  return(list("out" = out, "Q" = Q, "P" = P, "gamma" = r))
}


#' Fit PATHpro by optimizing `PATHpro_loss()`
#'
#' @param initial_parameter_guess Vector of paameter value guesses.
#' @param pair_freq_list List of cell state pair frequencies from phylogeny. 
#' @param state_freqs Vector of cell state frequencies. 
#' @param pair_pat_dist Vector of times/patristic distances corresponding to pair_freq_list.
#' @param total_prolif_rate Total proliferation rate. 
#' @param cell_prolif_rates Vector of cell state proliferation rates. 
#' @param alph Regulizer penalty weight to use when comparing pair_freq_list with those using guessed parameters.
#' @param xtol,rtol,atol Tolerances used for optimization using \code{nlminb()}. 
#' @return A list containing:\tabular{ll}{
#' \code{out} \tab Minimal sum of Euclidean distances and regulizer penalty for input and inference-based state pair frequencies.  \cr
#' \code{Q} \tab Inferred transition rate matrix. \cr
#' \code{P} \tab Inferred transition probability matrix. \cr
#' \code{gamma} \tab Inferred cell state specific proliferation rates (if not provided).
#' }
#' @importFrom stats nlminb

fitPATHpro <- function(initial_parameter_guess, pair_freq_list, state_freqs, 
                       pair_pat_dist = 2, total_prolif_rate = NULL, 
                       cell_prolif_rates = NULL, alph = 0,
                       xtol = 1e-15, rtol = 1e-15, atol = 1e-15) {
  
  op_search <- stats::nlminb(initial_parameter_guess,
                      function(g) PATHpro_loss(guess = g, pair_freqs = pair_freq_list, 
                                               times = pair_pat_dist, 
                                               state_freqs = state_freqs, 
                                               total_prolif = total_prolif_rate, 
                                               prolif_rates = cell_prolif_rates, 
                                               alpha = alph)$out,
                       
                       control = list(x.tol = xtol, rel.tol = rtol, 
                                      abs.tol = atol, 
                                      eval.max = 1000, iter.max = 1000, 
                                      step.min = 1e-3, step.max = 0.1))
  
  PATHpro_loss(guess = op_search$par, pair_freqs = pair_freq_list, 
               times = pair_pat_dist, 
               state_freqs = state_freqs, total_prolif = total_prolif_rate, 
               prolif_rates = cell_prolif_rates, 
               alpha = alph)
}

#' PATHpro inference using input cell state and weight matrices. 
#'
#' @param X Cell state matrix. 
#' @param W Phylogenetic weight matrix, or list of weight matrices.
#' @param t Mean pairwise patristic distances corresponding to \code{W}.
#' @param guess_list List of initial parameter guesses (optional).
#' @param total_prolif_rate_est Estimate of total proliferation rate (optional). 
#' @param cell_prolifs Vector of cell state-specific prolfieration rates, if known (optional). 
#' @param alpha0 Regulizer penalty weight. 
#' @return A list containing:\tabular{ll}{
#' \code{out} \tab Minimal sum of Euclidean distances and regulizer penalty for input and inference-based state pair frequencies.  \cr
#' \code{Q} \tab Inferred transition rate matrix, \eqn{\mathbf{Q}}. \cr
#' \code{P} \tab Inferred transition probability matrix, \eqn{\mathbf{P} = e^{\mathbf{Q}}}. \cr
#' \code{gamma} \tab Inferred cell state-specific proliferation rates (if not provided to \code{cell_prolifs}), \eqn{\mathbf{\gamma}}.
#' }
#' @export
#' @details \eqn{\sum_{1}^{\text{depth}} \lVert \mathbb{P}^\top \mathbf{A} \mathbb{P} - \mathbf{X}^\top \overline{\mathbf{W}} \mathbf{X} \rVert + \alpha \sum_{i \neq j} Q_{ij}^2}
#'
#' @importFrom Matrix colMeans
#' @importFrom stats runif
PATHpro.XW <- function(X, W, t, guess_list = NULL, total_prolif_rate_est = NULL, 
                       cell_prolifs = NULL, alpha0 = 0) {
  if(!is.list(W) & length(t) == 1) {
    W <- list(W)
  }
  Freqs <- list()
  for(i in 1:length(W)) {
    W[[i]] <- rowNorm(W[[i]])
    W[[i]] <- W[[i]]/sum(W[[i]])
    Freqs[[i]] <- as.matrix(Matrix::t(X)%*%W[[i]]%*%X)
  }
  if(is.null(guess_list)) {
    n <- ncol(X)
    guess_list <- list(stats::runif(n^2, 0, 0.1))
  }
  infs <- list()
  for(i in 1:length(guess_list)) {
    infs[[i]] <- fitPATHpro(initial_parameter_guess = guess_list[[i]],
                            pair_freq_list = Freqs, state_freqs = Matrix::colMeans(X), 
                            pair_pat_dist = t, 
                            total_prolif_rate = total_prolif_rate_est, 
                            cell_prolif_rates = cell_prolifs, 
                            alph = alpha0)
  }
  mg <- which.min(sapply(infs, '[[', "out"))
  inf0 <- infs[[mg]]
  #inf0$guess <- mg
  return(inf0)
}


#' PATHpro
#'
#' PATHpro (PATH proliferation) infers cell state transition and proliferation rates from phylogenetic correlations.
#'
#' @param tree Phylogeny.
#' @param cell_states A vector of cell states. 
#' @param depth Max node depth to use for computing phylogenetic correlations. 
#' @param prolif Vector of cell state-specific proliferation rates, if known (optional). 
#' @param total_prolif Total proliferation rate, if known (optional). 
#' @param guess_list List of initial parameter guess vectors (optional). 
#' @param alpha Regulizer penalty weight.
#' @param nstates Number of possible cell states (optional).
#' @param cell_state_order Cell state ordering (optional).
#' @return A list containing:\tabular{ll}{
#' \code{out} \tab Parameter fit; the sum of Euclidean distances between input and inference-based state pair frequencies, plus regulizer penalty.  \cr
#' \code{Q} \tab Inferred transition rate matrix, \eqn{\mathbf{Q}}. \cr
#' \code{P} \tab Inferred transition probability matrix, \eqn{\mathbf{P}(t=1) = e^{\mathbf{Q}t}}. \eqn{P_{ij}} is the probability a cell in state \eqn{i} transitions to state \eqn{j} after \eqn{t=1} time, where \eqn{i} is the row, \eqn{j} is the column number, and the units of \eqn{t} are determined by the branch lengths of \code{tree}. \cr
#' \code{gamma} \tab Inferred cell state-specific proliferation rates (if not provided to \code{prolif}), \eqn{\mathbf{\gamma}}. If \code{total_prolif} is not provided, \code{gamma} is not scaled and represents the differences in state-specific proliferation rates, with the lowest rate set to 0.
#' }
#' @export
#' @importFrom castor get_all_pairwise_distances
#' @importFrom data.table data.table
#'
PATHpro <- function(tree, cell_states, depth = 1, prolif = NULL, total_prolif = NULL,
                    guess_list = NULL, alpha = NULL, nstates = NULL, cell_state_order = NULL) {
	X <- catMat(cell_states, num_states = nstates, state_order = cell_state_order)
	N <- nrow(X)
  	if(is.null(alpha)) {
		n <- ncol(X)
  		alpha <- depth/(n^2 - n)
	}

	if(depth > 1) {
		wa <- castor::get_all_pairwise_distances(tree, only_clades = 1:N, 
                                             as_edge_counts = T, check_input = F)
		wp <- castor::get_all_pairwise_distances(tree, only_clades = 1:N,
                                             as_edge_counts = F, check_input = F)
		W <- data.table::data.table("Var1" = rep(1:N, times = N), "Var2" = rep(1:N, each = N), 
                                "edge" = c(wa), "pat" = c(wp))
		out <- lapply(seq(1, depth, 1), function(i) getWt(W, i, N))
    		wl <- lapply(out, '[[', "w")
    		t <- sapply(out, '[[', "t")} 
	else if (depth == 1) {
      		W <- one_node_depth(tree)
      		wl <- list(W$W)
      		t <- W$mean.pat
	}
	inf <- PATHpro.XW(X, wl, t, cell_prolifs = prolif, 
			  total_prolif_rate_est = total_prolif, 
                    	  guess_list = guess_list, 
                     	  alpha0 = alpha)
	return(inf)
}

