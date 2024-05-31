#' PATHpro loss function for optimization
#'
#' @param guess Vector of parameter value guesses.
#' @param xcors List of cell state pair frequencies from phylogeny.
#' @param times Vector of times/patristic distances corresponding to xcors.
#' @param freqs Vector of cell state frequencies.
#' @param total_prolif Total proliferation rate.
#' @param prolif_rates Vector of cell state proliferation rates. 
#' @param cost_func Cost function to use when comparing xcors with those using guessed parameters.  

PATHpro_loss <- function(guess, xcors, times, freqs, total_prolif, prolif_rates,
                         cost_func) {
  n <- length(freqs)
  Q <- makeQ(parameters = guess, dimensions = n)
  seq1 <- seq(1, (n-1)*n, n-1)
  seq2 <- seq(n-1, (n-1)*n, n-1)
  if(is.null(prolif_rates)) {
    r <- constrain_prolif(abs(guess[(rev(seq2)[1]+1):(rev(seq2)[1]+n)]), 
                          total_prolif, freqs)
    if(is.null(total_prolif)) {
      r <- abs(guess[(rev(seq2)[1]+1):(rev(seq2)[1]+n)])
      r <- r - min(r)
    }
  } else {
    r <- prolif_rates
  }
  A <- diag(r) + Q
  out <- 0
  for(i in 1:length(xcors)) {
    ti <- times[i]/2
    nu <- c(expm::expm(A*ti)%*%matrix(1, n))
    Delta <- diag(normalize(freqs*nu))
    Psi <- rowNorm(expm::expm(A*ti))
    out <- out + cost_func(xcors[[i]], t(Psi)%*%Delta%*%Psi, Q)
  }
  P <- expm::expm(Q)
  return(list("out"=out, "P"=P, "Q"=Q, "r"=r))
}


#' Fit PATHpro by optimizing `PATHpro_loss()`
#'
#' @param initial_parameter_guess Vector of paameter value guesses.
#' @param xcor_list List of cell state pair frequencies from phylogeny. 
#' @param state_freqs Vector of cell state frequencies. 
#' @param pair_pat_dist Vector of times/patristic distances corresponding to xcor_list.
#' @param total_prolif_rate Total proliferation rate. 
#' @param cell_prolif_rates Vector of cell state proliferation rates. 
#' @param cf Cost function to use when comparing xcor_list with those using guessed parameters. 
fitPATHpro <- function(initial_parameter_guess, xcor_list, state_freqs, 
                       pair_pat_dist = 2, total_prolif_rate = NULL, 
                       cell_prolif_rates = NULL, cf = euc_reg,
                       xtol = 1e-15, rtol = 1e-15, atol = 1e-15) {
  
  op_search <- stats::nlminb(initial_parameter_guess,
                      function(g) PATHpro_loss(guess = g, xcors = xcor_list, 
                                               times = pair_pat_dist, 
                                               freqs = state_freqs, 
                                               total_prolif = total_prolif_rate, 
                                               prolif_rates = cell_prolif_rates, 
                                               cost_func = cf)$out,
                       
                       control = list(x.tol = xtol, rel.tol = rtol, 
                                      abs.tol = atol, 
                                      eval.max = 1000, iter.max = 1000, 
                                      step.min = 1e-3, step.max = 0.1))
  
  PATHpro_loss(guess = op_search$par, xcors = xcor_list, 
               times = pair_pat_dist, 
               freqs = state_freqs, total_prolif = total_prolif_rate, 
               prolif_rates = cell_prolif_rates, 
               cost_func = cf)
}

#' PATHpro inference using input cell state and weight matrices. 
#'
#' @param z Cell state matrix. 
#' @param w Phylogenetic weight matrix. 
#' @param guess_list Vector of initial parameter guesses.
#' @param total_prolif_rate_est Estimate of total proliferation rate. 
#' @param cell_prolifs Vector of cell state specific prolfieration rates. 
#' @param cf0 Cost function. 
#' @export
PATHpro.ZW <- function(z, w, t, guess_list = NULL, total_prolif_rate_est = NULL, 
                       cell_prolifs = NULL, cf0 = euc_reg) {
  if(!is.list(w) & length(t) == 1) {
    w <- list(w)
  }
  Freqs <- list()
  for(i in 1:length(w)) {
    w[[i]] <- rowNorm(w[[i]])
    w[[i]] <- w[[i]]/sum(w[[i]])
    Freqs[[i]] <- as.matrix(Matrix::t(z)%*%w[[i]]%*%z)
  }
  if(is.null(guess_list)) {
    n <- ncol(z)
    guess_list <- list(runif(n^2, 0, 0.1))
  }
  infs <- list()
  for(i in 1:length(guess_list)) {
    infs[[i]] <- fitPATHpro(initial_parameter_guess = guess_list[[i]],
                            xcor_list = Freqs, state_freqs = Matrix::colMeans(z), 
                            pair_pat_dist = t, 
                            total_prolif_rate = total_prolif_rate_est, 
                            cell_prolif_rates = cell_prolifs, 
                            cf = cf0)
  }
  mg <- which.min(sapply(infs, '[[', "out"))
  inf0 <- infs[[mg]]
  inf0$guess <- mg
  return(inf0)
}


#' PATHpro
#' 
#' @param tree Phylogeny. 
#' @param z Cell state matrix. 
#' @param depth Max node depth to use for computing phylogenetic correlations. 
#' @param prolif Vector of cell state proliferation rates. 
#' @param total_prolif Total proliferation rate. 
#' @param guess_list List of initial parameter guess vectors. 
#' @param cf Cost function. 
#' @export
PATHpro <- function(tree, z, depth = 1, prolif = NULL, total_prolif = NULL,
                    guess_list = NULL, cf = euc_reg) {
  N <- nrow(z)
  if(depth > 1) {
    wa <- castor::get_all_pairwise_distances(tree, only_clades = 1:N, 
                                             as_edge_counts = T, check_input = F)
    wp <- castor::get_all_pairwise_distances(tree, only_clades = 1:N,
                                             as_edge_counts = F, check_input = F)
    W <- data.table::data.table("Var1"=rep(1:N, times=N), "Var2"=rep(1:N, each=N), 
                                "edge"=c(wa), "pat"=c(wp))
    out <- lapply(seq(1, depth, 1), function(i) getWt(W, i, N))
    wl <- lapply(out, '[[', "w")
    t <- sapply(out, '[[', "t")} else if (depth == 1) {
      W <- one_node_depth(tree)
      wl <- list(W$W)
      t <- W$mean.pat
    }
  inf <- PATHpro.ZW(z, wl, t, cell_prolifs = prolif, 
                    total_prolif_rate_est = total_prolif, 
                    guess_list = guess_list, 
                    cf0 = cf)
  return(inf)
}

