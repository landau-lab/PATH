ran_tree <- function(N=1000, rho=1, lambda=1, mu=0) {
  castor::generate_tree_hbd_reverse(Ntips=N, lambda = lambda, 
                                    mu = mu, rho = rho)$trees[[1]]
}
