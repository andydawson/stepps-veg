veg_build_inits <- function(K, N, N_knots, eta, rho, mu, d_knots, d_inter){
  
  bt = TRUE
  
  if (length(rho) == K){
    bt = FALSE
    W = K
  } else {
    W = K-1
  }
  
  alpha_init <- matrix(0, nrow=W, ncol=N_knots)
  g_init     <- matrix(0, nrow=W, ncol=N)
  
#   if (bt){
#     alpha_init <- matrix(0, nrow=K, ncol=N_knots)
#     g_init     <- matrix(0, nrow=K, ncol=N)
#     nfit = K
#   } else {
#     alpha_init <- matrix(0, nrow=K-1, ncol=N_knots)
#     g_init     <- matrix(0, nrow=K-1, ncol=N)
#     nfit = W
#   }
 
  if (length(eta) < (K-1)){
    eta = rep(eta, W)
  }
  
  for (k in 1:W){  
    
    C_star <- eta[k]*exp(-d_knots/rho[k]) # construct spatial covariance matrix

    alpha_init[k,] = rmvnorm(n=1, mean=rep(0,N_knots), sigma=C_star)
    
    C_star_inv = chol2inv(chol(C_star))
    c = eta[k] * exp(-d_inter/rho[k])
    
    H = c %*% C_star_inv
    g_init[k,] = mu[k] + H %*% alpha_init[k,]

  }
  
  return(list(alpha_init=alpha_init, g_init=g_init))
}