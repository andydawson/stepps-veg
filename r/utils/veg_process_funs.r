###########################################################################################
# computes proportion chains from a stan fit object
###########################################################################################
library(Rcpp)
library(inline)
# library(bayesSurv)

# source('r/utils/simDataFuns.r')
# source('prediction/r/utils/pred_plot_funs.r')

# # do matrix multiply
# cppFunction('
#   NumericVector scale_up(NumericMatrix c, NumericMatrix C_star_inv, NumericVector alpha) {
#     //int N=c.nrow();
#     //NumericVector X(N);
#     
#     //X = c * (C_star_inv * alpha);
# 
#     return c * (C_star_inv * alpha);
#   }
# ', verbose=TRUE)

# # do matrix multiply
# code<- 'arma::mat c = Rcpp::as<arma::mat>(mc); 
#         arma::mat C_star_inv = Rcpp::as<arma::mat>(mC_star_inv);
#         arma::vec alpha = Rcpp::as<arma::vec>(malpha);
#     //int N=c.nrow();
#     //NumericVector X(N);
#     
#     //X = c * (C_star_inv * alpha);
# 
#     return wrap(c * (C_star_inv * alpha));'
# 
# scale_up2 <- cxxfunction( signature(mc="numeric", mC_star_inv="numeric", malpha="numeric"), 
#                                      code, plugin="RcppArmadillo")


# code <- 'arma::mat x = Rcpp::as<arma::mat>(X);
# arma::mat mu = Rcpp::as<arma::mat>(Mu);
# arma::mat sigma = Rcpp::as<arma::mat>(Sigma);
# int n = x.n_rows;
# arma::vec md(n);
# for (int i=0; i<n; i++){
# arma::mat x_i = x.row(i) - mu;
# arma::mat Y = arma::solve( sigma, arma::trans(x_i) );
# md(i) = arma::as_scalar(x_i * Y);
# }
# return wrap(md);'
# mahalanobis_RcppArma <- cxxfunction( signature(X="numeric",
#                                                Mu="numeric", Sigma="numeric"), code, plugin="RcppArmadillo")


# # additive log_ratio transformation
# cppFunction('
#   NumericMatrix sum2one_constraint(int K, int N, NumericMatrix g, NumericVector sum_exp_g) {
#     NumericMatrix r(N, K);
#     for (int k = 0; k<(K-1); k++)
#       for (int j = 0; j<N; j++)
#         r(j,k) = exp(g(j,k))/(1+sum_exp_g(j));
# 
#     for (int j = 0; j<N; j++)
#       r(j,K) = 1.0 / (1 + sum_exp_g[j]);
#    
# return r;
#   }
# ', verbose=TRUE)

# # compute Halpha matrix mults; doesn't work, no matrix mults in Rcpp without a lin alg package
# cppFunction('
#   NumericVector Halpha(NumericMatrix M, NumericMatrix c, NumericMatrix C_star_inv, NumericVector alpha) {
#     int N = M.nrow();
#     NumericVector Halpha(N);
#     
#     std::cout << N << std::endl;
# 
#     std::cout << "Check alpha size :" << alpha.size() << std::endl;
#     std::cout << "Check C_star_inv size :" << C_star_inv.nrow() << std::endl;
# 
#     Halpha = M * (c * (C_star_inv * alpha));
#    
#   return Halpha;
#   }
# ', verbose=TRUE)
# 
# HalphaR <- function(M, c, C_star_inv, alpha){
#   return(M %*% c %*% (C_star_inv %*% alpha))
# }

compute_prop_chains <- function(fit, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp, nug=FALSE){
  
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  N = nrow(d_inter)
  N_knots = ncol(d_inter)
  niters   = nrow(post[,1,]) 
  
  # W is the number of fitted taxa
  if (bt){
    W = K-1
  } else {
    W = K
  }
  
  g = array(NA, c(W, N, niters))
  
  if (decomp){
    
    x = matrix(1, nrow=N, ncol=1)
    N_p = N
    
    temp = qr(x)
    Q = qr.Q(temp)
    R = qr.R(temp)
    
#     P = Q %*% t(Q)
#     # M = diag(N) - P
#     
#     if (all(P-P[1,1]<1.0e-12)){
#       P = P[1,1]
#       N_p = 1
#     }

#     P = matrix(P, nrow=N,ncol=N)
#     M = diag(N) - P
  }
  
  for (k in 1:W){
    print(k)
    time.start.k <- Sys.time()
    
    if (mpp){
      if (nug) {
        g_cols = seq((4*W + W*N_knots + k), (4*W + W*N_knots + W*N), by=W) 
      } else if (length(which(col_substr == 'et')) == 1) {
        g_cols = seq((2*W + 1 + k + W*N_knots), (2*W + 1 + W*N_knots + W*N), by=W)
      } else {
        g_cols = seq((3*W + k + W*N_knots), (3*W + W*N_knots + W*N), by=W) 
      }
      
      g[k,,] = t(post[,1,g_cols])
      
    } else {
    
      col_substr = substr(colnames(post[,1,]),1,2)
      
      if (length(which(col_substr == 'et')) == 1){
        eta   = post[,1,which(col_substr == 'et')]
      } else {
        eta   = post[,1,which(col_substr == 'et')[k]]
      }
      
      rho   = post[,1,which(col_substr == 'rh')[k]]
      mu    = post[,1,which(col_substr == 'mu')[k]]
      
      #     if (nug) sigma = post[,1,which(col_substr == 'si')[k]]
      
      #     eta   = post[,1,k]
      #     sigma = post[,1,k + W]
      #     rho   = post[,1,k + 2*W]
      #     mu    = post[,1,k + 3*W]
      
      #     if (nug) {
      #       knot_cols = seq((4*W + k + 1), (4*W + W*N_knots + 1), by=W) 
      #       } else {
      #         knot_cols = seq((3*W + k + 1), (3*W + W*N_knots + 1), by=W) 
      #     }
      
      if (nug) {
        knot_cols = seq((4*W + k), (4*W + W*N_knots), by=W) 
      } else if (length(which(col_substr == 'et')) == 1) {
        knot_cols = seq((2*W +1 + k), (2*W + 1 + W*N_knots), by=W)
      } else {
        knot_cols = seq((3*W + k), (3*W + k + W*N_knots), by=W) 
      }
      
      alpha     = post[,1, knot_cols]
      
      for (i in 1:niters){
        print('Start matrix mult!')
        time.start <- Sys.time()
        C_star = eta[i]^2 * exp(-1/rho[i] * d_knots)
        #       diag(C_star) = diag(C_star) + sigma_sq[i]
        C_star_inv = chol2inv(chol(C_star))
        
        c = eta[i]^2 * exp(-1/rho[i] * d_inter)
        
        #       H = c %*% C_star_inv
        #       H = scale_up2(c, C_star_inv, alpha[i,])
        if (decomp){
          g[k,,i] = mu[i] + Q %*% (t(Q) %*% (c %*% (C_star_inv %*% alpha[i,])))
#           g[k,,i] = mu[i] + M %*% c %*% (C_star_inv %*% alpha[i,])
        } else {
          g[k,,i] = mu[i] + c %*% (C_star_inv %*% alpha[i,])
        }    
        time.end <- Sys.time()
        print(time.end-time.start)
      }
    
    }
  time.end.k <- Sys.time()
  print(time.end.k - time.start.k)
  }
  



  sum_exp_g = matrix(NA, nrow=N, ncol=niters)
  exp_g     = exp(g)
  for (i in 1:niters){
    sum_exp_g[,i] = colSums(exp_g[,,i])
  }
  
  r = array(NA, c(K, N, niters))
  
  if (bt){
    for (i in 1:niters){
      for (k in 1:W){  
        r[k,,i] = exp_g[k,,i]/(1 + sum_exp_g[,i])
      }
      r[W+1,,i] = 1/(1 + sum_exp_g[,i])
    }
  } else {
    for (i in 1:niters){
      for (k in 1:W){  
        r[k,,i] = exp_g[k,,i]/sum_exp_g[,i]
      }
    }
  }
  
#   r[,,i] <- sum2one_constraint(K, N, g, sum_exp_g)
  
  return(list(r=r, g=g))
  
}


compute_prop_chains_Halpha <- function(post, d_knots, d_inter, K, P, decomp=decomp, bt=bt, mpp=mpp, nug=FALSE){
  
  N = nrow(d_inter)
  N_knots = ncol(d_inter)
  niters   = nrow(post[,1,]) 
  
  # W is the number of fitted taxa
  if (bt){
    W = K-1
  } else {
    W = K
  }
  
  g = array(NA, c(W, N, niters))
  Halpha = array(NA, c(W, N, niters))
  
  if (decomp){
    x = matrix(1, nrow=N, ncol=1)
    N_p = N
    
    temp = qr(x)
    Q = qr.Q(temp)
    R = qr.R(temp)
#     
#     P = matrix(P, nrow=N,ncol=N)
#     M = diag(N) - P
  }
  
  for (k in 1:W){
    start.time.k <-Sys.time()
    print(k)
    
    if (mpp){
      
      col_substr = substr(colnames(post[,1,]),1,2)
      
      if (length(which(col_substr == 'et')) == 1){
        eta   = post[,1,which(col_substr == 'et')]
      } else {
        eta   = post[,1,which(col_substr == 'et')[k]]
      }
      
      rho   = post[,1,which(col_substr == 'rh')[k]]
      mu    = post[,1,which(col_substr == 'mu')[k]]
      
      if (nug) {
        knot_cols = seq((4*W + k), (4*W + W*N_knots), by=W) 
      } else if (length(which(col_substr == 'et')) == 1) {
        knot_cols = seq((2*W +1 + k), (2*W + 1 + W*N_knots), by=W)
      } else {
        knot_cols = seq((3*W + k), (3*W + W*N_knots), by=W) 
      }
      
      
      alpha     = post[,1, knot_cols]
      
      for (i in 1:niters){
#         print('Start matrix mult!')
#         start.time <- Sys.time()
        C_star = eta[i]^2 * exp(-1/rho[i] * d_knots)
        #       diag(C_star) = diag(C_star) + sigma_sq[i]
        C_star_inv = chol2inv(chol(C_star))
        
        c = eta[i]^2 * exp(-1/rho[i] * d_inter)
        
        #       H = c %*% C_star_inv
        #       H = scale_up2(c, C_star_inv, alpha[i,])
        if (decomp){
          cCinvalpha = c %*% (C_star_inv %*% alpha[i,] )
          Halpha[k,,i] = cCinvalpha - Q %*% ( t(Q) %*% ( cCinvalpha ) )
        } else {
          Halpha[k,,i] = c %*% (C_star_inv %*% alpha[i,])
        }  
#         end.time <- Sys.time()
#         print(end.time-start.time)
      }
      
      if (nug) {
        g_cols = seq((4*W + W*N_knots + k), (4*W + W*N_knots + W*N), by=W) 
      } else if (length(which(col_substr == 'et')) == 1) {
        g_cols = seq((2*W + 1 + k + W*N_knots), (2*W + 1 + W*N_knots + W*N), by=W)
      } else {
        g_cols = seq((3*W + k + W*N_knots), (3*W + W*N_knots + W*N), by=W) 
      }
      
      g[k,,] = t(post[,1,g_cols])
      
    } else {
      
      col_substr = substr(colnames(post[,1,]),1,2)
      
      if (length(which(col_substr == 'et')) == 1){
        eta   = post[,1,which(col_substr == 'et')]
      } else {
        eta   = post[,1,which(col_substr == 'et')[k]]
      }
      
      rho   = post[,1,which(col_substr == 'rh')[k]]
      mu    = post[,1,which(col_substr == 'mu')[k]]
      
      #     if (nug) sigma = post[,1,which(col_substr == 'si')[k]]
      
      #     eta   = post[,1,k]
      #     sigma = post[,1,k + W]
      #     rho   = post[,1,k + 2*W]
      #     mu    = post[,1,k + 3*W]
      
      #     if (nug) {
      #       knot_cols = seq((4*W + k + 1), (4*W + W*N_knots + 1), by=W) 
      #       } else {
      #         knot_cols = seq((3*W + k + 1), (3*W + W*N_knots + 1), by=W) 
      #     }
      
      if (nug) {
        knot_cols = seq((4*W + k), (4*W + W*N_knots), by=W) 
      } else if (length(which(col_substr == 'et')) == 1) {
        knot_cols = seq((2*W +1 + k), (2*W + 1 + W*N_knots), by=W)
      } else {
        knot_cols = seq((3*W + k), (3*W + W*N_knots), by=W) 
      }
      
      alpha     = post[,1, knot_cols]
      
      for (i in 1:niters){
#         print('Start matrix mult!')
#         start.time <- Sys.time()
        C_star = eta[i]^2 * exp(-1/rho[i] * d_knots)
        #       diag(C_star) = diag(C_star) + sigma_sq[i]
        C_star_inv = chol2inv(chol(C_star))
        
        c = eta[i]^2 * exp(-1/rho[i] * d_inter)
        
        #       H = c %*% C_star_inv
        #       H = scale_up2(c, C_star_inv, alpha[i,])
        if (decomp){
          Halpha[k,,i] = Q %*% (t(Q) %*% (c %*% (C_star_inv %*% alpha[i,])))
          #           Halpha[k,,i] = M %*% c %*% (C_star_inv %*% alpha[i,])
          g[k,,i] = mu[i] + Halpha[k,,i]#M %*% c %*% (C_star_inv %*% alpha[i,])
        } else {
          Halpha[k,,i] = c %*% (C_star_inv %*% alpha[i,])
          g[k,,i] = mu[i] +  Halpha[k,,i]#c %*% (C_star_inv %*% alpha[i,])
        }
#         end.time <- Sys.time()
#         print(end.time-start.time)
#         
      }
      
    }
    
    end.time.k <- Sys.time()
    print(end.time.k-start.time.k)
    
  }
  
  sum_exp_g = matrix(NA, nrow=N, ncol=niters)
  exp_g     = exp(g)
  for (i in 1:niters){
    sum_exp_g[,i] = colSums(exp_g[,,i])
  }
  
  r = array(NA, c(K, N, niters))
  
  if (bt){
    for (i in 1:niters){
      for (k in 1:W){  
        r[k,,i] = exp_g[k,,i]/(1 + sum_exp_g[,i])
      }
      r[W+1,,i] = 1/(1 + sum_exp_g[,i])
    }
  } else {
    for (i in 1:niters){
      for (k in 1:W){  
        r[k,,i] = exp_g[k,,i]/sum_exp_g[,i]
      }
    }
  }
  
  #   r[,,i] <- sum2one_constraint(K, N, g, sum_exp_g)
  
  return(list(r=r, g=g, Halpha=Halpha))
  
}


get_mu <- function(post, W){
  
  col_substr = substr(colnames(post[,1,]),1,2)
  idx_mu = which(col_substr == 'mu')
  
  mu = t(post[,1,idx_mu])
  
  return(mu)
}

# ###########################################################################################
# # computes proportion chains from a stan fit object
# ###########################################################################################
# compute_prop_chains <- function(fit, d_knots, d_inter, dim_pars){
#   
#   post = extract(fit, permuted=FALSE, inc_warmup=FALSE)
#   
#   N = nrow(d_inter)
#   N_knots = ncol(d_inter)
#   niters   = nrow(post[,1,]) 
#   W        = K - 1
#   
#   g = array(NA, c(W, N, niters))
#   
#   for (k in 1:W){
#     print(k)
#     
#     eta      = post[,1,k]
#     rho      = post[,1,k + W]
#     sigma_sq = post[,1,k + 2*W]
#     mu       = post[,1,k + 3*W]
#     
# #     knot_cols = seq((4*W + k), (4*W + W*N_knots), by=W)
# #     alpha     = post[,1, knot_cols]
#     
#     for (i in 1:niters){
#       C_star = eta[i] * exp(-1/rho[i] * d_knots)
#       #       diag(C_star) = diag(C_star) + sigma_sq[i]
#       C_star_inv = chol2inv(chol(C_star))
#       
#       c = eta[i] * exp(-1/rho[i] * d_inter)
#       
#       H = c %*% C_star_inv
#       alpha[i,] = rmvnorm(n=1, mean=rep(0,N_knots*T), sigma=C_star)
#       g[k,,i] = mu[k] + H %*% alpha[i,]
#       
#     }
#   }
#   
#   sum_exp_g = matrix(NA, nrow=N, ncol=niters)
#   exp_g     = exp(g)
#   for (i in 1:niters){
#     sum_exp_g[,i] = colSums(exp_g[,,i])
#   }
#   
#   r = array(NA, c(W+1, N, niters))
#   for (i in 1:niters){
#     for (k in 1:W){  
#       r[k,,i] = exp_g[k,,i]/(1 + sum_exp_g[,i])
#     }
#     r[W+1,,i] = 1/(1 + sum_exp_g[,i])
#   }
#   
#   return(r)
#   
# }

###########################################################################################
# computes proportion chain means for a proportions object that is
# ntaxa x ncells x niters in dimension
###########################################################################################
compute_prop_means <- function(r){
  
  ntaxa  = dim(r)[1]
  ncells = dim(r)[2]
  niters = dim(r)[3]
  
  r_means = matrix(NA, ntaxa, ncells)
  for (k in 1:ntaxa){
    r_means[k, ] = rowMeans(r[k,,])  
  }
  
  return(r_means)
}


###########################################################################################
# compute effective sample size form a stanfit object
###########################################################################################
ess <- function(summary_fit, ntaxa){
  ess = summary_fit[,"n_eff"]
  return(c(ess[1:(3*ntaxa)], ess[length(ess)]))
}

###########################################################################################
# compute acceptance rate for a stanfit object
###########################################################################################
accept_rate <- function(fit){
  post = extract(fit, permuted=FALSE, inc_warmup=FALSE)
  ar = apply(post, "chains", FUN = function(x) nrow(unique(as.data.frame(x)))) / nrow(post) # acceptance rates
  return(ar)
}


# ###########################################################################################
# # plot predicted and observed proportion
# ###########################################################################################
# 
# prop_maps <- function(y, r, meta){
# 
#   props <<- data.frame(t(r_means), x = meta$x, y = meta$y)
#   foo <<- melt(props, c('x','y'))
# 
#   pdf(paste('figures/', N_knots, 'knots_', N, 'cells_', K, 'taxa_props.pdf', sep=''), width=8, height=5)
#   d <<- ggplot() + geom_raster(data=foo, aes(x=x, y=y, fill=value)) + scale_fill_gradient(low = "white", high = "red", limits=c(0,1))
#   # d <- d + geom_text(aes(label=prop), data=dat, size=4, x=x1, y=y1)
#   # d <- d + geom_point(aes(x=x0, y=y0), data=dat, colour='black', size=2, shape=1)
#   d <<- d + facet_wrap(~variable, ncol=3)
#   d <<- d + scale_y_continuous(name="") + scale_x_continuous(name="") + 
#        theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank()) +
#        ggtitle("Predictions")
#   d
#   dev.off()
# }
# 
# 
# colnames(y) = c(colnames(counts)[1:(K-1)], 'other')
# props_data = t(apply(y, 1, function(x) x/sum(x)))
# 
# 
# props_data = data.frame(props_data, x = meta$x, y = meta$y)
# foo = melt(props_data, c('x','y'))
# 
# pdf(paste('figures/', N_knots, 'knots_', N, 'cells_', K, 'taxa_props_obs.pdf', sep=''), width=8, height=5)
# d <- ggplot() + geom_raster(data=foo, aes(x=x, y=y, fill=value)) + scale_fill_gradient(low = "white", high = "red", limits=c(0,1))
# # d <- d + geom_text(aes(label=prop), data=dat, size=4, x=x1, y=y1)
# # d <- d + geom_point(aes(x=x0, y=y0), data=dat, colour='black', size=2, shape=1)
# d <- d + facet_wrap(~variable, ncol=3)
# d <- d + scale_y_continuous(name="") + scale_x_continuous(name="") + 
#   theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank()) +
#   ggtitle("Observed")
# d
# dev.off()
