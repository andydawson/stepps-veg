// Fit tree taxa count data where taxa composition is modelled using a
// predictive process model

// decompose g into complimentary pieces to help sampling

data {
  int<lower=0> K;                  // number of taxa
  int<lower=0> N;                  // number of cells
  int<lower=0> N_knots;            // number of knots
  int y[N,K];                      // pls count data
  matrix[N_knots,N_knots] d_knots; // knot distance matrix
  matrix[N,N_knots] d_inter;       // cells-knots distance matrix

  int<lower=0> N_p;                // size of P in data file; will be 1 or N
  real P;                          // FIXME: want this to be a matrix OR a real
}

transformed data {
  vector[N] ones;
  vector[N_knots] zeros; 
  int W; 
  matrix[N_knots,N_knots] eye_knots;
  matrix[N, N] M;

  for (j in 1:N_knots) {
    for (i in 1:N_knots) {
      eye_knots[i,j] <- 0.0;
    }
    eye_knots[j,j] <- 1.0;
  }

  for (j in 1:N) {
    for (i in 1:N) {
      M[i,j] <- -P;
    }
    M[j,j] <- 1.0-P;
  }

  for (i in 1:N) ones[i] <- 1.0;
  for (i in 1:N_knots) zeros[i] <- 0.0;
  W <- K-1;
}

parameters {
  //vector<lower=0, upper=50>[W] eta;  
  vector<lower=0.1, upper=sqrt(100)>[W] eta;
  //vector<lower=0.1>[W] sigma;  
  vector<lower=1.0e-6, upper=1.0>[W] rho;  

  vector[W] mu;  
  vector[N_knots] alpha[W];
  vector[N] g[W];
}


transformed parameters {
  //vector<lower=0, upper=50>[W] eta; 
  //vector<lower=0, upper=50>[W] sigma_sq;  
  
  //for (i in 1:W) eta[i] <- pow(sqrt_eta[i],2);
}

model {
  //declarations
  matrix[N_knots,N_knots] C_star[W];
  matrix[N_knots,N_knots] C_star_L[W];
  matrix[N_knots,N_knots] C_star_inv[W];
  
  print("eta = ", eta);

  // priors 
  eta      ~ uniform(0.1,sqrt(100));//cauchy(0,25);
  rho      ~ uniform(1.0e-6,1.0);
  mu       ~ normal(0,5);
  
  // construct knots covariance matrix C*
  for (k in 1:W){
    for (i in 1:N_knots)
      for (j in i:N_knots)
        C_star[k,i,j] <- pow(eta[k],2) * exp(-1/rho[k] * d_knots[i,j]);
          //+ if_else(i==j, sigma_sq[k], 0.0);
          
    for (i in 1:N_knots)
      for (j in (i+1):N_knots)
        C_star[k,j,i] <- C_star[k,i,j];  
    
    C_star_L[k]   <- cholesky_decompose(C_star[k]);
    C_star_inv[k] <-  mdivide_right_tri_low(mdivide_left_tri_low(C_star_L[k], eye_knots)', C_star_L[k])';

    // latent GPs; one for each taxa
    alpha[k] ~ multi_normal_prec(zeros, C_star_inv[k]);
    //alpha[k] ~ multi_normal_cholesky(zeros, C_star_L[k]);
  } // end k loop 
 
  { 
  matrix[N,N_knots] c[W];
  //vector[N] g[W];
  vector[N] sum_exp_g;
  vector[K] r[N];
  vector[N] sig;
  vector[N] sqrt_var;
  row_vector[N_knots] c_i; 
  
  vector[N] mu_2;   

  for (k in 1:W){	    
    for (i in 1:N){
      for (j in 1:N_knots){
        c[k][i,j] <- pow(eta[k],2) * exp(-1/rho[k] * d_inter[i,j]);   
      }
    }
  }
  
  for (k in 1:W){
      
    mu_2 <- M * (c[k] * (C_star_inv[k] * alpha[k]));			

    for (i in 1:N){

      c_i <- row(c[k],i);
      sig[i] <- c_i * C_star_inv[k] * c_i';
      sqrt_var[i]  <- sqrt(pow(eta[k],2) + sig[i]);


      g[k][i] ~ normal(mu[k] + mu_2[i], sqrt_var[i]);

    }

    //g[k] ~ normal(mu[k] * ones + M * (c * (C_star_inv[k] * alpha[k])), sqrt_var);
      
  } // end k loop
    
  // sum process vals for each i
  for (i in 1:N) {
    sum_exp_g[i] <- 0.0;
    for (k in 1:W)
      sum_exp_g[i] <- sum_exp_g[i] + exp(g[k,i]);
  }

  // additive log-ratio transformation
  for (k in 1:W)
    for (i in 1:N)
      r[i,k] <- exp(g[k,i]) / (1 + sum_exp_g[i]);
  for (i in 1:N)
      r[i,K] <- 1 / (1 + sum_exp_g[i]);
  
  // link composition vector to count data through multinomial
  for (i in 1:N) 
    y[i] ~ multinomial(r[i]);
  
  } // end scope
}
