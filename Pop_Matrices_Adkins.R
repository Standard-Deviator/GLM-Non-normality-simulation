set.seed(7272)

# create empty list
fixed_objects <- list()

fixed_objects[["sigma_sq"]] <- 1

# define population values for IV standard deviations
fixed_objects[["IV_sds"]] <- round(runif(max(num_preds),1,3),2)

# scale IV_sds so that their the sum of their variances are equal to sigma_sq
# rationale here is that the overall R^2 should be approximatey 0.5 or less to avoid over fitting
fixed_objects[["IV_sds"]] <- sqrt((fixed_objects[["sigma_sq"]]/ sum(fixed_objects[["IV_sds"]]^2)) * fixed_objects[["IV_sds"]]^2)

# define y intercept, or mean of Y in this case
y_int <- 10

# define population betas for IV (additive and multiplicative models)
fixed_objects[["beta"]] <- c(y_int,round(runif(max(num_preds),0,2),2))

fixed_objects[["alpha"]] <- 0.05

for (ss in sample_sizes){
  for (p in num_preds){
    for (kappa in mult_colls){
      # Use Jones & Waller 2013 algorithm for predictor correlation matrix
      # generate largest and smallest eigenvalues using the following simultaneous equations
      # 1. lambda_1 + lambda_2 = 0.6667*p    ----- This was derived such that the largest and smallest sum to 2/3s of the trace
      # 2a. kappa = sqrt(lambda_1/lambda_2)
      # 2b. lambda_1 - kappa^2 * lambda_2 = 0 ----- This was derived using the formula for condition number
      lambda <- solve(matrix(c(1,-kappa^2,1,1),2,2,byrow=TRUE),c(0, 0.6667*p))
    
      # randomly generate remainging eigenvalues bounded by the largest and smallest
      tmp.evals <- c(lambda[1], runif(p-2,lambda[2],lambda[1]), lambda[2])
    
      # Sort and scale the eigenvalues from largest to smallest
      (evals <- sort(tmp.evals*(p/sum(tmp.evals)),dec=TRUE))
    
      # Use Marsaglia and Olkin's (1984) Method III
      # to construct a correlation matrix with the above eigenvalues (from Jones, 2010).
      Rx <- fungible::genCorr(evals,seed = 7272)
      
      # convert to covariance matrix
      COVx <- psych::cor2cov(Rx, sigma= fixed_objects[["IV_sds"]][1:p])
      
      # name COVx based on numberof IVs and kappa
      COVx_name <- paste("COVx_p",p,"k",kappa,sep="")
      fixed_objects[[COVx_name]] <- COVx
      
      # name fixed X matrix
      X_name <- paste("X_ss",ss,"p",p,"k",kappa,sep="")
      X <- MASS::mvrnorm(n=ss,mu=rep(0,p), Sigma=COVx)
      
      # append unit vector to include intercept
      X <- cbind(rep(1,ss),X)
      fixed_objects[[X_name]] <- X
    }
  }
}

# reset seed so that data generation is random across conditions
set.seed(NULL)
