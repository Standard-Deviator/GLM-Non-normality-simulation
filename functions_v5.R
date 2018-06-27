# libraries needed
#require("plyr")
#require("MASS") # needed for multivariate data creation
#require("psych") # needed only for debugging
#require("car") # produces HCCMs: HC0, HC1, HC2, HC3
#require("boot") # creates bootstrap samples and CIs
#require("fungible") # needed for the genCorr function to find corr matrix using eigenvalues

# Generate data based on design matrix and constants defines within "Pop_Matricies.R"
Generate <- function(condition, fixed_objects = NULL) {
  
  # grab current condition values
  N <- condition$sample_size
  res_dist <- condition$res_dist
  kappa <- condition$mult_coll
  p <- condition$num_pred
  
  # grab betas, and only use the subset relevent to this condition
  beta <- fixed_objects[["beta"]][1:I(p+1)]
  
  sigma_sq <- fixed_objects[["sigma_sq"]]
  
  # retrieve fixed X matrix from fixed_objects
  X_name <- paste("X_ss",N,"p",p,"k",kappa,sep="")
  dat <- as.matrix(fixed_objects[[X_name]],nrow = N, ncol = p + 1)
  
  # generate residuals based on condition
  if(as.character(res_dist) == "norm"){
    resids <- rnorm(N, mean=0, sd=sqrt(sigma_sq))} else
      if(as.character(res_dist) == "cont_norm"){
        prob1 <- sample(1:2, prob=c(0.9, 0.1),size=N, replace=TRUE)
        mus1 <- c(0, 0)
        # standard deviations are weighted such that the second sd is 10 times in magnitude of the first
        sd1 <- c(sqrt(sigma_sq/10.9), 10 * sqrt(sigma_sq/10.9))
        resids <- rnorm(N, mean=mus1[prob1], sd=sd1[prob1])} else
          if(as.character(res_dist) == "high_kurt"){
            resids <- SimDesign::rValeMaurelli(N,mean=0.0, sigma = sigma_sq, skew = 0.0, kurt= 100.00)} else
              if(as.character(res_dist) == "log_norm"){
                resids <- stats::rlnorm(N,meanlog = 0, sdlog = sqrt(sigma_sq))}
  
  # create DV using general linear model formula
  y <- dat %*% beta + resids
  y_exp <- exp(y)
  
  # create DV using multipicative model
  #y_mult <- as.vector(rep(0,length =N))
  # row_index tracks the row number
  #for (row_index in 1:N){
  #  beta_index tracks which beta coefficient is currently being appended
  #  temp <- 0
  #  for (beta_index in 1:length(beta)){
      # start with the y intercept or beta 0
  #    if (beta_index == 1)
  #      {temp <- exp(beta[1])} 
  #    else # multiply the current value by the next coefficient to the power of its X value from "dat"
  #    {temp <- temp * exp(beta[beta_index]*dat[row_index,beta_index])}
  #  }
  #  y_mult[row_index] <- temp * exp(resids[row_index])
  #}
  
  
  # convert to dataframe with appropriate names
  dat <- as.data.frame(cbind(y,y_exp,dat))
  names <-c("y","y_exp","unit_vec","x1","x2","x3","x4","x5","x6")
  names <- names[1:(p+3)]
  colnames(dat) <- names
  
  list(dat=dat, parameters = beta)
}

# Analyse function returns confidence intervals using four approaches,
# as well as estimates of effect size and fit statistics
Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  # retrieve data and parameters from list
  data <- dat[[1]]
  parameters <- dat[[2]]
  
  ret <- list() # to collect output
  
  # set constants and condition values
  alpha <- fixed_objects[["alpha"]]
  N <- condition$sample_size
  
  # create formulas for linear model
  formula <- "y ~ x1"
  formula_exp <- "y_exp ~ x1"
  for (index in 2:condition$num_pred){
    temp <- paste(" + x",as.character(index), sep="")
    formula <- paste(formula, temp, sep = "")
    formula_exp <- paste(formula_exp, temp, sep = "")
  }
  
  mod_add <- lm(formula, data = data)
  mod_exp <- lm(formula_exp, data = data)
  
  # save coefficients from both models
  ret$estimates.add <- mod_add$coefficients
  ret$estimates.exp <- mod_exp$coefficients
  
  ret$r.squared.add <- summary(mod_add)$r.squared
  ret$r.squared.exp <- summary(mod_exp)$r.squared
  
  ret$adj.r.squared.add <- summary(mod_add)$adj.r.squared
  ret$adj.r.squared.exp <- summary(mod_exp)$adj.r.squared
  
  # record the first four moments of the residuals for follow-up descriptives
  ret$resid.mean.add <- mean(mod_add$residuals)
  ret$resid.var.add <- var(mod_add$residuals)
  ret$resid.skew.add <- psych::skew(mod_add$residuals)
  ret$resid.kurt.add <- psych::kurtosi(mod_add$residuals)

  # record the first four moments of the residuals for follow-up descriptives  
  ret$resid.mean.exp <- mean(mod_exp$residuals)
  ret$resid.var.exp <- var(mod_exp$residuals)
  ret$resid.skew.exp <- psych::skew(mod_exp$residuals)
  ret$resid.kurt.exp <- psych::kurtosi(mod_exp$residuals)
  
  ret$parameters <- parameters
  ret$bias <- mod_add$coefficients - parameters
  
  ###################
  # CIs relying on robustness of GLM
  wald.CIs <- as.data.frame(confint(mod_add))
  wald.CIs <- cbind(wald.CIs,wald.CIs[,2] - wald.CIs[,1])
  colnames(wald.CIs) <- c("lower","upper","efficiency")

  ret$wald <- wald.CIs

  ###################
  # CIs after log transform
  exp.CIs <- as.data.frame(confint(mod_exp))
  exp.CIs <- cbind(exp.CIs,exp.CIs[,2] - exp.CIs[,1])
  colnames(exp.CIs) <- c("lower","upper","efficiency")
  
  ret$exp <- exp.CIs  
  ###################
  # Bootstrap 95% CI for regression coefficients
  
  # function to obtain regression weights
  # http://www.statmethods.net/advstats/bootstrapping.html
  bs.fun <- function(data, indices, formula) {
    d <- data[indices,] # allows boot to select sample
    fit <- lm(formula, data=d)
    return(coef(fit))
  }
  
  bs <- boot(data=data, statistic=bs.fun,
             R=1000, formula=formula, stype="i")
  
  # initialize data frame to collect perc CI properties
  bs.CIs.perc <- data.frame(matrix(NA, nrow=length(mod_add$coefficients), ncol=2))
  colnames(bs.CIs.perc) <- c("lower","upper")
  rownames(bs.CIs.perc) <- names(mod_add$coefficients)
  
  # initialize data frame to collect bca CI properties
  bs.CIs.bca <- data.frame(matrix(NA, nrow=length(mod_add$coefficients), ncol=2))
  colnames(bs.CIs.bca) <- c("lower","upper")
  rownames(bs.CIs.bca) <- names(mod_add$coefficients)
  
  for (i in 1:length(mod_add$coefficients)) {
    # calculate bootstrap intervals
    # warnings about using extreme order statistics as endpoints are expected in n=10
    temp <- boot.ci(bs, type=c("perc","bca"),index=i)
    
    bs.CIs.perc[i,] <- c(temp$percent[4:5])
    bs.CIs.bca[i,] <- c(temp$bca[4:5])
  }
  
  ret$bs.perc <- bs.CIs.perc
  ret$bs.bca <- bs.CIs.bca
  
  ret_ <- matrix(unlist(ret),1,length(unlist(ret)),byrow=FALSE)
  nms <- names(unlist(ret))
  return <- as.numeric(data.frame(ret_))
  names(return) <- nms
  
  return
}

# Returns, upper and lower probs, as well as coverage for each estimate
Summarise <- function(condition, results, fixed_objects = NULL, parameters_list = NULL) {
  # wald ECR with tails for the additive model
  properties <- list()
  for (index in 1:I(condition$num_pred+1)){
    lowername <- paste("wald.lower",index,sep="")
    uppername <- paste("wald.upper",index,sep="")
    
    properties[[index]] <- ECR(CIs = cbind(results[,lowername], results[,uppername]),
                   parameter = fixed_objects[["beta"]][index], 
                   tails = TRUE, 
                   names=paste("wald.",index,sep=""))
    properties[[index]] <- c(properties[[index]], 1 - properties[[index]][1] - properties[[index]][2])
    names(properties[[index]]) <- c(paste("wald.Lower Prob",index,sep=""),
                                        paste("wald.Upper Prob",index,sep=""),
                                        paste("wald.Coverage",index,sep=""))
  }
  
  wald.prop <- matrix(unlist(properties),1,length(unlist(properties)),byrow=FALSE)
  nms <- names(unlist(properties))
  wald.prop <- as.numeric(data.frame(wald.prop))
  names(wald.prop) <- nms
  
  ###################
  
  # log ECR with tails for the additive model
  properties <- list()
  for (index in 1:I(condition$num_pred+1)){
    lowername <- paste("exp.lower",index,sep="")
    uppername <- paste("exp.upper",index,sep="")
    
    properties[[index]] <- ECR(CIs = cbind(results[,lowername], results[,uppername]),
                               parameter = fixed_objects[["beta"]][index], 
                               tails = TRUE, 
                               names=paste("exp.",index,sep=""))
    properties[[index]] <- c(properties[[index]], 1 - properties[[index]][1] - properties[[index]][2])
    names(properties[[index]]) <- c(paste("exp.Lower Prob",index,sep=""),
                                    paste("exp.Upper Prob",index,sep=""),
                                    paste("exp.Coverage",index,sep=""))
  }
  
  exp.prop <- matrix(unlist(properties),1,length(unlist(properties)),byrow=FALSE)
  nms <- names(unlist(properties))
  exp.prop <- as.numeric(data.frame(exp.prop))
  names(exp.prop) <- nms
  
  ###################
  
  # percentile ECR with tails
  properties <- list()
  for (index in 1:I(condition$num_pred+1)){
    lowername <- paste("bs.perc.lower",index,sep="")
    uppername <- paste("bs.perc.upper",index,sep="")
    
    properties[[index]] <- ECR(CIs = cbind(results[,lowername], results[,uppername]),
                               parameter = fixed_objects[["beta"]][index], 
                               tails = TRUE, 
                               names=paste("bs.perc",index,sep=""))
    properties[[index]] <- c(properties[[index]], 1 - properties[[index]][1] - properties[[index]][2])
    names(properties[[index]]) <- c(paste("bs.perc.Lower Prob",index,sep=""),
                                    paste("bs.perc.Upper Prob",index,sep=""),
                                    paste("bs.perc.Coverage",index,sep=""))
  }
  
  bs.perc.prop <- matrix(unlist(properties),1,length(unlist(properties)),byrow=FALSE)
  nms <- names(unlist(properties))
  bs.perc.prop <- as.numeric(data.frame(bs.perc.prop))
  names(bs.perc.prop) <- nms
  
  ###################
  
  # bca ECR with tails
  properties <- list()
  for (index in 1:I(condition$num_pred+1)){
    lowername <- paste("bs.bca.lower",index,sep="")
    uppername <- paste("bs.bca.upper",index,sep="")
    
    properties[[index]] <- ECR(CIs = cbind(results[,lowername], results[,uppername]),
                               parameter = fixed_objects[["beta"]][index], 
                               tails = TRUE, 
                               names=paste("bs.bca",index,sep=""))
    properties[[index]] <- c(properties[[index]], 1 - properties[[index]][1] - properties[[index]][2])
    names(properties[[index]]) <- c(paste("bs.bca.Lower Prob",index,sep=""),
                                    paste("bs.bca.Upper Prob",index,sep=""),
                                    paste("bs.bca.Coverage",index,sep=""))
  }
  
  bs.bca.prop <- matrix(unlist(properties),1,length(unlist(properties)),byrow=FALSE)
  nms <- names(unlist(properties))
  bs.bca.prop <- as.numeric(data.frame(bs.bca.prop))
  names(bs.bca.prop) <- nms
  
  return <- c(colMeans(results),wald.prop,exp.prop,bs.perc.prop,bs.bca.prop)
  (return <- return)
}

# code from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
# used to calculates length of vector after removing NA cases
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}

# code from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
# used to summerize based on a variable of interest, and a vector of grouping variables
get_summary <- function(var_to_summerize,grouping_vars){
    summary_return <- plyr::ddply(final_results_sub, .variables = grouping_vars, 
                                .fun = function(xx, col) {
                                  c(N    = length2(xx[[col]], na.rm=TRUE),
                                    Mean = mean   (xx[[col]], na.rm=TRUE),
                                    SD   = sd     (xx[[col]], na.rm=TRUE)
                                  )
                                },
                                var_to_summerize)
  #summary_return$se <- summary_return$sd / sqrt(summary_return$N)  # Calculate standard error of the mean
  summary_return
}
