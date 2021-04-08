# distribution-plots-helpers.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 14, 2020
#
# Helper functions for distribution-plots.R

library(tidyverse)

# Downsample counts matrix
# Source: https://github.com/willtownes/scrna2019/blob/028363f04139a58b19143f3058ca8fa4a3533b63/util/functions.R
Down_Sample_Matrix<-function(expr_mat,min_lib_size=NULL){
  min_sz<-min(colSums(expr_mat))
  if(is.null(min_lib_size)){
    min_lib_size<-min_sz
  } else {
    stopifnot(min_lib_size<=min_sz)
  }
  down_sample<-function(x){
    prob <- min_lib_size/sum(x)
    unlist(lapply(x,function(y){rbinom(1, y, prob)}))
  }
  apply(expr_mat, 2, down_sample)
}

# Plot P(X_i = 0) against average expression level mu_i
plot_prob = function(dat_sub){
  n = round(median(colSums(dat_sub)))
  means = rowMeans(dat_sub)
  vars = apply(dat_sub, MARGIN = 1, var)
  
  # empirical P(X_i = 0)
  emp_probs_0 = apply(dat_sub, MARGIN = 1, function(r) sum(r==0)/ncol(dat_sub))
  plot_dt = tibble(means, emp_probs_0)
  
  emp_props = rowSums(dat_sub)/sum(colSums(dat_sub))
  # Model P(X_i = 0) under Binomial
  binom_probs_0 = dbinom(x = 0, size = n, prob = emp_props)
  # Model P(X_i = 0) under Poisson
  poiss_probs_0 = dpois(x = 0, lambda = n*emp_props)
  # Model P(X_i = 0) under Negative Binomial
  # Estimate size/dispersion parameter
  model = lm(vars ~ 1*means + I(means^2) + 0, tibble(means, vars))
  phi = 1/coef(model)["I(means^2)"]
  nbinom_probs_0 = dnbinom(x = 0, size = phi, mu = n*emp_props) 
  
  # Tibble for plot
  plot_lines_dt = tibble(means = means,
                         binomial = binom_probs_0,
                         poisson = poiss_probs_0,
                         nbinomial = nbinom_probs_0) %>%
    pivot_longer(-means, names_to = "model", values_to = "probs_0")
  
  # Plot
  plt = plot_lines_dt %>%
    ggplot(aes(x = means, y = probs_0)) +
    geom_point(data = plot_dt, aes(x = means, y = emp_probs_0), alpha = 0.4) + # Add data points
    geom_line(aes(color = model)) + # Add lines for models
    scale_x_log10() +
    labs(x = "Log of mean expression",
         y = "Fraction of zero droplets") +
    theme_bw() +
    theme(text = element_text(size = 15))
  
  # Return object
  out = list(plot = plt,
             lines_dt = plot_lines_dt)
  return(out)
}

# BIC functions
# Multinomial
mult_bic<-function(m){
  #multinomial model
  n<-colSums(m)
  p<-rowSums(m)/sum(n)
  ll<-sum(apply(m,2,dmultinom,prob=p,log=TRUE))
  df<-nrow(m)-1
  print(paste0("Loglikelihood is ", ll, " df is ", df))
  -2*ll+df*log(prod(dim(m)))
}

# Dirichlet multinomial (for overdispersion)
dmn_bic<-function(m){
  fit<-DirichletMultinomial::dmn(t(m),1)
  alpha<-drop(fit@fit$Estimate)
  ll<-sum(extraDistr::ddirmnom(t(m), colSums(m), alpha, log = TRUE))
  df = nrow(m)
  print(paste0("Loglikelihood is ", ll, " df is ", df))
  -2*ll+nrow(m)*log(prod(dim(m)))
}

# Poisson
poi_fit<-function(m,X=NULL,sz=NULL,maxit=100){
  #poisson
  if(is.null(sz)){ sz<-log(colMeans(m)) }
  if(is.null(X)){ #no covariates => closed form solution
    sz<-exp(sz) # sz is library size / number of genes = colMeans(m)
    lam<-rowSums(m)/sum(sz) # each gene has an empirically estimated lambda (~ proportion of counts in that gene)
    mu<-outer(lam,sz) # outer product of lambda and sz (equivalent to prop of counts in gene * library size), each cell x gene has a mu
    ll<-matrix(dpois(m,mu,log=TRUE),nrow=nrow(m)) # Calculate loglikelihood in every cell x gene
    return(data.frame(ll=rowSums(ll),converged=TRUE)) # Sum loglikelihoods in each row
  } else { #covariates included
    stopifnot(all(X[,1]==1))
    k<-ncol(X)
    fam<-poisson()
    #default maxit is 25 for glm.fit, but some fits take longer to converge
    ctl<-list(maxit=maxit) 
    f<-function(y){
      fit<-glm.fit(X,y,offset=sz,control=ctl,family=fam)
      c(ll=k-fit$aic/2, converged=fit$converged)
    }
    res<-as.data.frame(t(apply(m,1,f)))
    res$converged<-as.logical(res$converged)
    return(res)
  }
}

poi_bic<-function(m,X=NULL,prefit=NULL,maxit=100){
  #poisson. prefit should be a data frame with column "ll" for log-likelihood
  k<-if(is.null(X)){ 1 } else { ncol(X) }
  df<-k*nrow(m)
  if(is.null(prefit)){
    prefit<-poi_fit(m,X,maxit=maxit)
  }
  ll<-sum(prefit$ll) # Sum the sums of loglikelihoods across genes
  #compute BIC: -2*loglik+df*log(n_obs)
  print(paste0("Loglikelihood is ", ll, " df is ", df))
  -2*ll+df*log(prod(dim(m)))
}

# Negative binomial (single overdispersion parameter for all genes)
nb_bic_1 = function(m){
  # estimate phi
  means = rowMeans(m)
  vars = apply(m, MARGIN = 1, var)
  model = lm(vars ~ 1*means + I(means^2) + 0, tibble(means, vars))
  phi = 1/coef(model)["I(means^2)"]
  
  # estimate mu
  sz <- colMeans(m) # sz is library size / number of genes = colMeans(m)
  lam <- rowSums(m)/sum(sz) # each gene has an empirically estimated lambda (~ proportion of counts in that gene)
  mu <- outer(lam, sz) # outer product of lambda and sz (equivalent to prop of counts in gene * library size), each cell x gene has a mu
  
  # get ll
  ll = matrix(dnbinom(m, size = phi, mu = mu, log = TRUE), nrow = nrow(m)) # get likelihood for every cell x gene
  ll = sum(ll) # sum likelihood
  df = nrow(m) + 1
  print(paste0("Loglikelihood is ", ll, " df is ", df))
  -2*ll+df*log(prod(dim(m)))
}

# Negative binomial (separate overdispersion parameter for each gene)
nb_fit<-function(m,X=NULL,sz=NULL){
  #neg binom
  if(is.null(sz)){ sz<-log(colMeans(m))} 
  if(is.null(X)){
    f<-function(y){
      mgcv::gam(y~1,family=mgcv::nb,offset=sz,method="ML") # Each cell is a different observation for a given gene (y = one row of m)
    }
  } else {
    f<-function(y){
      mgcv::gam(y~X-1,family=mgcv::nb,offset=sz,method="ML")
    }
  }
  g<-function(y){
    fit<-f(y)
    th<-fit$family$getTheta(TRUE) #FALSE (default) gives log(theta)
    ll<-as.numeric(logLik(fit))
    c(theta=th,ll=ll,converged=fit$converged) # different theta for every gene
  }
  res<-as.data.frame(t(apply(m,1,g)))
  res$converged<-as.logical(res$converged)
  res
}

nb_bic_2<-function(m,X=NULL,prefit=NULL){
  if(is.null(prefit)){ 
    prefit<-nb_fit(m,X)
  } else {
    stopifnot(nrow(m)==nrow(prefit))
  }
  ll<-sum(prefit$ll)
  k<-if(is.null(X)){ 2 } else { ncol(X)+1 }
  df = k*nrow(m)
  print(paste0("Loglikelihood is ", ll, " df is ", df))
  -2*ll+df*log(prod(dim(m)))
}

# Pearson's chi-squared test
library(MASS)

# Bin by unique value
p_chisq_test = function(m, distribution = "poisson"){
  ind = c() # rows that fail
  p_values = c() # p-values for every row
  chi_squares = c() # chi-squared value for every row
  dfs = c() # degrees of freedom for every row
  
  if(distribution == "poisson"){
    for(i in 1:nrow(m)){
      x = m[i, ]
      poisson_fit = fitdistr(x, "poisson")
      lambda_fit = poisson_fit$estimate["lambda"]
      n_samples = length(x)
      f_obs = table(x)
      
      bins = as.numeric(names(f_obs))
      probs = dpois(bins, lambda = lambda_fit)
      f_hyp = probs*n_samples
      other_value = (1 - sum(probs))*n_samples
      
      f_obs = c(f_obs, 0) # Observed value
      f_hyp = c(f_hyp, other_value) # Expected value
      
      # Pool counts if < 5
      if(any(f_obs < 5)){
        start_pool_bin = which(f_obs < 5)[1]
        f_obs_pool = sum(f_obs[start_pool_bin:length(f_obs)])
        f_hyp_pool = sum(f_hyp[start_pool_bin:length(f_hyp)])
        
        f_obs = c(f_obs[1:(start_pool_bin - 1)], f_obs_pool)
        f_hyp = c(f_hyp[1:(start_pool_bin - 1)], f_hyp_pool)
      }

      p = 1 
      df = (length(bins) + 1) - p
      chiSquare = sum((f_obs-f_hyp)^2/f_hyp)
      
      p_values = c(p_values, 1 - pchisq(chiSquare, df = df))
      chi_squares = c(chi_squares, chiSquare)
      dfs = c(dfs, df)
    }
  } else if (distribution == "nb"){
    for(i in 1:nrow(m)){
      t = tryCatch({
        x = m[i, ]
        
        nb_fit = suppressWarnings(fitdistr(x, "negative binomial"))
        r_fit = nb_fit$estimate['size']
        p_fit = r_fit / (nb_fit$estimate['mu'] + r_fit)
        n_samples = length(x)
        f_obs = table(x)
        
        bins = as.numeric(names(f_obs))
        probs = dnbinom(bins, size = r_fit, prob = p_fit)
        f_hyp = probs*n_samples
        other_value = (1 - sum(probs))*n_samples
        
        f_obs = c(f_obs, 0) # Observed value
        f_hyp = c(f_hyp, other_value) # Expected value
        
        # Pool counts if < 5
        if(any(f_obs < 5)){
          start_pool_bin = which(f_obs < 5)[1]
          f_obs_pool = sum(f_obs[start_pool_bin:length(f_obs)])
          f_hyp_pool = sum(f_hyp[start_pool_bin:length(f_hyp)])
          
          f_obs = c(f_obs[1:(start_pool_bin - 1)], f_obs_pool)
          f_hyp = c(f_hyp[1:(start_pool_bin - 1)], f_hyp_pool)
        }
        
        p = 2
        df = (length(bins) + 1) - p
        chiSquare = sum((f_obs-f_hyp)^2/f_hyp)

        p_values = c(p_values, 1 - pchisq(chiSquare, df = df))
        chi_squares = c(chi_squares, chiSquare)
        dfs = c(dfs, df)}, 
        error = function(e){})
      if("NULL" %in% class(t)){
        ind = c(ind, i)
      }
    }
  }
  
  return(list(p_values, ind, chi_squares, dfs))
}

# Bin by sample
p_chisq_test_2 = function(m, distribution = "poisson"){
  total_sum = sum(m)
  lambda_j = rowSums(m)/total_sum # proportion of counts in gene j
  c_i = colSums(m)/total_sum # proportion of counts in sample i (should be approximately same across i for downsampled matrix)
  mu_ij = outer(lambda_j, c_i)*total_sum
  
  f_obs = m
  f_hyp = mu_ij
  
  if(distribution == "poisson"){
    chi_square = rowSums((f_obs-f_hyp)^2/f_hyp)
  } else if(distribution == "nb 1"){ # single overdispersion parameter
    # estimate phi
    # Option 1 (regression)
    # means = rowMeans(m)
    # vars = apply(m, MARGIN = 1, var)
    # model = lm(vars ~ 1*means + I(means^2) + 0, tibble(means, vars))
    # phi = 1/coef(model)["I(means^2)"]
    
    # Option 2 (edgeR)
    phi = 1/edgeR::estimateCommonDisp(m)
    
    f_var = mu_ij + mu_ij^2/phi
    chi_square = rowSums((f_obs-f_hyp)^2/f_var)
  } else if(distribution == "nb 2"){
    # Option 1
    # nb_var = function(x){
    #   estimates = suppressWarnings(fitdistr(x, "negative binomial")$estimate)
    #   var = estimates['mu'] + estimates['mu']^2/estimates['size']
    #   return(var)
    # }
    # chi_square = rowSums((f_obs-f_hyp)^2)/apply(m, 1, var)
    
    # Option 2
    # nb_phi = function(x){
    #   suppressWarnings(try({estimates = fitdistr(x, "negative binomial")$estimate}, silent = TRUE))
    #   if(!exists("estimates")) return(NA)
    #   size = estimates['size']
    #   return(size)
    # }
    # f_phi = as.numeric(apply(m, 1, nb_phi))

    # Option 3 (edgeR)
    f_phi = 1/edgeR::estimateDisp(m)$tagwise.dispersion
    
    # Option 4 (glmGamPoi)
    # f_phi = 1/glmGamPoi::glm_gp(m, design = ~ 1, size_factors = "deconvolution", overdispersion = TRUE)$overdispersions
    
    remove_na_rows = which(!is.na(f_phi))
    f_var = mu_ij[remove_na_rows, ] + mu_ij[remove_na_rows, ]^2/f_phi[remove_na_rows]
    
    chi_square = rowSums((f_obs[remove_na_rows, ]-f_hyp[remove_na_rows, ])^2/f_var)
    return(list(chi_square, f_phi)) # diagnosing purposes
  } 
  
  return(chi_square)
}



