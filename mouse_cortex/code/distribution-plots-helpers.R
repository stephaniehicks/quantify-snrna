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
    ggplot(aes(x = log10(means), y = probs_0)) +
    geom_point(data = plot_dt, aes(x = log(means), y = emp_probs_0), alpha = 0.4) + # Add data points
    geom_line(aes(color = model),
              size = 1) + # Add lines for models
    labs(x = "Log of mean expression",
         y = "Fraction of zeros droplets") +
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

