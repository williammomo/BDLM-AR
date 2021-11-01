###################################################
#                                                 #
#              Methods - main function            #
#                                                 #
###################################################

source("BDLM-AR Methods supporting.R")

library(MASS)

##' Hybrid MH/Gibbs sampler function
##' Create function to perform hybrid MH/Gibbs sampler for each parameter
##' 
##' @slot x the design matrix of covariates including intercept (n*k, k=L+2)
##' @slot x the vector of the treatment indicator 
##' @slot y the outcome (n*1)
##' @slot L the lag length of distributed lag model
##' @slot p order of autoregressive error
##' @slot beta0 prior mean of regression coefficient (k*1)
##' @slot a0 prior shape parameter of IG distribution (sigma2)
##' @slot b0 prior scale parameter of IG distribution (sigma2)
##' @slot phi0 prior mean of Normal distribution of autoregressive coefficient (p*1)
##' @slot sigma20 initial value of sigma2
##' @slot gamma0 initial value of hyperparameter gamma in precision matrix 
##' @slot iter iteration number
##' @slot burnin burnin number

MCMC.sampler=function(x, y, L, p, beta0=rep(0, L+2), a0=0, b0=0, phi0=rep(0,p), sigma20=1000, gamma0=1, iter=100, burnin=round(iter/2)){
  x = BDLM_design_matrix(x, L)
  
  beta.matrix=matrix(NA, nrow=L+2, ncol=iter)
  beta.matrix[,1]=beta0
  
  phi.matrix=matrix(NA, nrow=p, ncol=iter)
  phi.matrix[,1]=phi0
  
  sigma2.vec=rep(NA,iter)
  sigma2.vec[1]=sigma20
  
  gamma1.vec=rep(NA,iter)
  gamma1.vec[1]=gamma0
  
  gamma2.vec=rep(NA,iter)
  gamma2.vec[1]=gamma0
  
  Lambda0=create.Lambda.matrix(L=L, gamma1=gamma0, gamma2=gamma0, var_mu=0.01)
  
  Psi0=1/200*diag(p) #Prior precision matrix of Normal distribution of polynomial in the lag operator (p*p)
  
  for(i in 2:iter) {
    if ((10*i) %% iter == 0){
      cat('\n',100*i/iter, "% of iterations completeted!")
    }
    
    #Backshift transformation of X and Y
    y.star=lag(phi=phi.matrix[,i-1], input=y)
    x.star=apply(x, 2, function(a) lag(phi=phi.matrix[,i-1], input=a))
    
    #Conditional posterior of beta
    beta.n=(chol2inv(chol(Lambda0 + t(x.star) %*% x.star))) %*% (Lambda0 %*% beta0 + t(x.star) %*% y.star)
    Lambda.n=(Lambda0 + t(x.star) %*% x.star)
    beta.matrix[,i]=mvrnorm(n=1, mu=beta.n, Sigma=sigma2.vec[i-1]*(chol2inv(chol(Lambda.n))))
    
    #Conditional posterior of sigma2
    q=(t(beta.matrix[,i]-beta0)) %*% Lambda0 %*% (beta.matrix[,i]-beta0)
    d=t(y.star-as.matrix(x.star)%*%beta.matrix[,i])%*%(y.star-as.matrix(x.star)%*%beta.matrix[,i])
    sigma2.vec[i]=1/rgamma(n=1, shape=(a0+NROW(x)-p+NCOL(x))/2, rate=(b0+d+q)/2)
    
    #Metropolis Hasting algorithm to sample Bayesian ridge parameter lambda
    logtarget=function(g){
      0.5*log(det(1/sigma2.vec[i]*(create.Lambda.matrix(L=L, gamma1=g[1], gamma2=g[2], var_mu=0.01)))) - 0.5*t(beta.matrix[,i]) %*% (create.Lambda.matrix(L=L, gamma1=g[1], gamma2=g[2], var_mu=0.01)/sigma2.vec[i]) %*% (beta.matrix[,i]) + log(dexp(g[1])) + log(dexp(g[2]))
    }
    
    current_gamma = c(gamma1.vec[i-1], gamma2.vec[i-1])
    proposed_gamma = current_gamma + rnorm(2, mean=0, sd=0.2)
    proposed_det = det(1/sigma2.vec[i]*(create.Lambda.matrix(L=L, gamma1=proposed_gamma[1], gamma2=proposed_gamma[2], var_mu=0.01)))
    if (proposed_det<=0 | checkTridiag(create.Lambda.matrix(L=L, gamma1=proposed_gamma[1], gamma2=proposed_gamma[2], var_mu=0.01)) == FALSE){
      gamma1.vec[i] = current_gamma[1]
      gamma2.vec[i] = current_gamma[2]      
    } else {
      A = logtarget(proposed_gamma)-logtarget(current_gamma)
      if (log(runif(1))<A & all(proposed_gamma>0) & all(proposed_gamma<5)){
        #if(log(runif(1))<A & all(proposed_gamma>0)){      
        gamma1.vec[i] = proposed_gamma[1]
        gamma2.vec[i] = proposed_gamma[2]
      } else {
        gamma1.vec[i] = current_gamma[1]
        gamma2.vec[i] = current_gamma[2]
      }      
    }
    
    Lambda0=create.Lambda.matrix(L=L, gamma1=gamma1.vec[i], gamma2=gamma2.vec[i], var_mu=0.01)
    
    #Conditional posterior of Phi
    e.hat=ematrix(x, y, beta.matrix[,i], p)
    epsilon=y-as.matrix(x)%*%beta.matrix[,i]
    Psi.n=Psi0+t(e.hat)%*%e.hat/sigma2.vec[i]
    
    #Note: In t(E)%*%epsilon, the first p elements in epsilon is removed
    phi.n=chol2inv(chol(Psi.n))%*%(Psi0%*%phi0+t(e.hat)%*%epsilon[(p+1):length(y)]/sigma2.vec[i])
    phi.matrix[,i]=mvrnorm(n=1, mu=phi.n, Sigma=chol2inv(chol(Psi.n)))
    
    #Add restriction indicator
    while (root.check(phi.matrix[,i])==FALSE) {phi.matrix[,i]=mvrnorm(n=1, mu=phi.n, Sigma=chol2inv(chol(Psi.n)))}
  }
  return(list(beta=beta.matrix[,(burnin+1):iter], sigma2=sigma2.vec[(burnin+1):iter], phi=phi.matrix[,(burnin+1):iter], gamma1=gamma1.vec[(burnin+1):iter], gamma2=gamma2.vec[(burnin+1):iter]))
}

# Test
set.seed(0)
x_test = rbinom(100, size = 1, prob = 0.5)
y_test = rnorm(100, mean = 10, sd = 2)
output_test = MCMC.sampler(x_test, y_test, L=7, beta0=rep(0, 9), p=7, phi0=rep(0,7), iter = 1e4)

##' Create function to generate summary statistics of samples from posterior distribution
##' 
##' @slot posterior.data the samples from posterior distribution
##' @slot quantile.prob the summary quantile of posterior samples

summary.posterior=function(posterior.data, quantile.prob = c(0.025, 0.5, 0.975)){
  res = NULL
  n_variable = length(posterior.data)
  temp = data.frame(t(posterior.data[[1]]))
  posterior.data[[n_variable+1]] = apply(temp, 1, sum) - temp[,1]
  posterior.data[[n_variable+2]] = apply(temp, 1, sum) - temp[,1] - temp[,2]
  temp = data.frame(t(posterior.data[[2]]))
  posterior.data[[n_variable+3]] = as.vector(apply(temp, 1, sqrt))
  names(posterior.data)[c(n_variable+1, n_variable+2, n_variable+3)] = c("sumBeta", "delay", "sigma")
  posterior.data = lapply(posterior.data, as.matrix)
  
  for (i in seq(n_variable+3)){
    temp = data.frame(t(posterior.data[[i]]))
    if (nrow(temp) == 1){
      temp = as.data.frame(t(temp))
    }
    mean = sapply(temp, mean)
    sd = sapply(temp, sd)
    q = sapply(temp, quantile, probs = quantile.prob)
    summary_table = data.frame(mean, sd, t(q))
    colnames(summary_table) = c("mean", "sd", paste0(quantile.prob*100, "%"))
    if (nrow(summary_table) > 1){
      if (names(posterior.data[i]) == "beta"){
        row.names(summary_table) = paste0(names(posterior.data[i]), seq(-1, nrow(posterior.data[[i]])-2))
        row.names(summary_table)[1] = "mu"
      } else {
        row.names(summary_table) = paste0(names(posterior.data[i]), seq(1, nrow(posterior.data[[i]])))
      }
    } else {
      row.names(summary_table) = names(posterior.data[i])
    }
    res = rbind(res, summary_table)
  }
  return (res)
}

# Test
summary_test = summary.posterior(output_test)

##' Create function to calculate Euclidean distance between true and posterior mean estimated DL curve
##' 
##' @slot summary.posterior.data the summarized results from posterior distribution
##' @slot fit_lag the lag length of working model
##' @slot truth_beta the true lag coefficients

posterior_mean_Euclidean_BDLMAR = function(summary.posterior.data, fit_lag, truth_beta){
  len_truth = length(truth_beta)
  mybeta = summary.posterior.data[,1][2:(fit_lag+2)]
  len_coef = length(mybeta)
  
  if (len_coef<len_truth) mybeta = append(mybeta, rep(0, len_truth-len_coef))
  if (len_coef>len_truth) mybeta = mybeta[1:len_truth]
  
  distance = sqrt(sum((truth_beta - mybeta)^2))
  return(distance)
}

# Test
posterior_mean_Euclidean_BDLMAR(summary_test, fit_lag = 7, truth_beta = rep(0,7))
