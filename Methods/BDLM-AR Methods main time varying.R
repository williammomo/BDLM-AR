###################################################
#                                                 #
#              Methods - main function            #
#                                                 #
###################################################

source("BDLM-AR Methods supporting.R")

library(MASS)
library(doParallel)

##' Hybrid MH/Gibbs sampler function
##' Create function to perform hybrid MH/Gibbs sampler for each parameter
##' 
##' @slot x the vector of the treatment indicator 
##' @slot y the outcome (n*1)
##' @slot z the matrix of the time varying covariates
##' @slot L the lag length of distributed lag model
##' @slot p order of autoregressive error
##' @slot z_dim the number of time varying covariates
##' @slot beta0 prior mean of regression coefficient (k*1)
##' @slot a0 prior shape parameter of IG distribution (sigma2)
##' @slot b0 prior scale parameter of IG distribution (sigma2)
##' @slot phi0 prior mean of Normal distribution of autoregressive coefficient (p*1)
##' @slot iter iteration number
##' @slot burnin burnin number
##' @slot seed set a given seed

MCMC.sampler=function(x, y, z=NULL, L, p, z_dim=0, beta0=rep(0, L+2+z_dim), a0=0, b0=0, iter=100, burnin=round(iter/2), seed=NULL){
  x = BDLM_design_matrix(x, z, L)
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  sigma20=runif(1,100,10000)
  # sigma20=1000
  gamma0=runif(1,0.1,1)
  # gamma0=1
  
  beta0=rep(0, L+2+z_dim)
  phi0=rep(0,p)
  
  beta.matrix=matrix(NA, nrow=L+2+z_dim, ncol=iter)
  beta.matrix[,1]=runif(L+2+z_dim, -2, 2)
  
  phi.matrix=matrix(NA, nrow=p, ncol=iter)
  phi.matrix[,1]=runif(p,-0.2,0.2)
  
  sigma2.vec=rep(NA,iter)
  sigma2.vec[1]=sigma20
  
  gamma1.vec=rep(NA,iter)
  gamma1.vec[1]=gamma0
  
  gamma2.vec=rep(NA,iter)
  gamma2.vec[1]=gamma0
  
  mh.vec=rep(0,iter)
  
  Lambda0=create.Lambda.matrix(L=L, q=z_dim, gamma1=gamma0, gamma2=gamma0, var_mu=0.01)
  
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
      0.5*log(det(1/sigma2.vec[i]*(create.Lambda.matrix(L=L, q=z_dim, gamma1=g[1], gamma2=g[2], var_mu=0.01)))) - 0.5*t(beta.matrix[,i]) %*% (create.Lambda.matrix(L=L, q=z_dim, gamma1=g[1], gamma2=g[2], var_mu=0.01)/sigma2.vec[i]) %*% (beta.matrix[,i]) + log(dexp(g[1])) + log(dexp(g[2]))
    }
    
    current_gamma = c(gamma1.vec[i-1], gamma2.vec[i-1])
    proposed_gamma = current_gamma + rnorm(2, mean=0, sd=0.2)
    proposed_det = det(1/sigma2.vec[i]*(create.Lambda.matrix(L=L, q=z_dim, gamma1=proposed_gamma[1], gamma2=proposed_gamma[2], var_mu=0.01)))
    if (proposed_det<=0){
      gamma1.vec[i] = current_gamma[1]
      gamma2.vec[i] = current_gamma[2]      
    } else {
      A = logtarget(proposed_gamma)-logtarget(current_gamma)
      if (log(runif(1))<A & all(proposed_gamma>0) & all(proposed_gamma<5)){
        #if(log(runif(1))<A & all(proposed_gamma>0)){      
        gamma1.vec[i] = proposed_gamma[1]
        gamma2.vec[i] = proposed_gamma[2]
        mh.vec[i] = 1
      } else {
        gamma1.vec[i] = current_gamma[1]
        gamma2.vec[i] = current_gamma[2]
      }      
    }
    
    Lambda0=create.Lambda.matrix(L=L, q=z_dim, gamma1=gamma1.vec[i], gamma2=gamma2.vec[i], var_mu=0.01)
    
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
  return(list(beta=beta.matrix[,(burnin+1):iter], sigma2=sigma2.vec[(burnin+1):iter], phi=phi.matrix[,(burnin+1):iter], gamma1=gamma1.vec[(burnin+1):iter], gamma2=gamma2.vec[(burnin+1):iter], acceptance=mh.vec[(burnin+1):iter]))
}

# Test
set.seed(0)
x_test = rbinom(100, size = 1, prob = 0.5)
y_test = rnorm(100, mean = 10, sd = 2)
z_test = matrix(rnorm(100*3), ncol=3)
output_test = MCMC.sampler(x_test, y_test, z=z_test, L=7, p=7, z_dim=3, iter = 1e4)

##' Create function to generate summary statistics of samples from posterior distribution
##' 
##' @slot posterior.data the samples from posterior distribution
##' @slot quantile.prob the summary quantile of posterior samples
##' @slot z_dim the number of time varying covariates

summary.posterior=function(posterior.data, quantile.prob = c(0.025, 0.5, 0.975), z_dim = 0, n_chain = 1){
  if (n_chain > 1){
    posterior.data = Reduce(function(x,y) Map(array.bind, x, y), posterior.data)
  }
  res = NULL
  n_variable = length(posterior.data)
  temp = data.frame(t(posterior.data[[1]]))
  posterior.data[[n_variable+1]] = apply(temp, 1, sum) - temp[,1]
  posterior.data[[n_variable+2]] = apply(temp, 1, sum) - temp[,1] - temp[,2]
  if (z_dim>0){
    for(i in 1:z_dim){
      posterior.data[[n_variable+1]] = posterior.data[[n_variable+1]] - temp[,(i+1)]
      posterior.data[[n_variable+2]] = posterior.data[[n_variable+2]] - temp[,(i+2)]
    }
  }
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
        row.names(summary_table) = paste0(names(posterior.data[i]), seq(-1-z_dim, nrow(posterior.data[[i]])-2-z_dim))
        row.names(summary_table)[1] = "mu"
        if (z_dim>0){
            row.names(summary_table)[2:(1+z_dim)] = paste0("b", seq(z_dim))
          }
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
output_test = foreach (i=1:5) %dopar% {
  MCMC.sampler(x_test, y_test, z=z_test, L=7, beta0=rep(0, 12), p=7, z_dim=3, iter = 1e4)
}
summary_test = summary.posterior(output_test, z_dim=3, n_chain=5)

##' Create function to calculate Euclidean distance between true and posterior mean estimated DL curve
##' 
##' @slot summary.posterior.data the summarized results from posterior distribution
##' @slot fit_lag the lag length of working model
##' @slot truth_beta the true lag coefficients

posterior_mean_Euclidean_BDLMAR = function(summary.posterior.data, fit_lag, truth_beta){
  len_truth = length(truth_beta)
  mybeta = summary.posterior.data[,1][startsWith(rownames(summary.posterior.data), "beta")]
  len_coef = length(mybeta)
  
  if (len_coef<len_truth) mybeta = append(mybeta, rep(0, len_truth-len_coef))
  if (len_coef>len_truth) mybeta = mybeta[1:len_truth]
  
  distance = sqrt(sum((truth_beta - mybeta)^2))
  return(distance)
}

# Test
posterior_mean_Euclidean_BDLMAR(summary_test, fit_lag = 7, truth_beta = rep(0,7))

##' Create function to calculate Euclidean distance between true and posterior median estimated DL curve
##' 
##' @slot summary.posterior.data the summarized results from posterior distribution
##' @slot fit_lag the lag length of working model
##' @slot truth_beta the true lag coefficients

posterior_median_Euclidean_BDLMAR = function(summary.posterior.data, fit_lag, truth_beta){
  len_truth = length(truth_beta)
  mybeta = summary.posterior.data[,4][startsWith(rownames(summary.posterior.data), "beta")]
  len_coef = length(mybeta)
  
  if (len_coef<len_truth) mybeta = append(mybeta, rep(0, len_truth-len_coef))
  if (len_coef>len_truth) mybeta = mybeta[1:len_truth]
  
  distance = sqrt(sum((truth_beta - mybeta)^2))
  return(distance)
}

# Test
posterior_median_Euclidean_BDLMAR(summary_test, fit_lag = 7, truth_beta = rep(0,7))

##' Create function to implement the Gelman–Rubin convergence diagnostic
##' 
##' @slot posterior.data the samples from posterior distribution
##' @slot z_dim the number of time varying covariates
##' 
gr.diag = function(posterior.data, z_dim = 0){
  J = length(posterior.data)
  L = min(lengths(posterior.data[[1]]))
  n_beta = NROW(posterior.data[[1]][[1]])
  var_names = names(posterior.data[1])
  mean_res = NULL
  var_res = NULL
  
  for(j in 1:J){
    temp = data.frame(t(posterior.data[[j]][[1]]))
    posterior.data[[j]][[1]] = rbind(posterior.data[[j]][[1]], apply(temp, 1, sum) - temp[,1], apply(temp, 1, sum) - temp[,1] - temp[,2])
    if (z_dim>0){
      for(i in 1:z_dim){
        posterior.data[[j]][[1]][n_beta+1,] = posterior.data[[j]][[1]][n_beta+1,] - temp[,(i+1)]
        posterior.data[[j]][[1]][n_beta+2,] = posterior.data[[j]][[1]][n_beta+2,] - temp[,(i+2)]
      }
    }
    temp = data.frame(t(posterior.data[[j]][[2]]))
    posterior.data[[j]] = append(posterior.data[[j]], list(sigma = as.vector(apply(temp, 1, sqrt))))
    mean_res = rbind(mean_res, unlist(lapply(posterior.data[[j]], array.mean)))
    var_res = rbind(var_res, unlist(lapply(posterior.data[[j]], array.var)))
  }
  B = apply(mean_res,2,var)
  W = apply(var_res,2,mean)
  R = sqrt((L-1)/L+B/W)
  if(z_dim>0){
    names(R)[1:(n_beta+2)] = c("mu", paste0("b", seq(1:z_dim)), paste0("beta", seq(0, n_beta-2-z_dim)), "sumBeta", "delay")
  } else{
    names(R)[1:(n_beta+2)] = c("mu", paste0("beta", seq(0, n_beta-2)), "sumBeta", "delay")
  }
  return(R)
}

# Test
gr_test = gr.diag(output_test, z_dim=3)

##' Create function to implement the stable Gelman–Rubin convergence diagnostic
##' 
##' @slot posterior.data the samples from posterior distribution
##' 
library(stableGR)
stable.GR.transformer = function(posterior.data){
  J = length(posterior.data)
  n_matrix = length(posterior.data[[1]])
  res = list()
  for(j in 1:J){
    temp = NULL
    for(i in 1:n_matrix){
      if(is.matrix(posterior.data[[j]][[i]])){
        temp = cbind(temp,t(posterior.data[[j]][[i]]))
      } else{
        temp = cbind(temp,posterior.data[[j]][[i]])
      }
    }
    res = append(res,list(temp))
  }
  return(res)
}

# Test
stable_GR_test = stable.GR(stable.GR.transformer(output_test), multivariate = TRUE)$psrf
