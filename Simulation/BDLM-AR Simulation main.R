###################################################
#                                                 #
#              Simulation - main function         #
#                                                 #
###################################################

source("BDLM-AR Methods supporting.R")
source("BDLM-AR Methods main.R")

library(doParallel)
nCores = detectCores(logical = F) # detect physical cores
# registerDoParallel(nCores-3)
# 
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)
# stopCluster(cl)

##--------------------------------------------------------------------
## Function to conduct simulation under specific simulation parameters
##--------------------------------------------------------------------
simulation_BDLMAR = function(fit_lag, fit_AR_order, truth_beta, mu, phi, sigma, design, sample_size, time_varying=NULL, index){
  set.seed(index)
  N = sample_size
  sigma2 = sigma^2
  w = rnorm(n = N, 0, sigma)
  
  error = rep(NA, N)
  error[1] = w[1]
  
  for (j in 2:N){
    error[j] = phi * error[j-1] + w[j]
  }
  
  if (design == 'ABBABAAB'){
    nblock = 8
    x = c(rep(1, N/nblock), rep(0, N/nblock), rep(0, N/nblock), rep(1, N/nblock), rep(0, N/nblock), rep(1, N/nblock), rep(1, N/nblock), rep(0, N/nblock))
  } else if (design == 'ABAB'){
    nblock = 4
    x = c(rep(1, N/nblock), rep(0, N/nblock), rep(1, N/nblock), rep(0, N/nblock))
  } else if (design == 'ABBA'){
    nblock = 4
    x = c(rep(1, N/nblock), rep(0, N/nblock), rep(0, N/nblock), rep(1, N/nblock))
  } else if (design == 'weekly'){
    x = c(rep(rep(c(1,0,0,1), each = 7), 4), rep(1, 7), 0)
  } else if(design == 'daily'){
    x = rep(rep(c(1,0,0,1), each = 1), N/4)
  } else if(design == '6day'){
    x = rep(rep(c(1,0,0,1), each = 6), N/4/6)
  }  else if(design == '10day'){
    x = rep(rep(c(1,0,0,1), each = 10), N/4/10)
  }  else if(design == '2day'){
    x = rep(rep(c(1,0,0,1), each = 2), N/4/2)
  }
  
  if (!is.null(time_varying)){
    if (time_varying == 'linear'){
      zt = 0.3*seq(N)
    } else if (time_varying == 'weekend'){
      zt = 3*rep(c(rep(0, 5), 1, 1), length.out = N)
    } else {
      zt = rep(0, N)
    }
  } 
  
  # Generate data
  y = rep(NA, N)
  
  # Generate data from true distributed lag coefficient curve
  for (i in 8:N){
    y[i] = mu + truth_beta[1] * x[i] + truth_beta[2] * x[i-1] + truth_beta[3] * x[i-2] + truth_beta[4] * x[i-3] + truth_beta[5] * x[i-4] + truth_beta[6] * x[i-5] + truth_beta[7] * x[i-6] + truth_beta[8] * x[i-7] + error[i]
  }
  y[1] = mu + truth_beta[1] * x[1] + error[1]
  y[2] = mu + truth_beta[1] * x[2] + truth_beta[2] * x[1] + error[2]
  y[3] = mu + truth_beta[1] * x[3] + truth_beta[2] * x[2] + truth_beta[3] * x[1] + error[3]
  y[4] = mu + truth_beta[1] * x[4] + truth_beta[2] * x[3] + truth_beta[3] * x[2] + truth_beta[4] * x[1] + error[4]
  y[5] = mu + truth_beta[1] * x[5] + truth_beta[2] * x[4] + truth_beta[3] * x[3] + truth_beta[4] * x[2] + truth_beta[5] * x[1] + error[5]
  y[6] = mu + truth_beta[1] * x[6] + truth_beta[2] * x[5] + truth_beta[3] * x[4] + truth_beta[4] * x[3] + truth_beta[5] * x[2] + truth_beta[6] * x[1] + error[6]
  y[7] = mu + truth_beta[1] * x[7] + truth_beta[2] * x[6] + truth_beta[3] * x[5] + truth_beta[4] * x[4] + truth_beta[5] * x[3] + truth_beta[6] * x[2] + truth_beta[7] * x[1] + error[7]
  if (!is.null(time_varying)){
    y = y+zt
  }
  
  # Fit BDLAR model    
  temp_model = foreach (i=1:5) %dopar% {
    MCMC.sampler(x=x, y=y, L=fit_lag, p=fit_AR_order, a0=0, b0=0, iter=50000)
  }
  
  summary_temp_model=summary.posterior(posterior.data = temp_model, quantile.prob = c(0.025, 0.5, 0.975), n_chain = 5)
  
  sumBeta_mean = summary_temp_model['sumBeta',1]
  delay_mean = summary_temp_model['delay',1]
  phi_mean = summary_temp_model['phi',1]
  sigma_mean = summary_temp_model['sigma',1]
  sce_mean = posterior_mean_Euclidean_BDLMAR(summary.posterior.data = summary_temp_model, fit_lag = fit_lag, truth_beta = truth_beta)
  mu_mean = summary_temp_model[1,1]
  beta_mean = summary_temp_model[2:(fit_lag+2),1]
  b_mean = summary_temp_model[2,1]
  file = "Bayesian_ARMA_Gibbs_TwoParameters"
  acceptance_mean = summary_temp_model['acceptance',1]
  result_mean = c(mu_mean, sumBeta_mean, beta_mean, phi_mean, sigma_mean, delay_mean, sce_mean, fit_lag, file, acceptance_mean)
  names(result_mean) = c("mu", "sumBeta", sprintf("beta%s", c(0:fit_lag)), "phi", "sigma", "delay", "SCE", "lag", "file", "acceptance")

  sumBeta_median = summary_temp_model['sumBeta',4]
  delay_median = summary_temp_model['delay',4]
  phi_median = summary_temp_model['phi',4]
  sigma_median = summary_temp_model['sigma',4] 
  sce_median = posterior_median_Euclidean_BDLMAR(summary.posterior.data = summary_temp_model, fit_lag = fit_lag, truth_beta = truth_beta)
  mu_median = summary_temp_model[1,4]
  beta_median = summary_temp_model[2:(fit_lag+2),4]
  b_median = summary_temp_model[2,4]
  file = "Bayesian_ARMA_Gibbs_TwoParameters"
  acceptance_median = summary_temp_model['acceptance',4]
  result_median = c(mu_median, sumBeta_median, beta_median, phi_median, sigma_median, delay_median, sce_median, fit_lag, file, acceptance_median)
  names(result_median) = c("mu", "sumBeta", sprintf("beta%s", c(0:fit_lag)), "phi", "sigma", "delay", "SCE", "lag", "file", "acceptance")
  
  sumBeta_CI = summary_temp_model['sumBeta', c(3,5)]
  delay_CI = summary_temp_model['delay', c(3,5)]
  phi_CI = summary_temp_model['phi', c(3,5)]
  sigma_CI = summary_temp_model['sigma', c(3,5)]
  mu_CI = summary_temp_model[1, c(3,5)]
  beta_CI = summary_temp_model[2:(fit_lag+2), c(3,5)]
  result_CI = rbind(mu_CI, sumBeta_CI, beta_CI, phi_CI, sigma_CI, delay_CI)
  
  result_range = result_CI[,2] - result_CI[,1]
  
  truth_all = c(mu, sum(truth_beta), truth_beta, phi, sigma, sum(truth_beta)-truth_beta[1])
  names(truth_all) = c("mu", "sumBeta", sprintf("beta%s", c(0:(length(truth_beta)-1))), "phi", "sigma", "delay")
  
  cvg_table = data.frame(result_CI) %>%
    mutate(coef_name = rownames(result_CI)) %>%
    mutate(coef_name = sub("_CI", "", coef_name)) %>%
    right_join(data.frame(truth_all, coef_name = names(truth_all), stringsAsFactors=FALSE), by = "coef_name")
  
  result_cvg = cvg_table[,1]<=cvg_table$truth_all & cvg_table[,2]>=cvg_table$truth_all
  result_covernull = cvg_table[,1] * cvg_table[,2] < 0
  
  psrf = gr.diag(temp_model)
  stable_output = stable.GR(stable.GR.transformer(temp_model), multivariate = TRUE)
  stable_psrf = stable_output$psrf
  names(stable_psrf) = c("mu", sprintf("beta%s", c(0:fit_lag)), "sigma2", "phi", "gamma1", "gamma2")
  stable_mpsrf = stable_output$mpsrf
  
  return(list(result_mean, result_median, result_cvg, result_covernull, result_range, psrf, stable_psrf, stable_mpsrf))
}

test = simulation_BDLMAR(fit_lag=7, fit_AR_order=1, truth_beta=beta1, mu=10, phi=0.5, sigma=10, design="ABBA", sample_size=120,  index=0)

#####################################
# Parallel computing for simulation #
#####################################
mycombine = function(list1, list2){
  l = length(list1)
  for (i in seq(l)){
    list1[[i]] = rbind(list1[[i]], list2[[i]])
  }
  return (list1)
}

simulation_summary_BDLMAR = function(fit_lag, fit_AR_order, truth_beta, mu, phi, sigma, design, sample_size, time_varying=NULL, iter = 100, start = 1){
  
  result_combine_mean = NULL #collect posterior mean
  result_combine_median = NULL #collect posterior median
  result_combine_cvg = NULL #collect whether 95% posterior credible interval cover the truth
  result_combine_range = NULL #collect the width of 95% posterior credible interval
  result_combine_covernull = NULL #collect whether 95% posterior credible interval cover 0
  result_combine_psrf = NULL
  result_combine_stable_psrf = NULL
  result_combine_stable_mpsrf = NULL
  
  # foreach (i = 1:iter, .combine = mycombine) %dopar%{
  #   output = simulation_BDLMAR(fit_lag, fit_AR_order, truth_beta, mu, phi, sigma, design, sample_size, index = 100+i)
  # }
  
  # Switch to regular loop
  for (i in 1:iter){
    cat('\n', '######################################')
    cat('\n','Simulation round', i, "is running")
    cl = makePSOCKcluster(5)
    doParallel::registerDoParallel(cl)
    # clusterExport(cl, list(as.vector(lsf.str())))
    clusterExport(cl, names(as.list(.GlobalEnv)))
    clusterEvalQ(cl, {
      library(forecast)
      library(ggplot2)
      library(lubridate)
      library(tseries)
      library(MASS)
      library(invgamma)
      library(dplyr)
      library(stableGR)
      library(R.utils)
    })
    # output = simulation_BDLMAR(fit_lag=fit_lag, fit_AR_order=fit_AR_order, truth_beta=truth_beta, mu=mu, sigma=sigma, phi=phi, design=design, sample_size=sample_size, index = i)
    output = R.utils::withTimeout(simulation_BDLMAR(fit_lag=fit_lag, fit_AR_order=fit_AR_order, truth_beta=truth_beta, mu=mu, sigma=sigma, phi=phi, design=design, sample_size=sample_size, time_varying=time_varying, index=i+start), timeout = 600, onTimeout = "warning")
    if(length(output)>1){
    result_combine_mean = rbind(result_combine_mean, output[[1]])
    result_combine_median = rbind(result_combine_median, output[[2]])
    result_combine_cvg = rbind(result_combine_cvg, output[[3]])
    result_combine_covernull = rbind(result_combine_covernull, output[[4]])
    result_combine_range = rbind(result_combine_range, output[[5]])
    result_combine_psrf = rbind(result_combine_psrf, output[[6]])
    result_combine_stable_psrf = rbind(result_combine_stable_psrf, output[[7]])
    result_combine_stable_mpsrf = rbind(result_combine_stable_mpsrf, output[[8]])
    }
    parallel::stopCluster(cl)
  }
  return(list(result_combine_mean, result_combine_median, result_combine_cvg, result_combine_covernull, result_combine_range, result_combine_psrf, result_combine_stable_psrf, result_combine_stable_mpsrf))
}

#############################
# End of parallel computing #
#############################
start = Sys.time()
print(start)
test = simulation_summary_BDLMAR(fit_lag=7, fit_AR_order=1, truth_beta=beta1, mu=10, phi=0.5, sigma=10, design="ABBA", sample_size=120, iter=2)
end = Sys.time()
print(end-start)

##-----------------------------------------------------------
## Utility function to calculate RMSE
##-----------------------------------------------------------
myRMSE = function(x, truth_all){
  MSE = rep(NA, ncol(x))
  name = colnames(x)
  for (i in name){
    MSE[which(i == name)] = mean((x[,i] - truth_all[i])^2)
  }
  names(MSE) = name
  return (sqrt(MSE))
}

##-----------------------------------------------------------
## Function to generate metrics for single simulation object
##-----------------------------------------------------------

simulation_metrics_BDLMAR = function(data, truth_all){
  truth_lag_length = tail(grep('alpha', names(truth_all), value = T), n=1)
  truth_lag_length = as.numeric(gsub("[^0-9.]", "",  truth_lag_length)) - 2
  truth_AR_length = tail(grep('phi', names(truth_all), value = T), n=1)
  truth_AR_length = as.numeric(gsub("[^0-9.]", "",  'truth_AR_length'))
  
  temp_name = names(truth_all)
  colnames(data[[1]])[c(1:(ncol(data[[1]])-4))] = temp_name
  colnames(data[[2]])[c(1:(ncol(data[[2]])-4))] = temp_name
  
  bias_mean = apply(data[[1]], 2, as.numeric) %>%
    colMeans()
  bias_mean = bias_mean[1:(length(bias_mean)-4)]
  
  bias_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(bias_mean)))
  bias_union[is.na(bias_union)] = 0
  bias_mean = unlist((bias_union[2,] - bias_union[1,])[1,])
  
  
  bias_median = apply(data[[2]], 2, as.numeric) %>%
    colMeans()
  bias_median = bias_median[1:(length(bias_median)-4)]
  
  bias_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(bias_median)))
  bias_union[is.na(bias_union)] = 0
  bias_median = unlist((bias_union[2,] - bias_union[1,])[1,])
  
  RMSE_mean = apply(data[[1]], 2, as.numeric)[,c(1:(ncol(data[[1]])-4))] %>%
    myRMSE(truth_all = truth_all)
  RMSE_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(RMSE_mean)))
  RMSE_union[is.na(RMSE_union)] = 0
  RMSE_Rmean = unlist((RMSE_union[2,])[1,])
  
  RMSE_median = apply(data[[2]], 2, as.numeric)[,c(1:(ncol(data[[2]])-4))] %>%
    myRMSE(truth = truth_all)
  RMSE_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(RMSE_median)))
  RMSE_union[is.na(RMSE_union)] = 0
  RMSE_median = unlist((RMSE_union[2,])[1,])
  
  CVG = colMeans(data[[3]])
  names(CVG) = truth_all

  CoverNull = colMeans(data[[4]])
  names(CoverNull) = truth_all

  width = colMeans(data[[5]])
  names(width) = temp_name
  width_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(width)))
  width_union[is.na(width_union)] = 0
  width = unlist((width_union[2,])[1,])

  return(data.frame(truth_all, bias_mean, RMSE_mean, bias_median, RMSE_median, CVG, CoverNull, width))
  # return(data.frame(truth_all, bias_mean, MSE_mean, bias_median, MSE_median))
}

# Test
mu = 10
truth_beta = c(5, 5/2, 5/4, 5/8, 5/16, 0, 0, 0)
phi = 0.5
sigma = 10
truth_all = c(mu, sum(truth_beta), truth_beta, phi, sigma, sum(truth_beta)-truth_beta[1])
names(truth_all) = c("alpha[1]", "sum_alpha", sprintf("alpha[%s]", c(2:(2+length(truth_beta)-1))), "phi[1]", "sigma", "delay")
simulation_metrics_BDLMAR(test, truth_all = truth_all)
