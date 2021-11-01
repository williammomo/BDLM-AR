###################################################
#                                                 #
#              Simulation - main function         #
#                                                 #
###################################################

source("BDLM-AR Methods supporting.R")
source("BDLM-AR Methods main.R")

library(doParallel)
nCores = detectCores(logical = F) # detect physical cores
registerDoParallel(nCores)

##--------------------------------------------------------------------
## Function to conduct simulation under specific simulation parameters
##--------------------------------------------------------------------
simulation_BDLMAR = function(fit_lag, truth_beta, mu, phi, sigma, design, sample_size, index){
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
  
  # Fit BDLAR model    
  temp_model = MCMC.sampler(x=x, y=y, p=1, a0=0, b0=0, Psi0=1e-4*diag(1), gamma0=0.1, sigma20=1, iter=50000)
  summary_temp_model=summary.posterior(posterior.data = temp_model, quantile.prob = c(0.025, 0.5, 0.975))
  
  sumBeta_mean = summary_temp_model['sumBeta',1]
  delay_mean = summary_temp_model['delay',1]
  phi_mean = summary_temp_model['phi',1]
  sigma_mean = summary_temp_model['sigma',1]
  sce_mean = posterior_mean_sce_BDLAR(summary.posterior.data = summary_temp_model, fit_lag = fit_lag, truth_beta = truth_beta)
  mu_mean = summary_temp_model[1,1]
  beta_mean = summary_temp_model[2:(fit_lag+2),1]
  file = "Bayesian_ARMA_Gibbs_TwoParameters"
  result_mean = c(mu_mean, sumBeta_mean, beta_mean, phi_mean, sigma_mean, delay_mean, sce_mean, fit_lag, file)
  names(result_mean) = c("mu", "sumBeta", sprintf("beta%s", c(0:fit_lag)), "phi", "sigma", "delay", "SCE", "lag", "file")
  
  sumBeta_median = summary_temp_model['sumBeta',4]
  delay_median = summary_temp_model['delay',4]
  phi_median = summary_temp_model['phi',4]
  sigma_median = summary_temp_model['sigma',4] 
  sce_median = posterior_median_sce_BDLAR(summary.posterior.data = summary_temp_model, fit_lag = fit_lag, truth_beta = truth_beta)
  mu_median = summary_temp_model[1,4]
  beta_median = summary_temp_model[2:(fit_lag+2),4]
  file = "Bayesian_ARMA_Gibbs_TwoParameters"
  result_median = c(mu_median, sumBeta_median, beta_median, phi_median, sigma_median, delay_median, sce_median, fit_lag, file)
  names(result_median) = c("mu", "sumBeta", sprintf("beta%s", c(0:fit_lag)), "phi", "sigma", "delay", "SCE", "lag", "file")
  
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
  
  return(list(result_mean, result_median, result_cvg, result_covernull, result_range))
}

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

simulation_summary_BDLMAR = function(fit_lag, truth_beta, mu, phi, sigma, design, sample_size, iter = 100){
  
  result_combine_mean = NULL #collect posterior mean
  result_combine_median = NULL #collect posterior median
  result_combine_cvg = NULL #collect whether 95% posterior credible interval cover the truth
  result_combine_range = NULL #collect the width of 95% posterior credible interval
  result_combine_covernull = NULL #collect whether 95% posterior credible interval cover 0
  
  foreach (i = 1:iter, .combine = mycombine) %dopar%{
    output = simulation_BDLAR(fit_lag, truth_beta, mu, phi, sigma, design, sample_size, index = 100+i)
  }
  
  # Switch to regular loop
  # for (i in 1:iter){
  #   cat('\n','Simulation round', i, "is running")
  #   output = simulation_BDLAR(fit_lag, truth_beta, mu, sigma, phi, design, sample_size, index = i)
  #   result_combine_mean = rbind(result_combine_mean, output[[1]])
  #   result_combine_median = rbind(result_combine_median, output[[2]])
  #   result_combine_cvg = rbind(result_combine_cvg, output[[3]])
  #   result_combine_covernull = rbind(result_combine_covernull, output[[4]])
  #   result_combine_range = rbind(result_combine_range, output[[5]])
  # }
  # return(list(result_combine_mean, result_combine_median, result_combine_cvg, result_combine_covernull, result_combine_range))
}

#############################
# End of parallel computing #
#############################



##-----------------------------------------------------------
## Function to generate metrics for single simulation object
##-----------------------------------------------------------

simulation_metrics_BDLMAR = function(data, truth_all){
  truth_lag_length = tail(grep('alpha', names(truth_all), value = T), n=1)
  truth_lag_length = as.numeric(gsub("[^0-9.]", "",  truth_lag_length)) - 2
  truth_AR_length = tail(grep('phi', names(truth_all), value = T), n=1)
  truth_AR_length = as.numeric(gsub("[^0-9.]", "",  'truth_AR_length'))
  
  temp_name = names(truth_all)
  colnames(data[[1]])[c(1:(ncol(data[[1]])-3))] = temp_name
  colnames(data[[2]])[c(1:(ncol(data[[2]])-3))] = temp_name
  
  bias_mean = apply(data[[1]], 2, as.numeric) %>%
    colMeans()
  bias_mean = bias_mean[1:(length(bias_mean)-3)]
  
  bias_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(bias_mean)))
  bias_union[is.na(bias_union)] = 0
  bias_mean = unlist((bias_union[2,] - bias_union[1,])[1,])
  
  
  bias_median = apply(data[[2]], 2, as.numeric) %>%
    colMeans()
  bias_median = bias_median[1:(length(bias_median)-3)]
  
  bias_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(bias_median)))
  bias_union[is.na(bias_union)] = 0
  bias_median = unlist((bias_union[2,] - bias_union[1,])[1,])
  
  MSE_mean = apply(data[[1]], 2, as.numeric)[,c(1:(ncol(data[[1]])-3))] %>%
    myMSE(truth_all = truth_all)
  MSE_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(MSE_mean)))
  MSE_union[is.na(MSE_union)] = 0
  MSE_mean = unlist((MSE_union[2,])[1,])
  
  MSE_median = apply(data[[2]], 2, as.numeric)[,c(1:(ncol(data[[2]])-3))] %>%
    myMSE(truth = truth_all)
  MSE_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(MSE_median)))
  MSE_union[is.na(MSE_union)] = 0
  MSE_median = unlist((MSE_union[2,])[1,])
  
  CVG = colMeans(data[[3]])
  names(CVG) = truth_all
  
  CoverNull = colMeans(data[[4]])
  names(CoverNull) = truth_all
  
  width = colMeans(data[[5]])
  names(width) = temp_name
  width_union = dplyr::bind_rows(as.data.frame(t(truth_all)), as.data.frame(t(width)))
  width_union[is.na(width_union)] = 0
  width = unlist((width_union[2,])[1,]) 
  
  return(data.frame(truth_all, bias_mean, MSE_mean, bias_median, MSE_median, CVG, CoverNull, width))
}
