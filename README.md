# BDLM-AR
 
Supporting information for the manuscript “Analysis of N-of-1 trial using Bayesian distributed lag model with autocorrelated errors” by Ziwei Liao, Min Qian, Ian M. Kronish and Ying Kuen Cheung. All code is written in R and was based on R version 3.5. Required packages are listed at the beginning of the R file.

## Overview

Code is organized into two folders, “Methods” and “Simulation”.
The “Methods” folder contains the code to run the proposed method using a hybrid MH/Gibbs algorithm. Other related functions are also included.

1. **Methods**: Script “BDLM-AR Methods main.R” contains the main function "MCMC.sampler" to run the proposed method and related function to generate summary statistics of samples from posterior distribution. Script “BDLM-AR Methods supporting.R” contains the supporting function used in main function "MCMC.sampler".

2. **Simulation**: The “Simulation” folder contains code to run the simulations in Section 4 and also plots of simulation results.



### Contact
Ziwei Liao - ziwei.liao.fdu@gmail.com
