###################################################
#                                                 #
#           Methods - supporting function         #
#                                                 #
###################################################


##' Create function to perform a change of variable: calculate phi(L)*x
##' 
##' @slot phi the coefficient of autoregressive error
##' @slot input the continuous vector of one covariate

lag = function(phi, input){
  p = length(phi)
  input_lagsum = rep(0, (NROW(input)-p))
  input_star = rep(0, (NROW(input)-p))
  for (i in 1:(NROW(input)-p)){
    s=1
    while ((p+i-s>0) & s<=length(phi)) {input_lagsum[i] = input_lagsum[i]+input[p+i-s]*phi[s]; s=s+1}
    input_star[i]=input[p+i]-input_lagsum[i]
  }
  return(input_star)
}

# Example to calculate phi(L)*x, when x is matrix
y = rnorm(100)
b = t(data.frame(rbind(y, y-1)))
b = as.matrix(b)
c = apply(b, 2, function(b) lag(phi=phi, input=b))

##' Create function to get E matrix, which is used in the conditional posterior distribution of Phi
##' 
##' @slot beta the coefficient of distributed lag curve
##' @slot x,y the original predictor and outcome
##' @slot p the order of autoregressive process

ematrix=function(x, y, beta, p){
  epsilon=y-as.matrix(x)%*%beta
  out_matrix=matrix(NA, nrow=NROW(x)-p, ncol=p)
  for (i in 1:(NROW(x)-p)){
    out_matrix[i,]=epsilon[(i+p-1):(i+p-p)]
  }
  return(out_matrix)
}

##' Create function to check the the roots of the characteristic polynomial, return True for accept sample (roots are all outside unit circle) and False for reject sample
##' 
##' @slot phi the autoregressive coefficients

root.check=function(phi){
  root=polyroot(c(1,-phi))
  return(all(abs(root)>1))
}

##' Create function to generate prior precision matrix on lag coefficient beta
##' 
##' @slot L the lag length of distributed lag model
##' @slot gamma1 the increasing rate of ridge penalty
##' @slot gamma2 the increasing rate of smoothness penalty
##' @slot var_mu the inverse is prior variance on the intercept mu

create.Lambda.matrix=function(L, gamma1, gamma2, var_mu){
  Lambda=matrix(0, nrow=L+1, ncol=L+1)
  lambda.diag = exp(gamma1*seq(L+1))-1
  lambda.tridiag = exp(gamma2*seq(L+1))-1
  Lambda[1, 1] = lambda.diag[1] + lambda.tridiag[1];
  Lambda[L+1,L+1] = (lambda.diag[L+1]+lambda.tridiag[L]+lambda.tridiag[L+1]);
  if (L>=2){
    for (i in 2:L){
      Lambda[i,i] = (lambda.diag[i]+lambda.tridiag[i-1]+lambda.tridiag[i]);
    }
  }
  for (j in 1:L){
    Lambda[j,j+1] = -lambda.tridiag[j];
    Lambda[j+1,j] = -lambda.tridiag[j];
  }
  Lambda0=cbind(c(var_mu, rep(0, L+1)), rbind(rep(0, L+1), Lambda))
  return(Lambda0)
}

##' Create function to implement DETGTRI algorithm to check whether the precision matrix is positive definite
##' 
##' @slot input the precision matrix

checkTridiag=function(input){
  n = ncol(input)
  d = diag(input)
  b = rep(NA,n-1)
  c = rep(NA,n)
  c[1] = d[1]
  for (i in 1:(n-1)){
    b[i] = input[i,i+1]
    c[i+1] = d[i+1] - b[i]^2/c[i]
  }
  return(all(c>0))
}

##' Create function to generate design matrix based on the treatment sequence
##' 
##' @slot x the vector of the treatment indicator
##' @slot L the lag length of distributed lag curve

BDLM_design_matrix = function(x, L){
  if (L<1) {stop("Lag length must be >= 1")}
  n = length(x)
  temp_matrix = cbind(rep(1,n), x)
  for(i in 1:L){
    lag_x = c(rep(0,i), x[1:(n-i)])
    temp_matrix = cbind(temp_matrix, lag_x)
    colnames(temp_matrix) = NULL
  }
  return(temp_matrix)
}
