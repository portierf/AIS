#################
# Rcode for the NIPS paper entitled "Asymptotic optimality of adaptive importance samling"
# by Bernard Delyon and François Portier
# This program gives the simulation results of Figure 1 and 2 for AIC and wAIS with : 
#      Dimension................................... p = 8 
#      Number of iterations........................ T = 20    
#      Total number of requests to the integrand... N = 100000  
# Other dimensions and number of iterations can be obtained
# by direcly changing these values.
# A point of the figures is actually obtained by computing the Euclidean distance between the estimate and 
#the true mean (5,...,5); repeating this 100 times and taking the average of the results

rm(list = ls())
library(mnormt)
library(mvtnorm)

# Aim is to compute the mean of a Gaussian with
# mean (5,...,5) and diagonal variance (1,...,1). 
# The dimension is 8.
T = 5                           # Number of iterations. 
p = 16                          # dimension
N = 100000                      # number of requests to the integrand
Nfix= floor(N / T)              # Fixed number for simulations at each step (nt in the paper) 
mu_true = rep(5,p)              # mean
cov_true = diag( rep(1,p))      # covariance matrix
ppi = function(x) {
  dmvnorm(x , mean = mu_true, sigma = cov_true)
}
################### params of AIS and wAIS #################### 
# The sampling policy is selected among the family of student
# distributions with fixed degree of freedom df and shape matrix Sig_ais :
df = 3
Sig_ais = diag( 5 , p, p) * (df - 2) / df
# The mean will be updated using the GMM method (as described in Section 3.1 of the article)
# The initial chosen mean value is
mu_current_ais = rep(0,p)
###########################################
num_actu = 0   
denom_actu = 0
num_actu_stab = 0
denom_actu_stab = 0
mu_num_actu = 0
mu_denom_actu = 0

for (t in 1:T)
{
  ##################################################
  # Exploit
  ##################################################
  Nt = Nfix  # number of simulations in the step. Here fixed. nt of the paper
  ##################################################
  betamc = rmt(Nt ,mu_current_ais, Sig_ais , df = df)
  dens = dmt(betamc, mean = mu_current_ais, Sig_ais, df = df)
  w = apply(betamc, 1, ppi) / dens
  X = betamc * w
  num =  colSums(X)
  denom = sum(w)
  num_actu = num_actu + num
  denom_actu = denom_actu + denom
  mat_var_1 =  (mean((w-1)^2))^{-1}                ## equation (12) in the paper
  num_stab =  mat_var_1 * num
  denom_stab = mat_var_1 * denom
  num_actu_stab = num_actu_stab + num_stab
  denom_actu_stab = denom_actu_stab + denom_stab
  ###############################################
  # Integral estimate
  ###############################################
  # Estimate for AIS and wAIS
  AIS_estimate = num_actu / denom_actu              
  wAIS_estimate = num_actu_stab / denom_actu_stab   ## equation (13) in the paper
  ###############################################
  # Update of the distribution
  ###############################################
  mu_num_actu = mu_num_actu + num
  mu_denom_actu = mu_denom_actu + denom
  mu_current_ais = mu_num_actu / mu_denom_actu      ## equation (7) in the paper
}
cat('AIS_estimate = ',round(AIS_estimate,2),'\n')
cat('wAIS_estimate = ',round(wAIS_estimate,2),'\n')

