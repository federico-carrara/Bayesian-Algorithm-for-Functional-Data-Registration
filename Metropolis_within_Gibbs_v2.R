library(invgamma)
library(matlib)
library(MASS)
library("readxl") 
library("readr")
library(fda)
library(coda)
library(tmvtnorm)
library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(pracma)
library(RColorBrewer)
library(knitr)

setwd("C:/Users/fede1/OneDrive - Politecnico di Milano/Progetto BAYESIANA/Codici/")
#setwd("/Users/chiarascrimieri/Desktop/OneDrive - Politecnico di Milano/Progetto BAYESIANA/Codici/")

main.dir <- getwd()
scripts.dir <- paste(main.dir, "/Scripts/", sep="")
data.dir <- paste(main.dir, "/Data/", sep="")

source(file = paste(scripts.dir, "import_dataset_Federico.R", sep=""), chdir = T)
source(file = paste(scripts.dir, "aux_fun_Matteo.R", sep=""), chdir = T)

################################################################################
################################################################################
# Fake data generator -> we want data to resamble a "two-peaks gaussian"
# Each curve is given by the mix of two gaussians which are described by 
# a couple of values of the mean and of the sd.
# The returned object is a dataset whose rows are the observations for each individual
# fake.data.generator <- function(n, means, stdevs, noise.sd, 
#                                 series.length, scale.factors)
# {
#   time.steps <- seq(0, 1, length.out = series.length)
#   df <- matrix(0, nrow = N, ncol = series.length)
#   for(i in 1:N)
#   {
#     # Generate 2-peaks gaussian
#     df[i,] <- scale.factors[i] * (dnorm(x = time.steps, mean = means[i,1], sd = stdevs[i,1]) +  
#                               dnorm(x = time.steps, mean = means[i,2], sd = stdevs[i,2]))
#     # Add random noise
#     df[i,] <- df[i,] + rnorm(n = series.length, mean = 0, sd = noise.sd)
#   }
#   
#   return(list(data = as.data.frame(df), abscissa = time.steps))
# }
# 
# # Generate values for the means and the sds
# set.seed(1234)
# N <- 30
# series_length <- 200
# scale.factors <- sample(x = 50:150, size = N)
# means <- runif(n = N, min = 0.25, max = 0.55)
# means <- cbind(means, runif(n = N, min = means+0.05, max = means+0.35))
# stdevs <- cbind(rnorm(n = N, mean = 0.08, sd = 0.01), 
#                 rnorm(n = N, mean = 0.08, sd = 0.01))
# data.obj <- fake.data.generator(n = N, means = means, stdevs = stdevs, 
#                                 noise.sd = 0.1, series.length = series_length, 
#                                 scale.factors = scale.factors)
# matplot(t(data.obj$data), type = "l")
# pat.norm.df <- data.obj$data
# time_steps <- data.obj$abscissa
# # for(i in 1:N)
# # {
# #   plot(time_steps, pat.norm.df[i,], lwd=2, col="red", type="l")
# #   Sys.sleep(1)
# # }
# 
# N_1 <- N_2 <- N_3 <- 10

#########################################################################
#########################################################################
# Load real dataset
setwd(data.dir)

healthy.obj <- load.patients.df(filename = "HOP_controlND_ThoraxPelvisHipKneeFoot_all.csv",
                                use.last.value = F)
therapy.obj <- load.patients.df(filename = "HOP_TPPI_ThoraxPelvisHipKneeFoot_all.csv",
                                use.last.value = F)
surgery.obj <- load.patients.df(filename = "HOP_RI_ThoraxPelvisHipKneeFoot_all.csv",
                                use.last.value = F)

# Save data frames
healthy.df <- healthy.obj$data
therapy.df <- therapy.obj$data
surgery.df <- surgery.obj$data

list.df <- list(healthy.df, therapy.df, surgery.df)

lapply(list.df, dim)

# Merge the 3 dataset
max.n.col <- max(ncol(healthy.df),
                 ncol(therapy.df),
                 ncol(surgery.df))

group.labels <- c("Control", "Therapy", "Surgery")

# Merge Data
patients.df <- merge.patients.df(list.df = list.df,
                                 max.n.col = max.n.col,
                                 group.labels = group.labels,
                                 replace.w.nan = T)

patients.df %>% dim()

patients.df <- patients.df[,-1]

# Remove columns from max(n_times)
n_times <- c(healthy.obj$n, therapy.obj$n, surgery.obj$n)
patients.df <- patients.df[,-((max(n_times)+1):ncol(patients.df))]

# Remove patients (59, 75) because they have NaN in the middle
patients.df <- patients.df[-c(59,75),]
n_times <- n_times[-c(59,75)]

################################################################################
# OPTIONAL: Reduce dimensions of the dataset to speed up computations
# patients.df <- patients.df[c(1:5, 33:37, 69:73),]
# n_times <- n_times[c(1:5, 33:37, 69:73)]
# N_1 <- N_2 <- N_3 <- 5
################################################################################
################################################################################
# Define all  the parameters which are needed

# Number of patients in the groups
# Legend: Group.1 = healthy, Group.2 = therapy, group.3 = surgery
N_1 <- nrow(healthy.df)
N_2 <- nrow(therapy.df)-1
N_3 <- nrow(surgery.df)-1

# Tot. patients
N <- sum(c(N_1, N_2, N_3))

# Indices that divide the groups
group_indices <- c(N_1, N_1+N_2, N_1+N_2+N_3)

# Number of groups
n_gr <- 3


################################################################################
### Rescale all data btw 0 and 1
# Step 1: Assign to each observation a vector of time instants between 0 and 1.
#         This vector of time steps redefine the tim ein which the measurement takes place.
# N.B. these values of abscissas must be stored in a list, as for each time series
# we have a different length. 
time_steps <- vector("list", N) #  one for each patient
for(i in 1:N)
{
  # Divide for the number of not NA objects and normalize
  time_steps[[i]] <- (1:n_times[i])/n_times[i]
}

# Check the lengths
lapply(time_steps, length) == n_times

# N.B. The interpolated functions are then evaluated over a grid of "series_length" points
# This means that in the end the time series will have a fixed length equal to "series_length".
series_length <- 500
pat.norm.df <- matrix(0, nrow = N, ncol = series_length)
for (i in 1:N) 
{
  # Interpolate with splines
  splint <- spline(x = time_steps[[i]], 
                   y = patients.df[i,1:n_times[i]], 
                   xout = seq(from = 0, to = 1, length.out = series_length))
  
  # Save the interpolated values into the df
  pat.norm.df[i,] <- splint$y
}

# Make a data.frame object
pat.norm.df <- data.frame(pat.norm.df, row.names = paste0("patient.", 1:N))
pat.norm.df <- setNames(pat.norm.df, paste0("T.", 1:series_length))
dim(pat.norm.df)

# Plot data to see if it worked 
x11()
par(mfrow=c(2,1))
for(i in 1:N)
{
  plot(1:n_times[i], patients.df[i,1:n_times[i]], type = "l", col="red")
  plot(1:series_length, pat.norm.df[i,], type = "l", col="blue")
  Sys.sleep(1)
}

# Set the new values of times
time_steps <- seq(from = 0, to = 1, length.out = series_length)

################################################################################
################################################################################
# Generate basis of splines 
m <- 4           # spline order 
degree <- m-1    # spline degree 

p1 <- 48 # number of knots for SMOOTHING
p2 <- 3 # number of knots for WARPING
K <- m + p1 - 2 # number of parameters = order + knots - 2 (constraints) 
Q <- m + p2 - 2
br1 <- seq(from = 0, to = 1, length.out = p1) # position of knots (equidistant)
br2 <- seq(from = 0, to = 1, length.out = p2)

# For the moment we consider the same basis of splines both for smoothing and warping
basis_smooth <- create.bspline.basis(rangeval = c(0, 1), 
                                     norder = m,  
                                     breaks = br1)
# dropind = c(1,K))
basis_warp <- create.bspline.basis(rangeval = c(0, 1), 
                                   norder = m,  
                                   breaks = br2)
# dropind = c(1,Q))

# Evaluate the basis on the grid of time points
basismat_beta <- eval.basis(time_steps, basis_smooth) # for smoothing
basismat_phi <- eval.basis(time_steps, basis_warp) # for warping

# N.B. basismat_phi is constant as it is B_mu(t), not B_mu(tau_ij)

# Basismat is a matrix in which:
# - each row is given by the evaluation of the B-spline basis in a given time point t
# - each column stores the evaluation of a given basis over the whole grid of time points

######################################################################
######################################################################
# Set the initial values of the hyperparameters
a_lambda <- 1
b_lambda <- 0.005
# b_lambda <- 5

a_eps <- 1
b_eps <- 0.005

a_a <- 0.1 
b_a <- 1

a_c <- 0.1
b_c <- 1

m_a0 <- 2
sigma2_a0 <- 1

m_c0 <- 20
sigma2_c0 <- 1

mean_phi <- get.mean.phi.coeff(n.knots = Q-m, spline.order = m) # length Q

m_0 <- rep(1, K)
SIGMA_0 <- diag(K) 

# Construct matrix Omega and P
omega <- diag(2,K)
for(i in 1:(K-1))
  omega[i,i+1] <- omega[i+1,i] <- -1
omega[K,K] <- 1

# Matrix  P
P <- diag(2,Q)
for(i in 1:(Q-1))
  P[i,i+1] <- P[i+1,i] <- -1
P[Q,Q] <- 1

P_inv <- solve(P)

# sigma2_phi (same role as lambda)
a_phi <- 0.1
b_phi <- 0.1
  
################################################################################
################################################################################
### Build matrices to store the values of the parameters produced by the MCMC; 
### Set the initial values of the chain

# Number of iterations for MCMC 
n_iter <- 40000

# One couple (a_ij, c_ij) of coefficients for each patient 
c_vect <- matrix(0, nrow = N, ncol = n_iter)
a_vect <- matrix(0, nrow = N, ncol = n_iter)
c_vect[,1] <- matrix(0, nrow = N)
a_vect[,1] <- matrix(1, nrow = N)

c0 <- rep(0, n_iter)
a0 <- rep(0, n_iter)
c0[1] <- 0.3
a0[1] <- 1

# One vector of Beta_j for each group
beta_1 <- matrix(0, nrow = K, ncol = n_iter)
beta_2 <- matrix(0, nrow = K, ncol = n_iter)
beta_3 <- matrix(0, nrow = K, ncol = n_iter)
beta_1[,1] <- 1 # Initialize first column to 1
beta_2[,1] <- 1 # Initialize first column to 1
beta_3[,1] <- 1 # Initialize first column to 1

# mu_beta is a k-dimensional vector common to all the groups
mu_beta <- matrix(0, nrow = K, ncol = n_iter)
mu_beta[,1] <- 1 # Initialize to one

# Initialize the B-spline basis parameters phi_ij
phi <- matrix(0, nrow = Q*N, ncol = n_iter)
phi[,1] <- rep(seq(0, 1, length.out =  Q), N)

################################################################################
# Following parameters are the same for each patient and group
sigma2_eps <- rep(0, n_iter)
sigma2_eps[1] <- 0.3

sigma2_a <- rep(0, n_iter)
sigma2_a[1] <- 0.3

sigma2_c <- rep(0, n_iter)
sigma2_c[1] <- 0.3

lambda <- rep(0, n_iter)
lambda[1] <- 0.3

sigma2_phi <- rep(0, n_iter)
sigma2_phi[1] <- 0.3

################################################################################
################################################################################
# Some variables useful for diagnostic 
accept_count <- 0
dot_prod_Y_m <- rep(0, N)
numeric.pat.norm.df <- as.matrix(pat.norm.df) 

phi_product <- rep(0, Q)

# Gibbs sampler 
for(i in 2:n_iter) 
{
  tic()
  ##############################################################################
  # Update of lambda
  a_lambda_star <- a_lambda + 3/2*K
  
  res_1 <- t(beta_1[,i-1]-mu_beta[,i-1]) %*% omega %*% (beta_1[,i-1]-mu_beta[,i-1])
  res_2 <- t(beta_2[,i-1]-mu_beta[,i-1]) %*% omega %*% (beta_2[,i-1]-mu_beta[,i-1])
  res_3 <- t(beta_3[,i-1]-mu_beta[,i-1]) %*% omega %*% (beta_3[,i-1]-mu_beta[,i-1])
  b_lambda_star <- b_lambda + 1/2*(res_1 + res_2 + res_3)
  
  lambda[i] <- 1/rgamma(1, shape=a_lambda_star, rate=b_lambda_star) 
  ##############################################################################
  # print("Lambda updated!!!") # OK DEBUG
  ##############################################################################
  # Update of sigma2_a and sigma2_c
  a_a_star <- a_a + N/2
  
  b_a_star <- b_a + 1/2 * sum((a_vect[,i-1] - a0[i-1])^2)
  
  sigma2_a[i] <- 1/rgamma(1, shape = a_a_star, rate = b_a_star)
  
  a_c_star <- a_c + N/2
  
  b_c_star <- b_c + 1/2 * sum((c_vect[,i-1] - c0[i-1])^2)
  
  aa <- runif(10)
  aa0 <- runif(1)
  
  sum((aa - aa0)^2)
  
  sigma2_c[i] <- 1/rgamma(1, shape = a_c_star, rate = b_c_star)
  ##############################################################################
  # print("Sigma2_a & Sigma2_c updated!!!") # DEBUG OK
  ##############################################################################
  # Update of a0 and c0
  sigma2_a0_star = 1/(N/sigma2_a[i] + (1/sigma2_a0))
  m_a0_star = sigma2_a0_star * (1/sigma2_a[i] * sum(a_vect[,i-1]) + m_a0/sigma2_a0)
  a0[i] = rnorm(1, mean = m_a0_star, sd = sqrt(sigma2_a0_star))
  
  sigma2_c0_star = 1/(N/sigma2_c[i] + (1/sigma2_c0))
  m_c0_star = sigma2_c0_star * (1/sigma2_c[i] * sum(c_vect[,i-1]) + m_c0/sigma2_c0)
  c0[i] = rnorm(1, mean = m_c0_star, sd = sqrt(sigma2_c0_star))
  ##############################################################################
  # print("a0 & c0 updated!!!") # DEBUG OK
  ##############################################################################
  ### Update tau_ij and m_ij values
  # Compute matrix of tau_ij at the current iteration
  # n.b. phi[i-1], since M-H step is at the end
  tau_ij <- get_tau(mat_basis = basismat_phi, 
                    coeff_vect = phi[,i-1], 
                    N_patients = N, 
                    N_coeff = Q) #OK
  
  # Compute the matrix of m_ij at the current iteration
  m_ij <- eval_model(basis_set = basis_smooth, tau_mat = tau_ij,  
                     curr_a = a_vect[,i-1], curr_c = c_vect[,i-1], 
                     curr_beta_list = list(beta_1[,i-1], beta_2[,i-1], beta_3[,i-1]), 
                     group_idxs = group_indices)
  
  # Update of sigma_eps
  a_eps_star <- a_eps + N*series_length / 2
  
  # Build  matrix with pointwise differences
  diff_Y_m <- numeric.pat.norm.df - m_ij
  # Do the dot product
  dot_prod_Y_m <- rowSums(diff_Y_m * diff_Y_m)
  
  b_eps_star <- b_eps + 1/2 * sum(dot_prod_Y_m)
  
  sigma2_eps[i] <- 1/rgamma(1, shape = a_eps_star, rate = b_eps_star)
  ##############################################################################
  # print("sigma_eps updated!!!") # DEBUG OK
  ##############################################################################
  # Update of mu_beta
  mu_n <- (n_gr*omega/lambda[i] + solve(SIGMA_0)) %*%
    ((omega/lambda[i]) %*% (beta_1[,i-1]+beta_2[,i-1]+beta_3[,i-1]) + solve(SIGMA_0) %*% m_0)
  
  SIGMA_n <- solve(n_gr*omega/lambda[i] + solve(SIGMA_0))
  
  mu_beta[,i] <- mvrnorm(n = 1, mu = mu_n, Sigma = SIGMA_n)
  ##############################################################################
  # print("mu_beta updated!!!") # DEBUG OK, PROBLEMS FOUND
  ##############################################################################
  # Update of Beta_j
  M_1 <- get_M_j(basis_set = basis_smooth, tau_mat = tau_ij, 
                 curr_a = a_vect[,i-1], j = 1, 
                 group_idxs = group_indices, 
                 n_coeff = K, series_len = series_length)
  M_2 <- get_M_j(basis_set = basis_smooth, tau_mat = tau_ij, 
                 curr_a = a_vect[,i-1], j = 2, 
                 group_idxs = group_indices,
                 n_coeff = K, series_len = series_length)
  M_3 <- get_M_j(basis_set = basis_smooth, tau_mat = tau_ij, 
                 curr_a = a_vect[,i-1], j = 3, 
                 group_idxs = group_indices, 
                 n_coeff = K, series_len = series_length)
  # Sembra OK
  
  # V_beta_1 <- solve((t(M_1) %*% M_1) / sigma2_eps[i] + omega / lambda[i])
  # V_beta_2 <- solve((t(M_2) %*% M_2) / sigma2_eps[i] + omega / lambda[i])
  # V_beta_3 <- solve((t(M_3) %*% M_3) / sigma2_eps[i] + omega / lambda[i])
  # 
  # Use pseudo-inverse to guarantee numeric feasibility
  V_beta_1 <- ginv((t(M_1) %*% M_1) / sigma2_eps[i] + omega / lambda[i])
  V_beta_2 <- ginv((t(M_2) %*% M_2) / sigma2_eps[i] + omega / lambda[i])
  V_beta_3 <- ginv((t(M_3) %*% M_3) / sigma2_eps[i] + omega / lambda[i])
  # Ok di conseguenza
  
  # Compute (Y_j - C_j)
  Y_minus_C <- numeric.pat.norm.df - c_vect[,i-1] # OK
  Y_minus_C_1 <- as.vector(t(as.matrix(Y_minus_C[1:group_indices[1],])))
  Y_minus_C_2 <- as.vector(t(as.matrix(Y_minus_C[(group_indices[1]+1):group_indices[2],])))
  Y_minus_C_3 <- as.vector(t(as.matrix(Y_minus_C[(group_indices[2]+1):group_indices[3],])))
  
  m_beta_1 <- (V_beta_1 %*% 
     as.vector((1/sigma2_eps[i]) * t(M_1) %*% Y_minus_C_1 + (omega/lambda[i]) %*% mu_beta[,i]))
  m_beta_2 <- (V_beta_2 %*% 
     as.vector((1/sigma2_eps[i]) * t(M_2) %*% Y_minus_C_2 + (omega/lambda[i]) %*% mu_beta[,i]))
  m_beta_3 <- (V_beta_3 %*% 
     as.vector((1/sigma2_eps[i]) * t(M_3) %*% Y_minus_C_3 + (omega/lambda[i]) %*% mu_beta[,i]))
  # NOT OK
  
  # Draw new vector of Beta_j's
  beta_1[,i] <- mvrnorm(n = 1, mu = as.vector(m_beta_1), Sigma = V_beta_1)
  beta_2[,i] <- mvrnorm(n = 1, mu = as.vector(m_beta_2), Sigma = V_beta_2)
  beta_3[,i] <- mvrnorm(n = 1, mu = as.vector(m_beta_3), Sigma = V_beta_3)
  ##############################################################################
  # print("Beta_j updated!!!")
  ##############################################################################
  # Update of a_ij & c_ij
  SIGMA_c_a_inv <- solve(diag(c(sigma2_c[i], sigma2_a[i])))
  
  # N.B. a_ij and c_ij are iid, therefore the extraction of the values can 
  # also happen in cascade. 
  # for loop to update a_ij, c_ij 
  for(p in 1:N)
  {
    # Choose Beta vector depending on the group
    beta_vect <- NULL 
    if(p <= group_indices[1])
      beta_vect <- beta_1[,i]
    else if (p > group_indices[1] && p <= group_indices[2])
      beta_vect <- beta_2[,i]
    else
      beta_vect <- beta_3[,i]
    
    W_ij <- get_W_ij(basis_set = basis_smooth,
                     tau_ij = tau_ij[p,], 
                     beta_j = beta_vect, 
                     length_series_ij = series_length)
    
    # 2x2 matrix
    # SIGMA_l_ij <- solve(SIGMA_c_a_inv + (1/sigma2_eps[i]) * t(W_ij) %*% W_ij)
    SIGMA_l_ij <- ginv(SIGMA_c_a_inv + (1/sigma2_eps[i]) * t(W_ij) %*% W_ij)
    
    Y_ij <- numeric.pat.norm.df[p,]
    
    # 2 elements vector
    m_l_ij <- as.numeric(SIGMA_l_ij %*% 
       (SIGMA_c_a_inv %*% c(c0[i], a0[i]) + (1/sigma2_eps[i]) * t(W_ij) %*% Y_ij))
    
    samples <- mvrnorm(n = 1, mu = m_l_ij, Sigma = SIGMA_l_ij)
    # # Generate samples from truncated multivariate normal
    # samples <- rtmvnorm(n = 1, mean = m_l_ij, sigma = SIGMA_l_ij,
    # lower = c(-Inf, 0), algorithm = "rejection")
    
    # Update the values
    c_vect[p,i] <- samples[1] 
    # a_vect[p,i] <- samples[2] 
    # If sample of a_ij is negative, keep the value of previous iteration
    a_vect[p,i] <- ifelse(samples[2]>=0, samples[2], a_vect[p,i-1])
  } 
  ##############################################################################
  # print("a_ij & c_ij updated")
  ##############################################################################
  # Update of sigma2_phi
  a_phi_star <- a_phi + Q*N/2

  # Store the current values of phi in a matrix to make easier the following operations
  curr_phi_mat <- matrix(phi[,i-1], nrow = N, ncol = Q, byrow = T)
  
  # Define the matrix with rows (phi_i - mean_phi)
  phi_minus_mean <- t(t(curr_phi_mat) - mean_phi)
  
  # Matrix product between (phi_i - mean_phi) and matrix P
  phi_product <- phi_minus_mean %*% P
  
  # Rowsums trick for the other product
  phi_product <- rowSums(phi_product * phi_minus_mean)
    
  b_phi_star <- b_phi + 0.5 * sum(phi_product)
  
  # Draw new value from full conditional distribution
  sigma2_phi[i] <- 1/rgamma(1, shape = a_phi_star, rate = b_phi_star)
  
  ##############################################################################
  ##############################################################################
  # Metropolis-Hastings step for PHI 
  # Update the S_l parameters  and Phi for each patient
  for(j in 1:N)
  {
    # Choose the right vector of beta
    # Take the vector of depending on the group
    curr_beta <- NULL
    if(j <= group_indices[1]) 
      curr_beta <- beta_1[,i]
    else if (j > group_indices[1] && j <= group_indices[2])
      curr_beta <- beta_2[,i]
    else
      curr_beta <- beta_3[,i]
    
    # Run MH step
    ### SIMPLE MH
    MH_res <- MH_step_phi(curr_phi = phi[((j-1)*Q+1):(j*Q), i-1], 
                          curr_dot_prod_vec = dot_prod_Y_m, 
                          mat_basis_phi = basismat_phi, 
                          basis_set_beta = basis_smooth,
                          curr_a = a_vect[j,i], 
                          curr_c = c_vect[j,i], 
                          curr_beta = curr_beta,  
                          curr_sigma2_eps = sigma2_eps[i],
                          curr_sigma2_phi = sigma2_phi[i],
                          phi_mean_vec = mean_phi, 
                          P_matrix_inv = P_inv,
                          Y_patient = numeric.pat.norm.df[j,], 
                          patient_index = j)
    
    phi[((j-1)*Q+1):(j*Q),i] <- MH_res$phi_update
    accept_count <- accept_count + MH_res$accept
    
  }
  
  ##############################################################################
  # print("MH updated!!!")
  print(paste("Iteration nr.", i-1, "has been completed!!!"))
  toc()
}

################################################################################
# Thinning and saving the data
thin_factor <- 10
thinning_mask <- seq(1, n_iter, by = thin_factor)

lambda <- lambda[thinning_mask]
sigma2_a <- sigma2_a[thinning_mask]
sigma2_c <- sigma2_c[thinning_mask]
sigma2_eps <- sigma2_eps[thinning_mask]
a0 <- a0[thinning_mask]
c0 <- c0[thinning_mask]
a_vect <- a_vect[,thinning_mask]
c_vect <- c_vect[,thinning_mask]
beta_1 <- beta_1[,thinning_mask]
beta_2 <- beta_2[,thinning_mask]
beta_3 <- beta_3[,thinning_mask]
mu_beta <- mu_beta[,thinning_mask]
# S <- S[,thinning_mask]
phi <- phi[,thinning_mask]

# Save workspace
r.data.dir <- paste(main.dir, "/RData/", sep="")
setwd(r.data.dir)
MCMC_results <- list(lambda, sigma2_a, sigma2_c, sigma2_eps, a0, c0, a_vect, c_vect,
                     beta_1, beta_2, beta_3, mu_beta, phi)
save(MCMC_results,
     file = "MCMC_samples_45k_completedata_500_SlVersion.RData")



