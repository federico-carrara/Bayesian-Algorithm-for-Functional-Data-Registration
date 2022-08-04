# Load Rdata
r.data.dir <- paste(main.dir, "/RData/", sep="")
load(paste(r.data.dir, "/MCMC_samples_45k_completedata_500_SlVersion.RData", sep=""))
n_iter <- 45000
series_length <- 500
Q <- 10

MCMC_results %>% str()

lambda <- MCMC_results[[1]]
sigma2_a <- MCMC_results[[2]]
sigma2_c <- MCMC_results[[3]]
sigma2_eps <- MCMC_results[[4]]
a0 <- MCMC_results[[5]]
c0 <- MCMC_results[[6]]
a_vect <- MCMC_results[[7]]
c_vect <- MCMC_results[[8]]
beta_1 <- MCMC_results[[9]]
beta_2 <- MCMC_results[[10]]
beta_3 <- MCMC_results[[11]]
mu_beta <- MCMC_results[[12]]
phi <- MCMC_results[[13]]


################################################################################
# Visualization of MCMC results
plot(log(lambda), type="l", main = "Traceplot of lambda")

plot(sigma2_a, type="l", main = "Traceplot of (sigma^2)_a")

plot(sigma2_c, type="l", main = "Traceplot of (sigma^2)_c")

plot(a0, type="l", main = "Traceplot of a_0")

plot(c0, type="l", main = "Traceplot of c_0")

x11()
for(i in 1:K)
{
  plot(beta_1[i,], type = "l", ylim=c(-50,50))
  Sys.sleep(1)
}

x11()
for(i in 1:K)
{
  plot(beta_2[i,], type = "l")
  Sys.sleep(5)
}

x11()
for(i in 1:K)
{
  plot(beta_3[i,], type = "l")
  Sys.sleep(1)
}

x11()
for(i in 1:K)
{
  plot(mu_beta[i,], type = "l")
  Sys.sleep(1)
}

x11()
for(i in 1:N)
{
  plot(a_vect[i,], type = "l")
  Sys.sleep(1)
}

x11()
for(i in 1:N)
{
  plot(c_vect[i,], type = "l")
  Sys.sleep(1)
}

cols <- rep(brewer.pal(Q, "Set1"), N)
idx <- 1

x11()
for(i in 1:(N*Q))
{
  plot(S[i,], type = "l", col = cols[idx])
  Sys.sleep(1)
  if(i %% Q == 0) 
    idx <- idx + 1
}

x11()
for(i in 1:(N*Q))
{
  plot(phi[i,], type = "l", col = cols[idx], ylim = c(0,1))
  Sys.sleep(1)
  if(i %% Q == 0) 
    idx <- idx + 1
}

x11()
par(mfrow=c(2,5))
for(i in 1:Q)
{
  plot(S[i,], type = "l", col = "red", main = paste("Posterior samples of S_", i))
}
for(i in 1:Q)
{
  plot(phi[i,], type = "l", col = "blue", main = paste("Posterior samples of PHI_", i))
}

# Diagnostic of MH
accept_count/(n_iter*N)

################################################################################
# Compute the posterior estimates of the warping and the models
# 1. POSTERIOR OF THE WARPING FUNCTION TAU_ij
MCMC_samples <- n_iter/10
burn_in <- 3000

time_identity <- matrix(seq(0,1,len=series_length), nrow=N, ncol=series_length, byrow = T)

# Compute the posterior of the alligned model m_ij for each single individual
# N.B. To compute the posterior samples we use the posterior samples ((tau_ij)^-1)(m) 
# as the abscissas in which each sample (m_ij)(m) is evaluated
model_samples_allign <- vector("list", MCMC_samples-burn_in)
for(i in 1:(MCMC_samples-burn_in))
{
  model_samples_allign[[i]] <- eval_model(basis_set = basis_smooth, 
                                          #tau_mat = tau_inv_samples[[i]],
                                          tau_mat = time_identity, 
                                          curr_a = a_vect[,(i+burn_in)], 
                                          curr_c = c_vect[,(i+burn_in)], 
                                          curr_beta_list = list(beta_1[,(i+burn_in)], 
                                                                beta_2[,(i+burn_in)], 
                                                                beta_3[,(i+burn_in)]),
                                          group_idxs = group_indices)
  print(paste("Iteration nr.", i))
}

# Compute the posterior means of the m_ij for each patient
model_est_allign <- matrix(0, nrow = N, ncol = series_length)
for(i in 1:length(model_samples_allign))
{
  model_est_allign <- model_est_allign + model_samples_allign[[i]]
}
model_est_allign <- model_est_allign / length(model_samples_allign)


################################################################################
# Compute posterior samples of tau_ij
tau_samples <- vector("list", MCMC_samples-burn_in)
for(i in 1:(MCMC_samples-burn_in))
{
  tau_samples[[i]] <- get_tau(mat_basis = basismat_phi,
                              coeff_vect = phi[,(i+burn_in)],
                              N_patients = N,
                              N_coeff = Q)
  print(paste("Iteration nr.", i))
}

# Compute the means for each patient
tau_est <- matrix(0, nrow = nrow(tau_samples[[1]]), ncol = ncol(tau_samples[[1]]))
for(i in 1:length(tau_samples))
{
  tau_est <- tau_est + tau_samples[[i]]
}
tau_est <- tau_est / length(tau_samples)

# Invert the posterior samples of tau to get the posterior samples of (tau)^-1
tau_inv_samples <- vector("list", MCMC_samples-burn_in)
for(i in 1:(MCMC_samples-burn_in))
{
  tau_inv_samples[[i]] <- get_tau_inverse_interp(tau_mat = tau_samples[[i]])
  print(paste("Iteration nr.", i))
}

# Compute the posterior mean of (tau)^-1
tau_inv_est <- matrix(0, nrow = N, ncol = series_length)
for(i in 1:length(tau_inv_samples))
{
  tau_inv_est <- tau_inv_est + tau_inv_samples[[i]]
}
tau_inv_est <- tau_inv_est / length(tau_inv_samples)

# 2. POSTERIOR OF THE MODEL FUNCTIONS
model_samples <- vector("list", MCMC_samples-burn_in)
for(i in 1:(MCMC_samples-burn_in))
{
  model_samples[[i]] <- eval_model(basis_set = basis_smooth, 
                                   tau_mat = tau_samples[[i]], 
                                   curr_a = a_vect[,(i+burn_in)], 
                                   curr_c = c_vect[,(i+burn_in)], 
                                   curr_beta_list = list(beta_1[,(i+burn_in)], 
                                                         beta_2[,(i+burn_in)], 
                                                         beta_3[,(i+burn_in)]),
                                   group_idxs = group_indices)
  print(paste("Iteration nr.", i))
}

# Compute the means for each patient
model_est <- matrix(0, nrow = nrow(model_samples[[1]]), ncol = ncol(model_samples[[1]]))
for(i in 1:length(model_samples))
{
  model_est <- model_est + model_samples[[i]]
}
model_est <- model_est / length(model_samples)

################################################################################
### 1ST SLIDE
# Original curve
# x11()
# plot(x = seq(0,1,len=series_length), y = pat.norm.df[14,], type="l", lwd=2, 
#      main="Observation", xlab = "Time [t]", ylab = "Knee Angle [y(t)]", col = "blue",
#      ylim = c(5,75), cex.lab=1.5, cex.main = 2)
# 
# # Original plus model
# x11()
# plot(x = seq(0,1,len=series_length), y = pat.norm.df[14,], type="l", lwd=2, 
#      main="Observation & Shape Function", xlab = "Time [t]", ylab = "Knee Angle [y(t)]", 
#      col = "lightblue", ylim = c(5,75), cex.lab=1.5, cex.main = 2)
# lines(x = seq(0,1,len=series_length), y = model_est[14,], type="l", lwd=2, 
#       xlab = "Time [t]", ylab = "Knee Angle [y(t)]", col = "red", ylim = c(5,75))
# 
# # for(i in 1:N)
# # {
# #   plot(x = seq(0,1,len=series_length), y = pat.norm.df[i,], type="l", lwd=2,
# #        main=paste("Observation & Shape Function nr.", i), xlab = "Time [t]", ylab = "Knee Angle [y(t)]",
# #        col = "lightblue")
# #   lines(x = seq(0,1,len=series_length), y = model_est[i,], type="l", lwd=2,
# #         xlab = "Time [t]", ylab = "Knee Angle [y(t)]", col = "red")
# #   Sys.sleep(1)
# # }
# 
# # Matplot of the warping
# x11()
# matplot(t(tau_est[seq(1,N,2),]), type = "l", main = "Warping functions", lwd=2, lty = 1,
#         ylab = "Time [t]", cex.lab=1.5, cex.main = 2)

################################################################################
### Plot the results

# 1. Plot of all the curves to check the overall alignment
x11()
matplot(t(pat.norm.df[1:group_indices[1],]), lwd=2, lty = 1,
        type = "l", main = "Original curves Control group", 
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", cex.main = 2, cex.lab=1.5)
x11()
matplot(t(model_est_allign[1:group_indices[1],]), lwd=2, lty = 1,
        type = "l", main = "Registered curves Control group",
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", cex.main = 2, cex.lab=1.5)
x11()
matplot(t(pat.norm.df[(group_indices[1]+1):group_indices[2],]), lwd=2, lty = 1,
        type = "l", main = "Original curves Therapy group",
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", cex.main = 2, cex.lab=1.5)
x11()
matplot(t(model_est_allign[(group_indices[1]+1):group_indices[2],]), lwd=2, lty = 1,
        type = "l", main = "Registered curves Therapy group",
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", cex.main = 2, cex.lab=1.5)
x11()
matplot(t(pat.norm.df[(group_indices[2]+1):group_indices[3],]), lwd=2, lty = 1,
        type = "l", main = "Original curves Surgery group",
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", cex.main = 2, cex.lab=1.5)
x11()
matplot(t(model_est_allign[(group_indices[2]+1):group_indices[3],]), lwd=2, lty = 1,
        type = "l", main = "Registered curves Surgery group", cex.main = 2, cex.lab=1.5,
        xlab = "Time instants", ylab = "Knee Angle [y(t)]")


# Curves with the means in red
# Group_means
model_means1 <- colMeans(model_est_allign[1:group_indices[1],])
model_means2 <- colMeans(model_est_allign[(group_indices[1]+1):group_indices[2],])
model_means3 <- colMeans(model_est_allign[(group_indices[2]+1):group_indices[3],])

x11()
matplot(t(model_est_allign[1:new_group_indices[1],]), lwd=2, lty = 1,
        type = "l", main = "Mean curve Control group", cex.main = 2, cex.lab=1.5,
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", col="grey")
lines(model_means1, lwd=3, lty = 1, type = "l", col ="red")
x11()
matplot(t(model_est_allign[(new_group_indices[1]+1):new_group_indices[2],]), lwd=2, lty = 1,
        type = "l", main = "Mean curve Therapy group", cex.main = 2, cex.lab=1.5,
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", col="grey")
lines(model_means2, lwd=3, lty = 1, type = "l", col ="red")
x11()
matplot(t(model_est_allign[(new_group_indices[2]+1):new_group_indices[3],]), lwd=2, lty = 1,
        type = "l", main = "Mean curve Surgery group", cex.main = 2, cex.lab=1.5,
        xlab = "Time instants", ylab = "Knee Angle [y(t)]", col="grey")
lines(model_means3, lwd=3, lty = 1, type = "l", col ="red")


# Traceplots of S and PHI
phi_ij_scritta <- expression(phi[ij])
S_ij_scritta <- expression(S[ij])
x11()
plot(phi[2,], lwd=2, type = "l", col = "blue", ylim = c(0,1), 
     main = paste("Traceplot of ", phi_ij_scritta, sep=""),
     xlab = "Iterations", ylab = expression(phi[ij]), cex.main = 2, cex.lab=1.5)
for(i in 3:(Q-1))
{
  lines(phi[i,], lwd=2, type = "l", col = "blue", ylim = c(0,1))
}
x11()
plot(S[2,], lwd=2, type = "l", col = "red", 
     main = paste("Traceplot of ", S_ij_scritta, sep=""),
     xlab = "Iterations", ylab = expression(S[ij]), cex.main = 2, cex.lab=1.5)
for(i in 3:(Q-1))
{
  lines(S[i,], lwd=2, type = "l", col = "red", ylim = c(0,1))
}


################################################################################
# Common shape function according to Telesca paper
model_samples_cs <- vector("list", MCMC_samples-burn_in)
for(i in 1:(MCMC_samples-burn_in))
{
  model_samples_cs[[i]] <- matrix(0, nrow = N, ncol = series_length)
  for(j in 1:N)
  {
    curr_beta <- NULL
    if(j <= group_indices[1]) 
      curr_beta <- beta_1[,i]
    else if (j > group_indices[1] && j <= group_indices[2])
      curr_beta <- beta_2[,i]
    else
      curr_beta <- beta_3[,i]
    
    curr_basismat_beta <- eval.basis(tau_samples[[i]][j,], basis_smooth)
    curr_model <- c0[burn_in+i] + a0[burn_in+i]*curr_basismat_beta%*%curr_beta
    model_samples_cs[[i]][j,] <- curr_model
  }
  print(paste("Iter nr.", i))
}
str(model_samples_cs)

# Compute the means for each patient
model_est_cs <- matrix(0, nrow = nrow(model_samples_cs[[1]]), ncol = ncol(model_samples_cs[[1]]))
for(i in 1:length(model_samples_cs))
{
  model_est_cs <- model_est_cs + model_samples_cs[[i]]
}
model_est_cs <- model_est_cs / length(model_samples_cs)

x11()
par(mfrow=c(2,1))
matplot(t(model_est_cs), type = "l", main = "Posterior estimate of the common shapes")
matplot(t(pat.norm.df), type = "l", main = "Real data")

# Common shape vs. Original curves
x11()
for(i in 1:N)
{
  plot(1:series_length, pat.norm.df[i,], col="red", type="l", 
       ylim = c(0,max(c(model_est_cs[i,], as.numeric(pat.norm.df[i,])))))
  # plot(1:series_length, model_est_cs[i,], col="blue", type="l")
  lines(model_est_cs[i,], col="blue")
  legend("topright", legend = c("Real data", "Model approx"), 
         col = c("red", "blue"), lty = 1, cex = 0.8)
  Sys.sleep(1)
}
################################################################################


