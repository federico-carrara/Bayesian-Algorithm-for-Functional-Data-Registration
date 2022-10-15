################################################################################
# SOME FUNCTIONS USED FOR THE COMPUTATIONS
# To compute Phi coefficients starting from S_l values
get_Phi <- function(S_vect)
{
  T_sum <- sum(S_vect)
  phi_vect <- rep(0, length(S_vect)) 
  for(i in 1:length(S_vect))
  {
    phi_vect[i] <- sum(S_vect[1:i])/T_sum
  }
  return(phi_vect)
}

# Alternative definition of phi using (S_l)^2
get_Phi_2 <- function(S_vect)
{
  T_sum <- sum(S_vect^2)
  phi_vect <- rep(0, length(S_vect)) 
  for(i in 1:length(S_vect))
  {
    phi_vect[i] <- sum(S_vect[1:i]^2)/T_sum
  }
  return(phi_vect)
}

# Function to compute tau_ij(t) for all the patients
get_tau <- function(mat_basis, coeff_vect, N_patients, N_coeff) # coeff_vect è il vettore delle phi di tutti i pazienti 
{
  tau_mat <- matrix(0, nrow = N_patients, ncol = nrow(mat_basis))
  for(i in 1:N_patients)
  {
    tau_mat[i,] <- mat_basis %*% coeff_vect[((i-1)*N_coeff+1):(i*N_coeff)]
  }
  
  # Transpose tau so that each row represent the time transformation for a single patient
  return(tau_mat)
}

# Get the evaluation of the model for all the patients in the new time steps
# Store it in each row of a matrix
### - curr_a, curr_c are the vectors containing the value of such
###   parameters for each patients at the current iteration of MCMC.
### - curr_beta is a list containing beta of the 3 groups. 
eval_model <- function(basis_set, tau_mat, curr_a, curr_c,
                       curr_beta_list, group_idxs)
{
  # Matrix to store all the models -> same size as the times tau
  m_mat <- matrix(0, nrow = nrow(tau_mat), ncol = ncol(tau_mat))
  
  # Number of patients
  n_pats <- nrow(tau_mat)
  
  # Get the model evaluation for each patient
  for(i in 1:n_pats)
  {
    # Take the vector of depending on the group
    curr_beta <- NULL
    if(i <= group_idxs[1]) # i-esimo paziente appartenente al primo gruppo
      curr_beta <- curr_beta_list[[1]]
    else if (i > group_idxs[1] && i <= group_idxs[2]) # i-esimo paziente appartenente al secondo gruppo
      curr_beta <- curr_beta_list[[2]]
    else
      curr_beta <- curr_beta_list[[3]]
    
    # Evaluate the set of basis in the transformed time steps tau
    curr_basismat <- eval.basis(tau_mat[i,], basis_set)
    m_mat[i,] <- curr_c[i] + curr_a[i] * curr_basismat %*% curr_beta
    # N.B. rows of basismat sums to 1!!
  }
  
  return(m_mat) # ritorna la matrice del modello
}

# Get the list of matrices for update of Beta
# - series_len è la lunghezza della time series di ciascun paz. (Nij)
# - j is the index of the group
# - group_idxs is a vector of the indices delimiting the groups
get_M_j <- function(basis_set, tau_mat, curr_a, j, group_idxs, 
                    n_coeff, series_len)
{
  # Indices for the current group j
  start <- ifelse(j>1, group_idxs[j-1]+1, 1)
  end <- group_idxs[j] # indica dove finisce ogni gruppo
  
  # Matrix to store the M_j
  M_j <- matrix(0, nrow = series_len*(end-start+1), ncol = n_coeff)
  
  curr_index <- 1
  
  # Get the model evaluation for each patient
  for(i in start:end)
  {
    # Evaluate the set of basis in the transformed time steps tau
    ### N.B. Siccome la time series di ogni paziente ha lunghezza N_ij diversa, 
    ### di conseguenza la basismat di ogni paziente avra dimensioni (N_ij x k), dove
    ### k rappresenta il numero di basi
    new_index <- curr_index + series_len - 1
    
    curr_basismat <- eval.basis(tau_mat[i,], basis_set) # dim: T_i x k
    
    M_j[curr_index:new_index,] <- curr_a[i] * curr_basismat # dim: sum(T_i) x k
    curr_index <- new_index + 1
  }
  
  return(M_j)
}

# Compute matrix W_ij for each patient: (n_times[i])x(2 columns)
get_W_ij <- function(basis_set, tau_ij, beta_j, length_series_ij)
{
  # Initialize the first column of the matrix with 1s
  W_ij <- matrix(1, nrow = length_series_ij, ncol = 1)
  
  # Evaluate the basis in tau_ij
  basismat <- eval.basis(tau_ij, basis_set)
  
  # Bind a second column with the dot prod btw basismat and Beta
  W_ij <- cbind(W_ij, basismat %*% beta_j)
  
  return(W_ij)
}


################################################################################

##############################################################################
# Step of Metropolis-Hastings 
### N.B. Random-Walk MH with N(0, 0.1) as innovation
### N.B. Update for all the parameters of the single patient
### N.B. Also ADAPTIVE version is implemented
MH_step_S <- function(curr_S, curr_dot_prod_vec,
                      innovation_covmat, 
                      gamma_l_hyp, mat_basis_phi, 
                      basis_set_beta, patient_index,
                      curr_a, curr_c, curr_beta, 
                      curr_sigma2_eps, Y_patient,
                      use_adaptive, curr_iter, 
                      curr_covmat_P = NULL, t_0) 
{
  # Record whether we accept or not
  accept <- F
  
  ##############################################################################
  # Set the standard deviation adaptively if needed
  # if(use_adaptive && curr_iter > t_0)
  if(use_adaptive && curr_iter >= ncol(curr_covmat_P) && curr_iter >= t_0)
  {
    # Compute the standard deviation of P = log(S)
    # N.B. We approximate P~N, hence S~logN, so we can use the formula for the 
    # realtionship between the sd of logN and N distributions
    innovation_covmat <- curr_covmat_P 
  }
  ##############################################################################
  
  # N.B. S_1=0 ALWAYS!!!! -> add innovation only to S_2, ..., S_Q 
  # Moreover as S_1 is fixed, it is not a r.v.
  # Therefore it can be exluded from the computations!!
  curr_P <- log(curr_S[-1]) # P = log(S) -> S = exp(P) SEMPRE POSITIVO!!!
  # n.b. curr_P has length 4
  
  # Draw new values for the vector of P
  innov_P <- mvrnorm(n = 1, mu = rep(0,length(curr_P)), Sigma = innovation_covmat)
  new_P <- curr_P + innov_P
  
  # Transform the P_l into S_l (to guarantee positivity, S_1 excluded)
  new_S <- exp(new_P)
  
  # Compute new coefficient phi for the current patient
  new_phi <- get_Phi_2(S_vect = c(0,new_S))
  
  # Compute the new vector of tau for the current patient
  new_tau <- as.numeric(mat_basis_phi %*% new_phi)
  
  # Given the new value of S evaluate the model for the current patient
  curr_mat_basis_beta <- eval.basis(new_tau, basis_set_beta)
  new_model <- curr_c + curr_a * curr_mat_basis_beta %*% curr_beta # aggiorno il modello 
  
  # Evaluate the new dot product and put it in the vector
  new_dot_prod_vec <- curr_dot_prod_vec
  new_dot_prod_vec[patient_index] <- t(Y_patient - new_model) %*% (Y_patient - new_model)
  
  # Compute acceptance probability (equal for every S_l)
  # 1. Evaluate the full conditional with the current and the new values
  new_log_density <- -(1/(2*curr_sigma2_eps))*sum(new_dot_prod_vec) + 
    sum(new_P * gamma_l_hyp - new_S)
  curr_log_density <- -(1/(2*curr_sigma2_eps))*sum(curr_dot_prod_vec) + 
    sum(curr_P * gamma_l_hyp - curr_S[-1])
  
  # 2. Compute density ratio
  log_density_ratio <- new_log_density - curr_log_density
  
  # 3. Compute accept prob
  accept_prob <- min(log_density_ratio, 0)
  
  # Check acceptance (for each)
  u_check <- log(runif(1, 0, 1))
  
  # This is the case for which we accept
  S_ij <- curr_S
  P_ij <- curr_P
  dot_prod_ij <- NULL
  model_ij <- NULL
  tau_ij <- NULL
  phi_ij <- NULL
  if (u_check < accept_prob) 
  {
    S_ij <- c(0,new_S)
    P_ij <- new_P # has length 4
    dot_prod_ij <- new_dot_prod_vec[patient_index]
    model_ij <- new_model
    tau_ij <- new_tau
    phi_ij <- new_phi
    accept <- T
  }
  
  return(list(S_update = S_ij, 
              P_update = P_ij,
              accept = accept, 
              dot_prod_update = dot_prod_ij,
              model_update = model_ij,
              tau_update = tau_ij, 
              phi_update = phi_ij))
}

MH_step_S_univar <- function (curr_S, curr_dot_prod_vec,
                              innovation_sd = 0.01, gamma_l_hyp,
                              mat_basis_phi, basis_set_beta, 
                              curr_a, curr_c, curr_beta, 
                              curr_sigma2_eps, Y_patient,
                              patient_index) 
{
  # Record whether we accept or not
  accept <- F
  
  # N.B. S_1=0 ALWAYS!!!! -> add innovation only to S_2, ..., S_Q 
  # Moreover as S_1 is fixed, it is not a r.v.
  # Therefore it can be exluded from the computations!!
  curr_P <- log(curr_S[-1]) # P = log(S) -> S = exp(P) SEMPRE POSITIVO!!!
  # n.b. curr_P has length 4
  
  # Draw new values for the vector of P
  innov_P <- rnorm(n = length(curr_P), mean = 0, sd = innovation_sd)
  new_P <- curr_P + innov_P
  
  # Transform the P_l into S_l (to guarantee positivity, S_1 excluded)
  new_S <- exp(new_P)
  
  # Compute new coefficient phi for the current patient
  new_phi <- get_Phi_2(S_vect = c(0,new_S))
  
  # Compute the new vector of tau for the current patient
  new_tau <- as.numeric(mat_basis_phi %*% new_phi)
  # if(sum(is.na(new_tau)))
  # {
  #   print(paste("P: ", new_P))
  #   print(paste("S: ", new_S))
  #   print(paste("PHI: ", new_phi))
  # }
  
  # Given the new value of S evaluate the model for the current patient
  curr_mat_basis_beta <- eval.basis(new_tau, basis_set_beta)
  new_model <- curr_c + curr_a * curr_mat_basis_beta %*% curr_beta # aggiorno il modello 
  
  # Evaluate the new dot product and put it in the vector
  new_dot_prod_vec <- curr_dot_prod_vec
  new_dot_prod_vec[patient_index] <- t(Y_patient - new_model) %*% (Y_patient - new_model)
  
  # Compute acceptance probability (equal for every S_l)
  # 1. Evaluate the full conditional with the current and the new values
  new_log_density <- -(1/(2*curr_sigma2_eps))*sum(new_dot_prod_vec) + 
    sum(new_P * gamma_l_hyp - new_S)
  curr_log_density <- -(1/(2*curr_sigma2_eps))*sum(curr_dot_prod_vec) + 
    sum(curr_P * gamma_l_hyp - curr_S[-1])
  
  # 2. Compute density ratio
  log_density_ratio <- new_log_density - curr_log_density
  
  # 3. Compute accept prob
  accept_prob <- min(log_density_ratio, 0)
  
  # Check acceptance (for each)
  u_check <- log(runif(1, 0, 1))
  
  # This is the case for which we accept
  S_ij <- curr_S
  P_ij <- curr_P
  dot_prod_ij <- NULL
  model_ij <- NULL
  tau_ij <- NULL
  phi_ij <- NULL
  if (u_check < accept_prob) 
  {
    S_ij <- c(0,new_S)
    P_ij <- new_P # has length 4
    dot_prod_ij <- new_dot_prod_vec[patient_index]
    model_ij <- new_model
    tau_ij <- new_tau
    phi_ij <- new_phi
    accept <- T
  }
  
  return(list(S_update = S_ij, 
              P_update = P_ij,
              accept = accept, 
              dot_prod_update = dot_prod_ij,
              model_update = model_ij,
              tau_update = tau_ij, 
              phi_update = phi_ij))
}


# Metropolis Hastings step for phi according to Telesca, Inoue (2008)
# i.e. assume phi_ij ~ N(Y,SIGMA_phi)
# N.B. MH step is done component-wise using the marginal distribution of 
# the single component phi_ij. The proposal is a uniform distribution with 
# its support defined to guarantee monotonicity.
MH_step_phi <- function (curr_phi, curr_dot_prod_vec,
                         mat_basis_phi, basis_set_beta, 
                         curr_a, curr_c, curr_beta, 
                         curr_sigma2_eps, curr_sigma2_phi,
                         phi_mean_vec, P_matrix_inv,
                         Y_patient, patient_index,
                         iteration) 
{
  # Record whether we accept or not
  accept <- 0
  
  # Define length of curr_phi
  Q <- length(curr_phi)
  
  # N.B. Phi_1=0 and Phi_Q=1 ALWAYS!!!! -> add innovation only to phi_2, ..., phi_(Q-1) 
  # Moreover as Phi_1 and Phi_Q are fixed, they're not a r.v.
  # Therefore they can be excluded from the computations!!
  
  # Draw new values for the vector of Phi
  for(j in 2:(Q-1))
  {
    new_phi <- curr_phi
    
    # Bounds wrt new values
    # lower_bound <- new_phi[j-1]
    # upper_bound <- curr_phi[j+1]
    # upper_bound <- ifelse(upper_bound >= 1, 1, upper_bound)
    
    # Bounds with respect previous value
    #epsilon = (iteration)^(1/3)/sqrt(iteration+1) 
    epsilon = 0.05
    
    lower_bound <- curr_phi[j] - epsilon*(curr_phi[j] - curr_phi[j-1])
    upper_bound <- curr_phi[j] + epsilon*(curr_phi[j+1] - curr_phi[j])
    curr_range <- upper_bound - lower_bound
    new_phi[j] <- runif(n = 1, min = lower_bound, max = upper_bound)
    
    # Compute the new range for acceptance prob
    new_range <- epsilon*(curr_phi[j+1] - curr_phi[j-1])
    
    # Compute the new vector of tau for the current patient
    new_tau <- as.numeric(mat_basis_phi %*% new_phi)
    
    # Given the new value of S evaluate the model for the current patient
    curr_mat_basis_beta <- eval.basis(new_tau, basis_set_beta)
    new_model <- curr_c + curr_a * curr_mat_basis_beta %*% curr_beta 
    
    # Evaluate the new dot product and put it in the vector
    new_dot_prod_vec <- curr_dot_prod_vec
    new_dot_prod_vec[patient_index] <- t(Y_patient - new_model) %*% (Y_patient - new_model)
    
    # # Evaluate the new phi product vector
    # new_phi_prod <- t(new_phi) %*% P_matrix %*% new_phi
    phi_sd <- diag(curr_sigma2_phi * P_matrix_inv)[j]
    phi_mean <- phi_mean_vec[j]
    
    # Compute acceptance probability (equal for every S_l)
    # 1. Evaluate the full conditional with the current and the new values
    new_log_density <- -(1/(2*curr_sigma2_eps))*sum(new_dot_prod_vec) - 
      ((new_phi[j]-phi_mean)^2/(2*(phi_sd)^2))
    
    curr_log_density <- -(1/(2*curr_sigma2_eps))*sum(curr_dot_prod_vec) -
      ((curr_phi[j]-phi_mean)^2/(2*(phi_sd)^2))
    
    # 2. Compute density ratio
    log_density_ratio <- new_log_density - curr_log_density + log(curr_range/new_range)
    
    # 3. Compute accept prob
    accept_prob <- min(log_density_ratio, 0)
    accept_prob <- exp(accept_prob)
    
    # Check acceptance (for each)
    u_check <- runif(1, 0, 1)
    
    # # Save current value of phi_j for future use
    # prev_phi_j <- curr_phi[j]
    
    if (u_check <= accept_prob) 
    {
      curr_phi <- new_phi
      accept <- accept + 1
    }
  }
  
  return(list(accept = accept, 
              phi_update = curr_phi))
}

################################################################################
################################################################################
# Function to compute the inverse of tau_ij using interpolation
get_tau_inverse_interp <- function(tau_mat)
{
  # Matrix to store the results
  tau_inv_mat <- matrix(0, nrow = nrow(tau_mat), ncol = ncol(tau_mat))
  
  # Vector of linear times
  times_vect <- seq(0, 1, length.out = ncol(tau_mat)) 
  
  for(i in 1:nrow(tau_mat))
  {
    # Consider each vector tau_ij, interpolate it exchanging role of x and y 
    # and finally evaluate on the grid of times
    tau_inv_mat[i,] <- spline(x = tau_mat[i,], y = times_vect,
                              xout = times_vect)$y
  }
  return(tau_inv_mat)
}

################################################################################
# To compute the coefficients vj of the mean of the phi vector (Telesca)
get.mean.phi.coeff <- function(t1 = 0, tn = 1, n.knots, spline.order)
{
  # Define the interior knots (we assume equispaced knots)
  knots.set <- seq(t1, tn, length.out = n.knots+2)
  
  # Add the first r-1 t1 and last r-1 tn
  knots.set <- c(rep(t1, spline.order-1), knots.set, rep(tn, spline.order-1))
  
  # Recursively define the coefficients
  coeff.vect <- rep(0, n.knots + spline.order)
  coeff.vect[1] <- t1
  for(i in 1:(n.knots+spline.order-1))
  {
    coeff.vect[i+1] <- ((knots.set[i+spline.order] - knots.set[i+1])/
                          (spline.order-1)) + coeff.vect[i] 
  }
  
  return(coeff.vect)
}
