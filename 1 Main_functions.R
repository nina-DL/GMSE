# Latest version: 05-11-24
# This script contains necessary functions to generate simulations, compute gmse, etc. 

load("Data_gen_process.RData")

require(nnet)
library(fastDummies)
library(MASS)
library(svMisc) #progress

# Main function for data simulation: Sim_data --------------------------------

simulate_ISTAT_edu = function(N_registry, sample_prop, props, fit = mod_fit_ISTAT, seed_cov = 123, seed = 123){
  set.seed(seed_cov)
  Sim_data = data.frame(
    # the following are our covariates
    cleta_19_new = factor(sample(10:14, N_registry, prob = props$age_prop, rep = T),
                          levels = 10:14, labels = age),
    SESSO = factor(sample(1:2, N_registry, prob = props$SESSO_prop, rep = T),
                   levels = 1:2, labels = c("M", "F")),
    FL_ITA = factor(sample(0:1, N_registry, prob = props$FL_prop, rep = T),
                    levels = 0:1, labels = c("No", "Si")),
    TITSTU8_CENS11 = factor(sample(1:8, N_registry, prob = props$TITSTU11_prop, rep = T),
                            levels = 1:8, labels = titstu),
    # the following is our domain
    COD_PROV_RES = factor(sample(1:9, N_registry, prob = props$PROVres_prop, rep = T),
                          levels = 1:9, labels = 1:9),
    row.names = NULL
  )
  
  set.seed(seed)
  # the following determines our sample
  Sim_data$prob_incl_v2 = sample_prop
  Sim_data$lambda = sample(0:1, N_registry, c(1-sample_prop, sample_prop), rep = T) #moved outside from fixed data.frame
  p.true = predict(fit, newdata = Sim_data, "probs")
  p.true = p.true[,titstu]
  y.true = apply(p.true, 1, function(x) sample(1:8, 1, prob=x))
  
  # ensure that we have each category observed at least once in the sample
  while(length(table(y.true[Sim_data$lambda == 1])) < 8){
    y.true = apply(p.true, 1, function(x) sample(1:8, 1, prob=x))
  }
  
  Sim_data$y.true = factor(y.true, levels = 1:8, labels = titstu)
  
  return(list(Sim_data = Sim_data, p.true = p.true))
}

# Create the design matrix: X_design --------------------------------

get_design = function(Sim_data, intercept = T){
  
  # Define design matrix with dummy variables
  dati_B12 = dummy_cols(Sim_data, select_columns = c("cleta_19_new", "SESSO", "FL_ITA", "TITSTU8_CENS11"))
  
  if(intercept == T){
    dati_B12$X_interc = c(1)
    
    X_design = as.matrix(dati_B12[,c("X_interc", "cleta_19_new_29-39", "cleta_19_new_40-49", "cleta_19_new_50-69", "cleta_19_new_70+", 
                                     "SESSO_F", "FL_ITA_Si", "TITSTU8_CENS11_2 Literate but no education",
                                     "TITSTU8_CENS11_3 Primary", "TITSTU8_CENS11_4 Lower secondary", 
                                     "TITSTU8_CENS11_5 Upper secondary", "TITSTU8_CENS11_6 Bachelor degree or equivalent",
                                     "TITSTU8_CENS11_7 Master degree or equivalent", "TITSTU8_CENS11_8 PhD level")])
  } else {
    X_design = as.matrix(dati_B12[,c("cleta_19_new_0-28", "cleta_19_new_29-39", "cleta_19_new_40-49", "cleta_19_new_50-69", "cleta_19_new_70+", 
                                     "SESSO_F", "FL_ITA_Si", "TITSTU8_CENS11_2 Literate but no education",
                                     "TITSTU8_CENS11_3 Primary", "TITSTU8_CENS11_4 Lower secondary", 
                                     "TITSTU8_CENS11_5 Upper secondary", "TITSTU8_CENS11_6 Bachelor degree or equivalent",
                                     "TITSTU8_CENS11_7 Master degree or equivalent", "TITSTU8_CENS11_8 PhD level")])
    
  }
  
  return(X_design = X_design)
}

# Compute GMSE --------------------------------

# |-- V1: long (4 linearizations) ---------------------------------------------------

GMSE_long = function(Sim_data, X_design, ref_k = 8, domain = "all", env){
  
  # estimate the coefficients and p.hat
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  Sim_sample = Sim_data[Sim_data$lambda == 1, ]
  mod_est <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
  beta.hat = coef(mod_est)
  p.hat = predict(mod_est, newdata = Sim_data, "probs")
  p.hat = p.hat[,titstu] # check if titstu is available
  
  pi_incl = Sim_data$prob_incl_v2
  
  # Sigma_matrix
  Sigma_matrix = list()
  for(i in 1:env$N_registry){
    Sigma_matrix[[i]] = -p.hat[i,]%*%t(p.hat[i,])
    diag(Sigma_matrix[[i]]) = p.hat[i,]*(1-p.hat[i,])
  }
  
  # X_design = get_design(Sim_data = Sim_data)
  
  F_terms = list()
  for(k in 1:env$K){
    F_term_1 = c()
    for(k1 in 1:env$K){
      F_k = -X_design*p.hat[,k]*p.hat[,k1]
      F_term_1 = cbind(F_term_1, F_k)
    }
    F_term_1[,((k-1)*env$J+1):(k*env$J)] = X_design[,]*p.hat[,k]*(1-p.hat[,k])
    F_terms[[k]] = F_term_1
  }
  
  # now let's compute the G terms of Eq. (16)
  G_matrix = list()
  for(i in 1:env$N_registry){
    Gi_term = NULL
    for(k in 1:env$K){
      der.gi_term = rep(0, env$H) # error adjusted 23/02/2024; previously this line was outside the for loop
      der.gi_term[(env$J*(k-1)+1):(env$J*k)] = X_design[i,]
      Gi_term = cbind(Gi_term, der.gi_term)
    }
    G_matrix[[i]] = Gi_term
  }
  
  # A with lambda_i = pi_i
  A_matrix = NULL
  for(k in 1:env$K){
    A_row_blocks = NULL
    for(j in 1:env$K){
      A_matrix_blocks = matrix(0, nrow = env$J, ncol = env$J)
      for(i in 1:env$N_registry){
        A_matrix_blocks = -pi_incl[i]*X_design[i,]%*%t(X_design[i,])*Sigma_matrix[[i]][j,k] + A_matrix_blocks
      }
      
      A_row_blocks = rbind(A_row_blocks, A_matrix_blocks)
    }
    
    A_matrix = cbind(A_matrix, A_row_blocks)
  }
  
  #invA = solve(A_matrix)
  invA = ginv(A_matrix)
  
  # evaluate the big sum with Ui & dervUi (it includes invA, G and Sigma)
  # Define the inner big sum of the GMSE
  big_sum = c(0)
  for(i in 1:env$N_registry){
    
    # this is the first addendum in Eq. (20), with U_i as defined in eq. (16)
    U_i = pi_incl[i]*invA%*%G_matrix[[i]]
    add1 = U_i%*%Sigma_matrix[[i]]%*%t(U_i)
    
    #compute the derivative of Ai
    derAi_matrix = NULL
    for(k in 1:env$K){
      derA_row_blocks = NULL
      for(j in 1:env$K){
        derA_matrix_blocks = -X_design[i,]%*%t(X_design[i,])*Sigma_matrix[[i]][j,k]
        derA_row_blocks = rbind(derA_row_blocks, derA_matrix_blocks)
      }
      derAi_matrix = cbind(derAi_matrix, derA_row_blocks)
    }
    
    # this is the second addendum in Eq. (20), with the derivative of U_i computed according to Eq. (18)
    derUi_matrix = (invA%*%derAi_matrix%*%invA*pi_incl[i] - invA)%*%G_matrix[[i]]
    add2 = (derUi_matrix%*%Sigma_matrix[[i]]%*%t(derUi_matrix))*(pi_incl[i])*(1-pi_incl[i])
    
    sum_i = add1 + add2
    big_sum = big_sum + sum_i 
  }
  
  if(domain == "all"){
    gamma = rep(1, env$N_registry)
    theta.hat = apply(p.hat[gamma==1,], 2, sum)
    GMSE = c()
    for(k in 1:env$K){
      GMSE[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
    }
    CV = sqrt(GMSE)/theta.hat
  } else {
    GMSE = theta.hat = NULL
    for(d in levels(Sim_data[[domain]])){
      gamma = (Sim_data[[domain]]==d)*1
      theta.hat = rbind(theta.hat, apply(p.hat[gamma==1,], 2, sum))
      GMSE_d = c()
      for(k in 1:env$K){
        GMSE_d[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
      }
      GMSE = rbind(GMSE, GMSE_d)
    }
    rownames(GMSE) = levels(Sim_data[[domain]])
    colnames(GMSE) = titstu
    CV = sqrt(GMSE)/theta.hat
  }
  
  return(list(GMSE = GMSE, CV = CV, theta.hat = theta.hat))
}

# |-- V2: short (2 linearizations) ---------------------------------------------------

GMSE_short = function(Sim_data, X_design, ref_k = 8, domain = "all", env){
  
  # estimate the coefficients and p.hat
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  Sim_sample = Sim_data[Sim_data$lambda == 1, ]
  mod_est <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
  beta.hat = coef(mod_est)
  p.hat = predict(mod_est, newdata = Sim_data, "probs")
  p.hat = p.hat[,titstu] # check if titstu is available
  
  pi_incl = Sim_data$prob_incl_v2
  
  # Sigma_matrix
  Sigma_matrix = list()
  for(i in 1:env$N_registry){
    Sigma_matrix[[i]] = -p.hat[i,]%*%t(p.hat[i,])
    diag(Sigma_matrix[[i]]) = p.hat[i,]*(1-p.hat[i,])
  }
  
  F_terms = list()
  for(k in 1:env$K){
    F_term_1 = c()
    for(k1 in 1:env$K){
      F_k = -X_design*p.hat[,k]*p.hat[,k1]
      F_term_1 = cbind(F_term_1, F_k)
    }
    F_term_1[,((k-1)*env$J+1):(k*env$J)] = X_design[,]*p.hat[,k]*(1-p.hat[,k])
    F_terms[[k]] = F_term_1
  }
  
  # now let's compute the G terms of Eq. (16)
  G_matrix = list()
  for(i in 1:env$N_registry){
    Gi_term = NULL
    for(k in 1:env$K){
      der.gi_term = rep(0, env$H) # error adjusted 23/02/2024; previously this line was outside the for loop
      der.gi_term[(env$J*(k-1)+1):(env$J*k)] = X_design[i,]
      Gi_term = cbind(Gi_term, der.gi_term)
    }
    G_matrix[[i]] = Gi_term
  }
  
  # A with lambda_i = pi_i
  A_matrix = NULL
  for(k in 1:env$K){
    A_row_blocks = NULL
    for(j in 1:env$K){
      A_matrix_blocks = matrix(0, nrow = env$J, ncol = env$J)
      for(i in 1:env$N_registry){
        A_matrix_blocks = -pi_incl[i]*X_design[i,]%*%t(X_design[i,])*Sigma_matrix[[i]][j,k] + A_matrix_blocks
      }
      
      A_row_blocks = rbind(A_row_blocks, A_matrix_blocks)
    }
    
    A_matrix = cbind(A_matrix, A_row_blocks)
  }
  
  #invA = solve(A_matrix)
  invA = ginv(A_matrix)
  
  # evaluate the big sum with Ui (it includes invA, G and Sigma)
  big_sum = c(0)
  for(i in 1:env$N_registry){
    
    # this is according to Piero version
    U_i = invA%*%G_matrix[[i]]
    sum_i = U_i%*%Sigma_matrix[[i]]%*%t(U_i)*pi_incl[i]
    
    big_sum = big_sum + sum_i 
  }
  
  if(domain == "all"){
    gamma = rep(1, env$N_registry)
    theta.hat = apply(p.hat[gamma==1,], 2, sum)
    GMSE = c()
    for(k in 1:env$K){
      GMSE[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
    }
    CV = sqrt(GMSE)/theta.hat
  } else {
    GMSE = theta.hat = NULL
    for(d in levels(Sim_data[[domain]])){
      gamma = (Sim_data[[domain]]==d)*1
      theta.hat = rbind(theta.hat, apply(p.hat[gamma==1,], 2, sum))
      GMSE_d = c()
      for(k in 1:env$K){
        GMSE_d[k] = gamma%*%F_terms[[k]]%*%big_sum%*%t(F_terms[[k]])%*%gamma
      }
      GMSE = rbind(GMSE, GMSE_d)
    }
    rownames(GMSE) = levels(Sim_data[[domain]])
    colnames(GMSE) = titstu
    CV = sqrt(GMSE)/theta.hat
  }
  
  return(list(GMSE = GMSE, CV = CV, theta.hat = theta.hat))
}

# Give it a try to GMSE V1 vs V2-----------------------------------------------------------

wannatry = F
if(wannatry==T){
  # Set hyperparameters and give it a try
  N_registry = 30000
  Sim_res = simulate_ISTAT_edu(N_registry = N_registry, sample_prop = sample_prop, props = props, fit = mod_fit_ISTAT, seed = 123)
  Sim_data = Sim_res$Sim_data
  X_design = get_design(Sim_data = Sim_data, intercept = T)
  J = length(beta_true[1,])
  K = length(beta_true[,1]) + 1
  H = J*K
  
  GMSE_V1_resS1$Time = system.time(GMSE_V1_resS1 <- GMSE_long(Sim_data = Sim_data, X_design = X_design, env = environment()))
  GMSE_V2_resS1$Time = system.time(GMSE_V2_resS1 <- GMSE_short(Sim_data = Sim_data, X_design = X_design, env = environment()))
  
  cbind(GMSE_V1_resS1$GMSE, GMSE_V2_resS1$GMSE)
  cbind(GMSE_V1_resS1$CV, GMSE_V2_resS1$CV)
  cbind(GMSE_V1_resS1$theta.hat, GMSE_V2_resS1$theta.hat)
  cbind(GMSE_V1_resS1$Time, GMSE_V2_resS1$Time)
}

# |-- MC ---------------------------------------------------

MC_Accuracy = function(Sim_data, X_design, ref_k = 8, gamma = "all", 
                       M_design = 10, M_model = 10, env){
  
  theta = table(Sim_data$y.true)[titstu]
  
  # initiate simulations
  Em_theta = Vm_theta = Em_theta.hat = Em_theta.hat_II = Vm_theta.hat = Vm_theta.hat_II = #Biasm_hat_tilde = Biasm_hat_tilde_II = 
    MSEm_hat_tilde = MSEm_hat_tilde_II = MSEm_hat_theta = MSEm_hat_theta_II = MSEm_hat_theta_III = matrix(NA, nrow = M_design, ncol = env$K)
    #Termine1_Em_Yhat_Ytilde2 = Termine2_Em_Ytilde_Yhat2 = Termine2_Vm_Y = Termine3_Em_Yhat_Ytilde_Ytilde_Y 
  
  for(i in 1:M_design){
    set.seed(i)
    
    # Step 1: variability due to the sampling/design process
    # Note: for now we assume same prob for all individuals
    idx_sample = sample(env$N_registry, env$n_sample) 
    
    # We can also switch to the inclusion probs indicated by the db
    # idx_sample = sample(1:env$N_registry, env$n_sample, prob = Sim_data$prob_incl_v2) 
    
    theta.hat_IA = theta.hat_IIA = theta.hat_IIIA = theta_mc = matrix(NA, nrow = M_model, ncol = env$K)
    for(m in 1:M_model){
      
      # Step 2. Variability due to model: regenerate y based on its assumed model (true beta)
      y.true_mc = apply(env$p.true, 1, function(x) sample(1:(env$K), 1, prob=x))
      # ensure that we have each category observed at least once in the sample, otherwise it returns an error when fitting the model
      while(length(table(y.true_mc[idx_sample])) < 8){
        y.true_mc = apply(env$p.true, 1, function(x) sample(1:8, 1, prob=x))
      }
      
      Sim_data$y.true_mc = factor(y.true_mc, levels = 1:(env$K), labels = levels(Sim_data$y.true))
      # again, we take as Y of reference the same as in the main model fit 
      Sim_data$y.true_mc_ref <- relevel(Sim_data$y.true_mc, ref = ref_k)
      
      # Now use the i-th run sample to get estimates and evaluate errors
      Sim_sample = Sim_data[idx_sample, ]
      
      # estimate the p.hat
      mod_est_MCsample <- multinom(y.true_mc_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
      
      # Option I. Y.hat = p.hat (deterministic) + We re-estimate Y = p for all the registry (included sample)
      p.hat = predict(mod_est_MCsample, newdata = Sim_data, "probs")
      p.hat = p.hat[,titstu] # check if titstu is available
      theta.hat_IA[m,] = apply(p.hat, 2, sum) # Estimator as defined in Eq. (2)
      
      # # Option IB. Y.hat \in {0,1} ~ multin(p) (stochastic). Note: In the paper we use Option I, so we ignore this
      # Y.hat = apply(p.hat, 1, function(x) sample(1:(env$K), 1, prob=x))
      # Y.hat = factor(Y.hat, levels = 1:(env$K), labels = levels(Sim_data$y.true_ref))
      # theta.hat_IB[m,] = table(Y.hat)
      
      # Option II. Y.hat = p.hat (deterministic) + We estimate Y = p only for the registry but sample (for which we assume to know the truth)
      p.hat = p.hat[-idx_sample,]
      p.hat = p.hat[,titstu] # check if titstu is available
      theta.hat_IIA[m,] = apply(p.hat, 2, sum) + apply(env$p.true[idx_sample, ], 2, sum)
      theta.hat_IIIA[m,] = apply(p.hat, 2, sum) + table(Sim_sample$y.true_mc_ref)[titstu]
      # # Option IB
      # Y.hat = Y.hat[-idx_sample]
      # theta.hat_IIB[m,] = table(Y.hat) + table(Sim_data$y.true_ref[idx_sample])
      
      # This helps to evaluate the variability of the population parameter Y_sum (minor order term in GMSE 13a)
      theta_mc[m,] = table(Sim_data$y.true_mc_ref)[titstu]
    }
    
    Em_theta[i,] = apply(theta_mc, 2, mean) # expectation of theta = sumY under M
    Vm_theta[i,] = apply(theta_mc, 2, var) # variance of theta = sumY under M
    
    Em_theta.hat[i,] = apply(theta.hat_IA, 2, mean) # expectation under M of theta.hat = sumYhat
    Em_theta.hat_II[i,] = apply(theta.hat_IIA, 2, mean) # expectation under M of theta.hat = sumYhat
    
    Vm_theta.hat[i,] = apply(theta.hat_IA, 2, var) # variance under M of theta.hat = sumYhat
    Vm_theta.hat_II[i,] = apply(theta.hat_IIA, 2, var) # variance under M of theta.hat = sumYhat
    
    MSEm_hat_theta[i,] = apply(apply(theta.hat_IA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat = sumYhat, theta = sumY) under M
    MSEm_hat_theta_II[i,] = apply(apply(theta.hat_IIA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat_II = sumYhat_II, theta = sumY) under M
    MSEm_hat_theta_III[i,] = apply(apply(theta.hat_IIIA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat_II = sumYhat_II, theta = sumY) under M
    
    # we should not need this
    MSEm_hat_tilde[i,] = apply(apply(theta.hat_IA, 1, function(x) (x-env$theta.tilde)^2), 1, mean) # MSE(theta.hat = sumYhat, theta.tilde = sump) under M
    MSEm_hat_tilde_II[i,] = apply(apply(theta.hat_IIA, 1, function(x) (x-env$theta.tilde)^2), 1, mean) # MSE(theta.hat_II = sumYhat, theta.tilde = sump) under M
  
    progress(i,M_design)
  }
  
  # "Global Variance" - Eq (8)
  GVar = apply(Vm_theta.hat, 2, mean)
  GVar_II = apply(Vm_theta.hat_II, 2, mean)
  
  # "Global bias"
  GBias_hat_theta = apply(Em_theta.hat, 2, mean) - apply(Em_theta, 2, mean)
  GBias_hat_theta_II = apply(Em_theta.hat_II, 2, mean) - apply(Em_theta, 2, mean)
  
  # Computation of GMSE as per definition - Eq (7)
  GMSE_hat_theta = apply(MSEm_hat_theta, 2, mean)
  GMSE_hat_theta_II = apply(MSEm_hat_theta_II, 2, mean)
  GMSE_hat_theta_III = apply(MSEm_hat_theta_III, 2, mean)
  
  # Now compare with formula in 13a
  # Note that we have the Var component that enters with a negative sign: true_Vm_theta
  GMSE_hat_theta_13a = GVar - apply(Vm_theta, 2, mean)
  GMSE_hat_theta_13a_II = GVar_II - apply(Vm_theta, 2, mean)
  # GMSE_hat_theta_13a_trueVm = GVar - env$true_Vm_theta # in practice, true_Vm_theta is unknown
  # GMSE_hat_theta_13a_trueVm_II = GVar_II - env$true_Vm_theta # in practice, true_Vm_theta is unknown
  
  # Riunione Piero 29 luglio: prova a mettere i p_ik veri e vedere la differenza
  CV_GVar = sqrt(GVar)/env$theta.tilde
  CV_GVar_II = sqrt(GVar_II)/env$theta.tilde
  CV_hat_theta = sqrt(GMSE_hat_theta)/env$theta.tilde
  CV_hat_theta_II = sqrt(GMSE_hat_theta_II)/env$theta.tilde
  CV_hat_theta_III = sqrt(GMSE_hat_theta_III)/env$theta.tilde
  CV_hat_theta_13a = sqrt(GMSE_hat_theta_13a)/env$theta.tilde
  CV_hat_theta_13a_II = sqrt(GMSE_hat_theta_13a_II)/env$theta.tilde
  # CV_hat_theta_13a_trueVm = sqrt(GMSE_hat_theta_13a_trueVm)/apply(Em_theta.hat, 2, mean)
  # CV_hat_theta_13a_trueVm_II = sqrt(GMSE_hat_theta_13a_trueVm_II)/apply(Em_theta.hat, 2, mean)
  
  # overall look
  Res = rbind(GVar, GVar_II, GBias_hat_theta, GBias_hat_theta_II, GMSE_hat_theta, GMSE_hat_theta_II, GMSE_hat_theta_III,
                 GMSE_hat_theta_13a, GMSE_hat_theta_13a_II, #GMSE_hat_theta_13a_trueVm, GMSE_hat_theta_13a_trueVm_II,
              CV_GVar, CV_GVar_II, CV_hat_theta, CV_hat_theta_II, CV_hat_theta_III, CV_hat_theta_13a, CV_hat_theta_13a_II#, CV_hat_theta_13a_trueVm, CV_hat_theta_13a_trueVm_II
              )
  
  return(Res)
}


# |-- Bootstrap ---------------------------------------------------

Boot_Accuracy = function(Sim_data, X_design, ref_k = 8, gamma = "all",
                         B = 100, env){
  
  if(is.numeric(gamma) == F){
    gamma = rep(1, env$N_registry)
    #gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0)
  }
  
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  Sim_data$gamma = gamma
  Sim_data = cbind(Sim_data, env$p.true)
  
  theta.hat_IA = theta.hat_IIA = theta.hat_IIIA = theta_boot = matrix(NA, nrow = B, ncol = env$K)
  for(b in 1:B){
    set.seed(b)
    
    # Step 0: Bootstrap the register to capture population variability
    B_idx_pop = sample(env$N_registry, env$N_registry, replace = T) 
    Boot_data = Sim_data[B_idx_pop,]
    #Boot_p.true = env$p.true[B_idx_pop,]
    
    # Step 1: variability due to the sampling/design process
    B_idx_sample = which(Boot_data$lambda==1)
    
    # Now use the b-th run to get estimates and evaluate errors
    Boot_sample = Boot_data[B_idx_sample, ]
    
    # ensure that we have each category observed at least once in the sample
    while(any(table(Boot_sample$y.true)==0)){
      B_idx_pop = sample(env$N_registry, env$N_registry, replace = T) 
      Boot_data = Sim_data[B_idx_pop,]
      B_idx_sample = which(Boot_data$lambda==1)
      Boot_sample = Boot_data[B_idx_sample, ]
    }
    
    # estimate the p.hat
    mod_est_Bsample <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Boot_sample)
    
    # Option I. We can re-estimate the parameters for all the registry (included sample)
    p.hat = predict(mod_est_Bsample, newdata = Boot_data, "probs") # nd. here newdata = Sim_data OR newdata = Boot_data
    p.hat = p.hat[,titstu]
    theta.hat_IA[b,] = apply(p.hat[Boot_data$gamma==1,], 2, sum)
    
    # Y.hat = apply(p.hat, 1, function(x) sample(1:(env$K), 1, prob=x))
    # Y.hat = factor(Y.hat, levels = 1:(env$K), labels = levels(Sim_data$y.true_ref))
    # theta.hat_IB[m,] = table(Y.hat)
    
    # Option II. Or we estimate the parameters only for the registry but sample (for which we assume to know the truth)
    p.hat = p.hat[-B_idx_sample,]
    p.hat = p.hat[,titstu]
    theta.hat_IIA[b,] = apply(p.hat[Boot_data$gamma[-B_idx_sample]==1,], 2, sum) + apply(Boot_data[(Boot_data$lambda==1)&(Boot_data$gamma==1), 11:18], 2, sum)
    theta.hat_IIIA[b,] = apply(p.hat[Boot_data$gamma[-B_idx_sample]==1,], 2, sum) + table(Boot_sample$y.true_ref[(Boot_data$lambda==1)&(Boot_data$gamma==1)])[titstu]
    #theta.hat_IIA[b,] = apply(p.hat, 2, sum) + apply(env$p.true[B_idx_sample, ], 2, sum) # old 7-10-24
    
    # Y.hat = Y.hat[-idx_sample]
    # theta.hat_IIB[m,] = table(Y.hat) + table(Sim_data$y.true_ref[idx_sample])
    
    # This helps to evaluate the variability of pure model variability of population parameter [Y (sum); minor order term in GMSE]
    theta_boot[b,] = table(Boot_data$y.true_ref[Boot_data$gamma==1])[titstu]
    
    progress(b,B)
  }
  
  Em_theta = apply(theta_boot, 2, mean) # expectation of theta = sumY under M
  Vm_theta = apply(theta_boot, 2, var) # variance of theta = sumY under M
  
  Em_theta.hat = apply(theta.hat_IA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  Em_theta.hat_II = apply(theta.hat_IIA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  Em_theta.hat_III = apply(theta.hat_IIIA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  
  # "Global Variance" - Eq (8)
  GVar = Vm_theta.hat = apply(theta.hat_IA, 2, var) # variance under both M and Pi of theta.hat = sumYhat
  GVar_II = Vm_theta.hat_II = apply(theta.hat_IIA, 2, var) # variance under both M and Pi of theta.hat = sumYhat
  
  # "Global bias"
  GBias_hat_theta = Em_theta.hat - Em_theta
  GBias_hat_theta_II = Em_theta.hat_II - Em_theta
  
  # Computation of GMSE as per definition - Eq (7)
  GMSE_hat_theta = apply(((theta.hat_IA-theta_boot)^2), 2, mean)
  GMSE_hat_theta_II = apply(((theta.hat_IIA-theta_boot)^2), 2, mean)
  GMSE_hat_theta_III = apply(((theta.hat_IIIA-theta_boot)^2), 2, mean)
  
  # Now compare with formula in 13a
  # Note that we have the Var component that enters with a negative sign: true_Vm_theta
  GMSE_hat_theta_13a = GVar - Vm_theta
  GMSE_hat_theta_13a_II = GVar_II - Vm_theta
  # GMSE_hat_theta_13a_trueVm = GVar - env$true_Vm_theta
  # GMSE_hat_theta_13a_trueVm_II = GVar_II - env$true_Vm_theta
  
  CV_GVar = sqrt(GVar)/Em_theta.hat
  CV_GVar_II = sqrt(GVar_II)/Em_theta.hat_II
  CV_hat_theta = sqrt(GMSE_hat_theta)/Em_theta.hat
  CV_hat_theta_II = sqrt(GMSE_hat_theta_II)/Em_theta.hat_II
  CV_hat_theta_III = sqrt(GMSE_hat_theta_III)/Em_theta.hat_III
  CV_hat_theta_13a = sqrt(GMSE_hat_theta_13a)/Em_theta.hat
  CV_hat_theta_13a_II = sqrt(GMSE_hat_theta_13a_II)/Em_theta.hat_II
  # CV_hat_theta_13a_trueVm = sqrt(GMSE_hat_theta_13a_trueVm)/Em_theta.hat
  # CV_hat_theta_13a_trueVm_II = sqrt(GMSE_hat_theta_13a_trueVm_II)/Em_theta.hat
  
  # overall look
  Res = rbind(GVar, GVar_II, GBias_hat_theta, GBias_hat_theta_II, GMSE_hat_theta, GMSE_hat_theta_II, GMSE_hat_theta_III,
              GMSE_hat_theta_13a, GMSE_hat_theta_13a_II, #GMSE_hat_theta_13a_trueVm, GMSE_hat_theta_13a_trueVm_II,
              CV_GVar, CV_GVar_II, CV_hat_theta, CV_hat_theta_II, CV_hat_theta_III, CV_hat_theta_13a, CV_hat_theta_13a_II#, CV_hat_theta_13a_trueVm, CV_hat_theta_13a_trueVm_II
              )
  
  return(Res)
}

# |-- Bootstrap v2 ---------------------------------------------------

Boot_Accuracy_v2 = function(Sim_data, X_design, ref_k = 8, gamma = "all",
                         B = 100, env){
  
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  
  theta.hat_IA = theta.hat_IIA = theta.hat_IIIA = theta_boot = matrix(NA, nrow = B, ncol = env$K)
  for(b in 1:B){
    set.seed(b)
    
    # Step 0: Bootstrap the register to capture population variability
    # Step 1: variability due to the sampling/design process
    idx_sample = which(Sim_data$lambda==1)
    B_idx_sample = sample(idx_sample, env$n_sample, replace = T) 
    
    # Now use the b-th run to get estimates and evaluate errors
    Boot_sample = Sim_data[B_idx_sample, ]
    
    # ensure that we have each category observed at least once in the sample
    while(any(table(Boot_sample$y.true)==0)){
      B_idx_sample = sample(idx_sample, env$n_sample, replace = T) 
      Boot_sample = Sim_data[B_idx_sample, ]
    }
    
    Boot_data = rbind(Sim_data[-idx_sample, ], Boot_sample)
    Boot_idx = c(dim(Sim_data[-idx_sample, ])[1], dim(Boot_sample)[1])
    
    # estimate the p.hat
    mod_est_MCsample <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Boot_sample)
    
    # Option I. We can re-estimate the parameters for all the registry (included sample)
    p.hat = predict(mod_est_MCsample, newdata = Boot_data, "probs")
    p.hat = p.hat[,titstu]
    theta.hat_IA[b,] = apply(p.hat, 2, sum)
    
    # Y.hat = apply(p.hat, 1, function(x) sample(1:(env$K), 1, prob=x))
    # Y.hat = factor(Y.hat, levels = 1:(env$K), labels = levels(Sim_data$y.true_ref))
    # theta.hat_IB[m,] = table(Y.hat)
    
    # Option II. Or we estimate the parameters only for the registry but sample (for which we assume to know the truth)
    p.hat = p.hat[1:Boot_idx[1],]
    p.hat = p.hat[,titstu]
    theta.hat_IIA[b,] = apply(p.hat, 2, sum) + apply(env$p.true[B_idx_sample, ], 2, sum)
    theta.hat_IIIA[b,] = apply(p.hat, 2, sum) + table(Boot_sample$y.true_ref)[titstu]
    # Y.hat = Y.hat[-idx_sample]
    # theta.hat_IIB[m,] = table(Y.hat) + table(Sim_data$y.true_ref[idx_sample])
    
    # This helps to evaluate the variability of pure model variability of population parameter [Y (sum); minor order term in GMSE]
    theta_boot[b,] = table(Boot_data$y.true_ref)[titstu]
    
    progress(b,B)
  }
  
  Em_theta = apply(theta_boot, 2, mean) # expectation of theta = sumY under M
  Vm_theta = apply(theta_boot, 2, var) # variance of theta = sumY under M
  
  Em_theta.hat = apply(theta.hat_IA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  Em_theta.hat_II = apply(theta.hat_IIA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  Em_theta.hat_III = apply(theta.hat_IIIA, 2, mean) # expectation under both M and Pi of theta.hat = sumYhat
  
  # "Global Variance" - Eq (8)
  GVar = Vm_theta.hat = apply(theta.hat_IA, 2, var) # variance under both M and Pi of theta.hat = sumYhat
  GVar_II = Vm_theta.hat_II = apply(theta.hat_IIA, 2, var) # variance under both M and Pi of theta.hat = sumYhat
  
  # "Global bias"
  GBias_hat_theta = Em_theta.hat - Em_theta
  GBias_hat_theta_II = Em_theta.hat_II - Em_theta
  
  # Computation of GMSE as per definition - Eq (7)
  GMSE_hat_theta = apply(((theta.hat_IA-theta_boot)^2), 2, mean)
  GMSE_hat_theta_II = apply(((theta.hat_IIA-theta_boot)^2), 2, mean)
  GMSE_hat_theta_III = apply(((theta.hat_IIIA-theta_boot)^2), 2, mean)
  
  # Now compare with formula in 13a
  # Note that we have the Var component that enters with a negative sign: true_Vm_theta
  GMSE_hat_theta_13a = GVar - Vm_theta
  GMSE_hat_theta_13a_II = GVar_II - Vm_theta
  # GMSE_hat_theta_13a_trueVm = GVar - env$true_Vm_theta
  # GMSE_hat_theta_13a_trueVm_II = GVar_II - env$true_Vm_theta
  
  CV_GVar = sqrt(GVar)/Em_theta.hat
  CV_GVar_II = sqrt(GVar_II)/Em_theta.hat_II
  CV_hat_theta = sqrt(GMSE_hat_theta)/Em_theta.hat
  CV_hat_theta_II = sqrt(GMSE_hat_theta_II)/Em_theta.hat_II
  CV_hat_theta_III = sqrt(GMSE_hat_theta_III)/Em_theta.hat_III
  CV_hat_theta_13a = sqrt(GMSE_hat_theta_13a)/Em_theta.hat
  CV_hat_theta_13a_II = sqrt(GMSE_hat_theta_13a_II)/Em_theta.hat_II
  # CV_hat_theta_13a_trueVm = sqrt(GMSE_hat_theta_13a_trueVm)/Em_theta.hat
  # CV_hat_theta_13a_trueVm_II = sqrt(GMSE_hat_theta_13a_trueVm_II)/Em_theta.hat
  
  # overall look
  Res = rbind(GVar, GVar_II, GBias_hat_theta, GBias_hat_theta_II, GMSE_hat_theta, GMSE_hat_theta_II, GMSE_hat_theta_III,
              GMSE_hat_theta_13a, GMSE_hat_theta_13a_II, #GMSE_hat_theta_13a_trueVm, GMSE_hat_theta_13a_trueVm_II,
              CV_GVar, CV_GVar_II, CV_hat_theta, CV_hat_theta_II, CV_hat_theta_III, CV_hat_theta_13a, CV_hat_theta_13a_II#, CV_hat_theta_13a_trueVm, CV_hat_theta_13a_trueVm_II
  )
  
  return(Res)
}

# |-- Bootstrap Parametric ---------------------------------------------------

PBoot_Accuracy = function(Sim_data, X_design, ref_k = 8, gamma = "all", 
                          M_design = 10, M_model = 10, env){
  
  theta = table(Sim_data$y.true)[titstu]
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  Sim_sample = Sim_data[Sim_data$lambda==1, ]
  mod_est <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
  p.est = predict(mod_est, newdata = Sim_data, "probs")
  p.est = p.est[,titstu]
  theta.est = apply(p.est, 2, sum)
  
  # initiate simulations
  Em_theta = Vm_theta = Em_theta.hat = Em_theta.hat_II = Vm_theta.hat = Vm_theta.hat_II = #Biasm_hat_tilde = Biasm_hat_tilde_II = 
    MSEm_hat_tilde = MSEm_hat_tilde_II = MSEm_hat_theta = MSEm_hat_theta_II = MSEm_hat_theta_III = matrix(NA, nrow = M_design, ncol = env$K)
  #Termine1_Em_Yhat_Ytilde2 = Termine2_Em_Ytilde_Yhat2 = Termine2_Vm_Y = Termine3_Em_Yhat_Ytilde_Ytilde_Y 
  
  for(i in 1:M_design){
    set.seed(i)
    
    # Step 1: variability due to the sampling/design process
    # Note: for now we assume same prob for all individuals
    idx_sample = sample(env$N_registry, env$n_sample) 
    
    # We can also switch to the inclusion probs indicated by the db
    # idx_sample = sample(1:env$N_registry, env$n_sample, prob = Sim_data$prob_incl_v2) 
    
    theta.hat_IA = theta.hat_IIA = theta.hat_IIIA = theta_mc = matrix(NA, nrow = M_model, ncol = env$K)
    for(m in 1:M_model){
      
      # Step 2. Variability due to model: regenerate y based on its assumed model (true beta)
      y.true_mc = apply(p.est, 1, function(x) sample(1:(env$K), 1, prob=x))
      # ensure that we have each category observed at least once in the sample, otherwise it returns an error when fitting the model
      while(length(table(y.true_mc[idx_sample])) < 8){
        y.true_mc = apply(p.est, 1, function(x) sample(1:8, 1, prob=x))
      }
      
      Sim_data$y.true_mc = factor(y.true_mc, levels = 1:(env$K), labels = levels(Sim_data$y.true))
      # again, we take as Y of reference the same as in the main model fit 
      Sim_data$y.true_mc_ref <- relevel(Sim_data$y.true_mc, ref = ref_k)
      
      # Now use the i-th run sample to get estimates and evaluate errors
      Sim_sample = Sim_data[idx_sample, ]
      
      # estimate the p.hat
      mod_est_MCsample <- multinom(y.true_mc_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
      
      # Option I. Y.hat = p.hat (deterministic) + We re-estimate Y = p for all the registry (included sample)
      p.hat = predict(mod_est_MCsample, newdata = Sim_data, "probs")
      p.hat = p.hat[,titstu] # check if titstu is available
      theta.hat_IA[m,] = apply(p.hat, 2, sum) # Estimator as defined in Eq. (2)
      
      # # Option IB. Y.hat \in {0,1} ~ multin(p) (stochastic). Note: In the paper we use Option I, so we ignore this
      # Y.hat = apply(p.hat, 1, function(x) sample(1:(env$K), 1, prob=x))
      # Y.hat = factor(Y.hat, levels = 1:(env$K), labels = levels(Sim_data$y.true_ref))
      # theta.hat_IB[m,] = table(Y.hat)
      
      # Option II. Y.hat = p.hat (deterministic) + We estimate Y = p only for the registry but sample (for which we assume to know the truth)
      p.hat = p.hat[-idx_sample,]
      p.hat = p.hat[,titstu] # check if titstu is available
      theta.hat_IIA[m,] = apply(p.hat, 2, sum) + apply(env$p.true[idx_sample, ], 2, sum)
      theta.hat_IIIA[m,] = apply(p.hat, 2, sum) + table(Sim_sample$y.true_mc_ref)[titstu]
      # # Option IB
      # Y.hat = Y.hat[-idx_sample]
      # theta.hat_IIB[m,] = table(Y.hat) + table(Sim_data$y.true_ref[idx_sample])
      
      # This helps to evaluate the variability of the population parameter Y_sum (minor order term in GMSE 13a)
      theta_mc[m,] = table(Sim_data$y.true_mc_ref)[titstu]
    }
    
    Em_theta[i,] = apply(theta_mc, 2, mean) # expectation of theta = sumY under M
    Vm_theta[i,] = apply(theta_mc, 2, var) # variance of theta = sumY under M
    
    Em_theta.hat[i,] = apply(theta.hat_IA, 2, mean) # expectation under M of theta.hat = sumYhat
    Em_theta.hat_II[i,] = apply(theta.hat_IIA, 2, mean) # expectation under M of theta.hat = sumYhat
    
    Vm_theta.hat[i,] = apply(theta.hat_IA, 2, var) # variance under M of theta.hat = sumYhat
    Vm_theta.hat_II[i,] = apply(theta.hat_IIA, 2, var) # variance under M of theta.hat = sumYhat
    
    MSEm_hat_theta[i,] = apply(apply(theta.hat_IA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat = sumYhat, theta = sumY) under M
    MSEm_hat_theta_II[i,] = apply(apply(theta.hat_IIA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat_II = sumYhat_II, theta = sumY) under M
    MSEm_hat_theta_III[i,] = apply(apply(theta.hat_IIIA, 1, function(x) (x-theta)^2), 1, mean) # MSE(theta.hat_II = sumYhat_II, theta = sumY) under M
    
    # we should not need this
    MSEm_hat_tilde[i,] = apply(apply(theta.hat_IA, 1, function(x) (x-theta.est)^2), 1, mean) # MSE(theta.hat = sumYhat, theta.tilde = sump) under M
    MSEm_hat_tilde_II[i,] = apply(apply(theta.hat_IIA, 1, function(x) (x-theta.est)^2), 1, mean) # MSE(theta.hat_II = sumYhat, theta.tilde = sump) under M
    
    progress(i,M_design)
  }
  
  # "Global Variance" - Eq (8)
  GVar = apply(Vm_theta.hat, 2, mean)
  GVar_II = apply(Vm_theta.hat_II, 2, mean)
  
  # "Global bias"
  GBias_hat_theta = apply(Em_theta.hat, 2, mean) - apply(Em_theta, 2, mean)
  GBias_hat_theta_II = apply(Em_theta.hat_II, 2, mean) - apply(Em_theta, 2, mean)
  
  # Computation of GMSE as per definition - Eq (7)
  GMSE_hat_theta = apply(MSEm_hat_theta, 2, mean)
  GMSE_hat_theta_II = apply(MSEm_hat_theta_II, 2, mean)
  GMSE_hat_theta_III = apply(MSEm_hat_theta_III, 2, mean)
  
  # Now compare with formula in 13a
  # Note that we have the Var component that enters with a negative sign: true_Vm_theta
  GMSE_hat_theta_13a = GVar - apply(Vm_theta, 2, mean)
  GMSE_hat_theta_13a_II = GVar_II - apply(Vm_theta, 2, mean)
  # GMSE_hat_theta_13a_trueVm = GVar - env$true_Vm_theta # in practice, true_Vm_theta is unknown
  # GMSE_hat_theta_13a_trueVm_II = GVar_II - env$true_Vm_theta # in practice, true_Vm_theta is unknown
  
  # Riunione Piero 29 luglio: prova a mettere i p_ik veri e vedere la differenza
  CV_GVar = sqrt(GVar)/theta.est
  CV_GVar_II = sqrt(GVar_II)/theta.est
  CV_hat_theta = sqrt(GMSE_hat_theta)/theta.est
  CV_hat_theta_II = sqrt(GMSE_hat_theta_II)/theta.est
  CV_hat_theta_III = sqrt(GMSE_hat_theta_III)/theta.est
  CV_hat_theta_13a = sqrt(GMSE_hat_theta_13a)/theta.est
  CV_hat_theta_13a_II = sqrt(GMSE_hat_theta_13a_II)/theta.est
  # CV_hat_theta_13a_trueVm = sqrt(GMSE_hat_theta_13a_trueVm)/apply(Em_theta.hat, 2, mean)
  # CV_hat_theta_13a_trueVm_II = sqrt(GMSE_hat_theta_13a_trueVm_II)/apply(Em_theta.hat, 2, mean)
  
  # overall look
  Res = rbind(GVar, GVar_II, GBias_hat_theta, GBias_hat_theta_II, GMSE_hat_theta, GMSE_hat_theta_II, GMSE_hat_theta_III,
              GMSE_hat_theta_13a, GMSE_hat_theta_13a_II, #GMSE_hat_theta_13a_trueVm, GMSE_hat_theta_13a_trueVm_II,
              CV_GVar, CV_GVar_II, CV_hat_theta, CV_hat_theta_II, CV_hat_theta_III, CV_hat_theta_13a, CV_hat_theta_13a_II#, CV_hat_theta_13a_trueVm, CV_hat_theta_13a_trueVm_II
  )
  
  return(Res)
}



# Give it a try to All -----------------------------------------------------------


# Note: check reproducibility
wannatry = F

if(wannatry==T){
  # Set hyperparameters and give it a try
  N_registry = 50000
  Sim_res = simulate_ISTAT_edu(N_registry = N_registry, sample_prop = sample_prop, props = props, fit = mod_fit_ISTAT, seed = 123)
  Sim_data = Sim_res$Sim_data
  n_sample = sum(Sim_data$lambda)
  X_design = get_design(Sim_data = Sim_data, intercept = T)
  J = length(beta_true[1,])
  K = length(beta_true[,1]) + 1
  H = J*K
  
  # true coef: p.true
  p.true = predict(mod_fit_ISTAT, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true, 2, sum)
  #true_Vm_theta = apply(p.true, 2, function(x) sum(x*(1-x)))
  
  resMC = MC_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                      M_design = 10, M_model = 10, env = environment())
  resBoot = Boot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                                                 B = 100, env = environment())
  resGMSE_long$Time = system.time(resGMSE_long <- GMSE_long(X_design = X_design, Sim_data = Sim_data, env = environment()))
  resGMSE_short$Time = system.time(resGMSE_short <- GMSE_short(X_design = X_design, Sim_data = Sim_data, env = environment()))
  
  
  compare_gmse = rbind(GMSEAnalytic_long = resGMSE_long$GMSE, GMSEAnalytic_short = resGMSE_short$GMSE, resMC[c(1,2,5:8),], resBoot[c(1,2,5:8),])
  compare_cv = rbind(CVAnalytic_long = resGMSE_long$CV, CVAnalytic_short = resGMSE_short$CV, resMC[c(9:14),], resBoot[c(9:14),])
}

