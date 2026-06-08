# Latest version: 14-05-26
# This script contains necessary functions to generate simulations, compute gmse, etc. 

load("Data_gen_process.RData")

require(nnet)
library(fastDummies)
library(MASS)
library(svMisc) #progress

# Main function for data simulation: Sim_data --------------------------------

simulate_ISTAT_edu = function(N_registry, sample_prop, props, fit = mod_fit_ISTAT, seed_cov = 2025, seed = 2025){
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
  # 10-07-25: Cambiato lo schema di campionamento: SRS-WOR 
  idx = sample(N_registry, sample_prop*N_registry, replace = F)
  membership = rep(0, N_registry)
  membership[idx] = 1
  Sim_data$lambda = membership
  #Sim_data$lambda = sample(0:1, N_registry, c(1-sample_prop, sample_prop), rep = T) # old version
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
  Sim_res = simulate_ISTAT_edu(N_registry = N_registry, sample_prop = sample_prop, props = props, fit = mod_fit_ISTAT, seed = NULL)
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

# |-- MC: OK ---------------------------------------------------

MC_Global = function(Sim_data, X_design, ref_k = 8, gamma = "all",
                     M_design = 10, M_model = 10, env,
                     true_model = mod_fit_ISTAT){
  # true coef: p.true
  p.true = predict(true_model, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true, 2, sum)
  
  theta = table(Sim_data$y.true)[titstu]
  
  # initiate simulations
  M = M_design*M_model
  theta.hat = matrix(NA, nrow = M, ncol = env$K)
  
  for(i in 1:M){
    set.seed(i)
    
    # Step 1: variability due to the sampling/design process
    # Note: for now we assume same prob for all individuals
    idx_sample = sample(env$N_registry, env$n_sample) 
    
    # We can also switch to the inclusion probs indicated by the db
    # idx_sample = sample(1:env$N_registry, env$n_sample, prob = Sim_data$prob_incl_v2) 
    
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
      
    p.hat = predict(mod_est_MCsample, newdata = Sim_data, "probs")
    p.hat = p.hat[,titstu] # check if titstu is available
    theta.hat[i,] = apply(p.hat, 2, sum) # Estimator as defined in Eq. (2)
    
    cat(sprintf("\r============================ MC iteration: %d/%d", i, M))
  }
  
  # Global Variance
  GVar = apply(theta.hat, 2, var)
  GVar2 = apply(theta.hat, 2, function(x) mean((x - mean(x))^2)) # this is unadjusted one (denominator B, not B-1)
  
  # Global bias
  GBias_true = apply(theta.hat, 2, mean) - theta.tilde
  
  # MSE 
  GMSE_true = GVar2 + GBias_true^2
  
  CV_GVar = sqrt(GVar)/theta.tilde
  CV_GVar2 = sqrt(GVar2)/theta.tilde
  CV_GMSE_true = sqrt(GMSE_true)/theta.tilde
  
  # overall look
  GMSE = rbind(GVar, GVar2, GMSE_true)
  CV = rbind(CV_GVar, CV_GVar2, CV_GMSE_true)
  
  return(list(MC_Bias = GBias_true, MC_GMSE = GMSE, MC_CV = CV))
}

# |-- MC Decomposition: OK ---------------------------------------------------

MC_Componentwise = function(Sim_data, X_design, ref_k = 8, gamma = "all", 
                       M_design = 5, M_model = 10, env,
                       true_model = mod_fit_ISTAT){
  # true coef: p.true
  p.true = predict(true_model, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true, 2, sum)
  true_Vm_theta = apply(p.true, 2, function(x) sum(x*(1-x)))

  # initiate simulations
  Em_theta = Vm_theta = Em_theta.hat = Vm_theta.hat = matrix(NA, nrow = M_design, ncol = env$K)
  
  for(i in 1:M_design){
    set.seed(i)
    
    # Step 1: variability due to the sampling/design process
    # Note: for now we assume same prob for all individuals
    idx_sample = sample(env$N_registry, env$n_sample) 

    # We can also switch to the inclusion probs indicated by the db
    # idx_sample = sample(1:env$N_registry, env$n_sample, prob = Sim_data$prob_incl_v2) 
    
    theta.hat = theta_mc = matrix(NA, nrow = M_model, ncol = env$K)
    for(m in 1:M_model){
      
      # Step 2. Variability due to model: regenerate y based on its assumed model (true beta: MC estimate)
      y.true_mc = apply(env$p.true, 1, function(x) sample(1:(env$K), 1, prob=x))
      # ensure that we have each category observed at least once in the sample, otherwise it returns an error when fitting the model
      while(length(table(y.true_mc[idx_sample])) < 8){
        y.true_mc = apply(env$p.true, 1, function(x) sample(1:8, 1, prob=x))
      }
      
      Sim_data$y.true_mc = factor(y.true_mc, levels = 1:(env$K), labels = levels(Sim_data$y.true))
      # again, we take as Y of reference the same as in the main model fit 
      Sim_data$y.true_mc_ref <- relevel(Sim_data$y.true_mc, ref = ref_k)
      
      Sim_sample = Sim_data[idx_sample, ]
      
      # estimate the p.hat
      mod_est_MCsample <- multinom(y.true_mc_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
      
      # Option I. Y.hat = p.hat (deterministic) + We re-estimate Y = p for all the registry (included sample)
      p.hat = predict(mod_est_MCsample, newdata = Sim_data, "probs")
      p.hat = p.hat[,titstu] # check if titstu is available
      theta.hat[m,] = apply(p.hat, 2, sum) # Estimator as defined in Eq. (2)
      
      # This helps to evaluate the variability of the population parameter Y_sum (minor order term in GMSE 13a)
      theta_mc[m,] = table(Sim_data$y.true_mc_ref)[titstu]
      
    }
    
    Em_theta[i,] = apply(theta_mc, 2, mean) # expectation of theta = sumY under M
    Vm_theta[i,] = apply(theta_mc, 2, var) # variance of theta = sumY under M
    
    Em_theta.hat[i,] = apply(theta.hat, 2, mean) # expectation under M of theta.hat = sumYhat
    Vm_theta.hat[i,] = apply(theta.hat, 2, var) # variance under M of theta.hat = sumYhat

    cat(sprintf("\r============================ MC iteration (design): %d/%d", i, M_design))
  }
  
  # Model variance (negligible)
  VarM_est = apply(Vm_theta, 2, mean)
  
  # Upper bound computed in analytic GMSE: EdVm
  GVar = apply(Vm_theta.hat, 2, mean)
  
  # Global bias
  GBias_true = apply(Em_theta.hat, 2, mean) - theta.tilde

  # Computation of GMSE - using PartI-II-III (internal domain)
  Approx_term_Th1 = apply(Em_theta.hat, 2, var)
  GMSE_3parts = GVar + Approx_term_Th1 - VarM_est
  GMSE_3parts_true = GVar + Approx_term_Th1 - true_Vm_theta
  
  CV_GVar = sqrt(GVar)/theta.tilde
  CV_3parts = sqrt(GMSE_3parts)/theta.tilde
  CV_3parts_true = sqrt(GMSE_3parts_true)/theta.tilde
    
  # overall look
  GMSE = rbind(VarM_est, true_Vm_theta, Approx_term_Th1, GVar, GMSE_3parts, GMSE_3parts_true)
  CV = rbind(CV_GVar, CV_3parts_true)
  
  return(list(MC_Bias = GBias_true, MC_GMSE = GMSE, MC_CV = CV))
}

# |-- Design-based estimator: OK ---------------------------------------------------

Design_based = function(Sim_data, X_design, ref_k = 8, gamma = "all", 
                        M_design = 10, env,
                        true_model = mod_fit_ISTAT){
  
  # 01-06-2026: restrict the domain of interest, if required
  if(is.numeric(gamma) == F){
    gamma = rep(1, env$N_registry)
    #gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0)
  }
  Sim_data$gamma = gamma
  
  # true coef: p.true
  p.true = predict(true_model, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true[gamma==1,], 2, sum)
  
  # initiate simulations
  theta.hat.design = matrix(NA, nrow = M_design, ncol = env$K)
  
  for(i in 1:M_design){
    set.seed(i)
    
    # Step 1: variability due to the sampling/design process
    # Note: for now we assume same prob for all individuals
    idx_sample = sample(env$N_registry, env$n_sample) 
    
    while(sum(table(Sim_data$y.true[idx_sample])>0) < 8){
      idx_sample = sample(env$N_registry, env$n_sample) 
    }
    
    Sim_sample = Sim_data[idx_sample, ]
    
    # Design-based estimator (Rev 1)
    mod_fit = multinom(y.true ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
    
    p.hat.design = predict(mod_fit, newdata = Sim_data, "probs")
    p.hat.design = p.hat.design[,titstu] # check if titstu is available
    theta.hat.design[i,] = apply(p.hat.design[gamma==1,], 2, sum) # Estimator as defined in Eq. (2)
    
    cat(sprintf("\r============================ Design-based iteration: %d/%d", i, M_design))
    # We can also switch to the inclusion probs indicated by the db
    # idx_sample = sample(1:env$N_registry, env$n_sample, prob = Sim_data$prob_incl_v2) 
  }
  
  # Design-based estimates of the variance (Rev1)
  Design_Bias_true = apply(theta.hat.design, 2, mean) - theta.tilde
  #Design_Var = apply(theta.hat.design, 2, var)
  Design_Var = apply(theta.hat.design, 2, function(x) mean((x - mean(x))^2)) # this is unadjusted one (denominator B, not B-1)
  
  # Design_MSE_true = Design_Bias_true^2 + Design_Var#*(M_design-1)/M_design
  # the discrepancy with the below is due to correction factor of variance M-1/M
  Design_MSE_true = colMeans((theta.hat.design - matrix(theta.tilde, nrow(theta.hat.design), 8, byrow=TRUE))^2)
  
  CV_Var = sqrt(Design_Var)/apply(theta.hat.design, 2, mean)
  #CV_Var2 = sqrt(Design_Var2)/apply(theta.hat.design, 2, mean)
  CV_MSE_true = sqrt(Design_MSE_true)/theta.tilde
  #CV_MSE_true2 = sqrt(Design_MSE_true2)/theta.tilde
  
  # overall look
  GMSE = rbind(Design_Var, Design_MSE_true)
  CV = rbind(CV_Var, CV_MSE_true)
  
  return(list(Design_Bias = Design_Bias_true, Design_GMSE = GMSE, Design_CV = CV))
}

# |-- Model-based Estimator: OK ---------------------------------------------------

Model_based = function(Sim_data, X_design, ref_k = 8, gamma = "all", 
                       M_model = 100, env,
                       true_model = mod_fit_ISTAT){
  
  # 01-06-2026: restrict the domain of interest, if required
  if(is.numeric(gamma) == F){
    gamma = rep(1, env$N_registry)
    #gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0)
  }
  Sim_data$gamma = gamma
  
  # true coef: p.true
  p.true = predict(true_model, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true[gamma==1,], 2, sum)
  
  # # initiate simulations
  # Em_theta = Vm_theta = Em_theta.hat = Vm_theta.hat = matrix(NA, nrow = M_design, ncol = env$K)
  idx_sample = which(Sim_data$lambda==1)
  Sim_sample = Sim_data[idx_sample, ]
  mod_fit = multinom(y.true ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
  
  p.hat.design = predict(mod_fit, newdata = Sim_data, "probs")
  p.hat.design = p.hat.design[,titstu] # check if titstu is available
  
  theta.model = theta.hat.model = matrix(NA, nrow = M_model, ncol = env$K)
  for(m in 1:M_model){
    
    # Model-based estimator for the variance (Rev 1) : revised 14-05-2026
    
    # Step 2. Variability due to model: regenerate y based on its assumed model (est beta based on the sample)
    y_model = apply(p.hat.design, 1, function(x) sample(1:(env$K), 1, prob=x))
    # ensure that we have each category observed at least once in the sample, otherwise it returns an error when fitting the model
    while(length(table(y_model[idx_sample])) < 8){
      y_model = apply(p.hat.design, 1, function(x) sample(1:8, 1, prob=x))
    }
    
    Sim_data$y_model = factor(y_model, levels = 1:(env$K), labels = levels(Sim_data$y.true))
    # again, we take as Y of reference the same as in the main model fit 
    Sim_data$y_model_ref <- relevel(Sim_data$y_model, ref = ref_k)
    
    Sim_sample = Sim_data[idx_sample, ]
    
    # estimate the p.hat
    mod_est <- multinom(y_model_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
    
    # Option I. Y.hat = p.hat (deterministic) + We re-estimate Y = p for all the registry (included sample)
    p.hat.model = predict(mod_est, newdata = Sim_data, "probs")
    p.hat.model = p.hat.model[,titstu] # check if titstu is available
    theta.hat.model[m,] = apply(p.hat.model[gamma==1,], 2, sum) # Estimator as defined in Eq. (2)
    
    cat(sprintf("\r============================ Model-based iteration: %d/%d", m, M_model))
  }
  
  # Model-based estimates of the variance (Rev1)
  Model_Bias_true = apply(theta.hat.model, 2, mean) - theta.tilde
  #Model_Var = apply(theta.hat.model, 2, var)
  Model_Var = apply(theta.hat.model, 2, function(x) mean((x - mean(x))^2)) # this is unadjusted one (denominator B, not B-1)
  
  # Model_MSE_true = Model_Bias_true^2 + Model_Var#*(M_design-1)/M_design
  # the discrepancy with the below is due to correction factor of variance M-1/M
  Model_MSE_true = colMeans((theta.hat.model - matrix(theta.tilde, nrow(theta.hat.model), 8, byrow=TRUE))^2)
  
  CV_Var = sqrt(Model_Var)/apply(theta.hat.model, 2, mean)
  #CV_Var2 = sqrt(Model_Var2)/apply(theta.hat.model, 2, mean)
  CV_MSE_true = sqrt(Model_MSE_true)/theta.tilde
  #CV_MSE_true2 = sqrt(Model_MSE_true2)/theta.tilde
  
  # overall look
  GMSE = rbind(Model_Var, Model_MSE_true)
  CV = rbind(CV_Var, CV_MSE_true)
  
  return(list(Model_Bias = Model_Bias_true, Model_GMSE = GMSE, Model_CV = CV))
}

# |-- NP Bootstrap: OK ---------------------------------------------------

# update: 24/09/2025 (Chambers pg 136)
# update: 16/05/2026 (MSE computation)
# update: 26/05/2026 (Only one version of GVar taken)

NPBoot_Accuracy = function(Sim_data, X_design, ref_k = 8, gamma = "all",
                           B = 100, env, true_model = mod_fit_ISTAT){
  
  # 01-06-2026: restrict the domain of interest, if required
  if(is.numeric(gamma) == F){
    gamma = rep(1, env$N_registry)
    #gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0)
  }
  Sim_data$gamma = gamma
  
  # true coef: p.true
  p.true = predict(true_model, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true[gamma==1,], 2, sum)
  Sim_data = cbind(Sim_data, p.true)
  
  Sim_data$y.true_ref <- relevel(Sim_data$y.true, ref = ref_k)
  Sim_sample = Sim_data[Sim_data$lambda==1, ]
  mod_est = multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_sample)
  p.est = predict(mod_est, newdata = Sim_data, "probs")
  p.est = p.est[,titstu]
  theta.est = apply(p.est[gamma==1,], 2, sum)
  
  theta.boot = matrix(NA, nrow = B, ncol = env$K)
  for(b in 1:B){
    set.seed(b)
    
    # Step 1: Bootstrap the n sample data from the register
    idx_sample = which(Sim_data$lambda==1)
    B_idx_sample = sample(idx_sample, env$n_sample, replace = T) 
    
    # Step 2: Bootstrap the N-n nonsample data from the register
    idx_nsample = which(Sim_data$lambda==0)
    B_idx_nsample = sample(idx_nsample, env$N_registry - env$n_sample, replace = T) 
    
    # Now use the b-th run to get estimates and evaluate errors
    Boot_sample = Sim_data[B_idx_sample, ]
    
    # ensure that we have each category observed at least once in the sample
    while(any(table(Boot_sample$y.true)==0)){
      B_idx_sample = sample(idx_sample, env$n_sample, replace = T) 
      Boot_sample = Sim_data[B_idx_sample, ]
    }
    
    Boot_data = rbind(Sim_data[B_idx_nsample, ], Boot_sample)
    
    # estimate the p.boot
    mod_est_MCsample <- multinom(y.true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Boot_sample)
    
    # Option I. We can re-estimate the parameters for all the registry (included sample)
    p.boot = predict(mod_est_MCsample, newdata = Boot_data, "probs")
    p.boot = p.boot[,titstu]
    theta.boot[b,] = apply(p.boot[Boot_data$gamma==1,], 2, sum)
    
    cat(sprintf("\r============================ Bootstrap iteration: %d/%d", b, B))
  }
  
  theta.boot.est = apply(theta.boot, 2, mean) # expectation under both M and Pi of theta.boot = sumYhat

  # "Global Variance" - Eq (8)
  #GVar = Vm_theta.boot = apply(theta.boot, 2, var) # variance under both M and Pi of theta.boot = sumYhat
  GVar = apply(theta.boot, 2, function(x) mean((x - mean(x))^2)) # this is unadjusted one (denominator B, not B-1)
  
  # "Global bias"
  GBias_est = theta.boot.est - theta.est[titstu]
  GBias_true = theta.boot.est - theta.tilde[titstu]

  # Computation of GMSE as per definition - Eq (7)
  GMSE_est = GVar + GBias_est^2
  GMSE_true = colMeans((theta.boot - matrix(theta.tilde, nrow = nrow(theta.boot), ncol = length(theta.tilde), byrow = TRUE))^2)
  
  #CV_GVar = sqrt(GVar)/theta.boot.est
  CV_GVar = sqrt(GVar)/theta.boot.est
  CV_est = sqrt(GMSE_est)/theta.boot.est
  CV_true = sqrt(GMSE_true)/theta.tilde

  # overall look
  Bias = rbind(GBias_est, GBias_true)
  GMSE = rbind(GVar, GMSE_est, GMSE_true)
  CV = rbind(CV_GVar, CV_est, CV_true)
  
  return(list(Boot_Bias = Bias, Boot_GMSE = GMSE, Boot_CV = CV))
}

# Give it a try to All -----------------------------------------------------------


# Note: check reproducibility
wannatry = F

if(wannatry==T){
  # Set hyperparameters and give it a try
  N_registry = 50000
  Sim_res = simulate_ISTAT_edu(N_registry = N_registry, sample_prop = sample_prop, props = props, fit = mod_fit_ISTAT, seed = NULL)
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
  resBoot = NPBoot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                                                 B = 100, env = environment())
  resBoot_V2 = PBoot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                                B = 100, env = environment())
  resGMSE_long$Time = system.time(resGMSE_long <- GMSE_long(X_design = X_design, Sim_data = Sim_data, env = environment()))
  resGMSE_short$Time = system.time(resGMSE_short <- GMSE_short(X_design = X_design, Sim_data = Sim_data, env = environment()))
  
  resDesign = Design_based(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                           M_design = 100, env = environment())
  resModel = Model_based(Sim_data = Sim_data, X_design = X_design, ref_k = 8, gamma = "all", 
                          M_model = 100, env = environment())
  compare_gmse = rbind(GMSEAnalytic_long = resGMSE_long$GMSE, GMSEAnalytic_short = resGMSE_short$GMSE, resMC[c(1,2,5:8),], 
                       resBoot[c(1,2,5:8),], resBoot_V2[c(1,2,5:8),],
                       resDesign = resDesign, resModel = resModel)
  compare_cv = rbind(CVAnalytic_long = resGMSE_long$CV, CVAnalytic_short = resGMSE_short$CV, resMC[c(9:14),], resBoot[c(9:14),], resBoot_V2[c(9:14),])
}

