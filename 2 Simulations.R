# Latest version: 05-11-24
# This scripts computes GMSE with the Analytic approach and MC approach
require(nnet)
library(fastDummies)
library(MASS)
library(svMisc)

#require(ggplot2)
#require(reshape2)
#require(MNLpred)
library(dplyr)
library(matlib)
require(foreign)
library(svMisc) #progress


# Load necessary data and functions --------------------------------

source("1 Main_functions.R")

# Full MC procedure for getting simulation performance: GMSE vs MC vs Boot  ---------------------------------------------------

# Set hyperparameters
N_registry = 500000 # population size (register)
N_registry*sample_prop # expected sample size

sim_Accuracy_res = function(N_registry, sample_prop, props, fit, beta_true, ref_k = 8,
                            sim_seed = 2024, seed_cov = 123, M_design = 10, M_model = 10, B = 100,
                            method = c("GMSE_long", "GMSE_short", "MC", "Boot", "Boot_v2", "PBoot")){
  
  # Set hyperparameters and simulate data mimicing ISTAT according to props
  J = length(beta_true[1,])
  K = length(beta_true[,1]) + 1
  H = J*K
  
  # Step 0: simulate data (principally covariates: fixed) according to real ISTAT data patterns
  Sim_res = simulate_ISTAT_edu(N_registry = N_registry, sample_prop = sample_prop, props = props, fit = fit, seed_cov = seed_cov, seed = sim_seed)
  Sim_data = Sim_res$Sim_data
  X_design = get_design(Sim_data = Sim_data, intercept = T)
  n_sample = sum(Sim_data$lambda)
  
  # true coef: p.true
  p.true = predict(fit, newdata = Sim_data, "probs")
  p.true = p.true[,titstu] # check if titstu is available
  theta.tilde = apply(p.true, 2, sum)
  true_Vm_theta = apply(p.true, 2, function(x) sum(x*(1-x)))
  
  if("GMSE_long" %in% method){
    # Analytic GMSE: v1 (long)
    GMSE_Analytic_Res_long = GMSE_long(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment())
    GMSE_Analytic_long = GMSE_Analytic_Res_long$GMSE
    CV_Analytic_long = GMSE_Analytic_Res_long$CV
    Ysum_Analytic_long = GMSE_Analytic_Res_long$theta.hat
  }else{
    GMSE_Analytic_long = CV_Analytic_long = NA
  }
  
  if("GMSE_short" %in% method){
    # Analytic GMSE: v2 (short)
    GMSE_Analytic_Res = GMSE_short(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment())
    GMSE_Analytic = GMSE_Analytic_Res$GMSE
    CV_Analytic = GMSE_Analytic_Res$CV
    Ysum_Analytic = GMSE_Analytic_Res$theta.hat
  }else{
    GMSE_Analytic = CV_Analytic = NA
  }
  
  if("MC" %in% method){
  resMC = MC_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all", 
                      M_design = M_design, M_model = M_model, env = environment())

  GMSE_MC = resMC[5:7,]
  CV_MC = resMC[12:14,]
  }else{
    resMC = GMSE_MC = CV_MC = NA
  }
  
  if("Boot" %in% method){
  resBoot = Boot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all", 
                          B = B, env = environment())
  GMSE_Boot = resBoot[5:7,]
  CV_Boot = resBoot[12:14,]
  }else{
    resBoot = GMSE_Boot = CV_Boot = NA
  }
  
  if("Boot_v2" %in% method){
    resBoot2 = Boot_Accuracy_v2(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all",
                               B = B, env = environment())
    GMSE_Boot2 = resBoot2[5:7,]
    CV_Boot2 = resBoot2[12:14,]
  }else{
    resBoot2 = GMSE_Boot2 = CV_Boot2 = NA
  }
  
  if("PBoot" %in% method){
    resPBoot = PBoot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all",
                             M_design = M_design, M_model = M_model, env = environment())
    
    GMSE_PBoot = resPBoot[5:7,]
    CV_PBoot = resPBoot[12:14,]
  }else{
    resPBoot = GMSE_PBoot = CV_PBoot = NA
  }
  
  return(list(GMSE_Analytic = GMSE_Analytic, CV_Analytic = CV_Analytic, 
              GMSE_Analytic_long = GMSE_Analytic_long, CV_Analytic_long = CV_Analytic_long,
              resMC = resMC, resBoot = resBoot, resBoot2 = resBoot2, resPBoot = resPBoot,
              GMSE_MC = GMSE_MC, CV_MC = CV_MC, 
              GMSE_Boot = GMSE_Boot, CV_Boot = CV_Boot, 
              GMSE_Boot2 = GMSE_Boot2, CV_Boot2 = CV_Boot2, 
              GMSE_PBoot = GMSE_PBoot, CV_PBoot = CV_PBoot))
}


N_registry = 100000
prova = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                         fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 10, M_model = 100, B = 1000,
                         method = c("GMSE_short", "MC", "Boot", "PBoot"))

GMSE = rbind(GMSE_Analytic = prova$GMSE_Analytic, prova$GMSE_Boot, prova$GMSE_MC)
CV = rbind(CV_Analytic = prova$CV_Analytic, prova$CV_Boot, prova$CV_MC)
noquote(paste0(round(t(CV*100),2)[1,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[2,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[3,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[4,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[5,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[6,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[7,c(1,3,6)], "%&"))
noquote(paste0(round(t(CV*100),2)[8,c(1,3,6)], "%&"))

save(GMSE, CV, file = "Res_all_100k.RData")

N_registry = 300000
prova = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                         fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 10, M_model = 100, B = 1000,
                         method = c("GMSE_short", "MC", "Boot", "PBoot"))

GMSE = rbind(GMSE_Analytic = prova$GMSE_Analytic, prova$GMSE_Boot, prova$GMSE_MC)
CV = rbind(CV_Analytic = prova$CV_Analytic, prova$CV_Boot, prova$CV_MC)

save(GMSE, CV, file = "Res_all_300k.RData")

N_registry = 500000
prova = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                         fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 10, M_model = 100, B = 1000,
                         method = c("GMSE_short", "MC", "Boot", "PBoot"))

GMSE = rbind(GMSE_Analytic = prova$GMSE_Analytic, prova$GMSE_Boot, prova$GMSE_MC)
CV = rbind(CV_Analytic = prova$CV_Analytic, prova$CV_Boot, prova$CV_MC)

save(GMSE, CV, file = "Res_all_500k.RData")


# Variability assessment of the GMSE estimators ---------------------------

N_registry = 500000
GMSE = CV = list()
for(t in 1:10){
prova = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                         fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = t, M_design = 10, M_model = 100, B = 1000,
                         method = c("GMSE_short", "MC", "Boot"))
GMSE[[t]] = rbind(GMSE_Analytic = prova$GMSE_Analytic, prova$GMSE_MC, prova$GMSE_Boot)
CV[[t]] = rbind(CV_Analytic = prova$CV_Analytic, prova$CV_MC, prova$CV_Boot)
}

save(GMSE, CV, file = "Res_all_10var_500k.RData")

CV_arr = array(unlist(CV), dim = c(13, 8, 100))
CV_avg = apply(CV_arr, c(1,2), mean)
CV_sd = apply(CV_arr, c(1,2), sd)
alpha = .1
CV_CI1 = apply(CV_arr, c(1,2), quantile, probs = alpha/2)
CV_CI2 = apply(CV_arr, c(1,2), quantile, probs = (1-alpha/2))

idx = c(1,3,6,12)
CV_avg[idx,]
CV_sd[idx,]

aa = NULL
for(i in 1:100){
  aa = rbind(aa, CV[[i]][1,])
}
hist(aa[,1])

GMSE_arr = array(unlist(GMSE), dim = c(5, 8, 50))
apply(GMSE_arr, c(1,2), mean)
