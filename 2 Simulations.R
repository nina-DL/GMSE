# Latest version: 16-05-26
# This scripts computes GMSE with the Analytic approach and MC approach
require(nnet)
library(fastDummies)
library(MASS)

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
N_registry = 50000 # population size (register)
N_registry*sample_prop # expected sample size

sim_Accuracy_res = function(N_registry, sample_prop, props, fit = mod_fit_ISTAT, beta_true, ref_k = 8,
                            sim_seed = 2025, seed_cov = 2025, M_design = 100, M_model = 100, B = 100,
                            method = c("GMSE_long", "GMSE_short", "MC_Global", "MC_Componentwise", 
                                       "NPBoot", "PBoot", "Design", "Model")){
  
  # Set hyperparameters and simulate data mimicking ISTAT according to props
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
  
  # estimate the p.hat
  m.hat = multinom(y.true ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, data = Sim_data)
  p.hat = predict(m.hat, newdata = Sim_data, "probs")
  p.hat = p.hat[,titstu]
  theta.hat = apply(p.hat, 2, sum) # Estimator as defined in Eq. (3)
  
  if("GMSE_long" %in% method){
    # Analytic GMSE: v1 (long)
    GMSE_Analytic_Res_long = GMSE_long(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment())
    GMSE_Analytic_long = GMSE_Analytic_Res_long$GMSE
    CV_Analytic_long = GMSE_Analytic_Res_long$CV
    print("GMSE_long completed")
  }else{
    GMSE_Analytic_long = CV_Analytic_long = NA
  }
  
  if("GMSE_short" %in% method){
    # Analytic GMSE: v2 (short)
    GMSE_Analytic_Res = GMSE_short(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment())
    GMSE_Analytic = GMSE_Analytic_Res$GMSE
    CV_Analytic = GMSE_Analytic_Res$CV
    print("GMSE_short completed")
  }else{
    GMSE_Analytic = CV_Analytic = NA
  }
  
  if("MC_Componentwise" %in% method){
    resMC = MC_Componentwise(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all", 
                             M_design = M_design/100, M_model = M_model/10, 
                             env = environment(), true_model = fit)
    
    MCdec_Bias = resMC$MC_Bias
    MCdec_GMSE = resMC$MC_GMSE
    MCdec_CV = resMC$MC_CV
    print("MC_Componentwise completed")
  }else{
    resMC = MCdec_Bias = MCdec_GMSE = MCdec_CV = NA
  }
  
  if("MC_Global" %in% method){
    resMC_global = MC_Global(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all", 
                             M_design = M_design/100, M_model = M_model/10, 
                             env = environment(), true_model = fit)
    
    MC_Bias = resMC_global$MC_Bias
    MC_GMSE = resMC_global$MC_GMSE
    MC_CV = resMC_global$MC_CV
    print("MC_global completed")
  }else{
    resMC_global = MC_Bias = MC_GMSE = MC_CV = NA
  }
  
  
  if("NPBoot" %in% method){
    resBoot = NPBoot_Accuracy(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, gamma = "all", 
                              B = B, env = environment(), true_model = fit)
    
    Boot_Bias = resBoot$Boot_Bias
    Boot_GMSE = resBoot$Boot_GMSE
    Boot_CV = resBoot$Boot_CV
    print("NPBoot completed")
  }else{
    resBoot = Boot_Bias = Boot_GMSE = Boot_CV = NA
  }
  
  if("Design" %in% method){
    resDesign = Design_based(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k,
                             gamma = "all", M_design = M_design, 
                             env = environment(), true_model = fit)
    
    Design_Bias = resDesign$Design_Bias
    Design_GMSE = resDesign$Design_GMSE
    Design_CV = resDesign$Design_CV
    print("Design-based completed")
  }else{
    resDesign = Design_Bias = Design_GMSE = Design_CV = NA
  }
  
  if("Model" %in% method){
    resModel = Model_based(Sim_data = Sim_data, X_design = X_design, ref_k = ref_k, 
                           gamma = "all", M_model = M_model, 
                           env = environment(), true_model = fit)
    
    Model_Bias = resModel$Model_Bias
    Model_GMSE = resModel$Model_GMSE
    Model_CV = resModel$Model_CV
    print("Model-based completed")
  }else{
    resModel = Model_Bias = Model_GMSE = Model_CV = NA
  }
  
  return(list(theta.tilde = theta.tilde, theta.hat = theta.hat,
              GMSE_Analytic = GMSE_Analytic, CV_Analytic = CV_Analytic, 
              GMSE_Analytic_long = GMSE_Analytic_long, CV_Analytic_long = CV_Analytic_long,
              resMC = resMC, resMC_global = resMC_global, resBoot = resBoot, 
              resDesign = resDesign, resModel = resModel,
              MCdec_Bias = MCdec_Bias, MCdec_GMSE = MCdec_GMSE, MCdec_CV = MCdec_CV,
              MC_Bias = MC_Bias, MC_GMSE = MC_GMSE, MC_CV = MC_CV,
              Boot_Bias = Boot_Bias, Boot_GMSE = Boot_GMSE, Boot_CV = Boot_CV, 
              Design_Bias = Design_Bias, Design_GMSE = Design_GMSE, Design_CV = Design_CV,
              Model_Bias = Model_Bias, Model_GMSE = Model_GMSE, Model_CV = Model_CV))
}

# give it a try
N_registry = 100000
prova1 = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                          fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 1000, M_model = 1000, B = 1000,
                          method = c("GMSE_short", "MC_Global", "MC_Componentwise",
                                     "NPBoot", "Design", "Model"))


GMSE_N1 = rbind(LinGMSE = prova1$GMSE_Analytic, MCGMSE = prova1$resMC_global$MC_GMSE[1,], 
                BootGMSE = prova1$resBoot$Boot_GMSE[1,], DesignGMSE = prova1$resDesign$Design_GMSE[1,], ModelGMSE = prova1$resModel$Model_GMSE[1,])

CV_N1 = rbind(LinCV = prova1$CV_Analytic, MCCV = prova1$resMC_global$MC_CV[1,], 
              BootCV = prova1$resBoot$Boot_CV[1,], DesignCV = prova1$resDesign$Design_CV[1,], ModelCV = prova1$resModel$Model_CV[1,])
t(CV_N1)

# Answer to Referee 3 on the GMSE decomposition and the minor-order components
t(rbind(prova1$GMSE_Analytic, prova1$resMC$MC_GMSE))

save.image(file = "Res_N1.RData")

# give it a try
N_registry = 300000
prova2 = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                          fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 1000, M_model = 1000, B = 1000,
                          method = c("GMSE_short", "MC_Global", "MC_Componentwise",
                                     "NPBoot", "Design", "Model"))

GMSE_N2 = rbind(LinGMSE = prova2$GMSE_Analytic, MCGMSE = prova2$resMC_global$MC_GMSE[1,], 
                BootGMSE = prova2$resBoot$Boot_GMSE[1,], 
                DesignGMSE = prova2$resDesign$Design_GMSE[1,], ModelGMSE = prova2$resModel$Model_GMSE[1,])

CV_N2 = rbind(LinCV = prova2$CV_Analytic, MCCV = prova2$resMC_global$MC_CV[1,], 
              BootCV = prova2$resBoot$Boot_CV[1,], 
              DesignCV = prova2$resDesign$Design_CV[1,], ModelCV = prova2$resModel$Model_CV[1,])
t(CV_N2)

# Answer to Referee 3 on the GMSE decomposition and the minor-order components
t(rbind(prova2$GMSE_Analytic, prova2$resMC$MC_GMSE))

save.image(file = "Res_N2.RData")

# give it a try
N_registry = 500000
prova3= sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                         fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 1000, M_model = 1000, B = 1000,
                         method = c("GMSE_short", "MC_Global", "MC_Componentwise",
                                    "NPBoot", "Design", "Model"))

GMSE_N3 = rbind(LinGMSE = prova3$GMSE_Analytic, MCGMSE = prova3$resMC_global$MC_GMSE[1,], 
                BootGMSE = prova3$resBoot$Boot_GMSE[1,], 
                DesignGMSE = prova3$resDesign$Design_GMSE[1,], ModelGMSE = prova3$resModel$Model_GMSE[1,])

CV_N3 = rbind(LinCV = prova3$CV_Analytic, MCCV = prova3$resMC_global$MC_CV[1,], 
              BootCV = prova3$resBoot$Boot_CV[1,], 
              DesignCV = prova3$resDesign$Design_CV[1,], ModelCV = prova3$resModel$Model_CV[1,])
t(CV_N3)

# Answer to Referee 3 on the GMSE decomposition and the minor-order components
t(rbind(prova3$GMSE_Analytic, prova3$resMC$MC_GMSE))

save.image(file = "Res_N3.RData")

# Variability assessment of the GMSE estimators ---------------------------

N_registry = 100000
GMSE = CV = list()
for(t in 1:10){
  prova = sim_Accuracy_res(N_registry = N_registry, sample_prop = sample_prop, props = props, ref_k = 8,
                           fit = mod_fit_ISTAT, beta_true = beta_true, sim_seed = 2023, M_design = 1000, M_model = 1000, B = 1000,
                           method = c("GMSE_short", "NPBoot", "MC_Global"))
  GMSE[[t]] = rbind(LinGMSE = prova$GMSE_Analytic, MCGMSE = prova$resMC_global$MC_GMSE[1,], BootGMSE = prova$resBoot$Boot_GMSE[1,])
  CV[[t]] = rbind(LinCV = prova$CV_Analytic, MCCV = prova$resMC_global$MC_CV[1,], BootCV = prova$resBoot$Boot_CV[1,])
}
