# Latest version: 01-06-26
# This script implements the GMSE in the real ISTAT dataset (variables, their distribution and association patterns with Y)

require(nnet)
library(fastDummies)
library(MASS)
library(svMisc) #progress

require(ggplot2)
library(dplyr)
library(matlib)
require(foreign)


# Load and inspect the pseudo-real dataset --------------------------------

# Use this code (uncomment) if your DB is in a local folder
# dati_ISTAT =  read.delim("/Users/blackmamba/Library/CloudStorage/GoogleDrive-nina.deliu@uniroma1.it/Il mio Drive/[ISTAT–SAPIENZA] Accuratezza Stime da Registro/ISTAT Nina-Piero/DatiSintetici_v2.txt")

# Use this code to import the DB directly from Gdrive
library(googledrive)

link <- "https://drive.google.com/file/d/13x3AQYuLZFT0W6mOfejon1Zz2yfXbz1Y/view"
drive_download(as_id(link), path = "DatiSintetici_v2.txt", overwrite = TRUE)
dati_ISTAT <- read.delim("DatiSintetici_v2.txt")

# Create age variable (in class format) + Unify some age classes
dati_ISTAT$cleta_19_new = ((dati_ISTAT$cleta_19==5)|(dati_ISTAT$cleta_19==6)|(dati_ISTAT$cleta_19==7)|(dati_ISTAT$cleta_19==8)|(dati_ISTAT$cleta_19==9)|(dati_ISTAT$cleta_19==10))*10 + (dati_ISTAT$cleta_19>10)*dati_ISTAT$cleta_19

age = c("0-28", "29-39", "40-49", "50-69", "70+")

titstu = c("1 Illiterate",
           "2 Literate but no education",
           "3 Primary", 
           "4 Lower secondary",
           "5 Upper secondary",
           "6 Bachelor degree or equivalent",
           "7 Master degree or equivalent",
           "8 PhD level")


dati_ISTAT$TITSTU8_MS <- factor(dati_ISTAT$TITSTU8_MS,
                                levels = 1:8, labels = titstu)

dati_ISTAT$TITSTU_IMP19 <- factor(dati_ISTAT$TITSTU_IMP19,
                                  levels = 1:8, labels = titstu)

dati_ISTAT$SESSO <- factor(dati_ISTAT$SESSO,
                           levels = 1:2, labels = c("M", "F"))

dati_ISTAT$FL_ITA <- factor(dati_ISTAT$FL_ITA,
                            levels = 0:1, labels = c("No", "Si"))

dati_ISTAT$cleta_19_new <- factor(dati_ISTAT$cleta_19_new,
                                  levels = 10:14, labels = age)

dati_ISTAT$TITSTU8_CENS11 <- factor(dati_ISTAT$TITSTU8_CENS11,
                                    levels = 1:8, labels = titstu)

prov_res = 1:9
dati_ISTAT$COD_PROV_RES <- factor(dati_ISTAT$COD_PROV_RES,
                                  levels = 1:9, labels = prov_res)

dati_ISTAT$y.true = dati_ISTAT$TITSTU_IMP19

# names(dati_ISTAT)
dati_ISTAT_red = dati_ISTAT[dati_ISTAT$FL_IMP19%in%c("B1.1", "BS.1"),c("SESSO", "FL_ITA", "COD_PROV_RES", "TITSTU8_CENS11", "cleta_19_new", "fl_MS", "prob_incl_v2", "y.true")]
N_registry = dim(dati_ISTAT_red)[1] # population size (register)
n_sample = sum(dati_ISTAT_red$fl_MS)

table(dati_ISTAT_red$SESSO)/N_registry
table(dati_ISTAT_red$FL_ITA)/N_registry
table(dati_ISTAT_red$TITSTU8_CENS11)/N_registry
table(dati_ISTAT_red$cleta_19_new)/N_registry
table(dati_ISTAT_red$COD_PROV_RES)/N_registry

# Estimate the coefficients in the ISTAT dataset and get population totals (estimates) --------------------------------

# we use as reference k = K = 8 (to reflect the paper)
ref_k = 8
dati_ISTAT_red$y_true_ref <- relevel(dati_ISTAT_red$y.true, ref = ref_k)
mod_fit_ISTAT <- multinom(y_true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, 
                          data = dati_ISTAT_red[dati_ISTAT_red$fl_MS==1, ])
beta_est = coef(mod_fit_ISTAT)
beta_est
# noquote(paste0(round(t(beta_est),3)[14,], "&"))

p_est = predict(mod_fit_ISTAT, newdata = dati_ISTAT_red, "probs")[, titstu]
apply(p_est, 2, mean)
round(apply(p_est, 2, sum), 0)

# Some data description --------------------------------

sample_data = dati_ISTAT_red[dati_ISTAT_red$fl_MS==1,]
sample_cat = table(sample_data$y_true, sample_data$SESSO)[titstu,]
myFrame <- as.data.frame(sample_cat)
myFrame$type = "Sample"

myFrame2 = myFrame
myFrame2$Freq <- round(c(apply(p_est[dati_ISTAT_red$SESSO=="M",], 2, sum), 
                  apply(p_est[dati_ISTAT_red$SESSO=="F",], 2, sum)), 0)
myFrame2$Freq[16] = 544 #issues with rounding
myFrame2$Var2 = c(rep("M", 8), rep("F", 8))
myFrame2$Var1 = myFrame$Var1
myFrame2$type = "Register"

full_data = rbind(myFrame, myFrame2)
names(full_data) = c("Education", "Gender", "Count", "type")

full_data$Gender = ifelse(full_data$Gender=="M", "Male", "Female")

#library(viridis)
library(hrbrthemes)
ggplot(full_data, aes(fill=Gender, y=Count, x=Education)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  #geom_text(aes(label = Count), colour = "black", position = position_stack(vjust = 0.3)) +
  #scale_fill_viridis(discrete = T, option = "E") +
  #ggtitle("Studying 4 species..") +
  facet_grid(~type, scales = "free") +
  theme_ipsum() +
  geom_label(aes(label = after_stat(y), group = Education), 
            stat = 'summary', fun = sum, hjust = -.5, position = position_stack(vjust = 0.7),
            alpha = .5, show.legend  = FALSE) +
  theme(legend.position="bottom") +
  xlab("")

round(table(sample_data$y_true)[titstu]/round(apply(p_est, 2, sum), 0), 3)*100

# Load necessary data and functions for Accuracy estimation--------------------------------

source("1 Main_functions.R")

# Getting accuracy estimates: GMSE vs Boot  ---------------------------------------------------

# Set hyperparameters
J = length(beta_est[1,])
K = length(beta_est[,1]) + 1
H = J*K

ISTAT_data = dati_ISTAT_red
X_design = get_design(ISTAT_data, intercept = T)

ISTAT_data$lambda = ISTAT_data$fl_MS

# Analytic GMSE: v2 (short)
resAnalytic = GMSE_short(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment())
Ysum_Analytic = round(resAnalytic$theta.hat, 0)
apply(sample_cat, 1, sum)
GMSE_Analytic = round(resAnalytic$GMSE, 0)
CV_Analytic = round(resAnalytic$CV, 4)*100

# Bootstrap GMSE
resBoot = NPBoot_Accuracy(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, 
                          gamma = "all", B = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_Boot = round(resBoot$Boot_GMSE[1,], 0)
CV_Boot = round(resBoot$Boot_CV[1,], 4)*100

# Design-based GMSE
resDesign = Design_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, 
                         gamma = "all", M_design = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_Design = round(resDesign$Design_GMSE[1,], 0)
CV_Design = round(resDesign$Design_CV[1,], 4)*100

# Model-based GMSE
resModel = Model_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, 
                       gamma = "all",  M_model = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_Model = round(resModel$Model_GMSE[1,], 0)
CV_Model = round(resModel$Model_CV[1,], 4)*100

GMSE = rbind(GMSE_Analytic = GMSE_Analytic, GMSE_Boot, GMSE_Design, GMSE_Model)
CV = rbind(CV_Analytic = CV_Analytic, CV_Boot, CV_Design, CV_Model)

t(GMSE)
t(CV)

# Accuracy estimates by internal domain M/F: GMSE vs Boot  ---------------------------------------------------

# Analytic GMSE: v2 (short)
GMSE_Analytic_ResM = GMSE_short(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, domain = "SESSO", env = environment())
XX = round(GMSE_Analytic_ResM$theta.hat, 0)
# XX = sample_cat
GMSE_AnalyticM = round(GMSE_Analytic_ResM$GMSE[1,], 4)
CV_AnalyticM = round(GMSE_Analytic_ResM$CV[1,], 4)*100
GMSE_AnalyticF = round(GMSE_Analytic_ResM$GMSE[2,], 4)
CV_AnalyticF = round(GMSE_Analytic_ResM$CV[2,], 4)*100
YY = log(round(GMSE_Analytic_ResM$CV, 4)*100)

df1 <- reshape2::melt(XX, c("Education", "Gender"), value.name = "Count")
df2 <- reshape2::melt(YY, c("Gender", "Education"), value.name = "CV")
df = merge(df1, df2)
#df$Gender = ifelse(df$Gender==1, "M", "F")

ggplot(df, aes(x = Count, y = CV, color = Education,
               shape = Gender, group=interaction(Education, Gender))) +
  geom_point(size = 2) +
  scale_shape_manual(values=rep(c(20, 8), 8)) +
  theme_classic() +
  #theme(legend.position="right") +
  theme(legend.position = "inside", legend.position.inside = c(.9, .65)) +
  ylab("CV (log scale)") +
  xlab("Counts (Sample)") +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Estimated CV by education class and gender domain with respect to their (sample) size")


# t(matrix(round(GMSE_Analytic_ResM$GMSE, 0)[2,]))
# t(matrix(round(GMSE_Analytic_ResM$CV[2,], 4)*100))

resBootM = NPBoot_Accuracy(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0), 
                           B = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_BootM = round(resBootM$Boot_GMSE[1,],0)
CV_BootM = round(resBootM$Boot_CV[1,], 4)*100

resBootF = NPBoot_Accuracy(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, gamma = ifelse(ISTAT_data$SESSO=="F", 1, 0), 
                           B = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_BootF = round(resBootF$Boot_GMSE[1,],0)
CV_BootF = round(resBootF$Boot_CV[1,], 4)*100

resDesignM = Design_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k,  gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0),
                          M_design = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_DesignM = round(resDesignM$Design_GMSE[1,], 0)
CV_DesignM = round(resDesignM$Design_CV[1,], 4)*100

resDesignF = Design_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k,  gamma = ifelse(ISTAT_data$SESSO=="F", 1, 0),
                          M_design = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_DesignF = round(resDesignF$Design_GMSE[1,], 0)
CV_DesignF = round(resDesignF$Design_CV[1,], 4)*100

resModelM = Model_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k,  gamma = ifelse(ISTAT_data$SESSO=="M", 1, 0),
                       M_model = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_ModelM = round(resModelM$Model_GMSE[1,], 0)
CV_ModelM = round(resModelM$Model_CV[1,], 4)*100

resModelF = Model_based(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k,  gamma = ifelse(ISTAT_data$SESSO=="F", 1, 0),
                        M_model = 1000, env = environment(), true_model = mod_fit_ISTAT)
GMSE_ModelF = round(resModelF$Model_GMSE[1,], 0)
CV_ModelF = round(resModelF$Model_CV[1,], 4)*100


GMSE_M = rbind(GMSE_Analytic = GMSE_AnalyticM, GMSE_BootM, GMSE_DesignM, GMSE_ModelM)
GMSE_F = rbind(GMSE_Analytic = GMSE_AnalyticF, GMSE_BootF, GMSE_DesignF, GMSE_ModelF)
CV_M = rbind(CV_Analytic = CV_AnalyticM, CV_BootM, CV_DesignM, CV_ModelM)
CV_F = rbind(CV_Analytic = CV_AnalyticF, CV_BootF, CV_DesignF, CV_ModelF)

t(GMSE_M)
t(GMSE_F)
t(CV_M)
t(CV_F)

# Accuracy estimates by external domain "Municipality" (Analytic only) ---------------------------------------------------


# Analytic GMSE: v2 (short)
GMSE_Analytic_ResP = GMSE_short(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, domain = "COD_PROV_RES", env = environment())
XX = round(GMSE_Analytic_ResP$theta.hat, 0)
table(ISTAT_data$COD_PROV_RES[ISTAT_data$lambda==1])
XX = table(ISTAT_data$COD_PROV_RES[ISTAT_data$lambda==1], ISTAT_data$y.true[ISTAT_data$lambda==1])
round(GMSE_Analytic_ResP$GMSE, 0)
YY = log(round(GMSE_Analytic_ResP$CV, 4))

df1 <- reshape2::melt(XX, c("Province", "Education"), value.name = "Count")
df2 <- reshape2::melt(YY, c("Province", "Education"), value.name = "CV")
df = merge(df1, df2)

ggplot(df, aes(x = Count, y = CV, group = Education, color = Education,
               shape = Education)) +
  geom_point(size = 2) +
  scale_shape_manual(values=c(3, 16:18, 8, 15, 4, 13)) +
  theme_classic() +
  #theme(legend.position="right") +
  theme(legend.position = "inside", legend.position.inside = c(.9, .65)) +
  ylab("CV (log scale)") +
  xlab("Counts (Sample)") +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Estimated CV by education class and municipality domain with respect to their (sample) size")


# Computational time ------------------------------------------------------

resGMSE_Analytic$Time = system.time(resGMSE_Analytic <- GMSE_short(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, domain = "all", env = environment()))
sim_seed = 2024; seed_cov = 123; p.true = p_est
resGMSE_Boot$Time = system.time(resGMSE_Boot <- Boot_Accuracy(Sim_data = ISTAT_data, X_design = X_design, ref_k = ref_k, gamma = "all", B = 1000, env = environment()))
