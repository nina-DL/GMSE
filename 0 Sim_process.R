# Latest version: 05-11-24
# This script estimates the model coefficients from the ISTAT dataset (variables, their distribution and association patterns with Y)
# and uses them as "true" coefficients for simulating the data. 

# Load and inspect the pseudo-real dataset --------------------------------

# Use this code (uncomment) if your DB is in a local folder
# dati_ISTAT =  read.delim("/Users/blackmamba/Library/CloudStorage/GoogleDrive-nina.deliu@uniroma1.it/Il mio Drive/[ISTAT–SAPIENZA] Accuratezza Stime da Registro/ISTAT Nina-Piero/DatiSintetici_v2.txt")

# Use this code to import the DB directly from Gdrive
library(googledrive)

# Authenticate and get the file ID
drive_auth()  # Follow the authentication process

# Replace 'YOUR_FILE_ID' with the actual file ID from the Google Drive link
file_id = "1VZHnYcv8ghyaaGMZcHS9a2HhlvPupHIg"
# Create a temporary file to download the dataset
temp_file = tempfile(fileext = ".xlsx")
# Download the file
drive_download(as_id(file_id), path = temp_file)
# Read the dataset into R
dati_ISTAT = read.delim(temp_file)
# Clean up: Delete the temporary file
unlink(temp_file)
rm(file_id, temp_file)


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


# Estimate the coefficients in the ISTAT dataset and get population totals (estimates) --------------------------------

require(nnet)
# we use as reference k = K = 8 (to reflect the paper)
ref_k = 8
dati_ISTAT_red$y_true_ref <- relevel(dati_ISTAT_red$y.true, ref = ref_k)
mod_fit_ISTAT <- multinom(y_true_ref ~ cleta_19_new + SESSO + FL_ITA + TITSTU8_CENS11, 
                          data = dati_ISTAT_red[dati_ISTAT_red$fl_MS==1, ])

# we use these as "true" coef: beta_true
beta_true = coef(mod_fit_ISTAT) 
# p_true = predict(mod_fit_ISTAT, newdata = dati_ISTAT_red, "probs")

# we will use the sample_prop value for defining the sample size in the simulated sample
sample_prop = round(table(dati_ISTAT_red$fl_MS)/length(dati_ISTAT_red$fl_MS), 2)[2]

# we will use the props values for defining the characteristics of the sampled covariates
props = list(age_prop = table(dati_ISTAT_red$cleta_19_new)/length(dati_ISTAT_red$cleta_19_new),
             SESSO_prop = table(dati_ISTAT_red$SESSO)/length(dati_ISTAT_red$SESSO),
             FL_prop = table(dati_ISTAT_red$FL_ITA)/length(dati_ISTAT_red$FL_ITA),
             TITSTU11_prop = table(dati_ISTAT_red$TITSTU8_CENS11)/length(dati_ISTAT_red$TITSTU8_CENS11),
             PROVres_prop = table(dati_ISTAT_red$COD_PROV_RES)/length(dati_ISTAT_red$COD_PROV_RES))

# save the characteristics of the data generation process
save(dati_ISTAT_red, mod_fit_ISTAT, beta_true, sample_prop, props, ref_k, age, prov_res, titstu, file = "Data_gen_process.RData")

