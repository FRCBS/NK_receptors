library(data.table)
library(tidyverse)
library(readxl)

# modify pheno and covariate files for association testing with plink

###########################################################################################################################################

# read in pheno
pheno <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/pheno/pheno_first_tx.txt", data.table=F, header=T) 

# have to have a header of #FID IID pheno1 pheno2 etc.

colnames(pheno)
# [1] "ID"                               
# [2] "GenotypingID"                     
# [3] "DonorGenotypingID"                
# [4] "DonorType"                        
# [5] "pheno_aGvHD_severe"               
# [6] "pheno_aGvHD_severe_vs_all"        
# [7] "pheno_aGvHD_all"                  
# [8] "pheno_cGvHD_severe"               
# [9] "pheno_cGvHD_severe_broader"       
# [10] "pheno_cGvHD_severe_vs_all"        
# [11] "pheno_cGvHD_severe_broader_vs_all"
# [12] "pheno_cGvHD_all"                  
# [13] "pheno_relapse"       

pheno_donors <- pheno[,c(3,3,5,7,8,9,12,13)]
colnames(pheno_donors)[1] <- "#FID"
colnames(pheno_donors)[2] <- "IID"

# remove rows where genotyping ID is NA
pheno_donors <- pheno_donors[!(is.na(pheno_donors$IID)),]

# save a separate file for relapse and aGvHD
pheno_cGvHD <- pheno_donors[,c(1,2,5:7)]
pheno_aGvHD <- pheno_donors[,c(1:4)]
pheno_relapse <- pheno_donors[,c(1,2,8)]

pheno_poland <- pheno_cGvHD[,c(1,2,5)] # pheno_cGvHD_all
pheno_poland_agvhd <- pheno_aGvHD[,c(1,2,4)] # pheno_aGvHD_severe_broader_vs_all & pheno_aGvHD_all

# these are what is used in the association analysis:
write.table(pheno_cGvHD, "./results/pheno_and_covars/pheno_separate_cGvHD_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(pheno_aGvHD, "./results/pheno_and_covars/pheno_separate_aGvHD_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(pheno_relapse, "./results/pheno_and_covars/pheno_separate_relapse_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(pheno_poland, "./results/pheno_and_covars/pheno_separate_cGvHD_POLAND_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(pheno_poland_agvhd, "./results/pheno_and_covars/pheno_separate_aGvHD_POLAND_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)

##########################################################################################################################################

# for the covariates from patient_information

# modify covariates: every class of each covariate as a 0/1 column in the covariate file
# encode two biggest levels of a class into one column and the rest separately

# from modified patient_information which includes all the modifying steps done separately in the NK project
covariates <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/covars/covars_first_tx_all_imputed_ages.txt")
# remove donors who do not have genotyping IDs
covariates <- covariates[!(is.na(covariates$DonorGenotypingID)),]

# for the files modified in the patient_information folder

covars_levels <- covariates[,c(3,3)] # gather all levels of all covars into this table, had FID and IID now
colnames(covars_levels) <- c("FID", "IID")

# donortype
unique(covariates$DonorType)
# [1] "sibling"  "haplo"    "register"
covars_levels$donor_sibling <- as.numeric(covariates$DonorType == "sibling") # sib = 1, reg = 0, haplo = 0
covars_levels$donor_haplo <- as.numeric(covariates$DonorType == "haplo") # sib = 0, reg = 0, haplo = 1

# graft
unique(covariates$Graft)
# [1] "PB"    "BM"    "PB+BM" NA 
covariates$Graft[is.na(covariates$Graft)] <- "NONE"
covars_levels$graft_PB <- as.numeric(covariates$Graft == "PB") # PB = 1, BM = 0, PB+BM = 0, NA = 0
covars_levels$graft_PB_BM <- as.numeric(covariates$Graft == "PB+BM") # PB = 0, BM = 0, PB+BM = 1, NA = 0
covars_levels$graft_NA <- as.numeric(covariates$Graft == "NONE") # PB = 0, BM = 0, PB+BM = 0, NA = 1

# D/R sex combo
unique(covariates$sex_combo_binary)
# [1]  0  1 NA
covariates$sex_combo_binary[is.na(covariates$sex_combo_binary)] <- "NONE"
covars_levels$sex_combo_risk <- as.numeric(covariates$sex_combo_binary == "1") # 1 = 1, 0 = 0, NA = 0
covars_levels$sex_combo_NA <- as.numeric(covariates$sex_combo_binary == "NONE") # 1 = 0, 0 = 0, NA = 1

# preconditioning
unique(covariates$preconditioning)
# [1] "MAC"  "RIC"  "sekv" NA     "NMA"
covariates$preconditioning[is.na(covariates$preconditioning)] <- "NONE"
covars_levels$preconditioning_MAC <- as.numeric(covariates$preconditioning == "MAC") # MAC = 1, RIC = 0, sekv = 0, NMA = 0, NA = 0
covars_levels$preconditioning_sekv <- as.numeric(covariates$preconditioning == "sekv") # MAC = 0, RIC = 0, sekv = 1, NMA = 0, NA = 0
covars_levels$preconditioning_NMA <- as.numeric(covariates$preconditioning == "NMA") # MAC = 0, RIC = 0, sekv = 0, NMA = 1, NA = 0
covars_levels$preconditioning_NMA <- as.numeric(covariates$preconditioning == "NONE") # MAC = 0, RIC = 0, sekv = 0, NMA = 0, NA = 1

# population
unique(covariates$population)
# [1] "finland"   "newcastle" "katalonia" "poland" 
covars_levels$population_fin <- as.numeric(covariates$population == "finland") # finland = 1, newcastle = 0, katalonia = 0, poland = 0
covars_levels$population_spain <- as.numeric(covariates$population == "katalonia") # finland = 0, newcastle = 0, katalonia = 1, poland = 0
covars_levels$population_poland <- as.numeric(covariates$population == "poland") # finland = 0, newcastle = 0, katalonia = 0, poland = 1

# add other columns
covars_levels <- cbind(covars_levels, covariates[,c(7,8,10,13,14)])

write.table(covars_levels[,c(1:17,19:20)], "./results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(covars_levels, "./results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt", sep = "\t", quote = F, row.names = F, col.names = T) 


#----------------------------------------------------------------------------------------

# for poland: graft only has levels PB and PB+BM -> these into one column and removing the rest of the graft covar columns

poland <- covars_levels[covars_levels$population_poland == 1,]
poland$graft <- poland$graft_PB

write.table(poland[,c(1:4,8:17,19:21)], "./results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_poland.txt", sep = "\t", quote = F, row.names = F, col.names = T) 
write.table(poland[,c(1:4,8:21)], "./results/pheno_and_covars/covars_donors_cGvHD_all_binary_poland.txt", sep = "\t", quote = F, row.names = F, col.names = T) 















