library(data.table)
library(tidyverse)
library(readxl)

# modify pheno and covariate files for association testing with plink

###########################################################################################################################################

# read in pheno
pheno <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/pheno/pheno_first_tx.txt", data.table=F, header=T) 

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
# covariates <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/covars/covars_first_tx_all_imputed_ages.txt")
covariates <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/covars/covars_first_tx_all_imputed_ages_more.txt")
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

# CMV
covariates$CMV <- paste0(covariates$CMVdonorAll, "_", covariates$CMVpatientAll)
unique(covariates$CMV)
# [1] "Positive_Positive" "NA_Positive"       "Negative_Positive"
# [4] "Negative_Negative" "Positive_Negative" "NA_NA"            
# [7] "NA_Negative"       "Positive_NA"       "Negative_NA"
covars_levels$CMV_risk <- as.numeric(covariates$CMV == "Positive_Negative") # not risk is 0 (c("Positive_Positive", "Negative_Positive", "Negative_Negative"))
covars_levels$CMV_missing <- as.numeric(covariates$CMV %in% c("NA_Positive", "NA_NA", "NA_Negative", "Positive_NA", "Negative_NA"))


# # aGvHD
# unique(covariates$aGvHD)
# # [1]  0  1 NA
# covariates$aGvHD[is.na(covariates$aGvHD)] <- "NONE"
# covars_levels$aGvHD_yes <- as.numeric(covariates$aGvHD == "1") # 1 = 1, 0 = 0, NA = 0
# covars_levels$aGvHD_NA <- as.numeric(covariates$aGvHD == "NONE") # 1 = 0, 0 = 0, NA = 1
# # in association, aGvHD_NA is being excluded in plink (no reason given) --> removnig it from here
# also correlates too much with CMV missing in UK train

#----------------------------------

# add HLA-C alleles

# using best alleles
# column for C1: at least one C1 allele
# column for C2: at least one C2 allele

files <- list.files("/home/nihtiju/work/HSCT_predictor/results/HLA_imputation/imputed_HLA/best_alleles_datasetwise", full.names = T)
files <- files[files %like% "_hla_c.txt"]

best <- data.frame("sample.id" = character(), "allele1" = character(), "allele2" = character(), "prob" = numeric(), "matching" = numeric())

for (i in 1:length(files)) {
  
  table <- fread(files[i])
  best <- rbind(best, table)
  
}


best$allele1 <- paste0("C*", best$allele1)
best$allele2 <- paste0("C*", best$allele2)

ligand.groups <- list(
  C1='C\\*01\\:02|C\\*03\\:02|C\\*03\\:03|C\\*03\\:04|C\\*07\\:01|C\\*07\\:02|C\\*07\\:04|C\\*08\\:02|C\\*12\\:02|C\\*12\\:03|C\\*14\\:02|C\\*16\\:01',
  C2='C\\*02\\:02|C\\*04\\:01|C\\*04\\:06|C\\*05\\:01|C\\*06\\:02|C\\*15\\:02|C\\*15\\:05|C\\*16\\:02|C\\*17\\:01|C\\*17\\:03',
  Bw4_80I='A\\*23\\:01|A\\*24\\:02|A\\*24\\:07|A\\*25\\:01|A\\*32\\:01|B\\*27\\:02|B\\*38\\:01|B\\*49\\:01|B\\*51\\:01|B\\*57\\:01',
  A3_A11='A\\*03\\:01|A\\*11\\:01'
)


best$C1 <- as.numeric(grepl(ligand.groups[[1]], best$allele1) | grepl(ligand.groups[[1]], best$allele2))
best$C2 <- as.numeric(grepl(ligand.groups[[2]], best$allele1) | grepl(ligand.groups[[1]], best$allele2))

best <- best[match(covars_levels$IID, best$sample.id),]
covars_levels$HLA_C1 <- best$C1
covars_levels$HLA_C2 <- best$C2

#----------------------------------

# add GvHD profylaxis
profylaxis <- read_xlsx("/home/nihtiju/work/HSCT_predictor/data/clinical_files/translations_etc/GvHD_profylaxis_classification/Kopio_GvHD_estolääkitys_luokittelu__urpu_062024_korjaukset.xlsx")


harmonized_processed <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt")
harmonized_processed <- harmonized_processed[match(covars_levels$IID, harmonized_processed$DonorGenotypingID),]
harmonized_processed <- harmonized_processed[,c(17,26)]

all(unique(harmonized_processed$GvHDpreventionOriginal) %in% profylaxis$`Nimi original`) # F

unique(harmonized_processed$GvHDpreventionOriginal)[!(unique(harmonized_processed$GvHDpreventionOriginal) %in% profylaxis$`Nimi original`)]
# [1] "" 
profylaxis$`Nimi original`[109] <- ""

all(unique(harmonized_processed$GvHDpreventionOriginal) %in% profylaxis$`Nimi original`) # T


profylaxis <- profylaxis[match(harmonized_processed$GvHDpreventionOriginal, profylaxis$`Nimi original`),]


harmonized_processed$GvHD_prevention <- profylaxis$`luokitus loppuun 2`
covars_levels$GvHD_prevention <- harmonized_processed$GvHD_prevention

#----------------------------------

# add other columns which are already numeric
covars_levels <- cbind(covars_levels, covariates[,c(7,8,10,13,14,15)])

write.table(covars_levels[,-"aGvHD"], "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt", sep = "\t", quote = F, row.names = F, col.names = T) # exclude aGvHD from this
write.table(covars_levels, "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt", sep = "\t", quote = F, row.names = F, col.names = T)


#----------------------------------------------------------------------------------------

# for poland: graft only has levels PB and PB+BM -> these into one column and removing the rest of the graft covar columns

poland <- covars_levels[covars_levels$population_poland == 1,]
poland$graft <- poland$graft_PB

# donor sibling and donor haplo are almost the same --> correlate too uch in association testing
# remove donor haplo 
# sibling = 1 = sibling
# sibling = 0 = haplo or register

# still too big correlation
# /home/nihtiju/work/NK_receptors/src/covariate_correlation_matrix.R 
# donor sibling correlates a lot with GvHD profylaxis and HLA match score
# -> remove col donor sibling

write.table(poland[,c(1:2,8:22,24:27)], "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_poland.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(poland[,c(1:2,8:27)], "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_poland.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#----------------------------------------------------------------------------------------

# for newcastle: graft NA and CMV NA correlate too much

newcastle <- covars_levels[covars_levels$population_fin == 0 & covars_levels$population_spain == 0 & covars_levels$population_poland == 0,]
testi <- newcastle[,c(2,7,17)]

# both are almost fully 0
# 7 have 1 in both, 4 more in CMV missing
# change CMV_missing to cmv or graft missing

newcastle <- newcastle[,-"graft_NA"]
colnames(newcastle)[16] <- "CMV_or_graft_missing"

write.table(newcastle[,-"aGvHD"], "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary_newcastle.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(newcastle, "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_newcastle.txt", sep = "\t", quote = F, row.names = F, col.names = T)






