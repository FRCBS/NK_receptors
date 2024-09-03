library(data.table)
library(tidyverse)
library(flextable)

setwd("/home/nihtiju/work/NK_receptors/")


# a table summarizing the cohorts available and their features/characteristics

# first for discovery dataset, then for the replication dataset

###############################################################################################################################################################

# read in the processed clinical information
harmonized_processed <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt") # 2719

# and the list of donor IDs in the training set for whom the table is to be assembled
train_donors <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/train_test_split/train_sample_ids_donors.txt", header = F) # 1071

harmonized_nk_donors <- harmonized_processed[harmonized_processed$DonorGenotypingID %in% train_donors$V2,] # 1071


# check with genotyping IDs
genotypes <- fread("./results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed.psam", header = T) # 4154, R and D

# remove the ones not in genotyping files
harmonized_nk_donors <- harmonized_nk_donors[harmonized_nk_donors$DonorGenotypingID %in% genotypes$IID,] # 1045

###############################################################################################################################################################

# create a results table
# row = population
# column = feature

# table <- data.table(finland=character(), spain=character(), uk=character(), poland=character())
table_all <- data.table(population=c("Finland", "Spain", "UK", "Poland"))

# to be transposed in the end

#---------------------------------------------------------------------------------------------------------

# the number of donors in each pop
unique(harmonized_nk_donors$population)
# [1] "ic1"       "ic3"       "mcgill"    "tyks"     
# [5] "hus"       "newcastle" "katalonia" "poland" 

table_all$n_donors <- c(sum(harmonized_nk_donors$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus")),
                        sum(harmonized_nk_donors$population == "katalonia"),
                        sum(harmonized_nk_donors$population == "newcastle"),
                        sum(harmonized_nk_donors$population == "poland")
)

# make tables for each population for the next bits
finland <- harmonized_nk_donors[harmonized_nk_donors$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus"),]
spain <- harmonized_nk_donors[harmonized_nk_donors$population  == "katalonia",]
uk <- harmonized_nk_donors[harmonized_nk_donors$population  == "newcastle",]
poland <- harmonized_nk_donors[harmonized_nk_donors$population  == "poland",]

#---------------------------------------------------------------------------------------------------------

# tx years

get_tx_years <- function(clin_table){
  
  return(paste0(min(clin_table$TxYear, na.rm = T), "–", max(clin_table$TxYear, na.rm = T)))
  
}


table_all$tx_years <- c(get_tx_years(finland),
                        get_tx_years(spain),
                        get_tx_years(uk),
                        get_tx_years(poland))

table_all$tx_years_NA <- c(sum(is.na(finland$TxYear)),
                           sum(is.na(spain$TxYear)),
                           sum(is.na(uk$TxYear)),
                           sum(is.na(poland$TxYear)))


#---------------------------------------------------------------------------------------------------------

# donor and recipient ages

R_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$AgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$AgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$AgeAll, na.rm = T), digits = 2), ")"))
  
}



table_all$recipient_age_median_range <- c(R_age(finland),
                                          R_age(spain),
                                          R_age(uk),
                                          R_age(poland))

table_all$recipient_age_NA <- c(sum(is.na(finland$AgeAll)),
                                sum(is.na(spain$AgeAll)),
                                sum(is.na(uk$AgeAll)),
                                sum(is.na(poland$AgeAll)))



D_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$DonorAgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$DonorAgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$DonorAgeAll, na.rm = T), digits = 2), ")"))
  
}


table_all$donor_age_median_range <- c(D_age(finland),
                                      D_age(spain),
                                      D_age(uk),
                                      NA) # no ages for donors in poland


table_all$donor_age_NA <- c(sum(is.na(finland$DonorAgeAll)),
                            sum(is.na(spain$DonorAgeAll)),
                            sum(is.na(uk$DonorAgeAll)),
                            sum(is.na(poland$DonorAgeAll)))

#---------------------------------------------------------------------------------------------------------

# donor and recipient sex match

sex_match <- function(table){
  
  table$sex_match <- paste0(table$GenderAll, table$DonorGenderAll)
  # coded as M or F 
  # many have NA as well -> remove them
  table$sex_match[table$sex_match %like% "NA"] <- NA
  
  # what is returned in which order
  # D_M+R_M   D_M+R_F   D_F+R_M   D_F+R_F   one or both sexes are NA
  
  return(c(sum(table$sex_match == "MM", na.rm = T), sum(table$sex_match == "FM", na.rm = T), sum(table$sex_match == "MF", na.rm = T), sum(table$sex_match == "FF", na.rm = T), sum(is.na(table$sex_match))))
  
  
}

# get all sex matches in each population into one table
sex_matches <- data.table()
sex_matches$fin <- sex_match(finland)
sex_matches$spain <- sex_match(spain)
sex_matches$uk <- sex_match(uk)
sex_matches$poland <- sex_match(poland)
sex_matches <- t(sex_matches)
colnames(sex_matches) <- c("donor M, recipient M", "donor M, recipient F", "donor F, recipient M", "donor F, recipient F", "either or both sexes NA")

table_all <- cbind(table_all, sex_matches)

#---------------------------------------------------------------------------------------------------------

# graft

get_grafts <- function(table){
  
  # what is returned in which order
  # PB, BM, PB+BM, NA
  
  return(c(sum(table$Graft == "PB", na.rm = T), sum(table$Graft == "BM", na.rm = T), sum(table$Graft == "PB+BM", na.rm = T), sum(is.na(table$Graft))))
  
  
}

# get all in each population into one table
graft <- data.table()
graft$fin <- get_grafts(finland)
graft$spain <- get_grafts(spain)
graft$uk <- get_grafts(uk)
graft$poland <- get_grafts(poland)
graft <- t(graft)
colnames(graft) <- c("graft PB", "graft BM", "graft PB+BM", "graft NA")

table_all <- cbind(table_all, graft)

#---------------------------------------------------------------------------------------------------------

# donor type

get_donor_type <- function(table){
  
  # what is returned in which order
  # sibling, register, haplo, NA
  
  return(c(sum(table$DonorType == "sibling", na.rm = T), sum(table$DonorType == "register", na.rm = T), sum(table$DonorType == "haplo", na.rm = T), sum(is.na(table$DonorType))))
  
  
}

# get all in each population into one table
donor_type <- data.table()
donor_type$fin <- get_donor_type(finland)
donor_type$spain <- get_donor_type(spain)
donor_type$uk <- get_donor_type(uk)
donor_type$poland <- get_donor_type(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
donor_type <- t(donor_type)
colnames(donor_type) <- c("sibling", "register", "haplo", "donor type NA")

table_all <- cbind(table_all, donor_type)

#---------------------------------------------------------------------------------------------------------

# conditioning

unique(harmonized_nk_donors$PreconditioningType)
# [1] "MAC"  "RIC"  "sekv" NA 
# does not contain other types in all data (NMA)

get_preconditioning <- function(table){
  
  # what is returned in which order
  # MAC, RIC, sekv, NA
  
  return(c(sum(table$PreconditioningType == "MAC", na.rm = T), sum(table$PreconditioningType == "RIC", na.rm = T), sum(table$PreconditioningType == "sekv", na.rm = T), sum(is.na(table$PreconditioningType))))
  
  
}

# get all in each population into one table
preconditioning <- data.table()
preconditioning$fin <- get_preconditioning(finland)
preconditioning$spain <- get_preconditioning(spain)
preconditioning$uk <- get_preconditioning(uk)
preconditioning$poland <- get_preconditioning(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
preconditioning <- t(preconditioning)
colnames(preconditioning) <- c("preconditioning MAC", "preconditioning RIC", "preconditioning Other", "preconditioning NA")

table_all <- cbind(table_all, preconditioning)

#---------------------------------------------------------------------------------------------------------

# GvHD profylaxis

covars <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")
covars_finland <- covars[match(finland$DonorGenotypingID, covars$IID),]
covars_uk <- covars[match(uk$DonorGenotypingID, covars$IID),]
covars_spain <- covars[match(spain$DonorGenotypingID, covars$IID),]
covars_poland <- covars[match(poland$DonorGenotypingID, covars$IID),]

gvhd_profylaxis_levels <- function(table){
  
  return(c(sum(table$GvHD_prevention == 1, na.rm = T), 
           sum(table$GvHD_prevention == 2, na.rm = T), 
           sum(table$GvHD_prevention == 3, na.rm = T),
           sum(table$GvHD_prevention == 4, na.rm = T),
           sum(table$GvHD_prevention == 5, na.rm = T),
           sum(is.na(table$GvHD_prevention))))
  
  
}

# get all sex matches in each population into one table
gvhd_profylaxis <- data.table()
gvhd_profylaxis$fin <- gvhd_profylaxis_levels(covars_finland)
gvhd_profylaxis$spain <- gvhd_profylaxis_levels(covars_spain)
gvhd_profylaxis$uk <- gvhd_profylaxis_levels(covars_uk)
gvhd_profylaxis$poland <- gvhd_profylaxis_levels(covars_poland)
gvhd_profylaxis <- t(gvhd_profylaxis)
colnames(gvhd_profylaxis) <- c("profylaxis_1", "profylaxis_2", "profylaxis_3", "profylaxis_4", "profylaxis_5", "profylaxis_missing")

table_all <- cbind(table_all, gvhd_profylaxis)

#---------------------------------------------------------------------------------------------------------

# CMV

cmv_match <- function(table){
  
  table$cmv_match <- paste0(table$CMVdonorAll, "_", table$CMVpatientAll)
  # many have NA as well -> remove them
  table$cmv_match[table$cmv_match %like% "NA"] <- NA
  
  # what is returned 
  # [1] "NA_NA"             "Positive_Positive"
  # [3] "Negative_Positive" "Positive_Negative"
  # [5] "Negative_Negative" "NA_Positive"      
  # [7] "NA_Negative"       "Positive_NA"  
  
  return(c(sum(table$cmv_match == "Positive_Negative", na.rm = T), sum(table$cmv_match %in% c("Positive_Positive", "Negative_Positive", "Negative_Negative"), na.rm = T), sum(is.na(table$cmv_match))))
  
  
}

# get all sex matches in each population into one table
cmv_matches <- data.table()
cmv_matches$fin <- cmv_match(finland)
cmv_matches$spain <- cmv_match(spain)
cmv_matches$uk <- cmv_match(uk)
cmv_matches$poland <- cmv_match(poland)
cmv_matches <- t(cmv_matches)
colnames(cmv_matches) <- c("donor pos, recipient neg", "other combinations", "cmv_missing")

table_all <- cbind(table_all, cmv_matches)

#---------------------------------------------------------------------------------------------------------

# HLA-C1

hla_c1 <- data.table()
hla_c1$fin <- c(sum(covars_finland$HLA_C1 == 1), sum(is.na(covars_finland$HLA_C1)))
hla_c1$spain <- c(sum(covars_spain$HLA_C1 == 1), sum(is.na(covars_spain$HLA_C1)))
hla_c1$uk <- c(sum(covars_uk$HLA_C1 == 1), sum(is.na(covars_uk$HLA_C1)))
hla_c1$poland <- c(sum(covars_poland$HLA_C1 == 1), sum(is.na(covars_poland$HLA_C1)))
hla_c1 <- t(hla_c1)
colnames(hla_c1) <- c("at least one HLA C1 allele", "hla_c1 missing")
table_all <- cbind(table_all, hla_c1)

# HLA-C2

hla_c2 <- data.table()
hla_c2$fin <- c(sum(covars_finland$HLA_C2 == 1), sum(is.na(covars_finland$HLA_C2)))
hla_c2$spain <- c(sum(covars_spain$HLA_C2 == 1), sum(is.na(covars_spain$HLA_C2)))
hla_c2$uk <- c(sum(covars_uk$HLA_C2 == 1), sum(is.na(covars_uk$HLA_C2)))
hla_c2$poland <- c(sum(covars_poland$HLA_C2 == 1), sum(is.na(covars_poland$HLA_C2)))
hla_c2 <- t(hla_c2)
colnames(hla_c2) <- c("at least one HLA C2 allele", "hla_c2 missing")
table_all <- cbind(table_all, hla_c2)

#---------------------------------------------------------------------------------------------------------

# HLA match score / 10

hla_match <- data.table()
hla_match$fin <- c(sum(covars_finland$match_score == 10, na.rm = T), sum(covars_finland$match_score != 10 & !is.na(covars_finland$match_score), na.rm = T), sum(is.na(covars_finland$match_score)))
hla_match$spain <- c(sum(covars_spain$match_score == 10, na.rm = T), sum(covars_spain$match_score != 10 & !is.na(covars_spain$match_score), na.rm = T), sum(is.na(covars_spain$match_score)))
hla_match$uk <- c(sum(covars_uk$match_score == 10, na.rm = T), sum(covars_uk$match_score != 10 & !is.na(covars_uk$match_score), na.rm = T), sum(is.na(covars_uk$match_score)))
hla_match$poland <- c(sum(covars_poland$match_score == 10, na.rm = T), sum(covars_poland$match_score != 10 & !is.na(covars_poland$match_score), na.rm = T), sum(is.na(covars_poland$match_score)))
hla_match <- t(hla_match)
colnames(hla_match) <- c("10/10", "other_match_score", "hla_match_missing")
table_all <- cbind(table_all, hla_match)

#---------------------------------------------------------------------------------------------------------

# endpoints: aGvHD
unique(harmonized_nk_donors$aGvHD_status)
# [1] "7" NA  "2" "4" "3" "1"

get_aGvHD <- function(table){
  
  return(c(sum(table$aGvHD_status == "7", na.rm = T), sum(table$aGvHD_status %in% c("1", "2"), na.rm = T), sum(table$aGvHD_status %in% c("3", "4"), na.rm = T), sum(is.na(table$aGvHD_status))))
  
  
}

# get all in each population into one table
aGvHD <- data.table()
aGvHD$fin <- get_aGvHD(finland)
aGvHD$spain <- get_aGvHD(spain)
aGvHD$uk <- get_aGvHD(uk)
aGvHD$poland <- get_aGvHD(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
aGvHD <- t(aGvHD)
colnames(aGvHD) <- c("aGvHD grade 0","aGvHD grade 1-2", "aGvHD grade 3-4", "aGvHD NA")

table_all <- cbind(table_all, aGvHD)

#---------------------------------------------------------------------------------------------------------

# endpoints: cGvHD
unique(harmonized_nk_donors$cGvHD_status)
# [1] "7"   NA    "2"   "1"   "Yes"

# yes found in 
# all of cGvHD in poland (they have no grade, juts a yes/no)
# some in newcastle
# one hus

get_cGvHD <- function(table){
  
  return(c(sum(table$cGvHD_status == "7", na.rm = T), sum(table$cGvHD_status == "Yes", na.rm = T), sum(table$cGvHD_status == "1", na.rm = T), sum(table$cGvHD_status == "2", na.rm = T), sum(is.na(table$cGvHD_status))))
  
  
}

# get all in each population into one table
cGvHD <- data.table()
cGvHD$fin <- get_cGvHD(finland)
cGvHD$spain <- get_cGvHD(spain)
cGvHD$uk <- get_cGvHD(uk)
cGvHD$poland <- get_cGvHD(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
cGvHD <- t(cGvHD)
colnames(cGvHD) <- c("cGvHD grade 0","cGvHD class not known", "cGvHD grade 1", "cGvHD grade 2", "cGvHD NA")

table_all <- cbind(table_all, cGvHD)


#---------------------------------------------------------------------------------------------------------

# endpoints: relapse

get_relapse <- function(table){
  
  return(c(sum(table$Relapse == "Y", na.rm = T), sum(table$Relapse == "N", na.rm = T), sum(is.na(table$Relapse))))
  
  
}

# get all in each population into one table
relapse <- data.table()
relapse$fin <- get_relapse(finland)
relapse$spain <- get_relapse(spain)
relapse$uk <- get_relapse(uk)
relapse$poland <- get_relapse(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
relapse <- t(relapse)
colnames(relapse) <- c("relapse yes", "relapse no", "relapse NA")

table_all <- cbind(table_all, relapse)

#---------------------------------------------------------------------------------------------------------

# diagnosis

# myeloid/lymphoid/nonmalignant

# endpoints: cGvHD
unique(harmonized_nk_donors$DiagnosisClass)
# [1] "lymphoid"          "myeloid"          
# [3] "nonmalignant"      "myeloid, lymphoid"
# [5] NA 

get_diagnosis_class <- function(table){
  
  return(c(sum(table$DiagnosisClass == "myeloid", na.rm = T), sum(table$DiagnosisClass == "lymphoid", na.rm = T), sum(table$DiagnosisClass == "myeloid, lymphoid", na.rm = T), sum(table$DiagnosisClass == "nonmalignant", na.rm = T), sum(is.na(table$DiagnosisClass))))
  
  
}

diagnosis_class <- data.table()
diagnosis_class$fin <- get_diagnosis_class(finland)
diagnosis_class$spain <- get_diagnosis_class(spain)
diagnosis_class$uk <- get_diagnosis_class(uk)
diagnosis_class$poland <- get_diagnosis_class(poland) 
diagnosis_class <- t(diagnosis_class)
colnames(diagnosis_class) <- c("myeloid","lymphoid", "myeloid and lymphoid", "nonmalignant", "disease NA")

table_all <- cbind(table_all, diagnosis_class)

# malignant/nonmalignant

get_diagnosis_class <- function(table){
  
  return(c(sum(table$DiagnosisClass %in% c("myeloid","lymphoid", "myeloid and lymphoid"), na.rm = T), sum(table$DiagnosisClass == "nonmalignant", na.rm = T), sum(is.na(table$DiagnosisClass))))
  
  
}

diagnosis_class <- data.table()
diagnosis_class$fin <- get_diagnosis_class(finland)
diagnosis_class$spain <- get_diagnosis_class(spain)
diagnosis_class$uk <- get_diagnosis_class(uk)
diagnosis_class$poland <- get_diagnosis_class(poland) 
diagnosis_class <- t(diagnosis_class)
colnames(diagnosis_class) <- c("malignant","nonmalignant", "disease NA")

table_all <- cbind(table_all, diagnosis_class)

# per disease

# get the harmonized disease names
diagnosis_names <- fread("/home/nihtiju/work/HSCT_predictor/data/clinical_files/translations_etc/diagnosis_classification_old/diseases_and_classifications.csv")

get_diagnosis_individually <- function(table){
  
  # get harmonized diagnosis names
  table$diagnosis_in_words <- diagnosis_names[match(table$DiagnosisOriginal, diagnosis_names$`diseases original`),2]
  
  # go over all diagnosis
  numbers <- c()
  diagnosis_all <- unique(diagnosis_names$disease_sanoin)
  diagnosis_all <- diagnosis_all[!(is.na(diagnosis_all))]
  for (i in diagnosis_all) {
    # print(i) # disease in words
    
    numbers <- c(numbers, sum(table$diagnosis_in_words == i, na.rm = T))
    
    
  }
  
  numbers <- c(numbers, sum(is.na(table$diagnosis_in_words)))
  
  return(numbers)
  
}

diagnosis_names_list <- data.table()
diagnosis_names_list$fin <- get_diagnosis_individually(finland)
diagnosis_names_list$spain <- get_diagnosis_individually(spain)
diagnosis_names_list$uk <- get_diagnosis_individually(uk)
diagnosis_names_list$poland <- get_diagnosis_individually(poland) 
diagnosis_names_list <- t(diagnosis_names_list)

diagnosis_all <- unique(diagnosis_names$disease_sanoin)
diagnosis_all <- diagnosis_all[!(is.na(diagnosis_all))]
colnames(diagnosis_names_list) <- c(diagnosis_all, "diagnosis NA")

# remove the diagnosis for whic there are no individuals
diagnosis_names_list_nomissing <- diagnosis_names_list[,colSums(diagnosis_names_list) != 0]
diagnosis_names_list_nomissing_ordered <- diagnosis_names_list_nomissing[,order(colSums(diagnosis_names_list_nomissing), decreasing = T)]
table_all <- cbind(table_all, diagnosis_names_list_nomissing_ordered)


#---------------------------------------------------------------------------------------------------------------------------------------

# transpose table to have populations as columns

table_all_T <- t(table_all)

write.table(table_all_T, file="./results/for_paper/table_1/table_1_train_text.txt", sep="\t", quote=F, row.names=T, col.names=T)


############################################################################################################################################################
############################################################################################################################################################

# table 1 for test set

# read in the processed clinical information
harmonized_processed <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt") # 2719

# and the list of donor IDs in the training set for whom the table is to be assembled
test_donors <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/train_test_split/test_sample_ids_donors.txt", header = F) # 448

harmonized_nk_donors <- harmonized_processed[harmonized_processed$DonorGenotypingID %in% test_donors$V2,] # 448

# check with the genotype file
genotypes <- fread("./results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed.psam", header = T) # 4154, R and D

harmonized_nk_donors <- harmonized_nk_donors[harmonized_nk_donors$DonorGenotypingID %in% genotypes$IID,] # 446

#---------------------------------------------------------------------------------------------------------

# make table 1
# otherwise the same as above for training set
# but preconditining changed: other includes NMA too
# aGvHD changed: now "Yes" too -> a row added for that

# create a results table
# row = population
# column = feature

# table <- data.table(finland=character(), spain=character(), uk=character(), poland=character())
table_all <- data.table(population=c("Finland", "Spain", "UK", "Poland"))

# to be transposed in the end

#---------------------------------------------------------------------------------------------------------

# the number of donors in each pop
unique(harmonized_nk_donors$population)
# [1] "ic3"       "mcgill"    "tyks"      "hus"       "newcastle"
# [6] "katalonia" "poland" 

table_all$n_donors <- c(sum(harmonized_nk_donors$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus")),
                        sum(harmonized_nk_donors$population == "katalonia"),
                        sum(harmonized_nk_donors$population == "newcastle"),
                        sum(harmonized_nk_donors$population == "poland")
)

# make tables for each population for the next bits
finland <- harmonized_nk_donors[harmonized_nk_donors$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus"),]
spain <- harmonized_nk_donors[harmonized_nk_donors$population  == "katalonia",]
uk <- harmonized_nk_donors[harmonized_nk_donors$population  == "newcastle",]
poland <- harmonized_nk_donors[harmonized_nk_donors$population  == "poland",]

#---------------------------------------------------------------------------------------------------------

# tx years

get_tx_years <- function(clin_table){
  
  return(paste0(min(clin_table$TxYear, na.rm = T), "–", max(clin_table$TxYear, na.rm = T)))
  
}


table_all$tx_years <- c(get_tx_years(finland),
                        get_tx_years(spain),
                        get_tx_years(uk),
                        get_tx_years(poland))

table_all$tx_years_NA <- c(sum(is.na(finland$TxYear)),
                           sum(is.na(spain$TxYear)),
                           sum(is.na(uk$TxYear)),
                           sum(is.na(poland$TxYear)))


#---------------------------------------------------------------------------------------------------------

# donor and recipient ages

R_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$AgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$AgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$AgeAll, na.rm = T), digits = 2), ")"))
  
}



table_all$recipient_age_median_range <- c(R_age(finland),
                                          R_age(spain),
                                          R_age(uk),
                                          R_age(poland))

table_all$recipient_age_NA <- c(sum(is.na(finland$AgeAll)),
                                sum(is.na(spain$AgeAll)),
                                sum(is.na(uk$AgeAll)),
                                sum(is.na(poland$AgeAll)))



D_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$DonorAgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$DonorAgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$DonorAgeAll, na.rm = T), digits = 2), ")"))
  
}


table_all$donor_age_median_range <- c(D_age(finland),
                                      D_age(spain),
                                      D_age(uk),
                                      NA) # no ages for donors in poland


table_all$donor_age_NA <- c(sum(is.na(finland$DonorAgeAll)),
                            sum(is.na(spain$DonorAgeAll)),
                            sum(is.na(uk$DonorAgeAll)),
                            sum(is.na(poland$DonorAgeAll)))

#---------------------------------------------------------------------------------------------------------

# donor and recipient sex match

sex_match <- function(table){
  
  table$sex_match <- paste0(table$GenderAll, table$DonorGenderAll)
  # coded as M or F 
  # many have NA as well -> remove them
  table$sex_match[table$sex_match %like% "NA"] <- NA
  
  # what is returned in which order
  # D_M+R_M   D_M+R_F   D_F+R_M   D_F+R_F   one or both sexes are NA
  
  return(c(sum(table$sex_match == "MM", na.rm = T), sum(table$sex_match == "FM", na.rm = T), sum(table$sex_match == "MF", na.rm = T), sum(table$sex_match == "FF", na.rm = T), sum(is.na(table$sex_match))))
  
  
}

# get all sex matches in each population into one table
sex_matches <- data.table()
sex_matches$fin <- sex_match(finland)
sex_matches$spain <- sex_match(spain)
sex_matches$uk <- sex_match(uk)
sex_matches$poland <- sex_match(poland)
sex_matches <- t(sex_matches)
colnames(sex_matches) <- c("donor M, recipient M", "donor M, recipient F", "donor F, recipient M", "donor F, recipient F", "either or both sexes NA")

table_all <- cbind(table_all, sex_matches)

#---------------------------------------------------------------------------------------------------------

# graft

get_grafts <- function(table){
  
  # what is returned in which order
  # PB, BM, PB+BM, NA
  
  return(c(sum(table$Graft == "PB", na.rm = T), sum(table$Graft == "BM", na.rm = T), sum(table$Graft == "PB+BM", na.rm = T), sum(is.na(table$Graft))))
  
  
}

# get all in each population into one table
graft <- data.table()
graft$fin <- get_grafts(finland)
graft$spain <- get_grafts(spain)
graft$uk <- get_grafts(uk)
graft$poland <- get_grafts(poland)
graft <- t(graft)
colnames(graft) <- c("graft PB", "graft BM", "graft PB+BM", "graft NA")

table_all <- cbind(table_all, graft)

#---------------------------------------------------------------------------------------------------------

# donor type

get_donor_type <- function(table){
  
  # what is returned in which order
  # sibling, register, haplo, NA
  
  return(c(sum(table$DonorType == "sibling", na.rm = T), sum(table$DonorType == "register", na.rm = T), sum(table$DonorType == "haplo", na.rm = T), sum(is.na(table$DonorType))))
  
  
}

# get all in each population into one table
donor_type <- data.table()
donor_type$fin <- get_donor_type(finland)
donor_type$spain <- get_donor_type(spain)
donor_type$uk <- get_donor_type(uk)
donor_type$poland <- get_donor_type(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
donor_type <- t(donor_type)
colnames(donor_type) <- c("sibling", "register", "haplo", "donor type NA")

table_all <- cbind(table_all, donor_type)


# add others?
# in poland: MUD,MMUD,MSD,MMSD,haplo
# in newcastle (all): MUD,SIB,other,haplo
# in mcgill: 5 haplos

#---------------------------------------------------------------------------------------------------------

# preconditioning

unique(harmonized_nk_donors$PreconditioningType)
# [1] "MAC"  "RIC"  "sekv" "NMA"

get_preconditioning <- function(table){
  
  # what is returned in which order
  # MAC, RIC, sekv, NA
  
  return(c(sum(table$PreconditioningType == "MAC", na.rm = T), sum(table$PreconditioningType == "RIC", na.rm = T), sum(table$PreconditioningType %in% c("sekv", "NMA"), na.rm = T), sum(is.na(table$PreconditioningType))))
  
  
}

# get all in each population into one table
preconditioning <- data.table()
preconditioning$fin <- get_preconditioning(finland)
preconditioning$spain <- get_preconditioning(spain)
preconditioning$uk <- get_preconditioning(uk)
preconditioning$poland <- get_preconditioning(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
preconditioning <- t(preconditioning)
colnames(preconditioning) <- c("preconditioning MAC", "preconditioning RIC", "preconditioning Other", "preconditioning NA")

table_all <- cbind(table_all, preconditioning)

#---------------------------------------------------------------------------------------------------------

# GvHD profylaxis

covars <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")
covars_finland <- covars[match(finland$DonorGenotypingID, covars$IID),]
covars_uk <- covars[match(uk$DonorGenotypingID, covars$IID),]
covars_spain <- covars[match(spain$DonorGenotypingID, covars$IID),]
covars_poland <- covars[match(poland$DonorGenotypingID, covars$IID),]

gvhd_profylaxis_levels <- function(table){
  
  return(c(sum(table$GvHD_prevention == 1, na.rm = T), 
           sum(table$GvHD_prevention == 2, na.rm = T), 
           sum(table$GvHD_prevention == 3, na.rm = T),
           sum(table$GvHD_prevention == 4, na.rm = T),
           sum(table$GvHD_prevention == 5, na.rm = T),
           sum(is.na(table$GvHD_prevention))))
  
  
}

# get all sex matches in each population into one table
gvhd_profylaxis <- data.table()
gvhd_profylaxis$fin <- gvhd_profylaxis_levels(covars_finland)
gvhd_profylaxis$spain <- gvhd_profylaxis_levels(covars_spain)
gvhd_profylaxis$uk <- gvhd_profylaxis_levels(covars_uk)
gvhd_profylaxis$poland <- gvhd_profylaxis_levels(covars_poland)
gvhd_profylaxis <- t(gvhd_profylaxis)
colnames(gvhd_profylaxis) <- c("profylaxis_1", "profylaxis_2", "profylaxis_3", "profylaxis_4", "profylaxis_5", "profylaxis_missing")

table_all <- cbind(table_all, gvhd_profylaxis)

#---------------------------------------------------------------------------------------------------------

# CMV

cmv_match <- function(table){
  
  table$cmv_match <- paste0(table$CMVdonorAll, "_", table$CMVpatientAll)
  # many have NA as well -> remove them
  table$cmv_match[table$cmv_match %like% "NA"] <- NA
  
  # what is returned 
  # [1] "NA_NA"             "Positive_Positive"
  # [3] "Negative_Positive" "Positive_Negative"
  # [5] "Negative_Negative" "NA_Positive"      
  # [7] "NA_Negative"       "Positive_NA"  
  
  return(c(sum(table$cmv_match == "Positive_Negative", na.rm = T), sum(table$cmv_match %in% c("Positive_Positive", "Negative_Positive", "Negative_Negative"), na.rm = T), sum(is.na(table$cmv_match))))
  
  
}

# get all sex matches in each population into one table
cmv_matches <- data.table()
cmv_matches$fin <- cmv_match(finland)
cmv_matches$spain <- cmv_match(spain)
cmv_matches$uk <- cmv_match(uk)
cmv_matches$poland <- cmv_match(poland)
cmv_matches <- t(cmv_matches)
colnames(cmv_matches) <- c("donor pos, recipient neg", "other combinations", "cmv_missing")

table_all <- cbind(table_all, cmv_matches)

#---------------------------------------------------------------------------------------------------------

# HLA-C1

hla_c1 <- data.table()
hla_c1$fin <- c(sum(covars_finland$HLA_C1 == 1), sum(is.na(covars_finland$HLA_C1)))
hla_c1$spain <- c(sum(covars_spain$HLA_C1 == 1), sum(is.na(covars_spain$HLA_C1)))
hla_c1$uk <- c(sum(covars_uk$HLA_C1 == 1), sum(is.na(covars_uk$HLA_C1)))
hla_c1$poland <- c(sum(covars_poland$HLA_C1 == 1), sum(is.na(covars_poland$HLA_C1)))
hla_c1 <- t(hla_c1)
colnames(hla_c1) <- c("at least one HLA C1 allele", "hla_c1 missing")
table_all <- cbind(table_all, hla_c1)

# HLA-C2

hla_c2 <- data.table()
hla_c2$fin <- c(sum(covars_finland$HLA_C2 == 1), sum(is.na(covars_finland$HLA_C2)))
hla_c2$spain <- c(sum(covars_spain$HLA_C2 == 1), sum(is.na(covars_spain$HLA_C2)))
hla_c2$uk <- c(sum(covars_uk$HLA_C2 == 1), sum(is.na(covars_uk$HLA_C2)))
hla_c2$poland <- c(sum(covars_poland$HLA_C2 == 1), sum(is.na(covars_poland$HLA_C2)))
hla_c2 <- t(hla_c2)
colnames(hla_c2) <- c("at least one HLA C2 allele", "hla_c2 missing")
table_all <- cbind(table_all, hla_c2)

#---------------------------------------------------------------------------------------------------------

# HLA match score / 10

hla_match <- data.table()
hla_match$fin <- c(sum(covars_finland$match_score == 10, na.rm = T), sum(covars_finland$match_score != 10 & !is.na(covars_finland$match_score), na.rm = T), sum(is.na(covars_finland$match_score)))
hla_match$spain <- c(sum(covars_spain$match_score == 10, na.rm = T), sum(covars_spain$match_score != 10 & !is.na(covars_spain$match_score), na.rm = T), sum(is.na(covars_spain$match_score)))
hla_match$uk <- c(sum(covars_uk$match_score == 10, na.rm = T), sum(covars_uk$match_score != 10 & !is.na(covars_uk$match_score), na.rm = T), sum(is.na(covars_uk$match_score)))
hla_match$poland <- c(sum(covars_poland$match_score == 10, na.rm = T), sum(covars_poland$match_score != 10 & !is.na(covars_poland$match_score), na.rm = T), sum(is.na(covars_poland$match_score)))
hla_match <- t(hla_match)
colnames(hla_match) <- c("10/10", "other_match_score", "hla_match_missing")
table_all <- cbind(table_all, hla_match)


#---------------------------------------------------------------------------------------------------------

# endpoints: aGvHD
unique(harmonized_nk_donors$aGvHD_status)
# [1] "7" NA  "2" "4" "3" "1"
# [1] "7"   "2"   "1"   "3"   "4"   "Yes" NA

get_aGvHD <- function(table){
  
  return(c(sum(table$aGvHD_status == "7", na.rm = T), sum(table$aGvHD_status %in% c("1", "2"), na.rm = T), sum(table$aGvHD_status %in% c("3", "4"), na.rm = T), sum(table$aGvHD_status == "Yes", na.rm = T), sum(is.na(table$aGvHD_status))))
  
  
}

# get all in each population into one table
aGvHD <- data.table()
aGvHD$fin <- get_aGvHD(finland)
aGvHD$spain <- get_aGvHD(spain)
aGvHD$uk <- get_aGvHD(uk)
aGvHD$poland <- get_aGvHD(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
aGvHD <- t(aGvHD)
colnames(aGvHD) <- c("aGvHD grade 0","aGvHD grade 1-2", "aGvHD grade 3-4", "aGvHD class not known", "aGvHD NA")

table_all <- cbind(table_all, aGvHD)



#---------------------------------------------------------------------------------------------------------

# endpoints: cGvHD
unique(harmonized_nk_donors$cGvHD_status)
# [1] "7"   NA    "2"   "1"   "Yes"

# yes found in 
# all of cGvHD in poland (they have no grade, juts a yes/no)
# some in newcastle
# one hus

get_cGvHD <- function(table){
  
  return(c(sum(table$cGvHD_status == "7", na.rm = T), sum(table$cGvHD_status == "Yes", na.rm = T), sum(table$cGvHD_status == "1", na.rm = T), sum(table$cGvHD_status == "2", na.rm = T), sum(is.na(table$cGvHD_status))))
  
  
}

# get all in each population into one table
cGvHD <- data.table()
cGvHD$fin <- get_cGvHD(finland)
cGvHD$spain <- get_cGvHD(spain)
cGvHD$uk <- get_cGvHD(uk)
cGvHD$poland <- get_cGvHD(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
cGvHD <- t(cGvHD)
colnames(cGvHD) <- c("cGvHD grade 0","cGvHD class not known", "cGvHD grade 1", "cGvHD grade 2", "cGvHD NA")

table_all <- cbind(table_all, cGvHD)


#---------------------------------------------------------------------------------------------------------

# endpoints: relapse

get_relapse <- function(table){
  
  return(c(sum(table$Relapse == "Y", na.rm = T), sum(table$Relapse == "N", na.rm = T), sum(is.na(table$Relapse))))
  
  
}

# get all in each population into one table
relapse <- data.table()
relapse$fin <- get_relapse(finland)
relapse$spain <- get_relapse(spain)
relapse$uk <- get_relapse(uk)
relapse$poland <- get_relapse(poland) # has some mismatched sib & reg -> do we want these more specifically listed?
relapse <- t(relapse)
colnames(relapse) <- c("relapse yes", "relapse no", "relapse NA")

table_all <- cbind(table_all, relapse)

#---------------------------------------------------------------------------------------------------------

# diagnosis

# myeloid/lymphoid/nonmalignant

# endpoints: cGvHD
unique(harmonized_nk_donors$DiagnosisClass)
# [1] "lymphoid"          "myeloid"          
# [3] "nonmalignant"      "myeloid, lymphoid"
# [5] NA 

get_diagnosis_class <- function(table){
  
  return(c(sum(table$DiagnosisClass == "myeloid", na.rm = T), sum(table$DiagnosisClass == "lymphoid", na.rm = T), sum(table$DiagnosisClass == "myeloid, lymphoid", na.rm = T), sum(table$DiagnosisClass == "nonmalignant", na.rm = T), sum(is.na(table$DiagnosisClass))))
  
  
}

diagnosis_class <- data.table()
diagnosis_class$fin <- get_diagnosis_class(finland)
diagnosis_class$spain <- get_diagnosis_class(spain)
diagnosis_class$uk <- get_diagnosis_class(uk)
diagnosis_class$poland <- get_diagnosis_class(poland) 
diagnosis_class <- t(diagnosis_class)
colnames(diagnosis_class) <- c("myeloid","lymphoid", "myeloid and lymphoid", "nonmalignant", "disease NA")

table_all <- cbind(table_all, diagnosis_class)

# malignant/nonmalignant

get_diagnosis_class <- function(table){
  
  return(c(sum(table$DiagnosisClass %in% c("myeloid","lymphoid", "myeloid and lymphoid"), na.rm = T), sum(table$DiagnosisClass == "nonmalignant", na.rm = T), sum(is.na(table$DiagnosisClass))))
  
  
}

diagnosis_class <- data.table()
diagnosis_class$fin <- get_diagnosis_class(finland)
diagnosis_class$spain <- get_diagnosis_class(spain)
diagnosis_class$uk <- get_diagnosis_class(uk)
diagnosis_class$poland <- get_diagnosis_class(poland) 
diagnosis_class <- t(diagnosis_class)
colnames(diagnosis_class) <- c("malignant","nonmalignant", "disease NA")

table_all <- cbind(table_all, diagnosis_class)

# per disease

# get the harmonized disease names
diagnosis_names <- fread("/home/nihtiju/work/HSCT_predictor/data/clinical_files/translations_etc/diagnosis_classification_old/diseases_and_classifications.csv")

get_diagnosis_individually <- function(table){
  
  # get harmonized diagnosis names
  table$diagnosis_in_words <- diagnosis_names[match(table$DiagnosisOriginal, diagnosis_names$`diseases original`),2]
  
  # go over all diagnosis
  numbers <- c()
  diagnosis_all <- unique(diagnosis_names$disease_sanoin)
  diagnosis_all <- diagnosis_all[!(is.na(diagnosis_all))]
  for (i in diagnosis_all) {
    # print(i) # disease in words
    
    numbers <- c(numbers, sum(table$diagnosis_in_words == i, na.rm = T))
    
    
  }
  
  numbers <- c(numbers, sum(is.na(table$diagnosis_in_words)))
  
  return(numbers)
  
}

diagnosis_names_list <- data.table()
diagnosis_names_list$fin <- get_diagnosis_individually(finland)
diagnosis_names_list$spain <- get_diagnosis_individually(spain)
diagnosis_names_list$uk <- get_diagnosis_individually(uk)
diagnosis_names_list$poland <- get_diagnosis_individually(poland) 
diagnosis_names_list <- t(diagnosis_names_list)

diagnosis_all <- unique(diagnosis_names$disease_sanoin)
diagnosis_all <- diagnosis_all[!(is.na(diagnosis_all))]
colnames(diagnosis_names_list) <- c(diagnosis_all, "diagnosis NA")

# remove the diagnosis for whic there are no individuals
diagnosis_names_list_nomissing <- diagnosis_names_list[,colSums(diagnosis_names_list) != 0]
diagnosis_names_list_nomissing_ordered <- diagnosis_names_list_nomissing[,order(colSums(diagnosis_names_list_nomissing), decreasing = T)]
table_all <- cbind(table_all, diagnosis_names_list_nomissing_ordered)

#---------------------------------------------------------------------------------------------------------------------------------------

# transpose table to have populations as columns

table_all_T <- t(table_all)

write.table(table_all_T, file="./results/for_paper/table_1/table_1_test_text.txt", sep="\t", quote=F, row.names=T, col.names=T)

########################################################################################################################################
########################################################################################################################################

# table 1 for train and test set together


# column names
table_train <- fread('./results/for_paper/table_1/table_1_train_text.txt', data.table=F, header = T)
table_test <- fread('./results/for_paper/table_1/table_1_test_text.txt', data.table=F, header = T)

# test table has one extra row for aGvHD grade unknown
# add that row to train
table_train <- rbind(table_train[1:43,], c("aGvHD class not known", 0, 0, 0, 0), table_train[44:nrow(table_train),]) # this messes up the row numbers
rownames(table_train) <- 1:nrow(table_train)

# get the part with no diseases
table_train <- table_train[1:53,]
table_test <- table_test[1:53,]

# join tables together
# also order the populations to be from the biggest to the smallest
all <- cbind(table_train[,1:2], table_test$Finland, table_train$UK, table_test$UK, table_train$Spain, table_test$Spain, table_train$Poland, table_test$Poland)

colnames(all) <- c("population", "Finland_train", "Finland_test", "UK_train", "UK_test", "Spain_train", "Spain_test", "Poland_train", "Poland_test")

all <- all[,c(1,1,2:9)] # to have indentations in the flex table

#---------------------------------------------------------------------------------------------------------------------------------------
# diseases

# disease names harmonized for the training set in
diagnosis_train <- fread("./results/for_paper/table_1/supple_table_diseases_for_translating_names_harmonized.txt", data.table=F, header = T, stringsAsFactors = F) # made by hand

diagnosis_test <- fread("./results/for_paper/table_1/supple_table_diseases_for_translating_names_test_harmonized.txt", data.table=F, header = T, stringsAsFactors = F)  # made by hand
# NAs as "diagnosis missing"

diagnosis_test_new <- diagnosis_test[!(diagnosis_test$harmonized %in% diagnosis_train$`translated disease name`),]
diagnosis_test_ordered <- diagnosis_test[match(diagnosis_train$`translated disease name`, diagnosis_test$harmonized),]

rows <- 1:nrow(diagnosis_test_ordered)
rows <- rows[is.na(diagnosis_test_ordered$harmonized)]

for (i in rows) {
  
  print(i)
  diagnosis_test_ordered[i,4:7] <- 0
  
}

diagnosis_train_all <- cbind(diagnosis_train[,c(3,3,4)], diagnosis_test_ordered$Finland, diagnosis_train$UK, diagnosis_test_ordered$UK, diagnosis_train$Spain, diagnosis_test_ordered$Spain, diagnosis_train$Poland, diagnosis_test_ordered$Poland)

diagnosis_test_all <- cbind(diagnosis_test_new[,c(3,3,4)], c(0,0,0,0), diagnosis_test_new$UK, c(0,0,0,0), diagnosis_test_new$Spain, c(0,0,0,0), diagnosis_test_new$Poland, c(0,0,0,0))


colnames(diagnosis_train_all) <- colnames(all)
colnames(diagnosis_test_all) <- colnames(all)

diagnosis_all <- rbind(diagnosis_train_all, diagnosis_test_all)

diagnosis_all$all_sum <- rowSums(diagnosis_all[,3:ncol(diagnosis_all)])
diagnosis_all <- diagnosis_all[order(diagnosis_all$all_sum, decreasing = T),]


# leave only  the top 5
table_diseases_5 <- diagnosis_all[1:5,]

# get the rest as "other"
other <- diagnosis_all[6:nrow(diagnosis_all),]
# add to top 5
table_diseases_5 <- rbind(table_diseases_5, table_diseases_5[1,]) # add one row, modify the name to be other next
table_diseases_5[6,1:2] <- "Other"
table_diseases_5[6,3:10] <- colSums(other[,3:10]) # for individual populations

# remove the sum of all populations & the column for harmonizing disease names
table_diseases_5 <- table_diseases_5[,c(1:10)]

table <- rbind(all, table_diseases_5)

#-------------------------------------------------------------------------------------------------------------------------------------

# make a flextable

first <- c("Number of HSCT donors, n", 
           "HSCT time, years", "HSCT time, years", 
           "Recipient age in years, median (range)", "Recipient age in years, median (range)", 
           "Donor age in years, median (range)", "Donor age in years, median (range)",
           "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)",
           "Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)",
           "Donor type, n (%)", "Donor type, n (%)", "Donor type, n (%)", "Donor type, n (%)", 
           "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", 
           "GvHD profylaxis, n (%)","GvHD profylaxis, n (%)","GvHD profylaxis, n (%)","GvHD profylaxis, n (%)","GvHD profylaxis, n (%)","GvHD profylaxis, n (%)",
           "Donor-recipient CMV combination, n (%)","Donor-recipient CMV combination, n (%)","Donor-recipient CMV combination, n (%)",
           "Donor with at least one HLA-C1 allele, n (%)","Donor with at least one HLA-C1 allele, n (%)",
           "Donor with at least one HLA-C2 allele, n (%)","Donor with at least one HLA-C2 allele, n (%)",
           "HLA-match, n (%)","HLA-match, n (%)","HLA-match, n (%)",
           
           "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", 
           "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", 
           "Relapse, n (%)", "Relapse, n (%)", "Relapse, n (%)", 
           "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", " ") 
# the very last one as " " to get the bottom border visible later in the flextable

second <- c("", 
            "", "Missing, n (%)", 
            "", "Missing, n (%)", 
            "", "Missing, n (%)",
            "Male-male", "Male-female", "Female-male", "Female-female", "Missing",
            "Peripheral blood", "Bone marrow", "Both", "Missing", 
            "Sibling","Register", "Haplo", "Missing", 
            "Myeloablative", "Reduced intensity", "Other", "Missing", 
            "1","2","3","4","5","Missing",
            "Donor positive, Recipient negative","Other combinations","Missing",
            "","Missing",
            "","Missing",
            "10/10", "Other", "Missing",
            "Grade 0", paste0("Grade ", as.roman(1), "-", as.roman(2)), paste0("Grade ", as.roman(3), "-", as.roman(4)), "Grade unknown", "Missing", 
            "grade 0", "Yes, classification unknown", "Limited", "Extensive", "Missing", 
            "Yes", "No", "Missing", 
            "Acute myeloid leukemia", "Acute lymphoblastic leukemia", "Myelodysplastic syndrome", "Chronic myeloid leukemia", "Multiple myeloma", "Other")

table[,1] <- first
table[,2] <- second

colnames(table) <- c(" ", "  ", "Finland", "    ", "UK", "     ", "Spain", "      ", "Poland", "       ")

# add %
rows <- c(3,5,7:nrow(table))

for (i in rows) {
  
  for (j in 3:10) {
    
    # i = row
    # j = col
    
    value <- (as.numeric(table[i,j]) / as.numeric(table[1,j])) * 100
    
    pasted <- paste0(table[i,j], " (", round(value), ")")
    
    table[i,j] <- pasted
    
    
    
  }
  
}

# headers
table_flex <- flextable(table)
table_flex <- add_header_row(table_flex, values = c("", "", rep(c("Discovery", "Replication"), 4)), top = F)

table_flex <- align(table_flex, part = "header", align = "left")

# merge vertical duplicated names
table_flex <- merge_v(table_flex, j = c(" ", "  "))
table_flex <- valign(table_flex, j = 1, valign = "top")

# add footnote
table_flex <- add_footer_lines(table_flex, "GvHD, graft-versus-host disease; aGvHD, acute GvHD; cGvHD, chronic GvHD; CMV, cytomagalovirus")
table_flex <- footnote(table_flex, i = c(5,7), j = 2,
                       value = as_paragraph(
                         c("Missing ages were imputed, see Materials and methods")
                       ),
                       ref_symbols = c("a"),
                       part = "body")



# table_flex <- footnote(table_flex, i = 25, j = 1,
#                        value = as_paragraph(
#                          c("1: CSA + MTX +/- ATG +/- MMF or CSA + MTX + steroid or CSA + MTX + ECP or evero + MTX + MMF \n
#                            2: CSA +/- ATG or CSA + ECP or CSA + steroidi or Tacro +/- ATG \n
#                            3: CSA + MMF +/- ATG or CSA + MMF +/- steroid or Tacro + MMF or Tacro + Sirolimus \n
#                            4: All combinations with PtCy \n
#                            5: Other combinations")
#                        ),
#                        ref_symbols = c("b"),
#                        part = "body")
# too big, creates problems when saving as an image

table_flex <- footnote(table_flex, i = 25, j = 1,
                       value = as_paragraph(
                         c("1: CSA + MTX +/- ATG +/- MMF or CSA + MTX + steroid or CSA + MTX + ECP or evero + MTX + MMF, 2: CSA +/- ATG or CSA + ECP or CSA + steroidi or Tacro +/- ATG, 3: CSA + MMF +/- ATG or CSA + MMF +/- steroid or Tacro + MMF or Tacro + Sirolimus, 4: All combinations with PtCy, 5: Other combinations")
                       ),
                       ref_symbols = c("b"),
                       part = "body")


table_flex <- footnote(table_flex, i = 54, j = 1,
                       value = as_paragraph(
                         c("Five most frequent diagnoses are presented here")
                       ),
                       ref_symbols = c("c"),
                       part = "body")

table_flex_onePage <- table_flex
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)

set_flextable_defaults(background.color = "white")

save_as_image(table_flex, path = "./results/for_paper/table_1/table_1_train_test.png", bg = "white")
save_as_docx(table_flex, path = "./results/for_paper/table_1/table_1_train_test_docx.docx", align = "left")
# and changing thd document to be in landscape format manually after saving the table
# -> wide enough to see all columns

save_as_docx(table_flex_onePage, path = "./results/for_paper/table_1/table_1_train_test_docx_onePage.docx", align = "left")


save_as_pptx(table_flex, path = "./results/for_paper/table_1/table_1_train_test_pptx.pptx")















