library(data.table)
library(tidyverse)

# a table summarizing the cohorts available and their features/characteristics

###############################################################################################################################################################

# read in the processed clinical information
harmonized_processed <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt") # 2719

# and the list of donor IDs in the training set for whom the table is to be assembled
train_donors <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/train_test_split/train_sample_ids_donors.txt", header = F) # 1071

harmonized_nk_donors <- harmonized_processed[harmonized_processed$DonorGenotypingID %in% train_donors$V2,] # 1071


# check with genotyping IDs
genotypes <- fread("./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed.psam", header = T) # 4154, R and D

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


# add others?
# in poland: MUD,MMUD,MSD,MMSD,haplo
# in newcastle (all): MUD,SIB,other,haplo
# in mcgill: 5 haplos

#---------------------------------------------------------------------------------------------------------

# preconditioning

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
diagnosis_names <- fread("/home/nithiju/work/HSCT_predictor/data/clinical_files/translations_etc/diseases_and_classifications.csv")

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

write.table(table_all_T, file="./results/for_paper/table_1.txt", sep="\t", quote=F, row.names=T, col.names=T)


############################################################################################################################################################

# table 1 for test set

# read in the processed clinical information
harmonized_processed <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt") # 2719

# and the list of donor IDs in the training set for whom the table is to be assembled
test_donors <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/train_test_split/test_sample_ids_donors.txt", header = F) # 448

harmonized_nk_donors <- harmonized_processed[harmonized_processed$DonorGenotypingID %in% test_donors$V2,] # 448

# check with the genotype file
genotypes <- fread("./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed.psam", header = T) # 4154, R and D

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
diagnosis_names <- fread("/home/nithiju/work/HSCT_predictor/data/clinical_files/translations_etc/diseases_and_classifications.csv")

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

write.table(table_all_T, file="./results/for_paper/table_1_test.txt", sep="\t", quote=F, row.names=T, col.names=T)


















