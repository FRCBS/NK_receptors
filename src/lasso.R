library(data.table)
library(tidyverse)
library(glmnet)

# install glmnet

# sudo apt install r-cran-gmlnet
# sudo apt install r-cran-RcppEigen

# sudo apt-get update
# sudo apt-get install -y r-cran-glmnet

# sudo apt-get update
# sudo apt-get install -y r-cran-glmnet

#########################################################################################################################################

setwd("/home/nihtiju/work/NK_receptors/")

dosage <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered_dosage.raw")
# decimals, not strictly 0/1/2 -> pgen dosage conversion has worked ok

donor_IDs <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
donor_IDs_relapse <- fread("./results/ID_lists/donors_train_relapse.txt", header = F)

covars_aGvHD_relapse <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")
covars_cGvHD <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt")

pheno_aGvHD <- fread("./results/pheno_and_covars/pheno_separate_aGvHD_plink.txt")
pheno_cGvHD <- fread("./results/pheno_and_covars/pheno_separate_cGvHD_plink.txt")
pheno_relapse <- fread("./results/pheno_and_covars/pheno_separate_relapse_plink.txt")

# merge these together

dosage_ordered <- dosage[match(pheno_aGvHD$IID, dosage$IID),]

merge_tables <- function(pheno, pheno_col, covars, dosages, ids){
  
  cols <- c(2,pheno_col)
  table <- pheno[,..cols] # = IID and pheno; this did not work without the ..
  # table <- pheno[,c(2,pheno_col)] # does not work
  
  table <- cbind(table, covars[,3:ncol(covars)])
  
  table <- cbind(table, dosages[,7:ncol(dosages)]) # all SNPs
  # table <- cbind(table, dosages[,7:600]) # a few SNP to test this out
  
  table <- table[table$IID %in% ids$V2,]  # leave in only the training set
  
  table <- table[,-1] # remove ID column
  
  table <- table[complete.cases(table),] # remove any rows with NAs
  table <- as.matrix(table) # needed for lasso
  
  return(table)
  
}

aGvHD_all <- merge_tables(pheno_aGvHD, 4, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
aGvHD_severe <- merge_tables(pheno_aGvHD, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 

relapse <- merge_tables(pheno_relapse, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs_relapse)

cGvHD_severe <- merge_tables(pheno_cGvHD, 3, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_severe_broader <- merge_tables(pheno_cGvHD, 4, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_all <- merge_tables(pheno_cGvHD, 5, covars_cGvHD, dosage_ordered, donor_IDs) 


#------------------------------------------------------------------------------------------------------

# lasso models



lasso <- function(data){
  
  y <- data[,1]
  x <- data[,-1]
  
  set.seed(1)
  cv_out <- cv.glmnet(x, y, alpha = 1, family = "binomial") # alpha = 1 -> lasso (alpha = 0 -> ridge)
  
  # the value of λ that results in the smallest cross-validation error is
  # bestlam <- cv_out$lambda.min
  
  # what variables were chosen for the model
  lasso_coef <- predict(cv_out, type = "coefficients", s = cv_out$lambda.min)
  
  coeffs <- data.table(variable = c("(Intercept)", colnames(x)), s1 = as.vector(lasso_coef))
  coeffs_nonZero <- coeffs[coeffs$s1 != 0,]
  
  return(list(coeffs_nonZero, cv_out))
  
}

lasso_aGvHD_all <- lasso(aGvHD_all)
lasso_aGvHD_severe <- lasso(aGvHD_severe)

lasso_relapse <- lasso(relapse)

lasso_cGvHD_severe <- lasso(cGvHD_severe)
lasso_cGvHD_severe_broader <- lasso(cGvHD_severe_broader)
lasso_cGvHD_all <- lasso(cGvHD_all)

# models into a list
train_models <- list(lasso_aGvHD_all[[2]], lasso_aGvHD_severe[[2]], lasso_relapse[[2]], lasso_cGvHD_severe[[2]], lasso_cGvHD_severe_broader[[2]], lasso_cGvHD_all[[2]])
# variables separately
lasso_aGvHD_all <- lasso_aGvHD_all[[1]]
lasso_aGvHD_severe <- lasso_aGvHD_severe[[1]]
lasso_relapse <- lasso_relapse[[1]]
lasso_cGvHD_severe <- lasso_cGvHD_severe[[1]]
lasso_cGvHD_severe_broader <- lasso_cGvHD_severe_broader[[1]]
lasso_cGvHD_all <- lasso_cGvHD_all[[1]]

write.table(lasso_aGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# add everything into one table

lasso_aGvHD_all$pheno <- "aGvHD_all"
lasso_aGvHD_severe$pheno <- "aGvHD_severe"
lasso_relapse$pheno <- "relapse"
lasso_cGvHD_severe$pheno <- "cGvHD_severe"
lasso_cGvHD_severe_broader$pheno <- "cGvHD_severe_broader"
lasso_cGvHD_all$pheno <- "cGvHD_all"

all_train <- rbind(lasso_aGvHD_all, lasso_aGvHD_severe, lasso_relapse, lasso_cGvHD_severe, lasso_cGvHD_severe_broader, lasso_cGvHD_all)
all_train$dataset <- "train"


#------------------------------------------------------------------------------------------------------

# the same for test set

donor_IDs <- fread("./results/ID_lists/test_sample_ids_donors.txt", header = F)
donor_IDs_relapse <- fread("./results/ID_lists/donors_test_relapse.txt", header = F)

# merge files together
aGvHD_all <- merge_tables(pheno_aGvHD, 4, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
aGvHD_severe <- merge_tables(pheno_aGvHD, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 

relapse <- merge_tables(pheno_relapse, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs_relapse)

cGvHD_severe <- merge_tables(pheno_cGvHD, 3, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_severe_broader <- merge_tables(pheno_cGvHD, 4, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_all <- merge_tables(pheno_cGvHD, 5, covars_cGvHD, dosage_ordered, donor_IDs) 



# lasso

lasso_aGvHD_all <- lasso(aGvHD_all)
lasso_aGvHD_severe <- lasso(aGvHD_severe)

lasso_relapse <- lasso(relapse)

lasso_cGvHD_severe <- lasso(cGvHD_severe)
lasso_cGvHD_severe_broader <- lasso(cGvHD_severe_broader)
lasso_cGvHD_all <- lasso(cGvHD_all)

# extract variables
lasso_aGvHD_all <- lasso_aGvHD_all[[1]]
lasso_aGvHD_severe <- lasso_aGvHD_severe[[1]]
lasso_relapse <- lasso_relapse[[1]]
lasso_cGvHD_severe <- lasso_cGvHD_severe[[1]]
lasso_cGvHD_severe_broader <- lasso_cGvHD_severe_broader[[1]]
lasso_cGvHD_all <- lasso_cGvHD_all[[1]]


# save tables

write.table(lasso_aGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# add everything into one table

lasso_aGvHD_all$pheno <- "aGvHD_all"
lasso_aGvHD_severe$pheno <- "aGvHD_severe"
lasso_relapse$pheno <- "relapse"
lasso_cGvHD_severe$pheno <- "cGvHD_severe"
lasso_cGvHD_severe_broader$pheno <- "cGvHD_severe_broader"
lasso_cGvHD_all$pheno <- "cGvHD_all"

all_test <- rbind(lasso_aGvHD_all, lasso_aGvHD_severe, lasso_relapse, lasso_cGvHD_severe, lasso_cGvHD_severe_broader, lasso_cGvHD_all)
all_test$dataset <- "test"

all <- rbind(all_train, all_test)
write.table(all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/all_variables_train_ans_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)


#------------------------------------------------------------------------------------------------------

# are there the same SNPs in the models form the training set & the test set?

common_snps <- function(train, test){
  
  train_results <- fread(train)
  test_results <- fread(test)
  
  common <- train_results[train_results$variable %in% test_results$variable,]
  
  return(common)
}

aGvHD_all <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_test.txt") # Intercept
aGvHD_severe <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_test.txt") # Intercept
relapse <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_test.txt") # Intercept
cGvHD_all <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_test.txt") # Intercept, aGvHD
cGvHD_severe <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_test.txt") # Intercept, aGvHD, population_poland
cGvHD_severe_broader <- common_snps("/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_train.txt", "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_test.txt") # Intercept

aGvHD_all$pheno <- "aGvHD_all"
aGvHD_severe$pheno <- "aGvHD_severe"
relapse$pheno <- "relapse"
cGvHD_all$pheno <- "cGvHD_all"
cGvHD_severe$pheno <- "cGvHD_severe"
cGvHD_severe_broader$pheno <- "cGvHD_severe_broader"
all <- rbind(aGvHD_all, aGvHD_severe, relapse, cGvHD_all, cGvHD_severe, cGvHD_severe_broader)
write.table(all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/all_common_variables_train_ans_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)

