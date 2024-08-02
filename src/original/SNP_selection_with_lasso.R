library(data.table)
library(tidyverse)
library(glmnet)


dosage <- fread("./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage.raw")
# decimals, not strictly 0/1/2 -> pgen dosage conversion has worked ok

donor_IDs <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
donor_IDs_relapse <- fread("./results/ID_lists/donors_train_relapse.txt", header = F)

covars_aGvHD_relapse <- fread("./results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")
covars_cGvHD <- fread("./results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt")

pheno_aGvHD <- fread("./results/pheno_and_covars/pheno_separate_aGvHD_plink.txt")
pheno_cGvHD <- fread("./results/pheno_and_covars/pheno_separate_cGvHD_plink.txt")
pheno_relapse <- fread("./results/pheno_and_covars/pheno_separate_relapse_plink.txt")

# merge these together

dosage_ordered <- dosage[match(pheno_aGvHD$IID, dosage$IID),]

merge_tables <- function(pheno, pheno_col, covars, dosages, ids){
  
  cols <- c(2,pheno_col)
  table <- pheno[,..cols] # = IID and pheno; did not work wihtout the ..
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
  cv_out <- cv.glmnet(x, y, alpha = 1, family = "binomial") # alpha = 1 -> lasso (alpha = 0 -> rodge)
  # plot(cv_out)
  
  # we see that the value of λ that results in the smallest cross-validation error is
  # bestlam <- cv_out$lambda.min

  # what variables were chosen for the model
  lasso_coef <- predict(cv_out, type = "coefficients", s = cv_out$lambda.min)
  # lasso_coef <- lasso_coef[lasso_coef != 0] # an error 
  
  coeffs <- data.table(variable = c("(Intercept)", colnames(x)), s1 = as.vector(lasso_coef))
  coeffs_nonZero <- coeffs[coeffs$s1 != 0,]
  
  return(coeffs_nonZero)
  
}

lasso_aGvHD_all <- lasso(aGvHD_all)
lasso_aGvHD_severe <- lasso(aGvHD_severe)

lasso_relapse <- lasso(relapse)

lasso_cGvHD_severe <- lasso(cGvHD_severe)
lasso_cGvHD_severe_broader <- lasso(cGvHD_severe_broader)
lasso_cGvHD_all <- lasso(cGvHD_all)


# kestä ajaa
# 400 SNPs: 1min40sec
# 600 SNPs: 3min50sec
# kaikki = ~800: 5min20sec
# oikeasti pelkkä train-setti: paljon nopeampi kaikilla SNP

write.table(lasso_aGvHD_all, "./results/select_SNPs_with_lasso/aGvHD_all.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "./results/select_SNPs_with_lasso/aGvHD_severe.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "./results/select_SNPs_with_lasso/relapse.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(lasso_cGvHD_severe, "./results/select_SNPs_with_lasso/cGvHD_severe.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "./results/select_SNPs_with_lasso/cGvHD_severe_broader.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "./results/select_SNPs_with_lasso/cGvHD_all.txt", sep = "\t", quote = F, row.names = F, col.names = T)

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

write.table(lasso_aGvHD_all, "./results/select_SNPs_with_lasso/aGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "./results/select_SNPs_with_lasso/aGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "./results/select_SNPs_with_lasso/relapse_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(lasso_cGvHD_severe, "./results/select_SNPs_with_lasso/cGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "./results/select_SNPs_with_lasso/cGvHD_severe_broader_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "./results/select_SNPs_with_lasso/cGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#------------------------------------------------------------------------------------------------------

# are there the same SNPs in the models form the training set & the test set?

common_snps <- function(train, test){
  
  train_results <- fread(train)
  test_results <- fread(test)
  
  common <- train_results[train_results$variable %in% test_results$variable,]
  
  return(common)
}

aGvHD_all <- common_snps("./results/select_SNPs_with_lasso/aGvHD_all.txt", "./results/select_SNPs_with_lasso/aGvHD_all_test.txt")
aGvHD_severe <- common_snps("./results/select_SNPs_with_lasso/aGvHD_severe.txt", "./results/select_SNPs_with_lasso/aGvHD_severe_test.txt")

relapse <- common_snps("./results/select_SNPs_with_lasso/relapse.txt", "./results/select_SNPs_with_lasso/relapse_test.txt")

cGvHD_all <- common_snps("./results/select_SNPs_with_lasso/cGvHD_all.txt", "./results/select_SNPs_with_lasso/cGvHD_all_test.txt")
cGvHD_severe <- common_snps("./results/select_SNPs_with_lasso/cGvHD_severe.txt", "./results/select_SNPs_with_lasso/cGvHD_severe_test.txt")
cGvHD_severe_broader <- common_snps("./results/select_SNPs_with_lasso/cGvHD_severe_broader.txt", "./results/select_SNPs_with_lasso/cGvHD_severe_broader_test.txt")

