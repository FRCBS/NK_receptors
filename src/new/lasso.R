library(data.table)
library(tidyverse)
library(glmnet)

# install glmnet


# install.packages("glmnet", dependencies = TRUE)
# # ERROR: dependency ‘RcppEigen’ is not available for package ‘glmnet’
# install.packages("RcppEigen")
# 
# # glmnet nettisivut: https://glmnet.stanford.edu/articles/glmnet.html
# install.packages("glmnet", repos = "https://cran.us.r-project.org")
# # Installing package into ‘/home/nihtiju/R/x86_64-pc-linux-gnu-library/4.3’
# # (as ‘lib’ is unspecified)
# # Warning in install.packages :
# #   unable to access index for repository https://cran.us.r-project.org/src/contrib:
# #   cannot open URL 'https://cran.us.r-project.org/src/contrib/PACKAGES'
# # Warning in install.packages :
# #   package ‘glmnet’ is not available for this version of R
# # 
# # A version of this package for your version of R might be available elsewhere,
# # see the ideas at
# # https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
# 
# # siellä teksti:
# 
# # The package is available, but not for the current version of R or for the type of OS (Unix/Windows). To retrieve the information on available versions of package pkg, use
# # av <- available.packages(filters=list())
# # av[av[, "Package"] == pkg, ]
# # in your R session, and look at the ‘Depends’ and ‘OS_type’ fields (there may be more than one matching entry). If the package depends on a version of R later than the one in use, it is possible that an earlier version is available which will work with your version of R: for CRAN look for ‘Old sources’ on the package’s CRAN landing page and manually retrieve an appropriate version (of comparable age to your version of R).
# 
# av <- available.packages(filters=list())
# av[av[, "Package"] == "glmnet", ]
# av[av[, "Package"] == "RcppEigen", ]
# 
# # Package 
# # "glmnet" 
# # Version 
# # "4.1-8" 
# # Priority 
# # NA 
# # Depends 
# # "R (>= 3.6.0), Matrix (>= 1.0-6)" 
# # Imports 
# # "methods, utils, foreach, shape, survival, Rcpp" 
# # LinkingTo 
# # "RcppEigen, Rcpp" 
# # Suggests 
# # "knitr, lars, testthat, xfun, rmarkdown" 
# # Enhances 
# # NA 
# # License 
# # "GPL-2" 
# # License_is_FOSS 
# # NA 
# # License_restricts_use 
# # NA 
# # OS_type 
# # NA 
# # Archs 
# # NA 
# # MD5sum 
# # "98a316b4857ea3f46b44285dc5764068" 
# # NeedsCompilation 
# # "yes" 
# # File 
# # NA 
# # Repository 
# # "https://cloud.r-project.org/src/contrib" 
# 
# R.version
# 
# # platform       x86_64-pc-linux-gnu         
# # arch           x86_64                      
# # os             linux-gnu                   
# # system         x86_64, linux-gnu           
# # status                                     
# # major          4                           
# # minor          3.3                         
# # year           2024                        
# # month          02                          
# # day            29                          
# # svn rev        86002                       
# # language       R                           
# # version.string R version 4.3.3 (2024-02-29)
# # nickname       Angel Food Cake   
# 
# 
# 
# # terminaalisssa
# # sudo apt install r-cran-gmlnet
# # sudo apt install r-cran-RcppEigen
# 
# 
# # sudo apt-get update
# # sudo apt-get install -y r-cran-glmnet
# # tämä toimi!

# sudo apt-get update
# sudo apt-get install -y r-cran-glmnet





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
  cv_out <- cv.glmnet(x, y, alpha = 1, family = "binomial") # alpha = 1 -> lasso (alpha = 0 -> ridge)
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

write.table(lasso_aGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)

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

# save tables

write.table(lasso_aGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_aGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/aGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_relapse, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/relapse_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_severe_broader, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_severe_broader_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(lasso_cGvHD_all, "/home/nihtiju/work/NK_receptors/results/lasso/SNPs_selected/cGvHD_all_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)


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

#########################################################################################################################

# lasso models on training data -> do they predict test data


# train data

donor_IDs <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
donor_IDs_relapse <- fread("./results/ID_lists/donors_train_relapse.txt", header = F)

aGvHD_all <- merge_tables(pheno_aGvHD, 4, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
aGvHD_severe <- merge_tables(pheno_aGvHD, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
relapse <- merge_tables(pheno_relapse, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs_relapse)
cGvHD_severe <- merge_tables(pheno_cGvHD, 3, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_severe_broader <- merge_tables(pheno_cGvHD, 4, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_all <- merge_tables(pheno_cGvHD, 5, covars_cGvHD, dosage_ordered, donor_IDs) 

# test data

donor_IDs_test <- fread("./results/ID_lists/test_sample_ids_donors.txt", header = F)
donor_IDs_relapse_test <- fread("./results/ID_lists/donors_test_relapse.txt", header = F)

aGvHD_all_test <- merge_tables(pheno_aGvHD, 4, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
aGvHD_severe_test <- merge_tables(pheno_aGvHD, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs) 
relapse_test <- merge_tables(pheno_relapse, 3, covars_aGvHD_relapse, dosage_ordered, donor_IDs_relapse)
cGvHD_severe_test <- merge_tables(pheno_cGvHD, 3, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_severe_broader_test <- merge_tables(pheno_cGvHD, 4, covars_cGvHD, dosage_ordered, donor_IDs) 
cGvHD_all_test <- merge_tables(pheno_cGvHD, 5, covars_cGvHD, dosage_ordered, donor_IDs) 


lasso_prediction <- function(train, test){
  
  # train model
  y <- train[,1]
  x <- train[,-1]
  
  set.seed(1)
  cv_out <- cv.glmnet(x, y, alpha = 1, family = "binomial") # alpha = 1 -> lasso (alpha = 0 -> ridge)
  
  # test model
  y_test <- test[,1]
  x_test <- test[,-1]
  
  # tests
  # lasso.pred <- predict(cv_out, s = cv_out$lambda.min, newx = x_test, type = "class") # 0/1, class label corresponding to the maximum probability
  # lasso.pred <- predict(cv_out, s = cv_out$lambda.min, newx = x_test, type = "response") # fitted probabilities
  # 
  # assess.glmnet(lasso.pred, newy = y_test, family = "binomial")
  # assess.glmnet(lasso.pred, newx = x_test, newy = y_test, s = "lambda.min")
  
  
  lasso.pred <- predict(cv_out, s = cv_out$lambda.min, newx = x_test, type = "class") # 0/1, class label corresponding to the maximum probability
  cnf <- confusion.glmnet(lasso.pred, newy = y_test, family = "binomial")
  
  print(cnf)
  
}

lasso_prediction(aGvHD_all, aGvHD_all_test)
#           True
# Predicted   0   1 Total
# 1     490 325   815
# 2       0   2     2
# Total 490 327   817
# 
# Percent Correct:  0.6022

lasso_prediction(aGvHD_severe, aGvHD_severe_test)
#           True
# Predicted   0  1 Total
# 1     490 95   585
# Total 490 95   585
# 
# Percent Correct:  0.8376 

lasso_prediction(relapse, relapse_test)
#           True
# Predicted   0   1 Total
# 1     532 266   798
# Total 532 266   798
# 
# Percent Correct:  0.6667 

lasso_prediction(cGvHD_severe, cGvHD_severe_test)
#           True
# Predicted   0   1 Total
# 1     340 162   502
# 2       9  33    42
# Total 349 195   544
# 
# Percent Correct:  0.6857 

lasso_prediction(cGvHD_severe_broader, cGvHD_severe_broader_test)
#           True
# Predicted   0   1 Total
# 1     338 195   533
# 2      11  30    41
# Total 349 225   574
# 
# Percent Correct:  0.6411

lasso_prediction(cGvHD_all, cGvHD_all_test)
#           True
# Predicted   0   1 Total
# 1     220 106   326
# 2     129 255   384
# Total 349 361   710
# 
# Percent Correct:  0.669





