library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

#----------------------------------------------------------------------------

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
assoc_results <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")

#----------------------------------------------------------------------------

# save the association results with additional info

SNP_info_ordered <- SNP_info[match(assoc_results$ID, SNP_info$ID),]
assoc_results <- cbind(assoc_results[,1:7], SNP_info_ordered[,7:ncol(SNP_info_ordered)], assoc_results[,8:ncol(assoc_results)])

write.table(assoc_results, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered_INFO.csv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(assoc_results, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered_INFO.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# add betas
betas <- assoc_results[,c(1,7,8,9,13:16)]
# add beta and its 95% CI
betas$beta <- log(betas$OR_all)
betas$beta_lower <- log(betas$L95_all) # to log odds format
betas$beta_upper <- log(betas$U95_all)

write.table(betas, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered_INFO_betas.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#----------------------------------------------------------------------------

# association results for all populations

files <- list.files("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/all_train/", pattern = ".glm.logistic.hybrid", full.names = T)
files <- files[!(files %like% "adjusted")]

assoc_results <- assoc_results[1:13,] # contains SNPs twice (train and test results) -> only one needed
assoc_results$pheno_file <- paste0("_pheno_", assoc_results$pheno, "_a")

test <- fread(files[1])
all <- test[1,]
all$rsID <- "test"
all$pheno <- "test"
all <- all[-1,]

for (i in 1:length(files)) {
  
  file <- fread(files[i])
  file_pheno <- str_split(files[i], "\\.")[[1]][2]
  res_sub <- assoc_results[assoc_results$pheno_file %like% paste0("_", file_pheno, "_a"),]
  
  file <- file[file$ID %in% res_sub$ID,]
  file <- file[file$P < 0.05,]
  file <- file[file$TEST != "ADD",]
  
  file$pheno <- res_sub$pheno[1]
  
  SNP_info_ordered <- SNP_info[match(file$ID, SNP_info$ID),]
  file$rsID <- SNP_info_ordered$rsID
  
  all <- rbind(all, file)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/covariate_results_chosen_SNPs_train_all_datasets.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#------------------------------------------------------------------------

# association results for all populations
# keep all covariates


SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
assoc_results <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")
SNP_info_ordered <- SNP_info[match(assoc_results$ID, SNP_info$ID),]
assoc_results <- cbind(assoc_results[,1:7], SNP_info_ordered[,7:ncol(SNP_info_ordered)], assoc_results[,8:ncol(assoc_results)])

files <- list.files("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/all_train/", pattern = ".glm.logistic.hybrid", full.names = T)
files <- files[!(files %like% "adjusted")]

assoc_results <- assoc_results[1:13,] # contains SNPs twice (train and test results) -> only one needed
assoc_results$pheno_file <- paste0("_pheno_", assoc_results$pheno, "_a")

test <- fread(files[1])
all <- test[1,]
all$rsID <- "test"
all$pheno <- "test"
all <- all[-1,]

for (i in 1:length(files)) {
  
  file <- fread(files[i])
  file_pheno <- str_split(files[i], "\\.")[[1]][2]
  # res_sub <- assoc_results[assoc_results$pheno_file %like% file_pheno,]
  res_sub <- assoc_results[assoc_results$pheno_file %like% paste0("_", file_pheno, "_a"),]
  
  file <- file[file$ID %in% res_sub$ID,]
  file <- file[file$TEST != "ADD",]
  
  file$pheno <- res_sub$pheno[1]
  
  SNP_info_ordered <- SNP_info[match(file$ID, SNP_info$ID),]
  file$rsID <- SNP_info_ordered$rsID
  
  all <- rbind(all, file)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/covariate_results_chosen_SNPs_train_all_datasets_all_covariates.txt", sep = "\t", quote = F, row.names = F, col.names = T)

############################################################################################################################

# get multiple testing p-values

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
assoc_results <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")
SNP_info_ordered <- SNP_info[match(assoc_results$ID, SNP_info$ID),]
assoc_results <- cbind(assoc_results[,1:7], SNP_info_ordered[,7:ncol(SNP_info_ordered)], assoc_results[,8:ncol(assoc_results)])

files <- list.files("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/all_train/", pattern = ".glm.logistic.hybrid.adjusted", full.names = T)

assoc_results <- assoc_results[1:13,] # contains SNPs twice (train and test results) -> only one needed
assoc_results$pheno_file <- paste0("_pheno_", assoc_results$pheno, "_a")

test <- fread(files[1])
all <- test[1,]
all$rsID <- "test"
all$pheno <- "test"
all <- all[-1,]

for (i in 1:length(files)) {
  
  file <- fread(files[i])
  file_pheno <- str_split(files[i], "\\.")[[1]][2]
  res_sub <- assoc_results[assoc_results$pheno_file %like% paste0("_", file_pheno, "_a"),]
  
  file <- file[file$ID %in% res_sub$ID,]
  file$rsID <- res_sub$rsID
  file$pheno <- res_sub$pheno[1]
  
  all <- rbind(all, file)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/multiple_testing_chosen_SNPs_train_all_datasets.txt", sep = "\t", quote = F, row.names = F, col.names = T)








