library(data.table)
library(tidyverse)

# correlation matrix for covariates

#----------------------------------------------------------------------------------------------------------------------------

# 2 biggest levels in one column, the smaller levels as separate columns

covars_big_levels_merged <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt")

cor_matrix <- cor(covars_big_levels_merged[,3:ncol(covars_big_levels_merged)], use = "complete.obs")
cor_matrix <- round(cor_matrix, 2)

write.table(cor_matrix, "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_correlation_matrix.txt", sep = "\t", quote = F, row.names = T, col.names = T) 

#-------------------------------------------------------

# poland test set


covars_big_levels_merged <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_poland.txt")
test_ids <- fread("./results/ID_lists/test_sample_ids_donors.txt", header = F)
covars_big_levels_merged <- covars_big_levels_merged[covars_big_levels_merged$IID %in% test_ids$V2,]
covars_big_levels_merged <- covars_big_levels_merged[,c(1:4,6,12,15:22)]


cor_matrix <- cor(covars_big_levels_merged[,3:ncol(covars_big_levels_merged)], use = "complete.obs")
cor_matrix <- round(cor_matrix, 2)

write.table(cor_matrix, "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_correlation_matrix_poland_testset.txt", sep = "\t", quote = F, row.names = T, col.names = T) 

# poland train set

covars_big_levels_merged <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_cGvHD_all_binary_poland.txt")
train_ids <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
covars_big_levels_merged <- covars_big_levels_merged[covars_big_levels_merged$IID %in% train_ids$V2,]
covars_big_levels_merged <- covars_big_levels_merged[,c(1:4,6,12,15:22)]


cor_matrix <- cor(covars_big_levels_merged[,3:ncol(covars_big_levels_merged)], use = "complete.obs")
cor_matrix <- round(cor_matrix, 2)

write.table(cor_matrix, "/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_correlation_matrix_poland_trainset.txt", sep = "\t", quote = F, row.names = T, col.names = T) 












