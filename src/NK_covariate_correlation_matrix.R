library(data.table)
library(tidyverse)

# correlation matrix for covariates

#----------------------------------------------------------------------------------------------------------------------------

# 2 biggest levels in one column, the smaller levels as separate columns

covars_big_levels_merged <- fread("./results/pheno_and_covars/covars_donors_cGvHD_all_binary.txt")

cor_matrix <- cor(covars_big_levels_merged[,3:20], use = "complete.obs")
cor_matrix <- round(cor_matrix, 2)

write.table(cor_matrix, "./results/covars_correlation_matrix.txt", sep = "\t", quote = F, row.names = T, col.names = T) 

