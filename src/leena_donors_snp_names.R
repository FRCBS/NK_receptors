library(data.table)
library(tidyverse)

results <- fread("./results/association_testing_train/chosen_snps_train.txt")
# leave in only the 4 selected SNPs
results_cGvHD <- results[results$pheno == "cGvHD_all",]
results_relapse <- results[results$pheno == "relapse",]
results_relapse <- results_relapse[1:2,]
results <- rbind(results_cGvHD, results_relapse)

donors_chr_1 <- fread("./results/leena_donors/leena_donors_chr1.bim")
donors_chr_18 <- fread("./results/leena_donors/leena_donors_chr18.bim")

sum(results$name_in_ours %in% donors_chr_1$V2) # 2
sum(results$name_in_ours %in% donors_chr_18$V2) # 2
# the same names in blood donors as in the HSCT dataseta

write.table(results$name_in_ours, "./results/leena_donors/snps_names.txt", quote = F, row.names = F, col.names = F)

