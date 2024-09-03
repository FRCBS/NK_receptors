library(data.table)
library(tidyverse)

results <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_chosen.txt", header = T)

donors_chr_1 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr1.bim")
donors_chr_6 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr6.bim")
donors_chr_7 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr7.bim")
donors_chr_12 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr12.bim")
donors_chr_18 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr18.bim")


sum(results$ID %in% donors_chr_1$V2) # 8
sum(results$ID %in% donors_chr_6$V2) # 1
sum(results$ID %in% donors_chr_7$V2) # 1
sum(results$ID %in% donors_chr_12$V2) # 1
sum(results$ID %in% donors_chr_18$V2) # 2
# the same names in blood donors as in the HSCT dataseta

write.table(results$ID, "/home/nihtiju/work/NK_receptors/results/blood_donors/snps_names.txt", quote = F, row.names = F, col.names = F)

