library(data.table)
library(tidyverse)

# merge SNP name lists for Leena and Katarzyna's NK receptor SNPs

katarzyna <- fread("./results/SNP_names_in_imputed_datasets_receptors.txt", header = F)
leena <- fread("./results/leena_SNPs/snp_names_receptors.txt", header = F)

all <- rbind(katarzyna, leena)
write.table(all, file="./results/SNP_names_receptors_all.txt", sep="\t", quote=F, row.names=F, col.names=F)
