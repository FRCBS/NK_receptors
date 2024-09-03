library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

# make range file for extracting SNPss
snps <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/extracted_ranges_ic1.bim")[,2]
snps_eqtl <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_eQTL.txt", header = F)

all <- unname(unlist(c(snps, snps_eqtl))) # 15297

# remove possible duplicates
sum(duplicated(all)) # 1166
all <- unique(all) # 14131


write.table(all, "/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names.txt", quote = F, col.names = F, row.names = F)








