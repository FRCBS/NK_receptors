library(data.table)
library(tidyverse)

# read in all SNP names in the data
all <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed.pvar")

old <- fread("./results/SNP_names_receptors_all.txt", header = F)
new <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names.txt", header = F)

sum(old$V1 %in% new$V1) # 749
sum(old$V1 %in% new$V1) / nrow(old) # 0.8679027

sum(new$V1 %in% old$V1) # 749
sum(new$V1 %in% old$V1) / nrow(new) # 0.05300403

new_filtered <- new[!(new$V1 %in% old$V1),]

write.table(new_filtered, "/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_filtered.txt", sep = "\t", quote = F, col.names = F, row.names = F)








