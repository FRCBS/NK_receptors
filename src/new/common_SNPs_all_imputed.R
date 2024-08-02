library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

# find common SNPs between all imputed NK SNPs

finns <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen.pvar") # 14187
uk <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen.pvar") # 11412
spain <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen.pvar") # 11412

# all(uk$V2 == spain$V2) # T -> as they should be

common_snps <- Reduce(intersect, list(finns$ID, uk$ID)) # 11358

write.table(as.data.table(common_snps), file="/home/nihtiju/work/NK_receptors/results/extract_range/common_SNPs_all_imputed.txt", sep="\t", quote=F, row.names=F, col.names=F)


# are the non-common ones because ref and alt allele are the other way around?
# if yes, then they would not have even been extracted
# all SNPs on the lists in total: 863
# out of these found and common: 797 / 863 = 0.93

