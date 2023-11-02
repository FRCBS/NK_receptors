library(data.table)
library(tidyverse)

# find common SNPs between all imputed NK SNPs

finns <- fread("./results/pgen_NK_SNPs/finns_plink_pgen.pvar") # 808
uk <- fread("./results/pgen_NK_SNPs/newcastle_plink_pgen.pvar") # 852
spain <- fread("./results/pgen_NK_SNPs/katalonia_plink_pgen.pvar") # 852

# all(uk$V2 == spain$V2) # T -> as they should be

common_snps <- Reduce(intersect, list(finns$ID, uk$ID)) # 797

write.table(as.data.table(common_snps), file="./results/common_SNPs_all_imputed.txt", sep="\t", quote=F, row.names=F, col.names=F)


# are the non-common ones because ref and alt allele are the other way around?
# if yes, then they would not have even been extracted
# all SNPs on the lists in total: 863
# out of these found and common: 797 / 863 = 0.93

