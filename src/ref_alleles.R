library(data.table)
library(tidyverse)

# make a file to force the same ref /alt alleles for all datasets

var_fin <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered.pvar")
var_eur <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered.pvar")

all(var_fin$ID == var_eur$ID) # T, as it should be
all(var_fin$ALT == var_eur$ALT) # the same already but will change for some populations in association testing if this step is skipped 

var_fin <- var_fin[,c(5,3)] # alt & name

write.table(var_fin, file="/home/nihtiju/work/NK_receptors/results/extract_range/change_ref_allele.txt", sep="\t", quote=F, row.names=F, col.names=F)
