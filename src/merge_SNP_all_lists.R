library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

old <- fread("./results/SNP_names_receptors_all.txt", header = F)
new <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names.txt", header = F)

all <- rbind(old, new)
write.table(all, file="/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_old_and_new.txt", sep="\t", quote=F, row.names=F, col.names=F)












