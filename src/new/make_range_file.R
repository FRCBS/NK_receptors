library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

# make range file for extracting SNPss
snps <- fread("/home/nihtiju/work/NK_receptors/results/gene_start_end.txt") # gene names
posititons <- fread("/home/nihtiju/work/NK_receptors/results/glist-hg38")

posititons <- posititons[posititons$V4 %in% snps$`checked gene`,]

# save table
write.table(snps, "/home/nihtiju/work/NK_receptors/results/range_file.txt", sep = " ", col.names = F, row.names = F, quote = F)









