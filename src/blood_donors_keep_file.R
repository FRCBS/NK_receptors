library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

# add dosage
dosage <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr1_extracted_dosage.raw")

# unify sample names

# # pseudonyms in genotype data do not match in vitro data's IDs 
# first change the IDs according to an ID file, then modify IDs to match the in vitro data's IDs
ids <- fread("/home/nihtiju/work/NK_receptors/data/blood_donors/blood_donors_IDs.txt", header = F)
ids <- ids[match(dosage$IID, ids$V1),]
dosage$NK_ID <- ids$V2

# this for when the IDs are correct & when fixing the code
dosage$sample <- paste0(str_sub(dosage$NK_ID, 1, 2), str_sub(dosage$NK_ID, 6, 7))
dosage$sample <- str_replace_all(dosage$sample, "NK0", "NK")

unique(dosage$sample)

# order dosages to be the same as other data
cytotoxicity_data_factors <- fread("./data/NK_assay_processed/NK_data_parsed.tsv", stringsAsFactors = T)
dosage <- dosage[match(cytotoxicity_data_factors$SAMPLE, dosage$sample),]

unique(dosage$sample)

write.table(unique(dosage[,1:2]), "/home/nihtiju/work/NK_receptors/results/frequencies/in_vitro_keep.txt", quote = F, col.names = F, row.names = F, sep = "\t")


