library(data.table)
library(tidyverse)
library(biomaRt)

# kokeilu biomartilla snippien sijaintia

# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")


files <- list.files(path = "./data/NK_cell_SNPs/leena", pattern = ".txt", all.files = T, full.names = T)

ensembl = useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp")

get_locations <- function(file){
   
   # read file
   snp_table <- read.table(file = file, sep="\t", header=F, stringsAsFactors=F)
   
   # get additional info with biomart
   info <- getBM(attributes=c("refsnp_id", "associated_gene", "chr_name", "chrom_start", "chrom_end", "allele"),
                      filters = 'snp_filter',
                      values = snp_table$V1,
                      mart = ensembl)
   return(info)
}


ADGRG1 <- get_locations(files[1])
CD2 <- get_locations(files[2])
CD226 <- get_locations(files[3])
CD244 <- get_locations(files[4])
CD48 <- get_locations(files[5])
CD58 <- get_locations(files[6])
FCGR3A <- get_locations(files[7])
NKG2A <- get_locations(files[8])

# save results
write.table(ADGRG1, "./results/leena_SNPs/locations/ADGRG1_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(CD2, "./results/leena_SNPs/locations/CD2_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(CD226, "./results/leena_SNPs/locations/CD226_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(CD244, "./results/leena_SNPs/locations/CD244_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(CD48, "./results/leena_SNPs/locations/CD48_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(CD58, "./results/leena_SNPs/locations/CD58_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(FCGR3A, "./results/leena_SNPs/locations/FCGR3A_locations.txt", quote = F, sep = "\t", row.names = F)
write.table(NKG2A, "./results/leena_SNPs/locations/NKG2A_locations.txt", quote = F, sep = "\t", row.names = F)












