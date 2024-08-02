library(data.table)
library(tidyverse)
library(biomaRt)

# read in locations for leena's SNPs
ADGRG1 <- fread('./results/leena_SNPs/locations/ADGRG1_locations.txt', data.table=F) 
CD2 <- fread('./results/leena_SNPs/locations/CD2_locations.txt', data.table=F) 
CD226 <- fread('./results/leena_SNPs/locations/CD226_locations.txt', data.table=F) 
CD244 <- fread('./results/leena_SNPs/locations/CD244_locations.txt', data.table=F) 
CD48 <- fread('./results/leena_SNPs/locations/CD48_locations.txt', data.table=F) 
CD58 <- fread('./results/leena_SNPs/locations/CD58_locations.txt', data.table=F) 
FCGR3A <- fread('./results/leena_SNPs/locations/FCGR3A_locations.txt', data.table=F) 
NKG2A <- fread('./results/leena_SNPs/locations/NKG2A_locations.txt', data.table=F) 

# are there any duplicates? and if so, are their locations the same?
sum(duplicated(ADGRG1$refsnp_id)) # 2
snps_dup <- ADGRG1[duplicated(ADGRG1$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- ADGRG1[ADGRG1$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
ADGRG1 <- ADGRG1[!(duplicated(ADGRG1$refsnp_id)),]

sum(duplicated(CD2$refsnp_id)) # 5
snps_dup <- CD2[duplicated(CD2$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- CD2[CD2$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
CD2 <- CD2[!(duplicated(CD2$refsnp_id)),]

sum(duplicated(CD226$refsnp_id)) # 7
snps_dup <- CD226[duplicated(CD226$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- CD226[CD226$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
CD226 <- CD226[!(duplicated(CD226$refsnp_id)),]

sum(duplicated(CD244$refsnp_id)) # 5
snps_dup <- CD244[duplicated(CD244$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- CD244[CD244$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
CD244 <- CD244[!(duplicated(CD244$refsnp_id)),]

sum(duplicated(CD48$refsnp_id)) # 1
snps_dup <- CD48[duplicated(CD48$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- CD48[CD48$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
CD48 <- CD48[!(duplicated(CD48$refsnp_id)),]

sum(duplicated(CD58$refsnp_id)) # 0
snps_dup <- CD58[duplicated(CD58$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- CD58[CD58$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
# CD58 <- CD58[!(duplicated(CD58$refsnp_id)),]

sum(duplicated(FCGR3A$refsnp_id)) # 14
snps_dup <- FCGR3A[duplicated(FCGR3A$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- FCGR3A[FCGR3A$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
FCGR3A <- FCGR3A[!(duplicated(FCGR3A$refsnp_id)),]

sum(duplicated(NKG2A$refsnp_id)) # 0
snps_dup <- NKG2A[duplicated(NKG2A$refsnp_id),1]
# all have the same chr
# is the pos the same
dup_all <- NKG2A[NKG2A$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
# NKG2A <- NKG2A[!(duplicated(NKG2A$refsnp_id)),]

# join together into one table
leena_all <- rbind(ADGRG1, CD2, CD226, CD244, CD48, CD58, FCGR3A, NKG2A)
leena_all <- leena_all[,c(1,3,4)]
colnames(leena_all) <- c("rs-code", "chr", "pos")

#-------------------------------------------------------------------------------------------------------------------

# add katarzyna's SNPs
# locations only written into './data/Suunnitelma - NK related SNP Katarzyna 2022_01 (ID 85710).docx'
# written into an excel table './data/katarzyna_snp_locations'
katarzyna <- fread('./data/NK_cell_SNPs/katarzyna_snp_locations.csv', data.table=F) 
# the locations fetched by hand -> get them from biomart too and compare to make sure they are correct

ensembl = useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp")
info <- getBM(attributes=c("refsnp_id", "associated_gene", "chr_name", "chrom_start", "chrom_end", "allele"),
              filters = 'snp_filter',
              values = katarzyna$`Rs-code`,
              mart = ensembl)
# remove duplicates
sum(duplicated(info$refsnp_id)) # 110
snps_dup <- info[duplicated(info$refsnp_id),1]
dup_all <- info[info$refsnp_id %in% snps_dup,] # same pos -> ok to remove duplicates
# many have weird chr names -> remove them
dup_all <- dup_all[dup_all$chr_name %in% 1:23,]
# all ok
info <- info[!(duplicated(info$refsnp_id)),]
# leave only rsID, chr and pos
info <- info[,c(1,3,4)]
colnames(info) <- c("rs-code", "chr", "pos")
# order info to be the same as katarzyna
info <- info[match(katarzyna$`Rs-code`, info$`rs-code`),]
# are the two the same?
all(info == katarzyna) # T -> all original locations correct


#-------------------------------------------------------------------------------------------------------------------

# join these to leena's snps
all <- rbind(info, leena_all)
write.table(all, "./results/locations_all_NK_SNPs.txt", sep = "\t", row.names=F, col.names = T, quote = F)
















