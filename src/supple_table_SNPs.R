library(data.table)
library(tidyverse)

# make an empty table to gather all SNPs

all <- data.table(gene=character(), rs_code=character(), name_in_ours=character(), chr=numeric(), pos_hg38=numeric(), in_finns=logical(), in_spain=logical(), in_uk=logical(), in_poland=logical())

# read in katarzyna's snps
katarzyna <- fread("./data/NK_cell_SNPs/Katarzyna_rs_codes_and_their_names_in_our_datasets.txt") # made by hand
# leave only receptors
katarzyna <- katarzyna[katarzyna$type == "receptor",]
# replace random text with NA
katarzyna[4,3:6] <- NA


# in which populations can these be found?
# read in population variant files

files_nk <- list.files(path = "./results/pgen_NK_SNPs", pattern = "plink_pgen.pvar", full.names = T)
files_nk
# [1] "./results/pgen_NK_SNPs/all_datasets_plink_pgen.pvar" # not filtered for common variants
# [2] "./results/pgen_NK_SNPs/finns_plink_pgen.pvar"       
# [3] "./results/pgen_NK_SNPs/katalonia_plink_pgen.pvar"   
# [4] "./results/pgen_NK_SNPs/newcastle_plink_pgen.pvar"   
# [5] "./results/pgen_NK_SNPs/poland_plink_pgen.pvar"  
# all imputed NK receptor SNPs, not filtered for common variants
files_nk <- files_nk[2:5]

katarzyna$in_finns <- NA
katarzyna$in_spain <- NA
katarzyna$in_uk <- NA
katarzyna$in_poland <- NA

for (i in 1:length(files_nk)) {
  
  var <- fread(files_nk[i])
  # order to be the same as snp list
  var <- var[match(katarzyna$name_in_ours, var$ID),]
  
  katarzyna[,i+6] <- !(is.na(var$ID))
  
}

# add chr and pos separtely

katarzyna$chr <- NA
katarzyna$pos <- NA

split <- str_split(katarzyna$location_in_hg38, ":")

for (i in 1:length(split)) {
  
  katarzyna$chr[i] <- split[[i]][1]
  katarzyna$pos[i] <- split[[i]][2]
  
}

katarzyna <- katarzyna[,c("Gene", "rs_code", "name_in_ours", "chr", "pos", "in_finns", "in_spain", "in_uk", "in_poland")]
colnames(katarzyna) <- colnames(all)

all <- rbind(all, katarzyna)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

# the same for leena's SNPs

files <- list.files(path = "./data/NK_cell_SNPs/leena", pattern = ".txt", all.files = T, full.names = T)
get_snps <- function(file){
  
  # read file
  snp_table <- read.table(file = file, sep="\t", header=F, stringsAsFactors=F)
  
  gene <- str_split(file, "/")[[1]][5]
  gene <- str_split(gene, "_")[[1]][1]
  
  snp_table$gene <- gene
  
  return(snp_table)
}


ADGRG1 <- get_snps(files[1])
CD2 <- get_snps(files[2])
CD226 <- get_snps(files[3])
CD244 <- get_snps(files[4])
CD48 <- get_snps(files[5])
CD58 <- get_snps(files[6])
FCGR3A <- get_snps(files[7])
NKG2A <- get_snps(files[8])

# merge all & only include receptors
all_leena <- rbind(ADGRG1, CD2, CD226, CD244, FCGR3A, NKG2A)
colnames(all_leena)[1] <- "rs_code"

all_leena_unique_dupGenes <- distinct(all_leena) # 915 -> some snps listed for more than one gene
all_leena_unique <- distinct(all_leena, rs_code, .keep_all=TRUE) # 879

# add multiple genes in the gene column where appropriate
for (i in 1:nrow(all_leena_unique)) {
  
  snp <- all_leena_unique$rs_code[i]
  
  genes <- all_leena_unique_dupGenes[all_leena_unique_dupGenes$rs_code == snp,2]
  
  all_leena_unique$gene[i] <- paste(genes, collapse=' ')
  
}


# get positions (not found for all)

# read in rs-code, chr, pos
ADGRG1 <- fread('./results/leena_SNPs/locations/ADGRG1_locations.txt', data.table=F) 
CD2 <- fread('./results/leena_SNPs/locations/CD2_locations.txt', data.table=F) 
CD226 <- fread('./results/leena_SNPs/locations/CD226_locations.txt', data.table=F) 
CD244 <- fread('./results/leena_SNPs/locations/CD244_locations.txt', data.table=F) 
CD48 <- fread('./results/leena_SNPs/locations/CD48_locations.txt', data.table=F) 
CD58 <- fread('./results/leena_SNPs/locations/CD58_locations.txt', data.table=F) 
FCGR3A <- fread('./results/leena_SNPs/locations/FCGR3A_locations.txt', data.table=F) 
NKG2A <- fread('./results/leena_SNPs/locations/NKG2A_locations.txt', data.table=F) 
# remove duplicates 
ADGRG1 <- ADGRG1[!(duplicated(ADGRG1$refsnp_id)),]
CD2 <- CD2[!(duplicated(CD2$refsnp_id)),]
CD226 <- CD226[!(duplicated(CD226$refsnp_id)),]
CD244 <- CD244[!(duplicated(CD244$refsnp_id)),]
CD48 <- CD48[!(duplicated(CD48$refsnp_id)),]
CD58 <- CD58[!(duplicated(CD58$refsnp_id)),]
FCGR3A <- FCGR3A[!(duplicated(FCGR3A$refsnp_id)),]
NKG2A <- NKG2A[!(duplicated(NKG2A$refsnp_id)),]

# join chr and pos together into one string
pos <- function(biomart_results){
  
  biomart_results$chr_pos <- paste0(biomart_results$chr_name, "_", biomart_results$chrom_start)
  
  return(biomart_results)
  
}

ADGRG1 <- pos(ADGRG1) 
CD2 <- pos(CD2) 
CD226 <- pos(CD226) 
CD244 <- pos(CD244) 
CD48 <- pos(CD48) 
CD58 <- pos(CD58) 
FCGR3A <- pos(FCGR3A) 
NKG2A <- pos(NKG2A) 

# add gene name cols
ADGRG1$gene <- "ADGRG1"
CD2$gene <- "CD2"
CD226$gene <- "CD226"
CD244$gene <- "CD244"
CD48$gene <- "CD48"
CD58$gene <- "CD58"
FCGR3A$gene <- "FCGR3A"
NKG2A$gene <- "NKG2A"

# join together into one table
leena_all <- rbind(ADGRG1, CD2, CD226, CD244, FCGR3A, NKG2A)

leena_all <- leena_all[,c(1,3,4,8,7)]
colnames(leena_all) <- c("rs_code", "chr", "pos", "gene", "chr_pos")

# order to be the same as the original unique list
# remove duplicates first
leena_all_unique_rs <- distinct(leena_all, rs_code, .keep_all=TRUE) # 861 
leena_all_unique_chrPos <- distinct(leena_all, chr_pos, .keep_all=TRUE) # 861 
# leena_all_unique_dupGenes <- distinct(leena_all) # 897 -> jollekin snplle geeni eri / useampi geeni, kuten yllä alkuperäisessä listassa
leena_all_matched <- leena_all[match(all_leena_unique$rs_code, leena_all$rs_code),] 

all_leena_unique$chr <- leena_all_matched$chr
all_leena_unique$pos <- leena_all_matched$pos
all_leena_unique$chr_pos <- leena_all_matched$chr_pos

# get SNP names in our datasets
snp_names <- fread("./results/leena_SNPs/snp_names_receptors.txt", data.table=F, header = F) 
# get only chr_pos for these
snp_names$chr_pos <- NA
split <- str_split(snp_names$V1, "_")

for (i in 1:length(split)) {
  
  name <- paste(split[[i]][1], split[[i]][2], sep = "_")
  name_noCHR <- str_sub(name, 4, nchar(name))
  snp_names$chr_pos[i] <- name_noCHR
  
}

sum(duplicated(snp_names$V1)) # 0
sum(duplicated(snp_names$chr_pos)) # 3

dup <- snp_names[duplicated(snp_names$chr_pos),2]

# ok to remove the second ones, now known why bims do not match the numbers
snp_names <- snp_names[match(all_leena_unique$chr_pos, snp_names$chr_pos),]

# add name in ours to leena's table/list
all_leena_unique$name_in_ours <- snp_names$V1

# in which populations can these be found?
# read in population bim files

files_nk <- list.files(path = "./results/pgen_NK_SNPs", pattern = "plink_pgen.pvar", full.names = T)
files_nk
# [1] "./results/pgen_NK_SNPs/all_datasets_plink_pgen.pvar"
# [2] "./results/pgen_NK_SNPs/finns_plink_pgen.pvar"       
# [3] "./results/pgen_NK_SNPs/katalonia_plink_pgen.pvar"   
# [4] "./results/pgen_NK_SNPs/newcastle_plink_pgen.pvar"   
# [5] "./results/pgen_NK_SNPs/poland_plink_pgen.pvar" 
files_nk <- files_nk[2:5]

all_leena_unique$in_finns <- NA
all_leena_unique$in_spain <- NA
all_leena_unique$in_uk <- NA
all_leena_unique$in_poland <- NA

for (i in 1:length(files_nk)) {
  
  var <- fread(files_nk[i])
  
  # order to be the same as snp list
  var <- var[match(all_leena_unique$name_in_ours, var$ID),]
  all_leena_unique[,i+6] <- !(is.na(var$ID)) # finns: 791
  
}

all_leena_unique <- all_leena_unique[,c("gene", "rs_code", "name_in_ours", "chr", "pos", "in_finns", "in_spain", "in_uk", "in_poland")]
colnames(all_leena_unique) <- colnames(all)
all <- rbind(all, all_leena_unique)

# all snps, their locations/positions, name in ours, if found from either finns or europeans or both, variants marked in more than one gene on one row
write.table(all, file="./results/for_paper/supple_table_all_snps.txt", sep="\t", quote=F, row.names=F, col.names=T)


