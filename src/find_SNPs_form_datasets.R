library(data.table)
library(tidyverse)

ADGRG1 <- read.table(file = "./results/leena_SNPs/locations/ADGRG1_locations.txt", sep="\t", header=T, stringsAsFactors=F)
CD2 <- read.table(file = "./results/leena_SNPs/locations/CD2_locations.txt", sep="\t", header=T, stringsAsFactors=F)
CD226 <- read.table(file = "./results/leena_SNPs/locations/CD226_locations.txt", sep="\t", header=T, stringsAsFactors=F)
CD244 <- read.table(file = "./results/leena_SNPs/locations/CD244_locations.txt", sep="\t", header=T, stringsAsFactors=F)
CD48 <- read.table(file = "./results/leena_SNPs/locations/CD48_locations.txt", sep="\t", header=T, stringsAsFactors=F)
CD58 <- read.table(file = "./results/leena_SNPs/locations/CD58_locations.txt", sep="\t", header=T, stringsAsFactors=F)
FCGR3A <- read.table(file = "./results/leena_SNPs/locations/FCGR3A_locations.txt", sep="\t", header=T, stringsAsFactors=F)
NKG2A <- read.table(file = "./results/leena_SNPs/locations/NKG2A_locations.txt", sep="\t", header=T, stringsAsFactors=F)

# which are ligands, which receptors:
# ligand	CD58/LFA-3
# ligand	CD48
# receptor	FCGR3A (CD16)
# receptor	ADGRG1 (GPR56)
# receptor	CD2
# receptor	CD244 (2B4)
# receptor	CD226 (DNAM-1)
# receptor	NKG2A rs1983526 (correlated to altered outcome of IL-2-based AML therapy)


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

all <- rbind(ADGRG1, CD2, CD226, CD244, FCGR3A, NKG2A)
receptors <- rbind(ADGRG1, CD2, CD226, CD244, FCGR3A, NKG2A)
ligands <- rbind(CD48, CD58)

#-------------------------------------------------------------------------------

# what are their names in imputed finnish datasets (the names should be the same in renamed european datasets too)
ic1 <- fread("/home/nithiju/work/HSCT_predictor/results/imputation/ic1/ic3_b38_imputed_info_all.bim")
ic1 <- ic1[ic1$V1 %in% unique(all$chr_name),]
ic1$location <- paste0(unlist(ic1[,1]), "_", unlist(ic1[,4]))

receptor_names <- ic1[ic1$location %in% receptors$chr_pos,2] # 794
ligand_names <- ic1[ic1$location %in% ligands$chr_pos,2] # 85

#------------------------
# europeans may have some more because their imputation ref was different and because not all found from finns
katalonia <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/katalonia/katalonia_all.bim")
katalonia <- katalonia[katalonia$V1 %in% unique(all$chr_name),]
katalonia$location <- paste0(unlist(katalonia[,1]), "_", unlist(katalonia[,4]))

receptor_names_eur <- katalonia[katalonia$location %in% receptors$chr_pos,] # 838
ligand_names_eur <- katalonia[katalonia$location %in% ligands$chr_pos,] # 95
all(receptor_names_eur$V2 %in% katalonia$V2) # T
all(ligand_names_eur$V2 %in% katalonia$V2) # T

# what are the names after renaming
renaming_file <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/eur_variants_renaming.txt", header = F)
# new value is read from column 2 and the (old) variant ID from column 1
# col 1 = old, col 2 = new name
all(katalonia$V2 %in% renaming_file$V1) # T
all(receptor_names_eur$V2 %in% renaming_file$V1) # T


renaming_file_res <- renaming_file[match(receptor_names_eur$V2, renaming_file$V1),]
renaming_file_lig <- renaming_file[match(ligand_names_eur$V2, renaming_file$V1),]

# are the names the same after renaming?
all(renaming_file_res$V1 == renaming_file_res$V2) # T
all(renaming_file_lig$V1 == renaming_file_lig$V2) # T

all_receptors <- c(receptor_names, renaming_file_res$V2)
all_receptors <- unlist(all_receptors) 
all_receptors <- unique(all_receptors) 

all_ligands <- c(ligand_names, renaming_file_lig$V2)
all_ligands <- unlist(all_ligands) 
all_ligands <- unique(all_ligands) 

write.table(all_receptors, "./results/leena_SNPs/snp_names_receptors.txt", quote = F, row.names = F, col.names = F)
write.table(all_ligands, "./results/leena_SNPs/snp_names_ligands.txt", quote = F, row.names = F, col.names = F)

