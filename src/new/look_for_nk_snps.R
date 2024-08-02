library(data.table)
library(tidyverse)

# look for NK cell SNPs from the hg19 datasets

# original data: bims
ic1 <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_IC1_IC3_sibling_cohort/IC1_cleaned.bim", stringsAsFactors=F)
ic1_2 <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_IC1_IC3_sibling_cohort/IC1_subject_filtered.bim", stringsAsFactors=F)
ic3 <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_IC1_IC3_sibling_cohort/IC3F.bim", stringsAsFactors=F)
tyks_hus <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_register/VPU_ILLUMINA_AUG_2018.bim", stringsAsFactors=F)
mcgill <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/McGill/McGill_mergedunion.bim", stringsAsFactors=F)
katalonia <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/spanish_sibling/IC2_3S_genotyped_cohort.bim", stringsAsFactors=F)
newcastle <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/english_newcastle/VPU_ILLUMINA_NOV_2019.bim", stringsAsFactors=F)
# newcastle and poland together
newcastle_poland <- read.table("/home/nithiju/work/HSCT_predictor/data/genotype/newcastle_donors_and_poland/VPU_ILLUMINA_OCT_2022_Plus_bed_bim_fam.bim", stringsAsFactors=F)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# look for the SNPs with their names
# rs_codes <- c("rs1051792", "rs1063635", "rs1065075", "rs3828903", "rs1264457", "rs1736927", "rs1063320", "rs2302489", "rs7301582", "rs2734414", "rs1049174", "rs1154831", "rs2255336", "rs1433097", "rs1052248", "rs11575840", "rs3179003", "rs9471577", "rs380573", "rs7246537", "rs8101605") # 21 

receptors <- c("rs2302489", "rs7301582", "rs2734414", "rs1049174", "rs1154831", "rs2255336", "rs1433097", "rs1052248", "rs11575840", "rs3179003", "rs9471577", "rs380573", "rs7246537", "rs8101605") # 14
ligands <- c("rs1051792", "rs1063635", "rs1065075", "rs3828903", "rs1264457", "rs1736927", "rs1063320") # 7

ic1[ic1$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
ic1[ic1$V2 %in% ligands,2] 
# [1] "rs3828903"

ic1_2[ic1_2$V2 %in% receptors,2] 
# 0
ic1_2[ic1_2$V2 %in% ligands,2] 
# 0

ic3[ic3$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
ic3[ic3$V2 %in% ligands,2] 
# [1] "rs1063320" "rs3828903"

tyks_hus[tyks_hus$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
# [4] "rs2302489"  "rs1049174"  "rs2734414" 
# [7] "rs7301582"  "rs380573"   "rs1433097" 
tyks_hus[tyks_hus$V2 %in% ligands,2] 
# [1] "rs1063320" "rs3828903"

mcgill[mcgill$V2 %in% receptors,2] 
# 0
mcgill[mcgill$V2 %in% ligands,2] 
# 0

katalonia[katalonia$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
katalonia[katalonia$V2 %in% ligands,2] 
# [1] "rs1063320" "rs3828903"

katalonia[katalonia$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
katalonia[katalonia$V2 %in% ligands,2] 
# [1] "rs1063320" "rs3828903"

newcastle[newcastle$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
# [4] "rs2302489"  "rs1049174"  "rs2734414" 
# [7] "rs7301582"  "rs380573"   "rs1433097" 
newcastle[newcastle$V2 %in% ligands,2] 
# [1] "rs1063320" "rs3828903"

newcastle_poland[newcastle_poland$V2 %in% receptors,2] 
# [1] "rs1052248"  "rs3179003"  "rs11575840"
# [4] "rs2302489"  "rs1049174"  "rs2734414" 
# [7] "rs7301582"  "rs380573"   "rs1433097" 
newcastle_poland[newcastle_poland$V2 %in% ligands,2] 
# [1] "rs1063320" "rs1264457" "rs3828903"

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# look for the SNPs with their locations in hg19

add_location <- function(bim_file) {
  bim_file$location <- paste0(bim_file[,1], ":", bim_file[,4])
  return(bim_file)
}
ic1 <- add_location(ic1)
ic1_2 <- add_location(ic1_2)
ic3 <- add_location(ic3)
tyks_hus <- add_location(tyks_hus)
mcgill <- add_location(mcgill)
katalonia <- add_location(katalonia)
newcastle <- add_location(newcastle)
newcastle_poland <- add_location(newcastle_poland)

snp_locations <- read.table("./data/NK_cell_SNPs/katarzyna_snp_locations_hg38_hg19.csv", stringsAsFactors=F, header = T, sep = "\t")

snp_locations_recepors <- snp_locations[snp_locations$ligand.receptor == "receptor",]
snp_locations_ligands <- snp_locations[snp_locations$ligand.receptor == "ligand",]

# receptors
ic1[ic1$location %in% snp_locations_recepors$pos_hg_19,2] 
# [1] "rs1052248"     "rs3179003"     "rs11575840"   
# [4] "SEQ-rs380573"  "SEQ-rs1433097"
ic1_2[ic1_2$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "RS1052248"     "RS3179003"     "RS11575840"   
# [4] "SEQ-RS380573"  "SEQ-RS1433097"
ic3[ic3$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "rs1052248"     "rs3179003"     "rs11575840"   
# [4] "SEQ-rs380573"  "SEQ-rs1433097"
tyks_hus[tyks_hus$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "rs1052248"     "rs3179003"    
# [3] "rs11575840"    "rs2302489"    
# [5] "rs1049174"     "rs2734414"    
# [7] "rs7301582"     "rs380573"     
# [9] "GSA-RS8101605" "rs1433097" 
mcgill[mcgill$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "." "." "." "." "." "." "." "." "." "." "."
# [12] "."
katalonia[katalonia$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "rs1052248"  "rs3179003"  "rs11575840"
newcastle[newcastle$location %in% snp_locations_recepors$pos_hg_19,2]
# [1] "rs1052248"     "rs3179003"    
# [3] "rs11575840"    "rs2302489"    
# [5] "rs1049174"     "rs2734414"    
# [7] "rs7301582"     "rs380573"     
# [9] "GSA-RS8101605" "rs1433097" 
newcastle_poland[newcastle_poland$location %in% snp_locations_recepors$pos_hg_19,2]
# 0 (koska on jo hg38)

# ligands
ic1[ic1$location %in% snp_locations_ligands$pos_hg_19,2] 
# [1] "rs3828903"
ic1_2[ic1_2$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "RS1063635" "RS3828903"
ic3[ic3$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "rs1063320" "rs3828903"
tyks_hus[tyks_hus$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "rs1063320"             
# [2] "ILMNSEQ_RS1264457_F2BT"
# [3] "rs3828903" 
mcgill[mcgill$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "." "." "." "." "." "."
katalonia[katalonia$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "rs1063320" "rs3828903"
newcastle[newcastle$location %in% snp_locations_ligands$pos_hg_19,2]
# [1] "rs1063320"             
# [2] "ILMNSEQ_RS1264457_F2BT"
# [3] "rs3828903" 
newcastle_poland[newcastle_poland$location %in% snp_locations_ligands$pos_hg_19,2]
# 0 (koska on jo hg38)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# look for the SNPs with their locations in hg38 (after imputation, before QC)
katalonia <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/katalonia/katalonia_all.bim")
newcastle <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/newcastle/newcastle_all.bim")
newcastle_donors <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/newcastle_donors/newcastle_donors_all.bim")
poland <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/poland/poland_all.bim")
ic1 <- fread("/home/nithiju/work/HSCT_predictor/results/imputation/ic1/ic3_b38_imputed_info_all.bim")
ic3 <- fread("/home/nithiju/work/HSCT_predictor/results/imputation/ic3_imputed_not_filtered/ic3_b38_imputed_info_all.bim")

all(ic1$V2 == ic3$V2) # T
all(katalonia$V2 == newcastle$V2) # T
all(katalonia$V2 == newcastle_donors$V2) # T
all(katalonia$V2 == poland$V2) # T

# riittää etsiä vain yhdestä eur ja yhdestä fin
ic1 <- ic1[ic1$V1 %in% c(6,12,19),]
katalonia <- katalonia[katalonia$V1 %in% c(6,12,19),]

# ic1 <- add_location(ic1)
# katalonia <- add_location(katalonia)
# turned out weird

ic1$location <- paste0(unlist(ic1[,1]), ":", unlist(ic1[,4]))
katalonia$location <- paste0(unlist(katalonia[,1]), ":", unlist(katalonia[,4]))

ic1[ic1$location %in% snp_locations_recepors$pos_hg38,2]
# V2
# 1:  chr6_31588804_T_A
# 2:  chr6_31589151_C_A
# 3:  chr6_31589863_C_T
# 4:  chr6_41336257_A_G
# 5: chr12_10307944_A_T
# 6: chr12_10372766_G_C
# 7: chr12_10379727_T_C
# 8: chr12_10406893_G_T
# 9: chr12_10446203_A_T
# 10: chr12_10449091_C_T
# 11: chr19_54281738_G_A
# 12: chr19_54617369_G_A
# 13: chr19_54637036_G_A
# 14: chr19_54907100_C_T

katalonia[katalonia$location %in% snp_locations_recepors$pos_hg38,2]
# V2
# 1:  chr6_31588804_T_A
# 2:  chr6_31589151_C_A
# 3:  chr6_31589863_C_T
# 4:  chr6_41336257_A_G
# 5: chr12_10307944_A_T
# 6: chr12_10372766_G_C
# 7: chr12_10379727_T_C
# 8: chr12_10406893_G_T
# 9: chr12_10446203_A_T
# 10: chr12_10449091_C_T
# 11: chr19_54281738_G_A
# 12: chr19_54617369_G_A
# 13: chr19_54637036_G_A
# 14: chr19_54907100_C_T
ic1[ic1$location %in% snp_locations_ligands$pos_hg38,2]
# V2
# 1: chr6_29828338_A_C
# 2: chr6_29830972_C_G
# 3: chr6_30490287_G_A
# 4: chr6_31411200_G_A
# 5: chr6_31412154_G_A
# 6: chr6_31496962_G_A
# 7: chr6_31505784_A_G
katalonia[katalonia$location %in% snp_locations_ligands$pos_hg38,2]
# V2
# 1: chr6_29828338_A_C
# 2: chr6_29830972_C_G
# 3: chr6_30490287_G_A
# 4: chr6_31411200_G_A
# 5: chr6_31412154_G_A
# 6: chr6_31496962_G_A
# 7: chr6_31505784_A_G

all(ic1[ic1$location %in% snp_locations_ligands$pos_hg38,2] == katalonia[katalonia$location %in% snp_locations_ligands$pos_hg38,2]) # T 
all(ic1[ic1$location %in% snp_locations_recepors$pos_hg38,2] == katalonia[katalonia$location %in% snp_locations_recepors$pos_hg38,2]) # T
# -> all potitions form from eur & fin

imputed_names_ligands <- ic1[ic1$location %in% snp_locations_ligands$pos_hg38,2]
imputed_names_receptors <- ic1[ic1$location %in% snp_locations_recepors$pos_hg38,2]

# will the names stay thee same when renaming the eur samples?
# they should stay the same because above the names are the same for fin & eur
# & renaming to be the same as fin
renaming_file <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/eur_variants_renaming.txt", header = F)

renaming_file_nk <- renaming_file[renaming_file$V1 %in% c("chr6_29828338_A_C", "chr6_29830972_C_G", "chr6_30490287_G_A", "chr6_31411200_G_A", "chr6_31412154_G_A", "chr6_31496962_G_A", "chr6_31505784_A_G", "chr6_31588804_T_A", "chr6_31589151_C_A", "chr6_31589863_C_T", "chr6_41336257_A_G", "chr12_10307944_A_T", "chr12_10372766_G_C", "chr12_10379727_T_C", "chr12_10406893_G_T", "chr12_10446203_A_T", "chr12_10449091_C_T", "chr19_54281738_G_A", "chr19_54617369_G_A", "chr19_54637036_G_A", "chr19_54907100_C_T"),]
all(renaming_file_nk$V1 == renaming_file_nk$V2) 
# T -> uudelleennimetään tismalleen samoiksi -> tulisi löytyä samalla nimellä siistityistä ja uudelleennimetyistä
all <- as.data.table(c("chr6_29828338_A_C", "chr6_29830972_C_G", "chr6_30490287_G_A", "chr6_31411200_G_A", "chr6_31412154_G_A", "chr6_31496962_G_A", "chr6_31505784_A_G", "chr6_31588804_T_A", "chr6_31589151_C_A", "chr6_31589863_C_T", "chr6_41336257_A_G", "chr12_10307944_A_T", "chr12_10372766_G_C", "chr12_10379727_T_C", "chr12_10406893_G_T", "chr12_10446203_A_T", "chr12_10449091_C_T", "chr19_54281738_G_A", "chr19_54617369_G_A", "chr19_54637036_G_A", "chr19_54907100_C_T"))
receptor_names <- all[all$V1 %in% imputed_names_receptors$V2]
ligand_names <- all[all$V1 %in% imputed_names_ligands$V2]
write.table(all, "./results/SNP_names_in_imputed_datasets.txt", quote = F, col.names = F, row.names = F)
write.table(receptor_names, "./results/SNP_names_in_imputed_datasets_receptors.txt", quote = F, col.names = F, row.names = F)
write.table(ligand_names, "./results/SNP_names_in_imputed_datasets_ligands.txt", quote = F, col.names = F, row.names = F)














