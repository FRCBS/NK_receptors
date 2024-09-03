library(data.table)
library(tidyverse)

# find eQTL SNP names from our datasets

snps <- fread("/home/nihtiju/work/NK_receptors/results/Blood_eQTL_NK_genes_significantSNPs.tsv")
snps$chr_pos <- paste0(snps$SNPChr, "_", snps$hg38_SNPpos)

# what are SNP names in imputed finnish datasets (the names should be the same in renamed european datasets too)
ic1 <- fread("/home/nihtiju/work/HSCT_predictor/results/imputation/ic1/ic3_b38_imputed_info_all.bim")
ic1$chr_pos <- paste0(ic1$V1, "_", ic1$V4)

# filter SNPs based on position
snps <- snps[snps$chr_pos %in% ic1$chr_pos,]

# some longer alleles in our data -> remove these
ic1 <- ic1[ic1$chr_pos %in% snps$chr_pos,]
unique(ic1$V5)
# [1] "A"                                                                                       
# [2] "T"                                                                                       
# [3] "G"                                                                                       
# [4] "C"                                                                                       
# [5] "CG"                                                                                      
# [6] "ACACAC"                                                                                  
# [7] "ACACACAC"                                                                                
# [8] "CACAATGTGCACATGTACCCTAAAACTTAAAGTATAATAAAAAATAAAAAAATAAAGATAAAACTAAAAAATAAAAAACAGAAAAAAA"
# [9] "GGTTTGGTTTGTAATGTGTACATTACAGAAAACTA"  

unique(ic1$V6)
# [1] "G"      "C"      "A"      "T"      "ATG"    "TAGA"   "CAAAAG"

# remove longer alleles
ic1 <- ic1[ic1$V5 %in% c("A", "T", "G", "C"),]
ic1 <- ic1[ic1$V6 %in% c("A", "T", "G", "C"),]
unique(ic1$V5)
# [1] "A" "T" "G" "C"
unique(ic1$V6)
# [1] "G" "C" "A" "T"


snps <- snps[snps$chr_pos %in% ic1$chr_pos,]
ic1 <- ic1[ic1$chr_pos %in% snps$chr_pos,]

# check alleles
ic1$alleles_found <- F
for (i in 1:nrow(ic1)) {
  
  print(i)
  
  ref <- snps[snps$chr_pos == ic1$chr_pos[i],]
  
  ic1$alleles_found[i] <- ic1$V5[i] %in% ref$AssessedAllele & ic1$V6[i] %in% ref$OtherAllele | 
    ic1$V5[i] %in% ref$OtherAllele & ic1$V6[i] %in% ref$AssessedAllele
  
}

sum(ic1$alleles_found) # 10131
ic1 <- ic1[ic1$alleles_found == T,]

write.table(ic1$V2, "/home/nihtiju/work/NK_receptors/results/extract_range/SNP_names_eQTL.txt", quote = F, col.names = F, row.names = F)
















