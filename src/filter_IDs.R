library(data.table)
library(tidyverse)

# filtering individuals - are ID lists for these still ok?

# IDs to be removed from newcastle and poland
remove <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/QC2_nc_donors_poland_remove_these_individuals.txt", header = F)
nc_IDs <- fread("./results/extract_range/newcastle_plink_pgen.psam")
pol_IDs <- fread("./results/extract_range/poland_plink_pgen.psam")

all(remove$V2 %in% c(nc_IDs$IID, pol_IDs$IID)) # F
remove$V2[!(remove$V2 %in% c(nc_IDs$IID, pol_IDs$IID))]

# compare pgen ID lists to older bed/bim/fam ID lists
# before any ID filtering or QC:
nc_IDs_old <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/newcastle_donors/newcastle_donors_chr1_plink.fam")
pol_IDs_old <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/final_results/poland/poland_chr1_plink.fam")

all(nc_IDs$IID %in% nc_IDs_old$V2) # T 
all(pol_IDs$IID %in% pol_IDs_old$V2) # T 

nc_IDs_old$V2[!(nc_IDs_old$V2 %in% nc_IDs$IID)]
# 0
pol_IDs_old$V2[!(pol_IDs_old$V2 %in% pol_IDs$IID)]
# 0

# remove-ID lista ok to filter still

remove_pol_nc <- remove

#----------------------------------------------------------------------------------------------------------------

# for finns

remove <- fread("/home/nithiju/work/HSCT_predictor/results/harmonize_imputed/non_HSCT_samples_to_remove.txt", header = F)
finns_IDs <- fread("/home/nithiju/work/HSCT_predictor/results/NK_cell_SNPs/pgen_NK_SNPs/finns_plink_pgen.psam")
all(remove$V2 %in% finns_IDs$IID) # T


ic1_fam <- fread("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_IC1_IC3_sibling_cohort/IC1_subject_filtered.fam")
ic3_fam <- fread("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_IC1_IC3_sibling_cohort/IC3F.fam")
mcgill_fam <- fread("/home/nithiju/work/HSCT_predictor/data/genotype/McGill/McGill_mergedunion.fam")
register_fam <- fread("/home/nithiju/work/HSCT_predictor/data/genotype/finnish_register/VPU_ILLUMINA_AUG_2018.fam")
finns_old <- c(ic1_fam$V2, ic3_fam$V2, mcgill_fam$V2, register_fam$V2)
all(finns_IDs$IID %in% finns_old) # F
sum(finns_IDs$IID %in% finns_old)
# 1344 

# get all duplicates:
finns_IDs$IID[finns_IDs$IID %like% ":"]
# these are duplicated genotyping IDs, one of the duplicates got file: as a prefix when merging finns
# all from the 4th fileset that was merged, this is mcgill

# only once in the phenotype file (duplicates have been removed from there)
# DNA/genotypes should be the same for both duplicates for they are from the same individual
# -> we can just remove these "4:" individuals
# are only once in finns_old

remove_fin <- finns_IDs[finns_IDs$IID %like% ":",1:2]

# join remove lists together
colnames(remove_pol_nc) <- c("#FID", "IID")
remove_all <- rbind(remove_pol_nc, remove_fin)
write.table(remove_all, file="./results/extract_range/individuals_to_remove.txt", sep="\t", quote=F, row.names=F, col.names=F)































