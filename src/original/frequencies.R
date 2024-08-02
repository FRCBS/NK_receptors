library(data.table)
library(tidyverse)

# calculate genotype frequencies and allele frequencies for the selected 2 SNPs

##########################################################################################################################################

# in vitro blood donors

# load data in
leena_data_factors <- fread("./data/NK_assay_processed/NK_data_parsed.tsv", stringsAsFactors = T) # the results are the same if factors or not
# add dosage
dosage_chr_1 <- fread("./results/leena_donors/leena_donors_chr1_dosage.raw")
dosage_chr_18 <- fread("./results/leena_donors/leena_donors_chr18_dosage.raw")
# unify sample names
dosage_chr_1$sample <- paste0(str_sub(dosage_chr_1$IID, 1, 2), str_sub(dosage_chr_1$IID, 6, 7))
dosage_chr_1$sample <- str_replace_all(dosage_chr_1$sample, "NK0", "NK")
dosage_chr_18$sample <- dosage_chr_1$sample
# order dosages to be the same as other data
dosage_chr_1 <- dosage_chr_1[match(leena_data_factors$SAMPLE, dosage_chr_1$sample),]
dosage_chr_18 <- dosage_chr_18[match(leena_data_factors$SAMPLE, dosage_chr_18$sample),]

# NK1 in leena's data but missing from the genotype data
leena_data_factors$chr1_161726992_G_C_G_0 <- as.factor(as.numeric(dosage_chr_1$chr1_161726992_G_C_G == 0))
leena_data_factors$chr1_161726992_G_C_G_1 <- as.factor(as.numeric(dosage_chr_1$chr1_161726992_G_C_G == 1))
leena_data_factors$chr1_161726992_G_C_G_2 <- as.factor(as.numeric(dosage_chr_1$chr1_161726992_G_C_G == 2))

leena_data_factors$chr1_161727208_C_G_C_0 <- as.factor(as.numeric(dosage_chr_1$chr1_161727208_C_G_C == 0))
leena_data_factors$chr1_161727208_C_G_C_1 <- as.factor(as.numeric(dosage_chr_1$chr1_161727208_C_G_C == 1))
leena_data_factors$chr1_161727208_C_G_C_2 <- as.factor(as.numeric(dosage_chr_1$chr1_161727208_C_G_C == 2))

leena_data_factors$chr18_70203271_C_A_C_0 <- as.factor(as.numeric(dosage_chr_18$chr18_70203271_C_A_C == 0))
leena_data_factors$chr18_70203271_C_A_C_1 <- as.factor(as.numeric(dosage_chr_18$chr18_70203271_C_A_C == 1))

leena_data_factors$chr18_70204107_A_C_A_0 <- as.factor(as.numeric(dosage_chr_18$chr18_70204107_A_C_A == 0))
leena_data_factors$chr18_70204107_A_C_A_1 <- as.factor(as.numeric(dosage_chr_18$chr18_70204107_A_C_A == 1))

# get cell lines separately
# using only K562 here
leena_data_factors_K562 <- leena_data_factors[leena_data_factors$CELL == "K562",]

duplicated(leena_data_factors_K562$SAMPLE)

data <- leena_data_factors_K562[!duplicated(leena_data_factors_K562$SAMPLE),]
data <- data[!(is.na(data$chr1_161726992_G_C_G_0)),] # one does not have genotype information, remove them

dosage_chr_1 <- dosage_chr_1[match(data$SAMPLE, dosage_chr_1$sample),]
dosage_chr_18 <- dosage_chr_18[match(data$SAMPLE, dosage_chr_18$sample),]

#-----------------------------------------------------------------------------

# calculate dosages / genotype frequencies

all(dosage_chr_1$chr1_161726992_G_C_G == dosage_chr_1$chr1_161727208_C_G_C) # T

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# GG
sum(dosage_chr_1$chr1_161726992_G_C_G == 2) # 4 --> 4/13 = 0.3076923
# GC
sum(dosage_chr_1$chr1_161726992_G_C_G == 1) # 4 --> 4/13 = 0.3076923
# CC
sum(dosage_chr_1$chr1_161726992_G_C_G == 0) # 5 --> 5/13 = 0.3846154


all(dosage_chr_18$chr18_70203271_C_A_C == dosage_chr_18$chr18_70204107_A_C_A) # T

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

# CC
sum(dosage_chr_18$chr18_70203271_C_A_C == 2) # 0 
# CA
sum(dosage_chr_18$chr18_70203271_C_A_C == 1) # 2 --> 2/13 = 0.1538462
# AA
sum(dosage_chr_18$chr18_70203271_C_A_C == 0) # 11 --> 11/13 = 0.8461538

#-------------------------------------------------------------------------------------------------

# calculated with plink

# a list of individuals to keep (not all have in votro data but have genotypes)
unique <- dosage_chr_18[!(duplicated(dosage_chr_18$IID)),]
unique <- unique[!(is.na(unique$chr18_70203271_C_A_C)),]
unique <- unique[,1:2]
colnames(unique)[1] <- "#FID"
write.table(unique, "./results/frequencies/in_vitro_keep.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# frequencies calculated with ./src/frequencies.sh

data_1 <- fread("./results/frequencies/blooddonors_chr_1_allelefreq.frq")
data_18 <- fread("./results/frequencies/blooddonors_chr_18_allelefreq.frq")
all <- rbind(data_1, data_18)
write.table(all, "./results/frequencies/blooddonors_maf.txt", quote = F, col.names = T, row.names = F, sep = "\t")


data_1 <- fread("./results/frequencies/blooddonors_chr_1_genotypefreq.frqx")
data_18 <- fread("./results/frequencies/blooddonors_chr_18_genotypefreq.frqx")
all <- rbind(data_1, data_18)
write.table(all, "./results/frequencies/blooddonors_genotypefrequencies.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##########################################################################################################################################

# HSCT - MAF

##########################################################################################################################################

# for HSCT data, discovery set

combined <- fread("./results/frequencies/combined_allelefreq.frq")
finland <- fread("./results/frequencies/finland_allelefreq.frq")
uk <- fread("./results/frequencies/uk_allelefreq.frq")
spain <- fread("./results/frequencies/spain_allelefreq.frq")
poland <- fread("./results/frequencies/poland_allelefreq.frq")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "discovery"
finland$set <- "discovery"
uk$set <- "discovery"
spain$set <- "discovery"
poland$set <- "discovery"

all <- rbind(combined, finland, uk, spain, poland)

#--------------------------------------------------------------------------------------------------

# for HSCT data, replication set

combined <- fread("./results/frequencies/combined_allelefreq_replication.frq")
finland <- fread("./results/frequencies/finland_allelefreq_replication.frq")
uk <- fread("./results/frequencies/uk_allelefreq_replication.frq")
spain <- fread("./results/frequencies/spain_allelefreq_replication.frq")
poland <- fread("./results/frequencies/poland_allelefreq_replication.frq")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "replication"
finland$set <- "replication"
uk$set <- "replication"
spain$set <- "replication"
poland$set <- "replication"

all <- rbind(all, combined, finland, uk, spain, poland)

write.table(all, "./results/frequencies/HSCT_maf.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##########################################################################################################################################

# HSCT - genotype frequencies

##########################################################################################################################################

# for HSCT data, discovery set

combined <- fread("./results/frequencies/combined_genotypefreq.frqx")
finland <- fread("./results/frequencies/finland_genotypefreq.frqx")
uk <- fread("./results/frequencies/uk_genotypefreq.frqx")
spain <- fread("./results/frequencies/spain_genotypefreq.frqx")
poland <- fread("./results/frequencies/poland_genotypefreq.frqx")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "discovery"
finland$set <- "discovery"
uk$set <- "discovery"
spain$set <- "discovery"
poland$set <- "discovery"

all <- rbind(combined, finland, uk, spain, poland)

#--------------------------------------------------------------------------------------------------

# for HSCT data, replication set

combined <- fread("./results/frequencies/combined_genotypefreq_replication.frqx")
finland <- fread("./results/frequencies/finland_genotypefreq_replication.frqx")
uk <- fread("./results/frequencies/uk_genotypefreq_replication.frqx")
spain <- fread("./results/frequencies/spain_genotypefreq_replication.frqx")
poland <- fread("./results/frequencies/poland_genotypefreq_replication.frqx")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$SNP %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "replication"
finland$set <- "replication"
uk$set <- "replication"
spain$set <- "replication"
poland$set <- "replication"

all <- rbind(all, combined, finland, uk, spain, poland)

write.table(all, "./results/frequencies/HSCT_genotypefrequencies.txt", quote = F, col.names = T, row.names = F, sep = "\t")

###################################################################################################################################################
###################################################################################################################################################

# calculated with plink2 

###################################################################################################################################################
###################################################################################################################################################

# HSCT - MAF

##########################################################################################################################################

# for HSCT data, discovery set

combined <- fread("./results/frequencies/combined_allelefreq_plink2.afreq")
finland <- fread("./results/frequencies/finland_allelefreq_plink2.afreq")
uk <- fread("./results/frequencies/uk_allelefreq_plink2.afreq")
spain <- fread("./results/frequencies/spain_allelefreq_plink2.afreq")
poland <- fread("./results/frequencies/poland_allelefreq_plink2.afreq")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "discovery"
finland$set <- "discovery"
uk$set <- "discovery"
spain$set <- "discovery"
poland$set <- "discovery"

all <- rbind(combined, finland, uk, spain, poland)

#--------------------------------------------------------------------------------------------------

# for HSCT data, replication set

combined <- fread("./results/frequencies/combined_allelefreq_replication_plink2.afreq")
finland <- fread("./results/frequencies/finland_allelefreq_replication_plink2.afreq")
uk <- fread("./results/frequencies/uk_allelefreq_replication_plink2.afreq")
spain <- fread("./results/frequencies/spain_allelefreq_replication_plink2.afreq")
poland <- fread("./results/frequencies/poland_allelefreq_replication_plink2.afreq")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "replication"
finland$set <- "replication"
uk$set <- "replication"
spain$set <- "replication"
poland$set <- "replication"

all <- rbind(all, combined, finland, uk, spain, poland)

write.table(all, "./results/frequencies/HSCT_maf_plink2.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##########################################################################################################################################

# HSCT - genotype frequencies

##########################################################################################################################################

# for HSCT data, discovery set

combined <- fread("./results/frequencies/combined_genotypefreq_plink2.gcount")
finland <- fread("./results/frequencies/finland_genotypefreq_plink2.gcount")
uk <- fread("./results/frequencies/uk_genotypefreq_plink2.gcount")
spain <- fread("./results/frequencies/spain_genotypefreq_plink2.gcount")
poland <- fread("./results/frequencies/poland_genotypefreq_plink2.gcount")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "discovery"
finland$set <- "discovery"
uk$set <- "discovery"
spain$set <- "discovery"
poland$set <- "discovery"

all <- rbind(combined, finland, uk, spain, poland)

#--------------------------------------------------------------------------------------------------

# for HSCT data, replication set

combined <- fread("./results/frequencies/combined_genotypefreq_replication_plink2.gcount")
finland <- fread("./results/frequencies/finland_genotypefreq_replication_plink2.gcount")
uk <- fread("./results/frequencies/uk_genotypefreq_replication_plink2.gcount")
spain <- fread("./results/frequencies/spain_genotypefreq_replication_plink2.gcount")
poland <- fread("./results/frequencies/poland_genotypefreq_replication_plink2.gcount")

# chr1_161726992_G_C_G
# chr1_161727208_C_G_C

# chr18_70203271_C_A_C
# chr18_70204107_A_C_A

combined <- combined[combined$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
finland <- finland[finland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
uk <- uk[uk$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
spain <- spain[spain$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]
poland <- poland[poland$ID %in% c("chr1_161726992_G_C", "chr1_161727208_C_G", "chr18_70203271_C_A", "chr18_70204107_A_C"),]

# add population
combined$population <- "combined"
finland$population <- "finland"
uk$population <- "uk"
spain$population <- "spain"
poland$population <- "poland"

# add discovery/replication
combined$set <- "replication"
finland$set <- "replication"
uk$set <- "replication"
spain$set <- "replication"
poland$set <- "replication"

all <- rbind(all, combined, finland, uk, spain, poland)

write.table(all, "./results/frequencies/HSCT_genotypefrequencies_plink2.txt", quote = F, col.names = T, row.names = F, sep = "\t")











