library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)

# sudo apt-get update
# sudo apt-get install -y r-cran-cowplot

setwd("/home/nihtiju/work/NK_receptors/")

#############################################################################################################################

# allele freq / MAF for HSCT
# calculated with plink

#############################################################################################################################

files_train <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies", pattern = "allelefreq_plink2_train.afreq", full.names = T)
files_test <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies", pattern = "allelefreq_plink2_test.afreq", full.names = T)

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[!duplicated(SNP_info$ID),]


all <- data.table("#CHROM" = numeric(), ID = character(), REF = character(), ALT = character(), ALT_FREQS = numeric(), OBS_CT = numeric(), population = character(), set = character())

for (i in 1:length(files_train)) {
  
  # read in file
  train <- fread(files_train[i])
  test <- fread(files_test[i])
  
  # filter wanted SNPs
  train <- train[train$ID %in% SNP_info$ID,]
  test <- test[test$ID %in% SNP_info$ID,]
  
  # add population
  if(files_train[i] %like% "all_datasets"){
    train$population <- "combined"
    test$population <- "combined"
  } else if (files_train[i] %like% "finns") {
    train$population <- "finland"
    test$population <- "finland"
  } else if (files_train[i] %like% "katalonia") {
    train$population <- "spain"
    test$population <- "spain"
  } else if (files_train[i] %like% "newcastle") {
    train$population <- "uk"
    test$population <- "uk"
  } else if (files_train[i] %like% "poland") {
    train$population <- "poland"
    test$population <- "poland"
  }
  
  # add dataset
  train$set <- "discovery"
  test$set <- "replication"
  
  all <- rbind(all, train, test)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/frequencies/HSCT_maf_plink2.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#############################################################################################################################

# allele freq / MAF for blood donors
# calculated with plink

#############################################################################################################################

files <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies/", pattern = "_allelefreq.frq", full.names = T)

all <- data.table(CHR = numeric(), SNP = character(), A1 = character(), A2 = character(), MAF = numeric(), NCHROBS = numeric())

for (i in 1:length(files)) {
 
  file <- fread(files[i])
  file <- file[file$SNP %in% SNP_info$ID,]
   
  all <- rbind(all, file)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/frequencies/blood_donors_maf.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#############################################################################################################################

# genotype freq for HSCT
# calculated with plink

#############################################################################################################################

files_train <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies", pattern = "_genotypefreq_plink2_train.gcount", full.names = T)
files_test <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies", pattern = "_genotypefreq_plink2_test.gcount", full.names = T)

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[!duplicated(SNP_info$ID),]


all <- data.table("#CHROM" = numeric(), ID = character(), REF = character(), ALT = character(), 
                  HOM_REF_CT = character(), HET_REF_ALT_CTS = character(), TWO_ALT_GENO_CTS = character(), 
                  HAP_REF_CT = character(), HAP_ALT_CTS = character(), MISSING_CT = character(),
                  population = character(), set = character())

for (i in 1:length(files_train)) {
  
  # read in file
  train <- fread(files_train[i])
  test <- fread(files_test[i])
  
  # filter wanted SNPs
  train <- train[train$ID %in% SNP_info$ID,]
  test <- test[test$ID %in% SNP_info$ID,]
  
  # add population
  if(files_train[i] %like% "all_datasets"){
    train$population <- "combined"
    test$population <- "combined"
  } else if (files_train[i] %like% "finns") {
    train$population <- "finland"
    test$population <- "finland"
  } else if (files_train[i] %like% "katalonia") {
    train$population <- "spain"
    test$population <- "spain"
  } else if (files_train[i] %like% "newcastle") {
    train$population <- "uk"
    test$population <- "uk"
  } else if (files_train[i] %like% "poland") {
    train$population <- "poland"
    test$population <- "poland"
  }
  
  # add dataset
  train$set <- "discovery"
  test$set <- "replication"
  
  all <- rbind(all, train, test)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/frequencies/HSCT_genotypefreq_plink2.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#############################################################################################################################

# genotype freq for HSCT
# calculated by hand from dosages 

#############################################################################################################################

# plot the distributions

#--------------------------------------------------------

# read in dosage
dosage <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered_dosage.raw")

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[!duplicated(SNP_info$ID),]

cols <- colnames(dosage)[7:ncol(dosage)]
cols <- str_sub(cols, 1, -3)
keep <- c(T,T,F,F,F,F,c(cols %in% SNP_info$ID))
dosage <- dosage[, ..keep]

train_IDs <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
test_IDs <- fread("./results/ID_lists/test_sample_ids_donors.txt", header = F)

# filter with donor IDs
all <- rbind(train_IDs, test_IDs)

dosage <- dosage[dosage$IID %in% all$V2,]

# make plots
plots <- list()
for (i in 3:ncol(dosage)) {
  
  data <- dosage[, ..i]
  colnames(data)[1] <- "snp"
  
  snp <- colnames(dosage)[i]
  
  p <- ggplot(data, aes(x = snp)) + geom_histogram(boundary = 0, binwidth = 0.1) +
    scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
    theme_minimal() +
    ggtitle(paste0("All datasets, ", snp)) +
    xlab("dosage")
  
  plots[[i-2]] <- p
  
}

p_all <- (plots[[1]] | plots[[2]] | plots[[3]]) / (plots[[4]] | plots[[5]] | plots[[6]]) / (plots[[7]] | plots[[8]] | plots[[9]]) / (plots[[10]] | plots[[11]])
ggsave("/home/nihtiju/work/NK_receptors/results/frequencies/dosage_histograms.png", bg = "white", width = 30, height = 14, dpi = 600)

#--------------------------------------------------------

# calculate genotype frequencies manually

# set thresholds to divide individuals into three groups
# use the plot /home/nihtiju/work/NK_receptors/results/frequencies/dosage_histograms.png to do this

# 0 - 0.65 homozygous
# 0.65 - 1.3 heterozygous
# 1.3 - 2 homozygous

# make results table
colnames_snps <- c()
for (i in 3:ncol(dosage)) {
  
  name <- colnames(dosage)[i]
  split <- str_split(name, "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  colnames_snps <- c(colnames_snps, paste0(name, "_", counted, counted), paste0(name, "_", counted, other), paste0(name, "_", other, other))
  
}

results <- data.frame(matrix(ncol = 1 + length(colnames_snps), nrow = 0))
colnames(results) <- c("subset", colnames_snps)
results <- as_tibble(results)

# go over train and test sets dosage
train <- dosage[dosage$IID %in% train_IDs$V2,]
test <- dosage[dosage$IID %in% test_IDs$V2,]

calculate_one_row <- function(subset, table){
  
  to_return <- c(subset)
  for (i in 3:ncol(table)) {
    to_return <- c(to_return, sum(table[, ..i] >= 1.3),
                   sum(table[, ..i] < 1.3 & table[, ..i] >= 0.65),
                   sum(table[, ..i] < 0.65))
  }
  
  return(to_return)
  
}

train_numbers <- calculate_one_row("Combination, discovery", train)
test_numbers <- calculate_one_row("Combination, validation", test)


# add results by population
fin_IDs <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/finns_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered.psam", header = T)
spain_IDs <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/katalonia_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered.psam", header = T)
uk_IDs <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/newcastle_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered.psam", header = T)
poland_IDs <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/poland_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered.psam", header = T)

filter_pop <- function(IDs){
  
  train_pop <- train[train$IID %in% IDs$IID,]
  test_pop <- test[test$IID %in% IDs$IID,]
  
  return(list(train = train_pop, test = test_pop))
  
}

finland <- filter_pop(fin_IDs)
spain <- filter_pop(spain_IDs)
uk <- filter_pop(uk_IDs)
poland <- filter_pop(poland_IDs)

results <- rbind(results, calculate_one_row("Finland, discovery", finland[["train"]]))
results <- rbind(results, calculate_one_row("Finland, validation", finland[["test"]]))
results <- rbind(results, calculate_one_row("UK, discovery", uk[["train"]]))
results <- rbind(results, calculate_one_row("UK, validation", uk[["test"]]))
results <- rbind(results, calculate_one_row("Spain, discovery", spain[["train"]]))
results <- rbind(results, calculate_one_row("Spain, validation", spain[["test"]]))
results <- rbind(results, calculate_one_row("Poland, discovery", poland[["train"]]))
results <- rbind(results, calculate_one_row("Poland, validation", poland[["test"]]))

colnames(results) <- c("subset", colnames_snps)

write.table(results, "/home/nihtiju/work/NK_receptors/results/frequencies/HSCT_genotypefrequencies_manual.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#############################################################################################################################

# genotype freq for blood donors
# calculated with plink

#############################################################################################################################

files <- list.files("/home/nihtiju/work/NK_receptors/results/frequencies/", pattern = "genotypefreq.frqx", full.names = T)

all <- data.table(CHR = numeric(), SNP = character(), A1 = character(), A2 = character(), "C(HOM A1)" = numeric(), "C(HET)" = numeric(), "C(HOM A2)" = numeric(), "C(HAP A1)" = numeric(), "C(HAP A2)" = numeric(), "C(MISSING)" = numeric())

for (i in 1:length(files)) {
  
  file <- fread(files[i])
  file <- file[file$SNP %in% SNP_info$ID,]
  
  all <- rbind(all, file)
  
}

write.table(all, "/home/nihtiju/work/NK_receptors/results/frequencies/blood_donors_genotype_freq.txt", quote = F, col.names = T, row.names = F, sep = "\t")

