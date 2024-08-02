library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)

setwd("/home/nihtiju/work/NK_receptors/")

# --------------------------------------------------------------------------------------------------------------------------

# for all datasets

# read in dosages
dosage <- fread("./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage.raw")


cols <- colnames(dosage)[3:ncol(dosage)]
cols_keep <- cols %in% c("chr1_161726992_G_C_C", "chr1_161727208_C_G_G", "chr18_70203271_C_A_A", "chr18_70204107_A_C_C")
cols_keep <- c(T,T,cols_keep)
dosage <- dosage[, ..cols_keep]

plot_1 <- ggplot(dosage, aes(x=chr1_161726992_G_C_C)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("All datasets, chr1_161726992_G_C_C") +
  xlab("dosage")

plot_2 <- ggplot(dosage, aes(x=chr1_161727208_C_G_G)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("All datasets, chr1_161727208_C_G_G") +
  xlab("dosage")

plot_3 <- ggplot(dosage, aes(x=chr18_70203271_C_A_A)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("All datasets, chr18_70203271_C_A_A") +
  xlab("dosage")

plot_4 <- ggplot(dosage, aes(x=chr18_70204107_A_C_C)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("All datasets, chr18_70204107_A_C_C") +
  xlab("dosage")

# --------------------------------------------------------------------------------------------------------------------------

# for katalonia

# read in dosages
dosage <- fread("./results/pgen_NK_SNPs/katalonia_dosage.raw")


cols <- colnames(dosage)[3:ncol(dosage)]
cols_keep <- cols %in% c("chr1_161726992_G_C_C", "chr1_161727208_C_G_G", "chr18_70203271_C_A_A", "chr18_70204107_A_C_C")
cols_keep <- c(T,T,cols_keep)
dosage <- dosage[, ..cols_keep]

plot_5 <- ggplot(dosage, aes(x=chr1_161726992_G_C_C)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("Spain, chr1_161726992_G_C_C") +
  xlab("dosage")

plot_6 <- ggplot(dosage, aes(x=chr1_161727208_C_G_G)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("Spain, chr1_161727208_C_G_G") +
  xlab("dosage")

plot_7 <- ggplot(dosage, aes(x=chr18_70203271_C_A_A)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("Spain, chr18_70203271_C_A_A") +
  xlab("dosage")

plot_8 <- ggplot(dosage, aes(x=chr18_70204107_A_C_C)) + geom_histogram(boundary = 0, binwidth=0.1) +
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.1)) + 
  theme_minimal() +
  ggtitle("Spain, chr18_70204107_A_C_C") +
  xlab("dosage")

#----------------------------------------------------------------------------------------------------

# combine all plots

p_all <- (plot_1 | plot_2 | plot_3 | plot_4) / (plot_5 | plot_6 | plot_7 | plot_8)
ggsave("./results/frequencies/dosage_histograms.png", bg = "white", width = 30, height = 14, dpi = 600)

#######################################################################################################################

# calculate genotype frequencies manually

# set thresholds to divide individuals into three groups
# use the plot ./results/frequencies/dosage_histograms.png to do this

# 0 - 0.65 homozygous
# 0.65 - 1.3 heterozygous
# 1.3 - 2 homozygous

# read in dosages
dosage <- fread("./results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage.raw")

cols <- colnames(dosage)[3:ncol(dosage)]
cols_keep <- cols %in% c("chr1_161726992_G_C_C", "chr1_161727208_C_G_G", "chr18_70203271_C_A_A", "chr18_70204107_A_C_C")
cols_keep <- c(T,T,cols_keep)
dosage <- dosage[, ..cols_keep]

results <- tibble("subset" = character(),
                      "chr1_161726992_G_C_C_CC" = numeric(), 
                      "chr1_161726992_G_C_C_CG" = numeric(), 
                      "chr1_161726992_G_C_C_GG" = numeric(),
                      "chr1_161727208_C_G_G_GG" = numeric(),
                      "chr1_161727208_C_G_G_GC" = numeric(),
                      "chr1_161727208_C_G_G_CC" = numeric(),
                      "chr18_70203271_C_A_A_AA" = numeric(),
                      "chr18_70203271_C_A_A_AC" = numeric(),
                      "chr18_70203271_C_A_A_CC" = numeric(),
                      "chr18_70204107_A_C_C_CC" = numeric(),
                      "chr18_70204107_A_C_C_CA" = numeric(),
                      "chr18_70204107_A_C_C_AA" = numeric())
# rbind does not work for data.tables, dataframes retired -> use tibble

#-------------------------------------------------------------------------

# filter the correct population

# combination -> no need to filter

# filter train and test sets
train_IDs <- fread("./results/ID_lists/train_sample_ids_donors.txt", header = F)
test_IDs <- fread("./results/ID_lists/test_sample_ids_donors.txt", header = F)

train <- dosage[dosage$IID %in% train_IDs$V2,]
test <- dosage[dosage$IID %in% test_IDs$V2,]


calculate_one_row <- function(subset, table){
  
  return(c(subset, 
           sum(table$chr1_161726992_G_C_C >= 1.3),
           sum(table$chr1_161726992_G_C_C < 1.3 & table$chr1_161726992_G_C_C >= 0.65),
           sum(table$chr1_161726992_G_C_C < 0.65),
           sum(table$chr1_161727208_C_G_G >= 1.3),
           sum(table$chr1_161727208_C_G_G < 1.3 & table$chr1_161727208_C_G_G >= 0.65),
           sum(table$chr1_161727208_C_G_G < 0.65),
           sum(table$chr18_70203271_C_A_A >= 1.3),
           sum(table$chr18_70203271_C_A_A < 1.3 & table$chr18_70203271_C_A_A >= 0.65),
           sum(table$chr18_70203271_C_A_A < 0.65),
           sum(table$chr18_70204107_A_C_C >= 1.3),
           sum(table$chr18_70204107_A_C_C < 1.3 & table$chr18_70204107_A_C_C >= 0.65),
           sum(table$chr18_70204107_A_C_C < 0.65)
  ))
  
}
    
results <- rbind(results, calculate_one_row("Combination, discovery", train))
results <- rbind(results, calculate_one_row("Combination, validation", test))

# add results by population
fin_IDs <- fread("./results/pgen_NK_SNPs/finns_plink_pgen_common_variants_filtered_ref_changed.psam", header = T)
spain_IDs <- fread("./results/pgen_NK_SNPs/katalonia_plink_pgen_common_variants_filtered_ref_changed.psam", header = T)
uk_IDs <- fread("./results/pgen_NK_SNPs/newcastle_plink_pgen_common_variants_filtered_ref_changed.psam", header = T)
poland_IDs <- fread("./results/pgen_NK_SNPs/poland_plink_pgen_common_variants_filtered_ref_changed.psam", header = T)

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

colnames(results) <- c("subset", "chr1_161726992_G_C_C_CC", "chr1_161726992_G_C_C_CG", "chr1_161726992_G_C_C_GG", "chr1_161727208_C_G_G_GG", "chr1_161727208_C_G_G_GC", "chr1_161727208_C_G_G_CC", "chr18_70203271_C_A_A_AA", "chr18_70203271_C_A_A_AC", "chr18_70203271_C_A_A_CC", "chr18_70204107_A_C_C_CC", "chr18_70204107_A_C_C_CA", "chr18_70204107_A_C_C_AA")

write.table(results, "./results/frequencies/HSCT_genotypefrequencies_manual.txt", quote = F, col.names = T, row.names = F, sep = "\t")




















