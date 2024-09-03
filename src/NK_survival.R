library(data.table)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(flextable)
library(cowplot)

# sudo apt-get update
# sudo apt-get install -y r-cran-gtsummary

# sudo apt-get update
# sudo apt-get install -y r-cran-ggsurvfit


# survival analysis for the four NK receptor SNPs

########################################################################################################################################

# add SNP dosage to survival data

########################################################################################################################################

setwd("/home/nihtiju/work/NK_receptors/")

# read in data
data <- fread("/home/nihtiju/work/HSCT_predictor/results/patient_information/survival/HSCT_survival_data.txt")
# read in SNP dosage 
dosage <- fread("/home/nihtiju/work/NK_receptors/results/survival/plink1_dosage.raw") # plink1 dosage, has values 2, 1, 0
# # pgen dosage
# dosage <- fread("/home/nihtiju/work/NK_receptors/results/extract_range/all_datasets_plink_pgen_common_variants_filtered_ref_changed_newSNPsLDfiltered_dosage.raw")


dosage_snps <- colnames(dosage)
dosage_snps <- str_sub(dosage_snps, 1, -3)

# filter dosage to only contain the wanted SNPs
SNPs <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_chosen.txt")
# only leave in relapse SNPs
SNPs <- SNPs[SNPs$pheno == "relapse",]

cols <- c(T,T, dosage_snps[3:length(dosage_snps)] %in% SNPs$ID)
dosage <- dosage[,..cols]

dosage <- dosage[match(data$DonorGenotypingID, dosage$IID),]

# train and test split
# train <- fread("/home/nihtiju/work/NK_receptors/results/ID_lists/train_sample_ids_donors.txt", header = F)
# test <- fread("/home/nihtiju/work/NK_receptors/results/ID_lists/test_sample_ids_donors.txt", header = F)
# only RFS -> relapse IDs:
train <- fread("./results/ID_lists/donors_train_relapse.txt", header = F)
test <- fread("./results/ID_lists/donors_test_relapse.txt", header = F)

# add time in years (in stead of days)
data$overall_survival_years <- data$overall_survival / 365.25
data$relapse_free_survival_years <- data$relapse_free_survival / 365.25

########################################################################################################################################

# univariate analysis 

########################################################################################################################################

survival_univariate <- function(survival_data, SNP_name, tested_set_ids, tested_set){
  
  SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  rsID <- SNP_info$rsID[SNP_info$ID == str_sub(SNP_name, 1, -3)]
  
  col <- colnames(dosage) %like% SNP_name
  
  split <- str_split(SNP_name, "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  survival_data$dosage <- dosage[, ..col]
  survival_data$dosage[survival_data$dosage == 2] <- paste0("Donor ", counted, counted)
  survival_data$dosage[survival_data$dosage == 1] <- paste0("Donor ", counted, other)
  survival_data$dosage[survival_data$dosage == 0] <- paste0("Donor ", other, other)
  survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted), paste0("Donor ", counted, other), paste0("Donor ", other, other)))
  
  RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ dosage, data = survival_data[survival_data$DonorGenotypingID %in% tested_set_ids$V2,]) %>% 
    ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
    labs(
      x = "Years after transplantation",
      y = "Relapse-free survival"
    ) + 
    add_confidence_interval() +
    scale_color_grey(start = 0, end = 0.4) +
    scale_fill_grey(start = 0, end = 0.4) +
    add_censor_mark(size = 5, shape=124) +
    add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
    add_pvalue("annotation", size = 7) +
    labs(title = paste0("Relapse-free survival ", rsID)) +
    theme_classic() +
    theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
  
  
  ggsave(file = paste0("/home/nihtiju/work/NK_receptors/results/survival/RFS_", SNP_name, "_", tested_set, "_univariate.png"), print(RFS_plot), dpi = 600, width = 10, height = 11)
  
}

# for all SNPs in a loop
for (i in 3:ncol(dosage)) {
  
  survival_univariate(data, colnames(dosage)[i], train, "train")
  survival_univariate(data, colnames(dosage)[i], test, "test")
}

#--------------------------------------------------------------------------------------------------------------------------

# univariate analysis
# merge 2 smallest genotype groups into one


survival_univariate_2groups <- function(survival_data, SNP_name, tested_set_ids, tested_set){
  
  SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  rsID <- SNP_info$rsID[SNP_info$ID == str_sub(SNP_name, 1, -3)]
  
  col <- colnames(dosage) %like% SNP_name
  
  split <- str_split(SNP_name, "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  survival_data$dosage <- dosage[, ..col]
  
  if(sum(survival_data$dosage == 2, na.rm = T) < sum(survival_data$dosage == 0, na.rm = T)){
    
    survival_data$dosage[survival_data$dosage == 2 | survival_data$dosage == 1] <- paste0("Donor ", counted, counted, "/", counted, other)
    survival_data$dosage[survival_data$dosage == 0] <- paste0("Donor ", other, other)
    
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted, "/", counted, other), paste0("Donor ", other, other)))
    
    
  } else if (sum(survival_data$dosage == 2, na.rm = T) > sum(survival_data$dosage == 0, na.rm = T)) {
    
    survival_data$dosage[survival_data$dosage == 2] <- paste0("Donor ", counted, counted)
    survival_data$dosage[survival_data$dosage == 0 | survival_data$dosage == 1] <- paste0("Donor ", counted, other, " / ", other, other)
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted), paste0("Donor ", counted, other, " / ", other, other)))
    
  }
  
  # plot data
  RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ dosage, data = survival_data[survival_data$DonorGenotypingID %in% tested_set_ids$V2,]) %>% 
    ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
    labs(
      x = "Years after transplantation",
      y = "Relapse-free survival"
    ) + 
    add_confidence_interval() +
    scale_color_grey(start = 0, end = 0.4) +
    scale_fill_grey(start = 0, end = 0.4) +
    add_censor_mark(size = 5, shape=124) +
    add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
    add_pvalue("annotation", size = 7) +
    labs(title = paste0("Relapse-free survival ", rsID)) +
    theme_classic() +
    theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
  
  
  ggsave(file = paste0("/home/nihtiju/work/NK_receptors/results/survival/RFS_", SNP_name, "_", tested_set, "_univariate_2groups.png"), print(RFS_plot), dpi = 600, width = 10, height = 11)
  
}

# for all SNPs in a loop
for (i in 3:ncol(dosage)) {
  
  survival_univariate_2groups(data, colnames(dosage)[i], train, "train")
  survival_univariate_2groups(data, colnames(dosage)[i], test, "test")
}

#--------------------------------------------------------------------------------------------------------------------------

# univariate analysis
# merge 2 smallest genotype groups into one
# SAVE PLOTS INTO A LIST AND MERGE INTO ONE PLOT


survival_univariate_2groups <- function(survival_data, SNP_name, tested_set_ids, tested_set){
  
  SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  rsID <- SNP_info$rsID[SNP_info$ID == str_sub(SNP_name, 1, -3)]
  
  col <- colnames(dosage) %like% SNP_name
  
  split <- str_split(SNP_name, "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  survival_data$dosage <- dosage[, ..col]
  
  if(sum(survival_data$dosage == 2, na.rm = T) < sum(survival_data$dosage == 0, na.rm = T)){
    
    survival_data$dosage[survival_data$dosage == 2 | survival_data$dosage == 1] <- paste0("Donor ", counted, counted, "/", counted, other)
    survival_data$dosage[survival_data$dosage == 0] <- paste0("Donor ", other, other)
    
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted, "/", counted, other), paste0("Donor ", other, other)))
    
    
  } else if (sum(survival_data$dosage == 2, na.rm = T) > sum(survival_data$dosage == 0, na.rm = T)) {
    
    survival_data$dosage[survival_data$dosage == 2] <- paste0("Donor ", counted, counted)
    survival_data$dosage[survival_data$dosage == 0 | survival_data$dosage == 1] <- paste0("Donor ", counted, other, " / ", other, other)
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted), paste0("Donor ", counted, other, " / ", other, other)))
    
  }
  
  # plot data
  RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ dosage, data = survival_data[survival_data$DonorGenotypingID %in% tested_set_ids$V2,]) %>% 
    ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
    labs(
      x = "Years after transplantation",
      y = "Relapse-free survival"
    ) + 
    add_confidence_interval() +
    scale_color_grey(start = 0, end = 0.4) +
    scale_fill_grey(start = 0, end = 0.4) +
    add_censor_mark(size = 5, shape=124) +
    add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
    add_pvalue("annotation", size = 7) +
    labs(title = paste0("Relapse-free survival ", rsID)) +
    theme_classic() +
    theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
  
  
  ggsave(file = paste0("/home/nihtiju/work/NK_receptors/results/survival/RFS_", SNP_name, "_", tested_set, "_univariate_2groups.png"), print(RFS_plot), dpi = 600, width = 10, height = 11)
  
  return(RFS_plot)
  
}

# for all SNPs in a loop
plot_list <- list()
for (i in 3:ncol(dosage)) {
  # make individual plots and save the train plots into a list
  plot_list[[i]] <-  survival_univariate_2groups(data, colnames(dosage)[i], train, "train")
  # make individual plots but do not save the test set plots into the list
  survival_univariate_2groups(data, colnames(dosage)[i], test, "test")
}

# individual numbers not in th eplots with this
# plot_grid(plot_list[[3]], plot_list[[4]], plot_list[[5]], ncol = 3)
# ggsave("./results/for_paper/survival/RFS_univariate.png", bg = "white", width = 20, height = 20, dpi = 600)

# this way the individual numbers are in the plots
built_plot_1 <- ggsurvfit_build(plot_list[[3]])
built_plot_2 <- ggsurvfit_build(plot_list[[4]])
built_plot_3 <- ggsurvfit_build(plot_list[[5]])

combined <- cowplot::plot_grid(built_plot_1, built_plot_2, built_plot_3, ncol = 3, labels = c("A", "B", "C"), label_size = 20)

ggsave(file = "./results/for_paper/survival/RFS_univariate.png", dpi = 600, width = 25, height = 10)

###########################################################################################################################

# multivariate analysis

###########################################################################################################################

# multivariate analysis
# merge 2 smallest genotype groups into one

# read in covariates
covars_binary <- fread("/home/nihtiju/work/NK_receptors/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")
# modify covariates - combine all levels into one column
covars <- covars_binary[,1:2]

covars$donortype <- NA
covars$donortype[covars_binary$donor_sibling == 1] <- "Sibling"
covars$donortype[covars_binary$donor_haplo == 1] <- "Haplo"
covars$donortype[covars_binary$donor_sibling == 0 & covars_binary$donor_haplo == 0] <- "Register"
covars$donortype <- as.factor(covars$donortype)

covars$graft <- NA
covars$graft[covars_binary$graft_PB == 1] <- "PB"
covars$graft[covars_binary$graft_PB_BM == 1] <- "Missing"
covars$graft[covars_binary$graft_NA == 1] <- "Missing"
covars$graft[covars_binary$graft_PB == 0 & covars_binary$graft_PB_BM == 0 & covars_binary$graft_NA == 0] <- "BM"
covars$graft <- as.factor(covars$graft)

covars$sex_combo <- NA
covars$sex_combo[covars_binary$sex_combo_risk == 1] <- "Donor F, recipient M"
covars$sex_combo[covars_binary$sex_combo_NA == 1] <- "Missing"
covars$sex_combo[covars_binary$sex_combo_risk == 0 & covars_binary$sex_combo_NA == 0] <- "Other combinations"
covars$sex_combo <- as.factor(covars$sex_combo)

covars$conditioning <- NA
covars$conditioning[covars_binary$preconditioning_MAC == 1] <- "MAC"
covars$conditioning[covars_binary$preconditioning_sekv == 1] <- "Other"
covars$conditioning[covars_binary$preconditioning_NMA == 1] <- "Other"
covars$conditioning[covars_binary$preconditioning_MAC == 0 & covars_binary$preconditioning_sekv == 0 & covars_binary$preconditioning_NMA == 0] <- "RIC"
covars$conditioning <- as.factor(covars$conditioning)

covars$population <- NA
covars$population[covars_binary$population_fin == 1] <- "Finland"
covars$population[covars_binary$population_spain == 1] <- "Spain"
covars$population[covars_binary$population_poland == 1] <- "Poland"
covars$population[covars_binary$population_fin == 0 & covars_binary$population_spain == 0 & covars_binary$population_poland == 0] <- "UK"
covars$population <- as.factor(covars$population)

covars$CMV <- NA
covars$CMV[covars_binary$CMV_risk == 1] <- "Donor neg, recipient pos"
covars$CMV[covars_binary$CMV_missing == 1] <- "Missing"
covars$CMV[covars_binary$CMV_risk == 0 & covars_binary$CMV_missing == 0] <- "Other combinations"
covars$CMV <- as.factor(covars$CMV)

covars$hla_c1 <- covars_binary$HLA_C1
covars$hla_c2 <- covars_binary$HLA_C2

covars$GvHD_profylaxis <- covars_binary$GvHD_prevention
covars$GvHD_profylaxis <- as.factor(covars$GvHD_profylaxis)

covars$age_r <- covars_binary$AgeAll
covars$age_d <- covars_binary$DonorAgeAll

covars$aml_not <- covars_binary$AML_not

covars$tx_year <- covars_binary$txyear

covars$hla_match <- covars_binary$match_score

# merge covars with the rest of the data
covars <- covars[match(data$DonorGenotypingID, covars$IID),]
data <- cbind(data, covars[,3:ncol(covars)])

# survival with covariates
survival_multivariate_2groups <- function(survival_data, SNP_name, tested_set_ids, tested_set){
  
  SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  rsID <- SNP_info$rsID[SNP_info$ID == str_sub(SNP_name, 1, -3)]
  
  col <- colnames(dosage) %like% SNP_name
  
  split <- str_split(SNP_name, "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  survival_data$dosage <- dosage[, ..col]
  
  # merge heterozygotes into the smaller homozygote group
  if(sum(survival_data$dosage == 2, na.rm = T) < sum(survival_data$dosage == 0, na.rm = T)){
    
    survival_data$dosage[survival_data$dosage == 2 | survival_data$dosage == 1] <- paste0("Donor ", counted, counted, "/", counted, other)
    survival_data$dosage[survival_data$dosage == 0] <- paste0("Donor ", other, other)
    
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted, "/", counted, other), paste0("Donor ", other, other)))
    
    
  } else if (sum(survival_data$dosage == 2, na.rm = T) > sum(survival_data$dosage == 0, na.rm = T)) {
    
    survival_data$dosage[survival_data$dosage == 2] <- paste0("Donor ", counted, counted)
    survival_data$dosage[survival_data$dosage == 0 | survival_data$dosage == 1] <- paste0("Donor ", counted, other, " / ", other, other)
    survival_data$dosage <- factor(survival_data$dosage, levels=c(paste0("Donor ", counted, counted), paste0("Donor ", counted, other, " / ", other, other)))
    
  }
  
  # survival
  
  
  res_RFS <- coxph(Surv(relapse_free_survival_years, RFS_status) ~ ., data = survival_data[survival_data$DonorGenotypingID %in% tested_set_ids$V2,-c(1:4,6,7)]) %>% 
    tbl_regression(exp = TRUE)  %>%
    as_flex_table() 
  res_RFS <- bg(res_RFS, bg = "white", part = "all")
  
  save_as_image(res_RFS, path = paste0("/home/nihtiju/work/NK_receptors/results/survival/RFS_", SNP_name, "_", tested_set, "_multivariate_2groups.png"))
  
}

# for all SNPs in a loop
for (i in 3:ncol(dosage)) {
  
  survival_multivariate_2groups(data, colnames(dosage)[i], train, "train")
  survival_multivariate_2groups(data, colnames(dosage)[i], test, "test")
  
}














