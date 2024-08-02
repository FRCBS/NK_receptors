library(data.table)
library(tidyverse)
library(lme4)
library(lmerTest) # uses lmer function from this package to get p-values 
# This function overloads lmer from the lme4-package (lme4::lmer) and adds a couple of slots needed for the computation of Satterthwaite denominator degrees of freedom. All arguments are the same as for lme4::lmer and all the usual lmer-methods work.

#----------------------------------------------------------------------------------------------------------------------------------------

# analyze in vitro NK cytotoxicity data

#----------------------------------------------------------------------------------------------------------------------------------------

# every level of dosage as a separate class

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

# run models with random intercept
# all levels in separate models

intercept_K562_1_0 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161726992_G_C_G_0 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_1_1 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161726992_G_C_G_1 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_1_2 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161726992_G_C_G_2 + (1|SAMPLE), data=leena_data_factors_K562)

intercept_K562_2_0 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161727208_C_G_C_0 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_2_1 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161727208_C_G_C_1 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_2_2 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr1_161727208_C_G_C_2 + (1|SAMPLE), data=leena_data_factors_K562)

intercept_K562_3_0 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr18_70203271_C_A_C_0 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_3_1 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr18_70203271_C_A_C_1 + (1|SAMPLE), data=leena_data_factors_K562)

intercept_K562_4_0 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr18_70204107_A_C_A_0 + (1|SAMPLE), data=leena_data_factors_K562)
intercept_K562_4_1 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + chr18_70204107_A_C_A_1 + (1|SAMPLE), data=leena_data_factors_K562)


get_tables_betas_levels_3 <- function(model1, model2, model3, name){
  
  coeffs_1 <- coef(summary(model1))
  coeffs_2 <- coef(summary(model2))
  coeffs_3 <- coef(summary(model3))
  
  snps <- rbind(coeffs_1[nrow(coeffs_1),], coeffs_2[nrow(coeffs_2),], coeffs_3[nrow(coeffs_3),])
  snps <- as.data.table(snps)
  
  snps$SNP <- paste0(name, c("_0", "_1", "_2"))
  
  #-------------------------------------
  
  # add another way to calculate CI
  # confint the same as confint.merMod in lme4
  
  CI_1 <- confint(model1)
  CI_1 <- CI_1[3:nrow(CI_1),]
  beta_and_CI_1 <- cbind(est=fixef(model1),CI_1)
  # OR_1 <- exp(beta_and_CI_1)            # in OR format
  OR_1 <- beta_and_CI_1                   # in log odds format
  
  CI_2 <- confint(model2)
  CI_2 <- CI_2[3:nrow(CI_2),]
  beta_and_CI_2 <- cbind(est=fixef(model2),CI_2)
  # OR_2 <- exp(beta_and_CI_2)
  OR_2 <- beta_and_CI_2
  
  CI_3 <- confint(model3)
  CI_3 <- CI_3[3:nrow(CI_3),]
  beta_and_CI_3 <- cbind(est=fixef(model3),CI_3)
  # OR_3 <- exp(beta_and_CI_3)
  OR_3 <- beta_and_CI_3
  
  OR_all <- rbind(OR_1[nrow(OR_1),], OR_2[nrow(OR_2),], OR_3[nrow(OR_3),])
  
  snps <- cbind(snps, OR_all)
  
  colnames(snps)[5] <- c("p-value")
  colnames(snps)[7:9] <- c("beta", "lower", "upper")  
  
  #-------------------------------------
  
  return(snps)
  
}

get_tables_betas_levels_2 <- function(model1, model2, name){
  
  coeffs_1 <- coef(summary(model1))
  coeffs_2 <- coef(summary(model2))
  
  snps <- rbind(coeffs_1[nrow(coeffs_1),], coeffs_2[nrow(coeffs_2),])
  snps <- as.data.table(snps)
  
  snps$SNP <- paste0(name, c("_0", "_1"))
  
  #-------------------------------------
  
  # add another way to calculate CI
  # confint the same as confint.merMod in lme4
  
  CI_1 <- confint(model1)
  CI_1 <- CI_1[3:nrow(CI_1),]
  beta_and_CI_1 <- cbind(est=fixef(model1),CI_1)
  # OR_1 <- exp(beta_and_CI_1)            # in OR format
  OR_1 <- beta_and_CI_1                   # in log odds format
  
  CI_2 <- confint(model2)
  CI_2 <- CI_2[3:nrow(CI_2),]
  beta_and_CI_2 <- cbind(est=fixef(model2),CI_2)
  # OR_2 <- exp(beta_and_CI_2)
  OR_2 <- beta_and_CI_2
  
  OR_all <- rbind(OR_1[nrow(OR_1),], OR_2[nrow(OR_2),])
  
  snps <- cbind(snps, OR_all)
  
  colnames(snps)[5] <- c("p-value")
  colnames(snps)[7:9] <- c("beta", "lower", "upper")  
  
  #-------------------------------------
  
  return(snps)
  
}

intercept_K562_betas_1 <- get_tables_betas_levels_3(intercept_K562_1_0, intercept_K562_1_1, intercept_K562_1_2, name = "chr1_161726992_G_C_G")
intercept_K562_betas_2 <- get_tables_betas_levels_3(intercept_K562_2_0, intercept_K562_2_1, intercept_K562_2_2, name = "chr1_161727208_C_G_C")
intercept_K562_betas_3 <- get_tables_betas_levels_2(intercept_K562_3_0, intercept_K562_3_1, name = "chr18_70203271_C_A_C")
intercept_K562_betas_4 <- get_tables_betas_levels_2(intercept_K562_4_0, intercept_K562_4_1, name = "chr18_70204107_A_C_A")
K562 <- rbind(intercept_K562_betas_1, intercept_K562_betas_2, intercept_K562_betas_3, intercept_K562_betas_4)


rsids <- c(paste(rep("rs11585450", 3), c("CC", "GC", "GG")), paste(rep("rs1875763", 3), c("GG", "CG", "CC")), paste(rep("rs8087187", 2), c("AA", "CA")), paste(rep("rs3911730", 2), c("CC", "AC")))

K562$rs_IDs <- rsids

colnames(K562)[5] <- "p_value"

colnames(K562) <- c("Estimate", "std_error", "df", "t value", "p_value", "SNP", "beta", "lower", "upper", "Rsid")

write.table(K562, file="./results/NK_killing_activity_results/results_K562.txt", sep="\t", quote=F, row.names=F, col.names=T)

