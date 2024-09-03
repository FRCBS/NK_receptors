library(data.table)
library(tidyverse)
library(lme4)
library(lmerTest) # uses lmer function from this package to get p-values 
# This function overloads lmer from the lme4-package (lme4::lmer) and adds a couple of slots needed for the computation of Satterthwaite denominator degrees of freedom. All arguments are the same as for lme4::lmer and all the usual lmer-methods work.

#-------------------------------------------------------------------------------------------------------------------------------------

# analyze in vitro NK cytotoxicity data

#-------------------------------------------------------------------------------------------------------------------------------------

# new cleaned code 
# with a loop

setwd("/home/nihtiju/work/NK_receptors/")

# every level of dosage as a separate class

# load data in
cytotoxicity_data_factors <- fread("./data/NK_assay_processed/NK_data_parsed.tsv", stringsAsFactors = T) # the results are the same if factors or not

# function for processing results from lme models
get_results <- function(model, dosage_num, SNP_col){
  
  # get result values
  coeffs <- coef(summary(model))
  coeffs <- cbind(variable = rownames(coeffs), coeffs)
  coeffs <- as.data.table(coeffs)
  
  cols <- SNP_col # an error without this
  coeffs$SNP <- paste0(colnames(dosage[, ..cols]), "_", dosage_num)
  
  print(dosage_num)

  # confint the same as confint.merMod in lme4
  CI <- confint(model)
  CI <- CI[3:nrow(CI),] # exclude .sig01 and .sigma
  beta_and_CI <- cbind(est=fixef(model),CI)
  OR <- exp(beta_and_CI)                                    # in OR format
  colnames(OR) <- c("OR", "OR_95CI_L", "OR_95CI_U")
  colnames(beta_and_CI) <- c("beta", "beta_95CI_L", "beta_95CI_U")  # in log odds format

  results <- cbind(coeffs, beta_and_CI, OR)
  results <- as.data.table(results)

  return(results)
  
}

results_0 <- get_results(intercept_K562_0, "0", j)
results_0 <- get_results(intercept_K562_0, "0", col)

# a table to save all results
all <- data.table("variable"=character(),
                  "Estimate"=numeric(),
                  "Std. Error"=numeric(),
                  "df"=numeric(),
                  "t value"=numeric(),
                  "Pr(>|t|)"=numeric(),
                  "SNP"=character(),
                  "beta"=numeric(),
                  "beta_95CI_L"=numeric(),
                  "beta_95CI_U"=numeric(),
                  "OR"=numeric(),
                  "OR_95CI_L"=numeric(),
                  "OR_95CI_U"=numeric())

# go over chrs with SNPs wanted
for (i in c(1,6,7,12,18)) {
  
  print(paste0("Now starting: CHR ", i))
  
  # add dosage
  dosage <- fread(paste0("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr",i,"_extracted_dosage.raw"))
  
  # unify sample names
  
  # # pseudonyms in genotype data do not match in vitro data's IDs 
  # first change the IDs according to an ID file, then modify IDs to match the in vitro data's IDs
  ids <- fread("/home/nihtiju/work/NK_receptors/data/blood_donors/blood_donors_IDs.txt", header = F)
  ids <- ids[match(dosage$IID, ids$V1),]
  dosage$FID <- ids[,2]
  dosage$IID <- ids[,2]
  
  # this for when the IDs are correct & when fixing the code
  dosage$sample <- paste0(str_sub(dosage$IID, 1, 2), str_sub(dosage$IID, 6, 7))
  dosage$sample <- str_replace_all(dosage$sample, "NK0", "NK")
  # order dosages to be the same as other data
  dosage <- dosage[match(cytotoxicity_data_factors$SAMPLE, dosage$sample),]
  
  # get cell lines separately
  # using only K562 here
  cytotoxicity_data_factors_K562 <- cytotoxicity_data_factors[cytotoxicity_data_factors$CELL == "K562",]
  dosage <- dosage[cytotoxicity_data_factors$CELL == "K562",] # this filtered too to have matching rows
  
  # go over the SNP dosage columns
  for (j in 7: (ncol(dosage) - 1)) {
    
    print(paste0("Column: ", j))
    
    cytotoxicity_data_factors_K562$dosage_0 <- as.factor(as.numeric(dosage[, ..j] == 0))
    cytotoxicity_data_factors_K562$dosage_1 <- as.factor(as.numeric(dosage[, ..j] == 1))
    cytotoxicity_data_factors_K562$dosage_2 <- as.factor(as.numeric(dosage[, ..j] == 2))
    
    # run models with random intercept
    # all levels in separate models
    # some do not have all dosages 0/1/2
    if(length(levels(cytotoxicity_data_factors_K562$dosage_0)) == 2) {
      # run model
      intercept_K562_0 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + dosage_0 + (1|SAMPLE), data = cytotoxicity_data_factors_K562)
      col <- j
      results_0 <- get_results(intercept_K562_0, "0", col)
      all <- rbind(all, results_0)
    } else {
      intercept_K562_0 <- NA
      results_0 <- NA
    }

    if(length(levels(cytotoxicity_data_factors_K562$dosage_1)) == 2) {
      intercept_K562_1 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + dosage_1 + (1|SAMPLE), data = cytotoxicity_data_factors_K562)
      col <- j
      results_1 <- get_results(intercept_K562_1, "1", col)
      all <- rbind(all, results_1)
    } else {
      intercept_K562_1 <- NA
      results_1 <- NA
    }

    if(length(levels(cytotoxicity_data_factors_K562$dosage_2)) == 2) {
      intercept_K562_2 <- lmer(VAL ~ DIL_NUM + HLAC + KIR + CMV + dosage_2 + (1|SAMPLE), data = cytotoxicity_data_factors_K562)
      col <- j
      results_2 <- get_results(intercept_K562_2, "2", col)
      all <- rbind(all, results_2)
    } else {
      intercept_K562_2 <- NA
      results_2 <- NA
    }
    
  }
  
}

# for some SNPs there is only one genotype in the blood donors
# for these, no model is ran & the SNPs are missing from the final result table

# add rsIDs?
# like this before
# rsids <- c(paste(rep("rs11585450", 3), c("CC", "GC", "GG")), paste(rep("rs1875763", 3), c("GG", "CG", "CC")), paste(rep("rs8087187", 2), c("AA", "CA")), paste(rep("rs3911730", 2), c("CC", "AC")))


# save table
write.table(all, file="/home/nihtiju/work/NK_receptors/results/cytotoxicity/results_K562.txt", sep="\t", quote=F, row.names=F, col.names=T)


















