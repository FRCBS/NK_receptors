library(data.table)
library(tidyverse)
library(flextable)

# clean tables which are going into the supplementary material

#-----------------------------------------------------------------------------------------------------------------------------------------

# association results

supple_tables <- function(paths, type){
  
  files <- list.files(path = paths, pattern = paste0("results_", type), full.names = T)
  files <- files[!(files %like% ".txt" | files %like% "numerical")]
  
  for (i in 1:length(files)) {
    
    print(i)
    
    results_file <- fread(files[i], data.table=F, header=T)
    
    # remove extra columns
    pheno <- results_file$pheno[1]
    
    if (ncol(results_file) == 21) {
      results_file <- results_file[,c(1,2,10,13,15,17,19,21)]
      colnames(results_file) <- c("Gene", "rsId", "Allele", "Finland", "Spain", "UK", "Poland", "Combined")
    } else {
      results_file <- results_file[,c(1,2,10,13,15,17,19)]
      colnames(results_file) <- c("Gene", "rsId", "Allele", "Finland", "Spain", "UK", "Combined")
    }
    
    write.table(results_file, file=paste0("./results/for_paper/supple_tables/association_", type, "/", pheno, "_", type), sep="\t", quote=F, row.names=F, col.names=T)
    
    # make a flextable
    
    table_flex <- flextable(results_file)
    table_flex <- align(table_flex, part = "header", align = "left")
    table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
    
    save_as_docx(table_flex, path = paste0("./results/for_paper/supple_tables/association_", type, "/", pheno, "_", type, "_flextable.docx"), align = "left")
    
  }
  
}


supple_tables("./results/association_testing_train", "train")
supple_tables("./results/association_testing_test", "test")


#-----------------------------------------------------------------------------------------------------------------------------------

# lasso

files <- list.files(path = "./results/select_SNPs_with_lasso", full.names = T)

for (i in 1:length(files)) {
  
  print(i)
  
  results_file <- fread(files[i], data.table=F, header=T)
  colnames(results_file) <- c("Variable", "Coefficient estimate")
  
  # change SNP names to be rsIds
  SNP_info_all <- fread("./results/for_paper/supple_table_all_snps.txt", data.table=F, header=T)
  
  
  change <- results_file$Variable %like% "chr"
  results_file$Variable[change] <- str_sub(results_file$Variable[change], 1, -3)
  
  SNP_info_all_matched <- SNP_info_all[match(results_file$Variable, SNP_info_all$name_in_ours),]
  results_file$Variable[change] <- SNP_info_all_matched$rs_code[change]
  
  
  write.table(results_file, file=paste0("./results/for_paper/supple_tables/lasso/", str_sub(files[i], 34, -5)), sep="\t", quote=F, row.names=F, col.names=T)
  
  # make a flextable
  
  table_flex <- flextable(results_file)
  table_flex <- align(table_flex, part = "header", align = "left")
  table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
  
  save_as_docx(table_flex, path = paste0("./results/for_paper/supple_tables/lasso/", str_sub(files[i], 34, -5), "_flextable.docx"), align = "left")
  
}

#-----------------------------------------------------------------------------------------------------------------------------------

# in vitro cytotoxicity



results_file <- fread("./results/NK_killing_activity_results/results_K562.txt", data.table=F, header=T)

# order columns better
results_file <- results_file[,c(10,10,7,8,9,5)]

# separate genotype into its own column
results_file$Rsid <- str_sub(results_file$Rsid, 1, -4)
results_file$Rsid.1 <- str_sub(results_file$Rsid.1, -2)
# colnames(results_file)[2] <- "Genotype"
colnames(results_file) <- c("rsId", "Genotype", "Beta", "Lower", "Upper", "P-value")

# make a flextable

table_flex <- flextable(results_file)
table_flex <- align(table_flex, part = "header", align = "left")
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)

save_as_docx(table_flex, path = paste0("./results/for_paper/supple_tables/cytotoxicity/K562_cytotoxicity_flextable.docx"), align = "left")

























