library(data.table)
library(tidyverse)

# look at new association results

# make resuls tables 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# paste variables together for making supple tables

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

SNP_info_all <- fread("./results/for_paper/supple_table_all_snps.txt", data.table=F, header=T) # duplicates (same rsID in different genes) combined

# add the tested allele
results <- fread("./results/association_testing_train/all_train/association_all_train_aGvHD.pheno_aGvHD_all.glm.logistic.hybrid", data.table=F, header=T) # duplicates (same rsID in different 
results <- results[results$TEST == "ADD",]
results <- results[match(SNP_info_all$name_in_ours, results$ID),]
SNP_info_all$allele <- results$A1

process_results <- function(pop_res_paths, save_to, test){
  
  
  pop_paths <- pop_res_paths
  
  test_type <- test
  
  # joka fenotyypille oma results table = kopio SNP_info_all:sta
  # sinne sarake jossa fenon nimi
  
  aGvHD_all <- cbind(SNP_info_all, pheno="aGvHD_all")
  aGvHD_severe <- cbind(SNP_info_all, pheno="aGvHD_severe")
  relapse <- cbind(SNP_info_all, pheno="relapse")
  cGvHD_all <- cbind(SNP_info_all, pheno="cGvHD_all")
  cGvHD_severe_broader <- cbind(SNP_info_all, pheno="cGvHD_severe_broader")
  cGvHD_severe <- cbind(SNP_info_all, pheno="cGvHD_severe")
  
  # go over populations
  for (i in 1:length(pop_paths)) {
    
    files <- list.files(path = pop_paths[i], pattern = ".glm", full.names = T)
    files <- files[!(files %like% "adjusted")]
    
    # go over phenotypes
    for (j in 1:length(files)) {
      
      results_file <- fread(files[j], data.table=F, header=T)
      results_file <- results_file[results_file$TEST == test_type,] # leave only SNPs
      
      # colnames(results_file)
      # [1] "#CHROM"     "POS"        "ID"         "REF"       
      # [5] "ALT"        "A1"         "FIRTH?"     "TEST"      
      # [9] "OBS_CT"     "OR"         "LOG(OR)_SE" "L95"       
      # [13] "U95"        "Z_STAT"     "P"          "ERRCODE"  
      
      # results_file$result <- paste0("P: ", round(as.numeric(results_file$P), 4), ", OR: ", round(as.numeric(results_file$OR), 2), " (", round(as.numeric(results_file$L95), 2), "–", round(as.numeric(results_file$U95), 2), ")")
      results_file$result <- paste0("P: ", signif(as.numeric(results_file$P), 3), ", OR: ", round(as.numeric(results_file$OR), 2), " (", round(as.numeric(results_file$L95), 2), "–", round(as.numeric(results_file$U95), 2), ")") # p-values with three significant digit
      # gives errors for poland if there is no as.numeric for L95, it is character for some reason
      
      results_file <- results_file[,c(3,15,17)]
      
      # if pheno == joku
      # käy kaikki fenot läpi
      # sen fenon tulostaulukkoon lisätään nämä tulokset: p-arvo ja result ja niiden sarakenimiin myös pop nimi
      # yhteen fenon tulostaulukkoon kaikkien populaatioiden tulokset ja niitä voi järjestää populaatiokohtaisilla p-arvoilla
      
      filename <- str_split(files[j], "/")[[1]][5]
      pheno <- str_split(filename, "\\.")[[1]][2]
      pop <- str_split(filename, "_")[[1]][2]
      
      # print(paste(pop, pheno))
      print(files[j])
      print("")
      
      
      # add pop label to colnames
      colnames(results_file)[2] <- paste0(colnames(results_file)[2], "_", pop)
      colnames(results_file)[3] <- paste0(colnames(results_file)[3], "_", pop)
      
      # ordr the results file to match the bigger pheno results file
      results_file <- results_file[match(aGvHD_all$name_in_ours, results_file$ID),] # all pheno results files have the same rows -> any is ok for matching
      
      if (pheno == "pheno_aGvHD_all") {
        aGvHD_all <- cbind(aGvHD_all, results_file[,2:3])
      } else if  (pheno == "pheno_aGvHD_severe") {
        aGvHD_severe <- cbind(aGvHD_severe, results_file[,2:3])
      } else if  (pheno == "pheno_relapse") {
        relapse <- cbind(relapse, results_file[,2:3])
      } else if  (pheno == "pheno_cGvHD_all") {
        cGvHD_all <- cbind(cGvHD_all, results_file[,2:3])
      } else if  (pheno == "pheno_cGvHD_severe_broader") {
        cGvHD_severe_broader <- cbind(cGvHD_severe_broader, results_file[,2:3])
      } else if  (pheno == "pheno_cGvHD_severe") {
        cGvHD_severe <- cbind(cGvHD_severe, results_file[,2:3])
      }
      
    }
    
  }
  
  # save the tables
  # order by population combination p-values 
  write.table(aGvHD_all[order(aGvHD_all$P_all),], file=save_to[1], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(aGvHD_severe[order(aGvHD_severe$P_all),], file=save_to[2], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(relapse[order(relapse$P_all),], file=save_to[3], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_all[order(cGvHD_all$P_all),], file=save_to[4], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_severe_broader[order(cGvHD_severe_broader$P_all),], file=save_to[5], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_severe[order(cGvHD_severe$P_all),], file=save_to[6], sep="\t", quote=F, row.names=F, col.names=T)
  
}

# training set, additive
process_results(paste0("./results/association_testing_train/", c("finns", "katalonia", "newcastle", "poland", "all"), "_train"),
                paste0("./results/association_testing_train/results_train_pheno_", c("aGvHD_all", "aGvHD_severe", "relapse", "cGvHD_all", "cGvHD_severe_broader", "cGvHD_severe")),
                "ADD")

# test set, additive
process_results(paste0("./results/association_testing_test/", c("finns", "katalonia", "newcastle", "poland", "all"), "_test"),
                paste0("./results/association_testing_test/results_test_pheno_", c("aGvHD_all", "aGvHD_severe", "relapse", "cGvHD_all", "cGvHD_severe_broader", "cGvHD_severe")),
                "ADD")

#-----------------------------------------------------------------------------------------------

# looking at the tables
# TRAIN set, ADDITIVE

aGvHD_all <- fread("./results/association_testing_train/results_train_pheno_aGvHD_all") 

aGvHD_severe <- fread("./results/association_testing_train/results_train_pheno_aGvHD_severe") 

relapse <- fread("./results/association_testing_train/results_train_pheno_relapse") 

cGvHD_all <- fread("./results/association_testing_train/results_train_pheno_cGvHD_all") 

cGvHD_severe <- fread("./results/association_testing_train/results_train_pheno_cGvHD_severe")

cGvHD_severe_broader <- fread("./results/association_testing_train/results_train_pheno_cGvHD_severe_broader") 

get_smallest <- function(table){
  
  table <- table[table$P_all < 0.01,]
  return(table)
  
}

aGvHD_all_filtered <- get_smallest(aGvHD_all) # 0
aGvHD_severe_filtered <- get_smallest(aGvHD_severe) # 1

relapse_filtered <- get_smallest(relapse) # 15

cGvHD_all_filtered <- get_smallest(cGvHD_all) # 2
cGvHD_severe_filtered <- get_smallest(cGvHD_severe) # 0
cGvHD_severe_broader_filtered <- get_smallest(cGvHD_severe_broader) # 1

# add  missing columns to merge the tables
aGvHD_severe_filtered <- cbind(aGvHD_severe_filtered[,1:17], P_poland = NA, result_poland = NA, aGvHD_severe_filtered[,18:19])
cGvHD_severe_broader_filtered <- cbind(cGvHD_severe_broader_filtered[,1:17], P_poland = NA, result_poland = NA, cGvHD_severe_broader_filtered[,18:19])

save <- rbind(aGvHD_severe_filtered, relapse_filtered, cGvHD_all_filtered, cGvHD_severe_broader_filtered)

write.table(save, "./results/association_testing_train/chosen_snps_train.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#-----------------------------------------------------------------------------------------------

# looking at the tables
# TEST set, ADDITIVE

aGvHD_all_test <- fread("./results/association_testing_test/results_train_pheno_aGvHD_all")

aGvHD_severe_test <- fread("./results/association_testing_test/results_train_pheno_aGvHD_severe")

relapse_test <- fread("./results/association_testing_test/results_train_pheno_relapse")

cGvHD_all_test <- fread("./results/association_testing_test/results_train_pheno_cGvHD_all")

cGvHD_severe_test <- fread("./results/association_testing_test/results_train_pheno_cGvHD_severe")

cGvHD_severe_broader_test <- fread("./results/association_testing_test/results_train_pheno_cGvHD_severe_broader")


# leave in the chosen SNPs
chosen <- fread("./results/association_testing_train/chosen_snps_train.txt") 

aGvHD_severe_test <- cbind(aGvHD_severe_test[,1:17], P_poland = NA, result_poland = NA, aGvHD_severe_test[,18:19])
cGvHD_severe_broader_test <- cbind(cGvHD_severe_broader_test[,1:17], P_poland = NA, result_poland = NA, cGvHD_severe_broader_test[,18:19])

aGvHD_severe_test <- aGvHD_severe_test[aGvHD_severe_test$rs_code %in% aGvHD_severe_filtered$rs_code,]
relapse_test <- relapse_test[relapse_test$rs_code %in% relapse_filtered$rs_code,]
cGvHD_all_test <- cGvHD_all_test[cGvHD_all_test$rs_code %in% cGvHD_all_filtered$rs_code,]
cGvHD_severe_broader_test <- cGvHD_severe_broader_test[cGvHD_severe_broader_test$rs_code %in% cGvHD_severe_broader_filtered$rs_code,]

save <- rbind(aGvHD_severe_test, relapse_test, cGvHD_all_test, cGvHD_severe_broader_test)
save <- save[match(chosen$rs_code, save$rs_code),]

write.table(save, "./results/association_testing_test/chosen_snps_test.txt", sep = "\t", quote = F, row.names = F, col.names = T)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# all individual variables for plotting them

#------------------------------------------------------------------------------------------------------------------------------------------------------------------


SNP_info_all <- fread("./results/for_paper/supple_table_all_snps.txt", data.table=F, header=T) 
# includes SNPs from different genes as one entry (whereas separately in supple_table_all_snps_locations_ifFound_all_imputed.txt)


# add the tested allele
results <- fread("./results/association_testing_train/all_train/association_all_train_aGvHD.pheno_aGvHD_all.glm.logistic.hybrid", data.table=F, header=T) # duplicates (same rsID in different 
results <- results[results$TEST == "ADD",]
results <- results[match(SNP_info_all$name_in_ours, results$ID),]
SNP_info_all$allele <- results$A1


process_results <- function(pop_res_paths, save_to, test){
  
  
  pop_paths <- pop_res_paths
  
  test_type <- test
  
  # joka fenotyypille oma results table = kopio SNP_info_all:sta
  # sinne sarake jossa fenon nimi
  
  aGvHD_all <- cbind(SNP_info_all, pheno="aGvHD_all")
  aGvHD_severe <- cbind(SNP_info_all, pheno="aGvHD_severe")
  relapse <- cbind(SNP_info_all, pheno="relapse")
  cGvHD_all <- cbind(SNP_info_all, pheno="cGvHD_all")
  cGvHD_severe_broader <- cbind(SNP_info_all, pheno="cGvHD_severe_broader")
  cGvHD_severe <- cbind(SNP_info_all, pheno="cGvHD_severe")
  
  # go over populations
  for (i in 1:length(pop_paths)) {
    
    files <- list.files(path = pop_paths[i], pattern = ".glm", full.names = T)
    files <- files[!(files %like% "adjusted")]
    
    # go over phenotypes
    for (j in 1:length(files)) {
      
      results_file <- fread(files[j], data.table=F, header=T)
      results_file <- results_file[results_file$TEST == test_type,] # leave only SNPs
      
      # colnames(results_file)
      # [1] "#CHROM"     "POS"        "ID"         "REF"       
      # [5] "ALT"        "A1"         "FIRTH?"     "TEST"      
      # [9] "OBS_CT"     "OR"         "LOG(OR)_SE" "L95"       
      # [13] "U95"        "Z_STAT"     "P"          "ERRCODE"  
      
      results_file <- results_file[,c(3,10,11,12,13,15)]
      # colnames(results_file)
      # [1] "ID"  "OR"  "LOG(OR)_SE"  "L95" "U95" "P"  
      
      # LOG(OR)_]SE	se	Standard error of log-odds (i.e. beta)
      
      colnames(results_file)[3] <- "se_beta"
      
      results_file$beta <- log(results_file$OR)
      
      # if pheno == joku
      # käy kaikki fenot läpi
      # sen fenon tulostaulukkoon lisätään nämä tulokset: p-arvo ja result ja niiden sarakenimiin myös pop nimi
      # yhteen fenon tulostaulukkoon kaikkien populaatioiden tulokset ja niitä voi järjestää populaatiokohtaisilla p-arvoilla
      
      filename <- str_split(files[j], "/")[[1]][5]
      pheno <- str_split(filename, "\\.")[[1]][2]
      pop <- str_split(filename, "_")[[1]][2]
      
      # print(paste(pop, pheno))
      print(files[j])
      print("")
      
      
      # add pop label to colnames
      colnames(results_file) <- paste0(colnames(results_file), "_", pop)
      
      # ordr the results file to match the bigger pheno results file
      results_file <- results_file[match(aGvHD_all$name_in_ours, results_file$ID),] # all pheno results files have the same rows -> any is ok for matching
      
      if (pheno == "pheno_aGvHD_all") {
        aGvHD_all <- cbind(aGvHD_all, results_file)
      } else if  (pheno == "pheno_aGvHD_severe") {
        aGvHD_severe <- cbind(aGvHD_severe, results_file)
      } else if  (pheno == "pheno_relapse") {
        relapse <- cbind(relapse, results_file)
      } else if  (pheno == "pheno_cGvHD_all") {
        cGvHD_all <- cbind(cGvHD_all, results_file)
      } else if  (pheno == "pheno_cGvHD_severe_broader") {
        cGvHD_severe_broader <- cbind(cGvHD_severe_broader, results_file)
      } else if  (pheno == "pheno_cGvHD_severe") {
        cGvHD_severe <- cbind(cGvHD_severe, results_file)
      }
      
    }
    
    
  }
  
  # save the tables
  write.table(aGvHD_all[order(aGvHD_all$P_all),], file=save_to[1], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(aGvHD_severe[order(aGvHD_severe$P_all),], file=save_to[2], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(relapse[order(relapse$P_all),], file=save_to[3], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_all[order(cGvHD_all$P_all),], file=save_to[4], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_severe_broader[order(cGvHD_severe_broader$P_all),], file=save_to[5], sep="\t", quote=F, row.names=F, col.names=T)
  write.table(cGvHD_severe[order(cGvHD_severe$P_all),], file=save_to[6], sep="\t", quote=F, row.names=F, col.names=T)
  
  
}

# additive train
process_results(paste0("./results/association_testing_train/", c("finns", "katalonia", "newcastle", "poland", "all"), "_train"),
                paste0("./results/association_testing_train/results_train_pheno_resultsPlots_numerical_", c("aGvHD_all", "aGvHD_severe", "relapse", "cGvHD_all", "cGvHD_severe_broader", "cGvHD_severe")),
                "ADD")
# additive test
process_results(paste0("./results/association_testing_test/", c("finns", "katalonia", "newcastle", "poland", "all"), "_test"),
                paste0("./results/association_testing_test/results_test_pheno_resultsPlots_numerical_", c("aGvHD_all", "aGvHD_severe", "relapse", "cGvHD_all", "cGvHD_severe_broader", "cGvHD_severe")),
                "ADD")

#-------------------------------------------------

relapse <- fread("./results/association_testing_train/results_train_pheno_resultsPlots_numerical_relapse") 
cGvHD_all <- fread("./results/association_testing_train/results_train_pheno_resultsPlots_numerical_cGvHD_all")
relapse_test <- fread("./results/association_testing_test/results_test_pheno_resultsPlots_numerical_relapse") 
cGvHD_all_test <- fread("./results/association_testing_test/results_test_pheno_resultsPlots_numerical_cGvHD_all")

# results for chosen SNPs
# more variants chosen above, here fewer (change this?)
all <- rbind(relapse[1:2,], cGvHD_all[1:2,])
all_test <- rbind(relapse_test[relapse_test$name_in_ours %in% unname(unlist(as.vector(relapse[1:2,3]))),], cGvHD_all_test[cGvHD_all_test$name_in_ours %in% unname(unlist(as.vector(cGvHD_all[1:2,3]))),])

all_test <- all_test[match(all$name_in_ours, all_test$name_in_ours),]

write.table(all, "./results/association_testing_train/chosen_snps_resultsPlots_numerical.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(all_test, "./results/association_testing_test/chosen_snps_resultsPlots_numerical_testset.txt", sep = "\t", quote = F, row.names = F, col.names = T)





