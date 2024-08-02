library(data.table)
library(tidyverse)

# look at new association results

# make resuls tables 

#-----------------------------------------------------------------------------------------------------------------------------

# extract OR and P-values from all populations for deciding on SNPs to continue with

#-----------------------------------------------------------------------------------------------------------------------------

setwd("/home/nihtiju/work/NK_receptors/")


make_association_tables <- function(results_path, tested_set){
  
  # all populations, all phenotypes
  
  directories <- list.dirs(path = results_path, full.names = T, recursive = F)
  
  files <- list.files(path = directories, pattern = ".glm", full.names = T)
  files <- files[!(files %like% "adjusted")]
  
  phenos <- c("pheno_aGvHD_all", "pheno_aGvHD_severe", "pheno_cGvHD_all", "pheno_cGvHD_severe_broader", "pheno_cGvHD_severe", "pheno_relapse")
  
  for (i in 1:length(phenos)) {
    
    print(phenos[i])
    
    files_pheno <- files[files %like% paste0(phenos[i], ".glm.logistic.hybrid")]
    
    # go over files
    for (j in 1:length(files_pheno)) {
      
      # first file --> make results table and save SNP names to match the other results to this
      if (j == 1) {
        
        file <- fread(files_pheno[j])
        file <- file[file$TEST == "ADD",]
        
        results <- file[,c(1:6,10,12,13,15,16)]
        
        # get population name for the column names
        split <- str_split(files_pheno[j], "/")[[1]][10]
        split <- str_sub(split, 1, -(nchar(tested_set) + 2))
        colnames(results)[7:11] <- c(paste0("OR_", split), paste0("L95_", split), paste0("U95_", split), paste0("P_", split), paste0("errorcode_", split))
        
        
      } else {
        
        # other files --> read and match results
        file <- fread(files_pheno[j])
        file <- file[file$TEST == "ADD",]
        file <- file[match(results$ID, file$ID),]
        
        # get population name for the column names
        split <- str_split(files_pheno[j], "/")[[1]][10]
        split <- str_sub(split, 1, -(nchar(tested_set) + 2))
        
        # check that alleles are ok (forced to be the same --> should be the same)
        # all(file[,1:6] == results[,1:6]) # T
        print(paste0("alleles ok, ", split, " : ", all(file[,1:6] == results[,1:6])))
        
        file <- file[,c(10,12,13,15,16)]
        
        # rename columns
        colnames(file) <- c(paste0("OR_", split), paste0("L95_", split), paste0("U95_", split), paste0("P_", split), paste0("errorcode_", split))
        
        results <- cbind(results, file)
        
      }
      
      
    }
    
    # poland does not have all phenos --> files missing for poland
    if (any(files_pheno %like% "poland") == F) {
      
      # add NAs
      cols <- c((ncol(results) - 4) : (ncol(results))) # 5 last cols = the last population
      poland <- results[, ..cols] 
      poland[,1:5] <- NA
      colnames(poland) <- c("OR_poland", "L95_poland", "U95_poland", "P_poland", "errorcode_poland")
      results <- cbind(results, poland)
      
      
    }
    
    # add column for indicating the pheno
    results <- cbind(rep(str_sub(phenos[i], 7), nrow(results)), results)
    colnames(results)[1] <- "pheno"
    
    # save results table
    results <- results[order(results$P_all),]
    write.table(results, paste0(results_path, "/results_", phenos[i], ".txt"), quote = F, col.names = T, row.names = F, sep = "\t")
    
  }

  
}


make_association_tables("/home/nihtiju/work/NK_receptors/results/association_LD/train", "train")
make_association_tables("/home/nihtiju/work/NK_receptors/results/association_LD/test", "test")

make_association_tables("/home/nihtiju/work/NK_receptors/results/association/train", "train")
make_association_tables("/home/nihtiju/work/NK_receptors/results/association/test", "test")

make_association_tables("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train", "train")
make_association_tables("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test", "test")


#------------------------------------------------------------------------------------------------------------

# look at the results

results_aGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_aGvHD_all.txt")
results_aGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_aGvHD_severe.txt")
results_cGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_cGvHD_all.txt")
results_cGvHD_severe_broader <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_cGvHD_severe_broader.txt")
results_cGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_cGvHD_severe.txt")
results_relapse <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train/results_pheno_relapse.txt")


# results_aGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_aGvHD_all.txt")
# results_aGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_aGvHD_severe.txt")
# results_cGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_cGvHD_all.txt")
# results_cGvHD_severe_broader <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_cGvHD_severe_broader.txt")
# results_cGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_cGvHD_severe.txt")
# results_relapse <- fread("/home/nihtiju/work/NK_receptors/results/association_LD/train/results_pheno_relapse.txt")


# results_aGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_aGvHD_all.txt")
# results_aGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_aGvHD_severe.txt")
# results_cGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_cGvHD_all.txt")
# results_cGvHD_severe_broader <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_cGvHD_severe_broader.txt")
# results_cGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_cGvHD_severe.txt")
# results_relapse <- fread("/home/nihtiju/work/NK_receptors/results/association/train/results_pheno_relapse.txt")

filter_results <- function(table, p_cutoff){
  
  return(table[table$P_all < p_cutoff,])
  
}

# all LD-rṕruned together

# # results_aGvHD_all <- filter_results(results_aGvHD_all, 0.01) # 4
# # results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.01) # 2
# # results_cGvHD_all <- filter_results(results_cGvHD_all, 0.01) # 9
# # results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.01) # 2
# # results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.01) # 4
# # results_relapse <- filter_results(results_relapse, 0.01) # 3
# 
# results_aGvHD_all <- filter_results(results_aGvHD_all, 0.005) # 2 kpl
# results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.005) # 0 kpl
# results_cGvHD_all <- filter_results(results_cGvHD_all, 0.005) # 3 kpl
# results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.005) # 1 kpl
# results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.005) # 0 kpl
# results_relapse <- filter_results(results_relapse, 0.005) # 1 kpl
# 
# # results_aGvHD_all <- filter_results(results_aGvHD_all, 0.001) # 0 kpl
# # results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.001) # 0 kpl
# # results_cGvHD_all <- filter_results(results_cGvHD_all, 0.001) # 2 kpl
# # results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.001) # 0 kpl
# # results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.001) # 0 kpl
# # results_relapse <- filter_results(results_relapse, 0.001) # 0 kpl

# only new ones LD pruned together

# results_aGvHD_all <- filter_results(results_aGvHD_all, 0.01) # 4
# results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.01) # 64
# results_cGvHD_all <- filter_results(results_cGvHD_all, 0.01) # 11
# results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.01) # 5
# results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.01) # 3
# results_relapse <- filter_results(results_relapse, 0.01) # 23

results_aGvHD_all <- filter_results(results_aGvHD_all, 0.005) # 2 kpl
results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.005) # 1 kpl
results_cGvHD_all <- filter_results(results_cGvHD_all, 0.005) # 4 kpl
results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.005) # 1 kpl
results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.005) # 2 kpl
results_relapse <- filter_results(results_relapse, 0.005) # 3 kpl

# results_aGvHD_all <- filter_results(results_aGvHD_all, 0.001) # 0 kpl
# results_aGvHD_severe <- filter_results(results_aGvHD_severe, 0.001) # 0 kpl
# results_cGvHD_all <- filter_results(results_cGvHD_all, 0.001) # 3 kpl
# results_cGvHD_severe <- filter_results(results_cGvHD_severe, 0.001) # 0 kpl
# results_cGvHD_severe_broader <- filter_results(results_cGvHD_severe_broader, 0.001) # 0 kpl
# results_relapse <- filter_results(results_relapse, 0.001) # 0 kpl


#-------------------------------------------------------------------

# get the names of SNPs to be checked in test data & in vitro data
# all <- rbind(results_aGvHD_all[1,], results_cGvHD_all[1,], results_relapse[1,])
# write.table(all$ID, "/home/nihtiju/work/NK_receptors/results/association_LD/SNPs_chosen.txt", quote = F, sep = "\t", col.names = F, row.names = F)


all <- rbind(results_aGvHD_all, results_aGvHD_severe, results_cGvHD_all, results_cGvHD_severe, results_cGvHD_severe_broader, results_relapse)
write.table(all, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_chosen.txt", quote = F, sep = "\t", col.names = T, row.names = F)


##################################################################################################################

# get the association results in the test set

results_aGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_aGvHD_all.txt")
results_aGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_aGvHD_severe.txt")
results_cGvHD_all <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_cGvHD_all.txt")
results_cGvHD_severe_broader <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_cGvHD_severe_broader.txt")
results_cGvHD_severe <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_cGvHD_severe.txt")
results_relapse <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test/results_pheno_relapse.txt")


results_aGvHD_all <- results_aGvHD_all[results_aGvHD_all$ID %in% all$ID[1:2],]
results_aGvHD_severe <- results_aGvHD_severe[results_aGvHD_severe$ID == all$ID[3],]
results_cGvHD_all <- results_cGvHD_all[results_cGvHD_all$ID %in% all$ID[4:7],]
results_cGvHD_severe <- results_cGvHD_severe[results_cGvHD_severe$ID == all$ID[8],]
results_cGvHD_severe_broader <- results_cGvHD_severe_broader[results_cGvHD_severe_broader$ID %in% all$ID[9:10],]
results_relapse <- results_relapse[results_relapse$ID %in% all$ID[11:13],]

all_test <- rbind(all, results_aGvHD_all, results_aGvHD_severe, results_cGvHD_all, results_cGvHD_severe, results_cGvHD_severe_broader, results_relapse)

all_test$tested_set <- c(rep("train", nrow(all)), rep("test", nrow(all)))


write.table(all_test, "/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test.txt", quote = F, col.names = T, row.names = F, sep = "\t")










