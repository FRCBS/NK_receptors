library(data.table)
library(tidyverse)
library(cowplot)
library(ggthemes)

setwd("/home/nihtiju/work/NK_receptors/")

# plots of NK in vitro data

# load data in
cytotoxicity_data_factors <- fread("/home/nihtiju/work/NK_receptors/data/NK_assay_processed/NK_data_parsed.tsv", stringsAsFactors = T) 
# add dosage
dosage_chr_1 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr1_extracted_dosage.raw")
dosage_chr_6 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr6_extracted_dosage.raw")
dosage_chr_7 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr7_extracted_dosage.raw")
dosage_chr_12 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr12_extracted_dosage.raw")
dosage_chr_18 <- fread("/home/nihtiju/work/NK_receptors/results/blood_donors/blood_donors_chr18_extracted_dosage.raw")

dosage <- cbind(dosage_chr_1, dosage_chr_6[,7:ncol(dosage_chr_6)], dosage_chr_7[,7:ncol(dosage_chr_7)], dosage_chr_12[,7:ncol(dosage_chr_12)], dosage_chr_18[,7:ncol(dosage_chr_18)])

# unify sample names

# # pseudonyms in genotype data do not match in vitro data's IDs 
# first change the IDs according to an ID file, then modify IDs to match the in vitro data's IDs
ids <- fread("/home/nihtiju/work/NK_receptors/data/blood_donors/blood_donors_IDs.txt", header = F)
ids <- ids[match(dosage$IID, ids$V1),]
dosage$IID <- ids[,2]

# this for when the IDs are correct & when fixing the code
dosage$sample <- paste0(str_sub(dosage$IID, 1, 2), str_sub(dosage$IID, 6, 7))
dosage$sample <- str_replace_all(dosage$sample, "NK0", "NK")

# order dosages to be the same as other data
dosage <- dosage[match(cytotoxicity_data_factors$SAMPLE, dosage$sample),]



for (i in 7:17) {
  
  cytotoxicity_data_factors[,as.character(as.name(colnames(dosage[,..i])))] <- dosage[,..i]
  
}

######################################################################################################################################

# line graphs

# add mean of 3 replicates
val_mean <- c()
for (i in seq(1, nrow(cytotoxicity_data_factors), by=3)) {
  # seq(1, nrow(cytotoxicity_data_factors), by=3) is every 3rd row
  
  # print(i)
  
  rep_1 <- cytotoxicity_data_factors$VAL[i]
  rep_2 <- cytotoxicity_data_factors$VAL[i+1]
  rep_3 <- cytotoxicity_data_factors$VAL[i+2]
  
  means <- mean(c(rep_1, rep_2, rep_3), na.rm=T) 
  val_mean <- c(val_mean, means, means, means)
  
}
cytotoxicity_data_factors$VAL_mean <- val_mean

cytotoxicity_data_factors_K562 <- cytotoxicity_data_factors[cytotoxicity_data_factors$CELL == "K562",]



# plots

# one line = mean for the replicates of one individual

plot_dilutions <- function(data, celline){
  
  plot_list <-  list()
  snps <- colnames(cytotoxicity_data_factors_K562)[9:19]
  snps_short <- str_sub(snps, 1, -3)
  
  # get rsIDs
  rsIDs <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  rsIDs <- rsIDs[match(snps_short, rsIDs$ID),]
  
  for (i in 1:length(snps)) {
    
    rsID <- rsIDs$rsID[i]
    
    snp <- snps[i]
    groups <- as.vector(data[,..snp])
    groups <- unname(unlist(groups))
    groups <- as.factor(groups) # dosage 0-2 
    
    split <- str_split(snp, "_")
    allele_1 <- split[[1]][5] # the one counted in dosage
    allele_2 <- c(split[[1]][3], split[[1]][4])
    allele_2 <- allele_2[allele_2 != allele_1] # the one not counted in dosage
    
    if (length(unique(groups)) == 3) {
      allele_all <- c(paste0(allele_2, allele_2), paste0(allele_1, allele_2), paste0(allele_1, allele_1))
    } else {
      allele_all <- c(paste0(allele_2, allele_2), paste0(allele_1, allele_2))
    }
    
    data$groups_variable <- groups # needed to get colours by genotype groups
    
    plot <- ggplot(data, aes(x=DIL_NUM, y=VAL_mean, color=groups_variable, group=SAMPLE)) +
      geom_point() +
      geom_line() + # geom_line(linewidth = 1) +
      ylab("Cytotoxicity") +
      xlab("Effector-target ratio") +
      ggtitle(paste0(celline, ", ", snp, ", ", rsID)) +
      theme_minimal() +
      theme(legend.title=element_blank(), legend.position = "bottom", legend.direction = "horizontal", text = element_text(size = 10)) +
      scale_color_discrete(name = "", labels = allele_all)
    
    plot_list[[i]] <-  plot
    
  }
  
  return(plot_list)
  
}

plots_K562 <- plot_dilutions(cytotoxicity_data_factors_K562, "K562")

plot_grid(plots_K562[[1]], plots_K562[[2]], plots_K562[[3]], 
          plots_K562[[4]], plots_K562[[5]], plots_K562[[6]], 
          plots_K562[[7]], plots_K562[[8]], plots_K562[[9]], 
          plots_K562[[10]], plots_K562[[11]],
          ncol = 3)

ggsave("./results/for_paper/in_vitro/killing_activity_by_dilutions_all_samples_K562.png", bg = "white", width = 15, height = 10, dpi = 600)

