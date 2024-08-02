library(data.table)
library(tidyverse)
library(cowplot)
library(ggthemes)

# plots of NK in vitro data

# load data in
leena_data_factors <- fread("/home/nihtiju/work/NK_receptors/data/NK_assay_processed/NK_data_parsed.tsv", stringsAsFactors = T) # the results are the same if factors or not
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
leena_data_factors$chr1_161726992_G_C_G <- as.factor(dosage_chr_1$chr1_161726992_G_C_G)
leena_data_factors$chr1_161727208_C_G_C <- as.factor(dosage_chr_1$chr1_161727208_C_G_C)
leena_data_factors$chr18_70203271_C_A_C <- as.factor(dosage_chr_18$chr18_70203271_C_A_C)
leena_data_factors$chr18_70204107_A_C_A <- as.factor(dosage_chr_18$chr18_70204107_A_C_A)

leena_data_factors <- leena_data_factors[leena_data_factors$SAMPLE != "NK1",]

# leena_data_factors_K562 <- leena_data_factors[leena_data_factors$CELL == "K562",]

######################################################################################################################################

# line graphs

# add mean of 3 replicates
val_mean <- c()
for (i in seq(1, nrow(leena_data_factors), by=3)) {
  # seq(1, nrow(leena_data_factors), by=3) is every 3rd row
  
  # print(i)
  
  rep_1 <- leena_data_factors$VAL[i]
  rep_2 <- leena_data_factors$VAL[i+1]
  rep_3 <- leena_data_factors$VAL[i+2]
  
  means <- mean(rep_1, rep_2, rep_3, na.rm=T) # THIS IS WRONG! SHOULD BE means <- mean(c(rep_1, rep_2, rep_3), na.rm=T)
  val_mean <- c(val_mean, means, means, means)
}
leena_data_factors$VAL_mean <- val_mean

leena_data_factors_K562 <- leena_data_factors[leena_data_factors$CELL == "K562",]

# plots

# one line = mean for the replicates of one individual

plot_dilutions <- function(data, celline){
  
  plot_list <-  list()
  snps <- c("chr1_161726992_G_C_G", "chr1_161727208_C_G_C", "chr18_70203271_C_A_C", "chr18_70204107_A_C_A")
  rsids <- c("rs11585450", "rs1875763", "rs8087187", "rs3911730")
  alleles <- c("GC", "CG", "CA", "AC") # the one tested first here
  
  for (i in 1:4) {
    
    snp <- snps[i]
    groups <- as.vector(data[,..snp])
    groups <- unname(unlist(groups)) # dosage 0-2
    
    rsid <- rsids[i]
    allele <- alleles[i]
    allele_1 <- substring(allele, 1, 1)
    allele_2 <- substring(allele, 2, 2)
    
    if (length(unique(groups)) == 3) {
      allele_all <- c(paste0(allele_2, allele_2), paste0(allele_1, allele_2), paste0(allele_1, allele_1))
    } else {
      allele_all <- c(paste0(allele_2, allele_2), paste0(allele_1, allele_2))
    }
    
    # print(snp)
    # print(unique(groups))
    
    data$groups_variable <- groups # needed to get colours by genotype groups
    
    plot <- ggplot(data, aes(x=DIL_NUM, y=VAL_mean, color=groups_variable, group=SAMPLE)) +
      geom_point() +
      geom_line() + # geom_line(linewidth = 1) +
      ylab("Cytotoxicity") +
      xlab("Effector-target ratio") +
      ggtitle(paste0(celline, ": ", rsid)) +
      theme_minimal() +
      theme(legend.title=element_blank(), legend.position = "bottom", legend.direction = "horizontal", text = element_text(size = 20)) +
      scale_color_discrete(name = "", labels = allele_all)
    
    plot_list[[i]] <-  plot
    
  }
  
  return(plot_list)
  
}

plots_K562 <- plot_dilutions(leena_data_factors_K562, "K562")

plot_grid(plots_K562[[1]], plots_K562[[2]], plots_K562[[3]], plots_K562[[4]], ncol = 2)

ggsave("./results/NK_killing_activity_results/killing_activity_by_dilutions_all_samples_K562.png", bg = "white", width = 20, height = 15, dpi = 600)

















