library(data.table)
library(tidyverse)

setwd("/home/nihtiju/work/NK_receptors/")

results <- fread("/home/nihtiju/work/NK_receptors/results/cytotoxicity/results_K562.txt")
results <- results[results$variable %like% "dosage",]

# round results
results$Estimate <- signif(results$Estimate, 3)
results$beta_95CI_L <- signif(results$beta_95CI_L, 3)
results$beta_95CI_U <- signif(results$beta_95CI_U, 3)
results$`Pr(>|t|)` <- signif(results$`Pr(>|t|)`, 3)


supple <- results[,7]
add <- data.table("GENOTYPE"=character(),
                  "BETA (95% CI)"=character(),
                  "P-VALUE"=character())
supple <- cbind(supple, add)


for (i in 1:nrow(supple)) {
  
  # add genotype
  split <- str_split(supple[i,1], "_")
  
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[!(other %in% counted)]
  
  if(split[[1]][6] == 0){
    supple$GENOTYPE[i] <- paste0(other, other)
  } else if (split[[1]][6] == 1) {
    supple$GENOTYPE[i] <- paste0(other, counted)
  } else if (split[[1]][6] == 2) {
    supple$GENOTYPE[i] <- paste0(counted, counted)
  }
  
  # modify SNP name
  supple[i,1] <- paste0(split[[1]][1], "_", split[[1]][2], "_", split[[1]][3], "_", split[[1]][4])
  
  # add results
  supple$`BETA (95% CI)`[i] <- paste0(results$Estimate[i], " (", results$beta_95CI_L[i], " - ", results$beta_95CI_U[i], ")")
  supple$`P-VALUE`[i] <- results$`Pr(>|t|)`[i]
  
}

write.table(supple, "/home/nihtiju/work/NK_receptors/results/for_paper/in_vitro/in_vitro_results_supple.txt", quote = F, col.names = T, row.names = F, sep = "\t")






