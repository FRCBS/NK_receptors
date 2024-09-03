library(data.table)
library(tidyverse)


modify_results <- function(results_path, dataset){
  
  files <- list.files(path = results_path, pattern = "results_pheno", full.names = T)
  
  for (i in 1:length(files)) {
    
    print(files[i])
    
    results <- fread(files[i])
    
    supple <- results[,1:7]
    colnames(supple)[1] <- "PHENOTYPE"
    colnames(supple)[2] <- "CHR"
    
    add <- data.table("COMBINED"=character(),
                      "FINLAND"=character(),
                      "UK"=character(),
                      "SPAIN"=character(),
                      "POLAND"=character())
    supple <- cbind(supple, add)
    
    for (j in 1:nrow(supple)) {
      
      supple[j,8] <- paste0("OR: ", results[j,8], " (", results[j,9], " - ", results[j,10], "), P: ", results[j,11]) # combined
      supple[j,9] <- paste0("OR: ", results[j,13], " (", results[j,14], " - ", results[j,15], "), P: ", results[j,16]) # finland
      supple[j,10] <- paste0("OR: ", results[j,23], " (", results[j,24], " - ", results[j,25], "), P: ", results[j,26]) # uk
      supple[j,11] <- paste0("OR: ", results[j,18], " (", results[j,19], " - ", results[j,20], "), P: ", results[j,21]) # spain
      supple[j,12] <- paste0("OR: ", results[j,28], " (", results[j,29], " - ", results[j,30], "), P: ", results[j,31]) # poland
      
    }
    
    write.table(supple, paste0("/home/nihtiju/work/NK_receptors/results/for_paper/association/", dataset, "_", supple$PHENOTYPE[1], ".txt"), quote = F, col.names = T, row.names = F, sep = "\t")
    
  }
  
}

modify_results("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/train", "train")
modify_results("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/test", "test")

