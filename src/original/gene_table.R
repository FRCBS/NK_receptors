library(data.table)
library(tidyverse)
library(flextable)

# a list of all the SNPs already made in ./src/supple_table_SNPs.R
# -> use that to make a gene-wise table

snp <- fread("../NK_receptors/results/for_paper/supple_table_all_snps.txt")
unique(snp$gene)
# [1] "CD94 (KLRD1)" 
# [2] "NKG2A (KLRC1)"
# [3] "NKG2C (KLRC2)"
# [4] "NKG2D (KLRK1)"
# [5] "NKG2D(KLRK1)" 
# [6] "NCR1"         
# [7] "NCR3"         
# [8] "NCR2"         
# [9] "ILT4"         
# [10] "ILT2"         
# [11] "ADGRG1"       
# [12] "CD2"          
# [13] "CD226"        
# [14] "CD244"        
# [15] "CD244 FCGR3A" 
# [16] "FCGR3A"       
# [17] "NKG2A"  
snp$gene[6] <- "NKG2D (KLRK1)" 
snp$gene[894] <- "NKG2A (KLRC1)"
unique(snp$gene)
# [1] "CD94 (KLRD1)" 
# [2] "NKG2A (KLRC1)"
# [3] "NKG2C (KLRC2)"
# [4] "NKG2D (KLRK1)"
# [5] "NCR1"         
# [6] "NCR3"         
# [7] "NCR2"         
# [8] "ILT4"         
# [9] "ILT2"         
# [10] "ADGRG1"       
# [11] "CD2"          
# [12] "CD226"        
# [13] "CD244"        
# [14] "CD244 FCGR3A" 
# [15] "FCGR3A"       
# [16] "NKG2A"  

genes <- data.table(gene=unique(snp$gene))
genes$type <- "receptor"
genes$SNPs <- 0

# add the number of SNPs for each gene
for (i in 1:nrow(genes)) {
  
  genes$SNPs[i] <- sum(snp$gene == genes$gene[i])
  
}

genes <- genes[order(genes$gene),]

colnames(genes) <- c("Gene", "Type", "Number of SNPs")

# headers
table_flex <- flextable(genes)

table_flex <- align(table_flex, part = "header", align = "left")
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
table_flex <- footnote(
  table_flex,
  i = 1,
  j = 1,
  value = as_paragraph("Gene(s) in which the polymorphisms are located in or which they regulate"),
  ref_symbols = "*",
  part = "header"
)


set_flextable_defaults(background.color = "white")

save_as_image(table_flex, path = "../NK_receptors/results/for_paper/flextables/genes.png", bg = "white")
save_as_docx(table_flex, path = "../NK_receptors/results/for_paper/flextables/genes.docx", align = "left")












