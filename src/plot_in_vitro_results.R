library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)

############################################################################################################################################################

## -------------------------------------------------------------
## Cytotoxicity results data
## -------------------------------------------------------------

K562 <- fread("/home/nihtiju/work/NK_receptors/results/cytotoxicity/results_K562.txt")
K562$group <- 'Cytotoxicity'
K562$SNP_short <- str_sub(K562$SNP, 1, -5)

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[match(K562$SNP_short, SNP_info$ID),]

K562$rsID <- SNP_info$rsID
K562$gene <- SNP_info$gene

# change rsIDs to include the tested genotypes
K562$geno <- NA
for (i in 1:nrow(K562)) {
  
  split <- str_split(K562$SNP[i], "_")
  counted <- split[[1]][5]
  other <- c(split[[1]][3], split[[1]][4])
  other <- other[other != counted]
  
  if(split[[1]][6] == 0){
    K562$geno[i] <- paste0(K562$rsID[i], " ", other, other)
  } else if (split[[1]][6] == 1) {
    K562$geno[i] <- paste0(K562$rsID[i], " ", counted, other)
  } else if (split[[1]][6] == 2) {
    K562$geno[i] <- paste0(K562$rsID[i], " ", counted, counted)
  }
  
}

# remove covariates
K562 <- K562[K562$variable %like% "dosage",]


## -------------------------------------------------------------
## Plots
## -------------------------------------------------------------

# common theme elements
th <- theme(panel.spacing = unit(.5, "lines"), 
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
            strip.background = element_rect(color = "black", linewidth = 1), 
            text = element_text(size = 15), 
            legend.position = "bottom", 
            legend.text = element_text(size = 14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.text = element_text(size=15), 
            axis.title.x = element_text(size=15), 
            axis.text.x = element_text(size=12))

# NK / all SNPs
p1 <- ggplot(K562, aes(beta, geno))  + 
  geom_point(size = 4) +
  geom_linerange(aes(xmin = beta_95CI_L, xmax = beta_95CI_U), linewidth = 1) + 
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  # the SNPs and genotypes in the same order as in the HSCT association plot
  # listed from top to bottom -> rev to make it from bottom to top
  scale_y_discrete(limits = rev(c("rs117914097 CC", "rs117914097 TC", 
                                  "rs116433596 GG", "rs116433596 AG", 
                                  "rs62483646 GG", "rs62483646 AG", 
                                  "rs4656994 GG", "rs4656994 AG", "rs4656994 AA", 
                                  "rs1875763 CC", "rs1875763 GC", "rs1875763 GG", 
                                  "rs11585450 GG", "rs11585450 CG", "rs11585450 CC", 
                                  "rs576627 CT", "rs576627 TT", 
                                  "rs3911730 AC", "rs3911730 CC", 
                                  "rs140677956 CC", "rs140677956 TC", "rs140677956 TT"))) +
  ylab("") + xlab("beta") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 5.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7.5, ymax = 10.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 16.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 18.5, ymax = 20.5, alpha = .25, fill = 'grey') +
  facet_grid(cols = vars(group), scales = 'free') +
  theme_minimal() + th

ggsave("./results/for_paper/in_vitro/in_vitro_results_plot_beta.png", bg = "white", width = 5, height = 10, dpi = 600)

# 2 SNPs missing (not in blood donor genotype data)
# aGvHD_2 SNP missing 4790
# the first relapse SNP missing 6312

# 2 SNPs do not have homozygotes for the allele in HSCT association results
# rs576627 CT
# rs3911730 AC

# SNPs twice in the HSCT association results - only once here, no need to duplicate
# rs4656994 & rs576627 



#-------------------------------------------------------

# plot with OR (in stead of beta)

K562$OR_95CI_U <- as.numeric(K562$OR_95CI_U) # this was character for some reason -> error when plotting

p2 <- ggplot(K562, aes(OR, geno))  + 
  geom_point(size = 4) +
  # geom_linerange(aes(xmin = OR_95CI_L, xmax = OR_95CI_U), linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  # the SNPs and genotypes in the same order as in the HSCT association plot
  # listed from top to bottom -> rev to make it from bottom to top
  scale_y_discrete(limits = rev(c("rs117914097 CC", "rs117914097 TC", 
                                  "rs116433596 GG", "rs116433596 AG", 
                                  "rs62483646 GG", "rs62483646 AG", 
                                  "rs4656994 GG", "rs4656994 AG", "rs4656994 AA", 
                                  "rs1875763 CC", "rs1875763 GC", "rs1875763 GG", 
                                  "rs11585450 GG", "rs11585450 CG", "rs11585450 CC", 
                                  "rs576627 CT", "rs576627 TT", 
                                  "rs3911730 AC", "rs3911730 CC", 
                                  "rs140677956 CC", "rs140677956 TC", "rs140677956 TT"))) +
  ylab("") + xlab("OR") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 5.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7.5, ymax = 10.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 16.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 18.5, ymax = 20.5, alpha = .25, fill = 'grey') +
  facet_grid(cols = vars(group), scales = 'free') +
  theme_minimal() + th

ggsave("./results/for_paper/in_vitro/in_vitro_results_plot_OR_noCIs.png", bg = "white", width = 5, height = 10, dpi = 600)


p2 <- ggplot(K562, aes(OR, geno))  + 
  geom_point(size = 4) +
  geom_linerange(aes(xmin = OR_95CI_L, xmax = OR_95CI_U), linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  # the SNPs and genotypes in the same order as in the HSCT association plot
  # listed from top to bottom -> rev to make it from bottom to top
  scale_y_discrete(limits = rev(c("rs117914097 CC", "rs117914097 TC", 
                                  "rs116433596 GG", "rs116433596 AG", 
                                  "rs62483646 GG", "rs62483646 AG", 
                                  "rs4656994 GG", "rs4656994 AG", "rs4656994 AA", 
                                  "rs1875763 CC", "rs1875763 GC", "rs1875763 GG", 
                                  "rs11585450 GG", "rs11585450 CG", "rs11585450 CC", 
                                  "rs576627 CT", "rs576627 TT", 
                                  "rs3911730 AC", "rs3911730 CC", 
                                  "rs140677956 CC", "rs140677956 TC", "rs140677956 TT"))) +
  ylab("") + xlab("OR") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 5.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7.5, ymax = 10.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 16.5, alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 18.5, ymax = 20.5, alpha = .25, fill = 'grey') +
  facet_grid(cols = vars(group), scales = 'free') +
  theme_minimal() + th

ggsave("./results/for_paper/in_vitro/in_vitro_results_plot_OR_yesCIs.png", bg = "white", width = 5, height = 10, dpi = 600)


















