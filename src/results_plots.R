library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)

############################################################################################################################################################

## -------------------------------------------------------------
## HSCT results data
## -------------------------------------------------------------

# read discovery and test set association results

all_disc <- fread("./results/association_testing_train/chosen_snps_resultsPlots_numerical.txt")[,-10]
all_test <- fread("./results/association_testing_test/chosen_snps_resultsPlots_numerical_testset.txt")[,-10]


get_tables <- function(x) {

  results <- data.table(gene=character(),
                        rs_code=character(),
                        name_in_ours=character(),
                        chr=numeric(), pos_hg38=numeric(),
                        pheno=character(), beta=numeric(),
                        se_beta=numeric(), lower=numeric(),
                        upper=numeric(), p_value=numeric(),
                        population=character())
  cols_beginning <- c(17, 13, 14, 15, 16)
  add_to_cols <- c(0, 7, 14, 21, 28)

  for (i in 1:nrow(x)) {

    add <- x[i,c(1:5, 10)] # basic info

    for (j in 1:5) {

      cols <- as.vector(cols_beginning+add_to_cols[j])
      values <- x[i,..cols] # association results
      values <- cbind(add, values)
      values$pop <- str_split(colnames(values)[ncol(values)], "_")[[1]][2]

      results <- rbind(results, values, use.names=FALSE)

    }

  }

  # add CIs in log odds format
  results$lower <- log(results$lower)
  results$upper <- log(results$upper)

  # add allele names to rscodes
  results$rs_code[results$rs_code == "rs3911730"] <- "rs3911730 A"
  results$rs_code[results$rs_code == "rs8087187"] <- "rs8087187 C"
  results$rs_code[results$rs_code == "rs11585450"] <- "rs11585450 G"
  results$rs_code[results$rs_code == "rs1875763"] <- "rs1875763 C"

  return(results)
}


# separate discovery and test sets
results_disc <- get_tables(all_disc)
results_test <- get_tables(all_test)

# add column to stagger/colour train and test
results_disc$type <- "train"
results_test$type <- "test"

results <- rbind(results_disc, results_test)

# population names for plotting
pop_names <- c(
  `finns` = "Finland",
  `newcastle` = "UK",
  `katalonia` = "Spain",
  `poland` = "Poland",
  `all` = "Combined"
)
# create population name variable
results$pop_name <- plyr::mapvalues(results$population, from = names(pop_names), to = pop_names) %>% 
  factor(., levels = pop_names)


# HSCT assoc plot for relapse and cGvHD
results_relapse <- results[results$pheno == "relapse",]
results_cGvHD   <- results[results$pheno == "cGvHD_all",]
results_cGvHD$pheno <- 'cGvHD'



## -------------------------------------------------------------
## Cytotoxicity results data
## -------------------------------------------------------------

K562 <- fread("./results/NK_killing_activity_results/results_K562.txt")
K562$group <- 'Cytotoxicity'

K562_relapse <- K562[K562$SNP %like% "chr18_",]
K562_relapse$pheno <- 'relapse' 

K562_cGvHD <- K562[K562$SNP %like% "chr1_",]
K562_cGvHD$pheno <- 'cGvHD' 


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


# HSCT / relapse
p1 <- ggplot(results_relapse, aes(beta, rs_code))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = lower, xmax = upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 1.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th


# HSCT / cGvHD
p2 <- ggplot(results_cGvHD, aes(beta, rs_code))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = lower, xmax = upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 1.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th


# NK / relapse
p3 <- ggplot(K562_relapse, aes(beta, Rsid))  + 
  geom_point(size = 4) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1) + 
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  scale_y_discrete(limits = c("rs3911730 CC", "rs3911730 AC", "rs8087187 AA", "rs8087187 CA")) +
  ylab("") + xlab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 2.5,
           alpha = .25, fill = 'grey') +
  facet_grid(cols = vars(group), scales = 'free') +
  theme_minimal() + th

# NK / cGvHD
p4 <- ggplot(K562_cGvHD, aes(beta, Rsid))  + 
  geom_point(size = 4) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1) + 
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  scale_y_discrete(limits = c("rs11585450 CC", "rs11585450 GC", "rs11585450 GG", "rs1875763 GG", "rs1875763 CG", "rs1875763 CC")) +
  ylab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 3.5,
           alpha = .25, fill = 'grey') +
  facet_grid(cols = vars(group), scales = 'free') +
  theme_minimal() + th


# patch all plots together, save as png
p_all <- 
  (
    ( 
      (p1 / p2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom') 
    ) |
      ( (p3 / p4) )
  ) + 
  plot_layout(widths = c(1, 0.26)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))

ggsave("./results/for_paper/all_results.png", bg = "white", width = 14, height = 8, dpi = 600)





