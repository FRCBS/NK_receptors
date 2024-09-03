library(data.table)
library(tidyverse)
library(cowplot)
library(patchwork)

setwd("/home/nihtiju/work/NK_receptors/")

############################################################################################################################################################

## -------------------------------------------------------------
## HSCT results data
## -------------------------------------------------------------

# read discovery and test set association results

all_disc <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")[1:13,]
all_test <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")[14:26,] # ordered to match training set results already

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[match(all_disc$ID, SNP_info$ID),]

get_tables <- function(x) {
  
  results <- data.table(gene=character(),
                        rs_code=character(),
                        name_in_ours=character(),
                        chr=numeric(), pos_hg38=numeric(),
                        pheno=character(), or=numeric(), 
                        lower=numeric(), upper=numeric(), 
                        p_value=numeric(),
                        population=character())
  cols_beginning <- c(8, 9, 10, 11) # or, l95, u95, p
  add_to_cols <- c(0, 5, 10, 15, 20) # how many to add to get the same values for all populations
  
  for (i in 1:nrow(x)) {
    
    add <- cbind(SNP_info[i,c(7,8,3,1,2)], x[i,1]) # basic info
    colnames(add) <- c("gene", "rs_code", "name_in_ours", "chr", "pos_hg38", "pheno")
    
    for (j in 1:5) {
      
      cols <- as.vector(cols_beginning + add_to_cols[j])
      values <- x[i,..cols] # association results
      values <- cbind(add, values)
      values$pop <- str_split(colnames(values)[ncol(values)], "_")[[1]][2]
      
      results <- rbind(results, values, use.names=FALSE)
      
    }
    
  }
  
  # add allele names to rscodes
  SNP_info_2 <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  SNP_info_2 <- SNP_info_2[match(results$name_in_ours, SNP_info_2$ID),]
  
  results$rs_code_allele <- paste0(results$rs_code, " ", SNP_info_2$A1)
  
  return(results)
}


# separate discovery and test sets
results_disc <- get_tables(all_disc)
results_test <- get_tables(all_test)

# add column to stagger/colour train and test
results_disc$type <- "train"
results_test$type <- "test"

results <- rbind(results_disc, results_test)

# add beta and its 95% CI
results$beta <- log(results$or)
results$beta_lower <- log(results$lower) # to log odds format
results$beta_upper <- log(results$upper)




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

# modify pheno names (these will be visible in the plots)
results$pheno[results$pheno == "aGvHD_all"] <- "aGvHD 1"
results$pheno[results$pheno == "aGvHD_severe"] <- "aGvHD 2"
results$pheno[results$pheno == "cGvHD_all"] <- "cGvHD 1"
results$pheno[results$pheno == "cGvHD_severe"] <- "cGvHD 2"
results$pheno[results$pheno == "cGvHD_severe_broader"] <- "cGvHD 3"
results$pheno[results$pheno == "relapse"] <- "Relapse"

# HSCT assoc plot for relapse and cGvHD
results_aGvHD_all <- results[results$pheno == "aGvHD 1",]
results_aGvHD_severe <- results[results$pheno == "aGvHD 2",]
results_cGvHD_all   <- results[results$pheno == "cGvHD 1",]
results_cGvHD_severe   <- results[results$pheno == "cGvHD 2",]
results_cGvHD_severe_broader   <- results[results$pheno == "cGvHD 3",]
results_relapse <- results[results$pheno == "Relapse",]

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

# HSCT / aGvHD_all
p_aGvHD_all <- ggplot(results_aGvHD_all, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
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

# HSCT / aGvHD_severe
p_aGvHD_severe <- ggplot(results_aGvHD_severe, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_all
p_cGvHD_all <- ggplot(results_cGvHD_all, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 1.5,
           alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5 , ymax = 3.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_severe
p_cGvHD_severe <- ggplot(results_cGvHD_severe, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_severe_broader
p_cGvHD_severe_broader <- ggplot(results_cGvHD_severe_broader, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
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


# HSCT / relapse
p_relapse <- ggplot(results_relapse, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("beta") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th



# patch all plots together, save as png
p_all <- 
  (p_aGvHD_all / p_aGvHD_severe / p_cGvHD_all / p_cGvHD_severe / p_cGvHD_severe_broader / p_relapse) + plot_layout(guides = 'collect', axis_titles = "collect") & theme(legend.position = 'bottom')

ggsave("./results/for_paper/association/association_results_plot_beta.png", bg = "white", width = 10, height = 15, dpi = 600)

############################################################################################################################################################

## -------------------------------------------------------------
## the same but errorcode == unfinished results removed
# 'UNFINISHED': Logistic/Firth regression didn't fail in an obvious manner, but the result didn't satisfy the usual convergence criteria when the iteration limit was hit. (The result is still reported in this case, but it's less accurate than usual.)
## -------------------------------------------------------------

# read discovery and test set association results

all_disc <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")[1:13,]
all_test <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/results_chosen_SNPs_train_and_test_ordered.csv")[14:26,] # ordered to match training set results already

# replace the "unfinished" results wit NA
replace_na <- function(table){
  
  cols <- 1:ncol(table)
  for (i in cols[colnames(table) %like% "error"]) {
    
    for (j in 1:nrow(table)) {
      
      if(is.na(table[j, ..i])){
        # table[j, ..i] <- "." # did not work with this
        table[j, i] <- "."
      }
      
      if(table[j, ..i] == "UNFINISHED"){

        for (a in 1:4) {

          col <- i - a
          table[j, col] <- NA
          
        }
        
      }
      
    }
    
  }
  
  return(table)
  
}

all_disc <- replace_na(all_disc)
all_test <- replace_na(all_test)

SNP_info <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
SNP_info <- SNP_info[match(all_disc$ID, SNP_info$ID),]

get_tables <- function(x) {
  
  results <- data.table(gene=character(),
                        rs_code=character(),
                        name_in_ours=character(),
                        chr=numeric(), pos_hg38=numeric(),
                        pheno=character(), or=numeric(), 
                        lower=numeric(), upper=numeric(), 
                        p_value=numeric(),
                        population=character())
  cols_beginning <- c(8, 9, 10, 11) # or, l95, u95, p
  add_to_cols <- c(0, 5, 10, 15, 20) # how many to add to get the same values for all populations
  
  for (i in 1:nrow(x)) {
    
    add <- cbind(SNP_info[i,c(7,8,3,1,2)], x[i,1]) # basic info
    colnames(add) <- c("gene", "rs_code", "name_in_ours", "chr", "pos_hg38", "pheno")
    
    for (j in 1:5) {
      
      cols <- as.vector(cols_beginning + add_to_cols[j])
      values <- x[i,..cols] # association results
      values <- cbind(add, values)
      values$pop <- str_split(colnames(values)[ncol(values)], "_")[[1]][2]
      
      results <- rbind(results, values, use.names=FALSE)
      
    }
    
  }
  
  # add allele names to rscodes
  SNP_info_2 <- fread("/home/nihtiju/work/NK_receptors/results/association_LD_newOnly/SNPs_wanted_info.csv")
  SNP_info_2 <- SNP_info_2[match(results$name_in_ours, SNP_info_2$ID),]
  
  results$rs_code_allele <- paste0(results$rs_code, " ", SNP_info_2$A1)
  
  return(results)
}


# separate discovery and test sets
results_disc <- get_tables(all_disc)
results_test <- get_tables(all_test)

# add column to stagger/colour train and test
results_disc$type <- "train"
results_test$type <- "test"

results <- rbind(results_disc, results_test)

# add beta and its 95% CI
results$beta <- log(results$or)
results$beta_lower <- log(results$lower) # to log odds format
results$beta_upper <- log(results$upper)




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

# modify pheno names (these will be visible in the plots)
results$pheno[results$pheno == "aGvHD_all"] <- "aGvHD 1"
results$pheno[results$pheno == "aGvHD_severe"] <- "aGvHD 2"
results$pheno[results$pheno == "cGvHD_all"] <- "cGvHD 1"
results$pheno[results$pheno == "cGvHD_severe"] <- "cGvHD 2"
results$pheno[results$pheno == "cGvHD_severe_broader"] <- "cGvHD 3"
results$pheno[results$pheno == "relapse"] <- "Relapse"

# HSCT assoc plot for relapse and cGvHD
results_aGvHD_all <- results[results$pheno == "aGvHD 1",]
results_aGvHD_severe <- results[results$pheno == "aGvHD 2",]
results_cGvHD_all   <- results[results$pheno == "cGvHD 1",]
results_cGvHD_severe   <- results[results$pheno == "cGvHD 2",]
results_cGvHD_severe_broader   <- results[results$pheno == "cGvHD 3",]
results_relapse <- results[results$pheno == "Relapse",]

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

# HSCT / aGvHD_all
p_aGvHD_all <- ggplot(results_aGvHD_all, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
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

# HSCT / aGvHD_severe
p_aGvHD_severe <- ggplot(results_aGvHD_severe, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_all
p_cGvHD_all <- ggplot(results_cGvHD_all, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 1.5,
           alpha = .25, fill = 'grey') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5 , ymax = 3.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_severe
p_cGvHD_severe <- ggplot(results_cGvHD_severe, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("") +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th

# HSCT / cGvHD_severe_broader
p_cGvHD_severe_broader <- ggplot(results_cGvHD_severe_broader, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
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


# HSCT / relapse
p_relapse <- ggplot(results_relapse, aes(beta, rs_code_allele))  +
  geom_point(aes(col = type), size = 4, position = position_dodge(width = 0.4)) +
  geom_linerange(aes(xmin = beta_lower, xmax = beta_upper, color = type), 
                 linewidth = 1, position = position_dodge(width = 0.4)) + 
  scale_color_manual(values=c("#CB2314", "#046C9A"), labels=c("Discovery", "Replication"), 
                     name = "", breaks = c("train", "test")) +
  scale_x_continuous(guide = guide_axis(check.overlap = T), 
                     expand = c(.1, .1)) +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("beta") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5,
           alpha = .25, fill = 'grey') +
  facet_grid(pheno ~ pop_name, scales = 'free') +
  theme_minimal() + th



# patch all plots together, save as png
p_all <- 
  (p_aGvHD_all / p_aGvHD_severe / p_cGvHD_all / p_cGvHD_severe / p_cGvHD_severe_broader / p_relapse) + plot_layout(guides = 'collect', axis_titles = "collect") & theme(legend.position = 'bottom')

ggsave("./results/for_paper/association/association_results_plot_beta_unfinished_removed.png", bg = "white", width = 10, height = 15, dpi = 600)





