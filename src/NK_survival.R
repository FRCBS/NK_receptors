library(data.table)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(flextable)
library(cowplot)

# survival analysis for the four NK receptor SNPs

########################################################################################################################################

# add SNP dosage to survival data

########################################################################################################################################

# read in data
data <- fread("./results/patient_information/survival/HSCT_survival_data.txt")
# read in SNP dosage (made in ./src/pgen_dosage.sh)
# dosage <- fread("/home/nithiju/work/NK_receptors/results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage.raw") # pgen dosage, has values 0-2
dosage <- fread("/home/nithiju/work/NK_receptors/results/survival/dosage.raw") # plink1 dosage, has values 2, 1, 0

dosage_snps <- colnames(dosage)
dosage_snps <- str_sub(dosage_snps, 1, -3)

# filter dosage to only contain the wanted four SNPs
SNPs <- fread("/home/nithiju/work/NK_receptors/results/leena_donors/snps_names.txt", header = F)
# SNPs$V1 %in% colnames(dosage) # 0
# SNPs$V1 %in% dosage_snps # 4

cols <- c(T,T, dosage_snps[3:length(dosage_snps)] %in% SNPs$V1)
dosage <- dosage[,..cols]

dosage <- dosage[match(data$DonorGenotypingID, dosage$IID),]

# only leave in relapse SNPs
cols <- c(T,T,F,F,T,T)
dosage <- dosage[,..cols]

# train and test split
train <- fread("/home/nithiju/work/NK_receptors/results/ID_lists/train_sample_ids_donors.txt", header = F)
test <- fread("/home/nithiju/work/NK_receptors/results/ID_lists/test_sample_ids_donors.txt", header = F)

# add time in years (in stead of days)
data$overall_survival_years <- data$overall_survival / 365.25
data$relapse_free_survival_years <- data$relapse_free_survival / 365.25

########################################################################################################################################

# univariate analysis in the training set

########################################################################################################################################

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC", "Donor CC"))



RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/RFS_chr18_70203271_C_A_A.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA", "Donor AA"))

RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/RFS_chr18_70204107_A_C_C.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


###########################################################################################################################

# univariate analysis

# combine the small homozygous group with the heterozygous group

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC/CC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor AC/CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC/CC"))



RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/RFS_chr18_70203271_C_A_A_two_groups.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA/AA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor CA/AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA/AA"))



RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/RFS_chr18_70204107_A_C_C_two_groups.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

#######################################################################################################################################

# univariate analysis in the training set, two groups

# figures for the paper

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC/CC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor AC/CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC/CC"))

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA/AA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor CA/AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA/AA"))


# data$overall_survival_years <- data$overall_survival / 365.25
# data$relapse_free_survival_years <- data$relapse_free_survival / 365.25
# already done at the begining now

plot_1 <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

plot_2 <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

built_plot_1 <- ggsurvfit_build(plot_1)
built_plot_2 <- ggsurvfit_build(plot_2)

combined <- cowplot::plot_grid(built_plot_1, built_plot_2, ncol = 2, labels = c("A", "B"), label_size = 20)

ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/fig_3.png", dpi = 600, width = 20, height = 11)


########################################################################################################################################

# univariate analysis in the testing set

########################################################################################################################################

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC", "Donor CC"))



RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% test$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/test_set/RFS_chr18_70203271_C_A_A.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA", "Donor AA"))



RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% test$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Relapse-free survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/test_set/RFS_chr18_70204107_A_C_C.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


###########################################################################################################################

# univariate analysis

# combine the small homozygous group with the heterozygous group

# no individuals in the small heterozygous group in the test set plots above -> cannot combine them with anything

###########################################################################################################################

# multivariate analysis

###########################################################################################################################

# read in covariates
# using the ones made in the thymic SNP project
covars <- fread("./results/thymic_SNP/survival/survival_data_thymus_factors.txt")
# all(covars$DonorGenotypingID == data$DonorGenotypingID, na.rm = T) # T
data <- cbind(data, covars[,7:16])

# use only two groups for the SNPs
data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC/CC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor AC/CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC/CC"))

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA/AA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor CA/AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA/AA"))

colnames(data)[20] <- "HLA-match score"
colnames(data)[9] <- "rs8087187"
colnames(data)[10] <- "rs3911730"

data_train <- data[data$DonorGenotypingID %in% train$V2,]
data_test <- data[data$DonorGenotypingID %in% test$V2,]

# RFS for chr18_70203271_C_A_A
res_RFS <- coxph(Surv(relapse_free_survival, RFS_status) ~ ., data = data_train[,-c(1,2,3,6,7,8,10)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/training_set/cox_chr18_70203271_C_A_A_train.png")

# RFS for chr18_70204107_A_C_C
res_RFS <- coxph(Surv(relapse_free_survival, RFS_status) ~ ., data = data_train[,-c(1,2,3,6,7,8,9)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/training_set/cox_chr18_70204107_A_C_C_train.png")

# RFS for chr18_70203271_C_A_A
res_RFS <- coxph(Surv(relapse_free_survival, RFS_status) ~ ., data = data_test[,-c(1,2,3,6,7,8,10)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/test_set/cox_chr18_70203271_C_A_A_test.png")

# RFS for chr18_70204107_A_C_C
res_RFS <- coxph(Surv(relapse_free_survival, RFS_status) ~ ., data = data_test[,-c(1,2,3,6,7,8,9)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/test_set/cox_chr18_70204107_A_C_C_test.png")

# for al:
# Warning message:
#   In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#                  Loglik converged before variable  8 ; coefficient may be infinite. 


###############################################################################################################################
###############################################################################################################################

# for overall survival

###############################################################################################################################
###############################################################################################################################

# read in data
data <- fread("./results/patient_information/survival/HSCT_survival_data.txt")
# read in SNP dosage (made in ./src/pgen_dosage.sh)
# dosage <- fread("/home/nithiju/work/NK_receptors/results/pgen_NK_SNPs/all_datasets_plink_pgen_common_variants_filtered_ref_changed_dosage.raw") # pgen dosage, has values 0-2
dosage <- fread("/home/nithiju/work/NK_receptors/results/survival/dosage.raw") # plink1 dosage, has values 2, 1, 0

dosage_snps <- colnames(dosage)
dosage_snps <- str_sub(dosage_snps, 1, -3)

# filter dosage to only contain the wanted four SNPs
SNPs <- fread("/home/nithiju/work/NK_receptors/results/leena_donors/snps_names.txt", header = F)
# SNPs$V1 %in% colnames(dosage) # 0
# SNPs$V1 %in% dosage_snps # 4

cols <- c(T,T, dosage_snps[3:length(dosage_snps)] %in% SNPs$V1)
dosage <- dosage[,..cols]

dosage <- dosage[match(data$DonorGenotypingID, dosage$IID),]

# only leave in relapse SNPs
cols <- c(T,T,F,F,T,T)
dosage <- dosage[,..cols]

# train and test split
train <- fread("/home/nithiju/work/NK_receptors/results/ID_lists/train_sample_ids_donors.txt", header = F)
test <- fread("/home/nithiju/work/NK_receptors/results/ID_lists/test_sample_ids_donors.txt", header = F)

# add time in years (in stead of days)
data$overall_survival_years <- data$overall_survival / 365.25
data$relapse_free_survival_years <- data$relapse_free_survival / 365.25

########################################################################################################################################

# univariate analysis in the training set

########################################################################################################################################

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC", "Donor CC"))



RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/OS_chr18_70203271_C_A_A.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA", "Donor AA"))

RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/OS_chr18_70204107_A_C_C.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


###########################################################################################################################

# univariate analysis

# combine the small homozygous group with the heterozygous group

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC/CC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor AC/CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC/CC"))



RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/OS_chr18_70203271_C_A_A_two_groups.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA/AA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor CA/AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA/AA"))



RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% train$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/training_set/OS_chr18_70204107_A_C_C_two_groups.png", print(RFS_plot), dpi = 600, width = 10, height = 11)

########################################################################################################################################

# univariate analysis in the testing set

########################################################################################################################################

# for chr18_70203271_C_A_A

data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC", "Donor CC"))



RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70203271_C_A_A, data = data[data$DonorGenotypingID %in% test$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs8087187") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/test_set/OS_chr18_70203271_C_A_A.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


#--------------------------------------------------------------------------------------------------------

# the same for the other SNP, chr18_70204107_A_C_C

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA", "Donor AA"))



RFS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr18_70204107_A_C_C, data = data[data$DonorGenotypingID %in% test$V2,]) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 5, theme = theme_risktable_default(axis.text.y.size = 15, plot.title.size = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = "Overall survival rs3911730") +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))
ggsave(file = "/home/nithiju/work/NK_receptors/results/survival/test_set/OS_chr18_70204107_A_C_C.png", print(RFS_plot), dpi = 600, width = 10, height = 11)


###########################################################################################################################

# univariate analysis

# combine the small homozygous group with the heterozygous group

# no individuals in the small heterozygous group in the test set plots above -> cannot combine them with anything

###########################################################################################################################

# multivariate analysis

###########################################################################################################################

# read in covariates
# using the ones made in the thymic SNP project
covars <- fread("./results/thymic_SNP/survival/survival_data_thymus_factors.txt")
# all(covars$DonorGenotypingID == data$DonorGenotypingID, na.rm = T) # T
data <- cbind(data, covars[,7:16])

# use only two groups for the SNPs
data$chr18_70203271_C_A_A <- dosage$chr18_70203271_C_A_A
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 2] <- "Donor AA"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 1] <- "Donor AC/CC"
data$chr18_70203271_C_A_A[data$chr18_70203271_C_A_A == 0] <- "Donor AC/CC"
data$chr18_70203271_C_A_A <- factor(data$chr18_70203271_C_A_A, levels=c( "Donor AA", "Donor AC/CC"))

data$chr18_70204107_A_C_C <- dosage$chr18_70204107_A_C_C
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 2] <- "Donor CC"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 1] <- "Donor CA/AA"
data$chr18_70204107_A_C_C[data$chr18_70204107_A_C_C == 0] <- "Donor CA/AA"
data$chr18_70204107_A_C_C <- factor(data$chr18_70204107_A_C_C, levels=c( "Donor CC", "Donor CA/AA"))

colnames(data)[20] <- "HLA-match score"
colnames(data)[9] <- "rs8087187"
colnames(data)[10] <- "rs3911730"

data_train <- data[data$DonorGenotypingID %in% train$V2,]
data_test <- data[data$DonorGenotypingID %in% test$V2,]

# RFS for chr18_70203271_C_A_A
res_RFS <- coxph(Surv(overall_survival, OS_status) ~ ., data = data_train[,-c(1,2,4,5,7,8,10)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/training_set/cox_OS_chr18_70203271_C_A_A_train.png")

# RFS for chr18_70204107_A_C_C
res_RFS <- coxph(Surv(overall_survival, OS_status) ~ ., data = data_train[,-c(1,2,4,5,7,8,9)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/training_set/cox_OS_chr18_70204107_A_C_C_train.png")

# RFS for chr18_70203271_C_A_A
res_RFS <- coxph(Surv(overall_survival, OS_status) ~ ., data = data_test[,-c(1,2,4,5,7,8,10)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/test_set/cox_OS_chr18_70203271_C_A_A_test.png")

# RFS for chr18_70204107_A_C_C
res_RFS <- coxph(Surv(overall_survival, OS_status) ~ ., data = data_test[,-c(1,2,4,5,7,8,9)]) %>% 
  tbl_regression(exp = TRUE)  %>%
  as_flex_table() 
res_RFS <- bg(res_RFS, bg = "white", part = "all")
save_as_image(res_RFS, path = "/home/nithiju/work/NK_receptors/results/survival/test_set/cox_OS_chr18_70204107_A_C_C_test.png")

# for al:
# Warning message:
#   In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#                  Loglik converged before variable  8 ; coefficient may be infinite. 




