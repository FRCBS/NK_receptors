library(data.table)
library(tidyverse)

# read in association results

relapse_train_all <- fread("/home/nihtiju/work/NK_receptors/results/association_testing_train/all_train/association_all_train_relapse.pheno_relapse.glm.logistic.hybrid")
relapse_train_all_chr18_70203271_C_A <- relapse_train_all[relapse_train_all$ID == "chr18_70203271_C_A",]
relapse_train_all_chr18_70204107_A_C <- relapse_train_all[relapse_train_all$ID == "chr18_70204107_A_C",]

cGvHD_train_all <- fread("/home/nihtiju/work/NK_receptors/results/association_testing_train/all_train/association_all_train_cGvHD.pheno_cGvHD_all.glm.logistic.hybrid")
cGvHD_train_all_chr1_161726992_G_C <- cGvHD_train_all[cGvHD_train_all$ID == "chr1_161726992_G_C",]
cGvHD_train_all_chr1_161727208_C_G <- cGvHD_train_all[cGvHD_train_all$ID == "chr1_161727208_C_G",]

write.table(relapse_train_all_chr18_70203271_C_A, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/relapse_train_all_populations_chr18_70203271_C_A.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(relapse_train_all_chr18_70204107_A_C, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/relapse_train_all_populations_cr18_70204107_A_C.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(cGvHD_train_all_chr1_161726992_G_C, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/cGvHD_train_all_populations_chr1_161726992_G_C.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(cGvHD_train_all_chr1_161727208_C_G, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/cGvHD_train_all_populations_chr1_161727208_C_G.txt", sep = "\t", quote = F, row.names = F, col.names = T)




relapse_test_all <- fread("/home/nihtiju/work/NK_receptors/results/association_testing_test/all_test/association_all_test_relapse.pheno_relapse.glm.logistic.hybrid")
relapse_test_all_chr18_70203271_C_A <- relapse_test_all[relapse_test_all$ID == "chr18_70203271_C_A",]
relapse_test_all_chr18_70204107_A_C <- relapse_test_all[relapse_test_all$ID == "chr18_70204107_A_C",]

cGvHD_test_all <- fread("/home/nihtiju/work/NK_receptors/results/association_testing_test/all_test/association_all_test_cGvHD.pheno_cGvHD_all.glm.logistic.hybrid")
cGvHD_test_all_chr1_161726992_G_C <- cGvHD_test_all[cGvHD_test_all$ID == "chr1_161726992_G_C",]
cGvHD_test_all_chr1_161727208_C_G <- cGvHD_test_all[cGvHD_test_all$ID == "chr1_161727208_C_G",]


write.table(relapse_test_all_chr18_70203271_C_A, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/relapse_test_all_populations_chr18_70203271_C_A.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(relapse_test_all_chr18_70204107_A_C, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/relapse_test_all_populations_chr18_70204107_A_C.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(cGvHD_test_all_chr1_161726992_G_C, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/cGvHD_test_all_populations_chr1_161726992_G_C.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(cGvHD_test_all_chr1_161727208_C_G, "/home/nihtiju/work/NK_receptors/results/for_paper/covariate_results/cGvHD_test_all_populations_chr1_161727208_C_G.txt", sep = "\t", quote = F, row.names = F, col.names = T)



















