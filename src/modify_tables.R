library(data.table)
library(tidyverse)
library(flextable)

#--------------------------------------------------------------------------------------------------------------------------------------

# table 1 for train and test set together


# column names
table_train <- fread('./results/for_paper/table_1.txt', data.table=F, header = T)
table_test <- fread('./results/for_paper/table_1_test.txt', data.table=F, header = T)

# test table has one extra row for aGvHD grade unknown
# add that row to train
table_train <- rbind(table_train[1:27,], c("aGvHD class not known", 0, 0, 0, 0), table_train[28:nrow(table_train),])

# get the part with no diseases
table_train <- table_train[1:37,]
table_test <- table_test[1:37,]

# join tables together
# also order the populations to be from the biggest to the smallest
all <- cbind(table_train[,1:2], table_test$Finland, table_train$UK, table_test$UK, table_train$Spain, table_test$Spain, table_train$Poland, table_test$Poland)

colnames(all) <- c("population", "Finland_train", "Finland_test", "UK_train", "UK_test", "Spain_train", "Spain_test", "Poland_train", "Poland_test")

all <- all[,c(1,1,2:9)]

#------------------
# diseases

# disease names harmonized for the training set in
diagnosis_train <- fread("./results/for_paper/supple_table_diseases_for_translating_names_harmonized.txt", data.table=F, header = T, stringsAsFactors = F) # made by hand

diagnosis_test <- fread("./results/for_paper/supple_table_diseases_for_translating_names_test_harmonized.txt", data.table=F, header = T, stringsAsFactors = F)  # made by hand
# NAs as "diagnosis missing"

diagnosis_test_new <- diagnosis_test[!(diagnosis_test$harmonized %in% diagnosis_train$`translated disease name`),]
diagnosis_test_ordered <- diagnosis_test[match(diagnosis_train$`translated disease name`, diagnosis_test$harmonized),]

rows <- 1:nrow(diagnosis_test_ordered)
rows <- rows[is.na(diagnosis_test_ordered$harmonized)]

for (i in rows) {
  
  print(i)
  diagnosis_test_ordered[i,4:7] <- 0
  
}

diagnosis_train_all <- cbind(diagnosis_train[,c(3,3,4)], diagnosis_test_ordered$Finland, diagnosis_train$UK, diagnosis_test_ordered$UK, diagnosis_train$Spain, diagnosis_test_ordered$Spain, diagnosis_train$Poland, diagnosis_test_ordered$Poland)

diagnosis_test_all <- cbind(diagnosis_test_new[,c(3,3,4)], c(0,0,0,0), diagnosis_test_new$UK, c(0,0,0,0), diagnosis_test_new$Spain, c(0,0,0,0), diagnosis_test_new$Poland, c(0,0,0,0))


colnames(diagnosis_train_all) <- colnames(all)
colnames(diagnosis_test_all) <- colnames(all)

diagnosis_all <- rbind(diagnosis_train_all, diagnosis_test_all)

diagnosis_all$all_sum <- rowSums(diagnosis_all[,3:ncol(diagnosis_all)])
diagnosis_all <- diagnosis_all[order(diagnosis_all$all_sum, decreasing = T),]


# leave only  the top 5
table_diseases_5 <- diagnosis_all[1:5,]

# get the rest as "other"
other <- diagnosis_all[6:nrow(diagnosis_all),]
# add to top 5
table_diseases_5 <- rbind(table_diseases_5, table_diseases_5[1,]) # add one row, modify the name to be other next
table_diseases_5[6,1:2] <- "Other"
table_diseases_5[6,3:10] <- colSums(other[,3:10]) # for individual populations

# remove the sum of all populations & the column for harmonizing disease names
table_diseases_5 <- table_diseases_5[,c(1:10)]


#------------------
table <- rbind(all, table_diseases_5)
#------------------


first <- c("Number of HSCT donors, n", "HSCT time, years", "HSCT time, years", "Recipient age in years, median (range)", "Recipient age in years, median (range)", "Donor age in years, median (range)", "Donor age in years, median (range)","Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)","Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)","Donor type, n (%)", "Donor type, n (%)", "Donor type, n (%)", "Donor type, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "Relapse, n (%)", "Relapse, n (%)", "Relapse, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", "Diagnosis, n (%)", " ") # the very last one as " " to get the bottom border visible ok
second <- c("", "", "Missing, n (%)", "", "Missing, n (%)", "", "Missing, n (%)","Male-male", "Male-female", "Female-male", "Female-female", "Missing","Peripheral blood", "Bone marrow", "Both", "Missing", "Sibling","Register", "Haplo", "Missing", "Myeloablative", "Reduced intensity", "Other", "Missing", "grade 0", paste0("grade ", as.roman(1), "-", as.roman(2)), paste0("grade ", as.roman(3), "-", as.roman(4)), "Grade unknown", "Missing", "grade 0", "Yes, classification unknown", "Limited", "Extensive", "Missing", "Yes", "No", "Missing", "Acute myeloid leukemia", "Acute lymphoblastic leukemia", "Myelodysplastic syndrome", "Chronic myeloid leukemia", "Multiple myeloma", "Other")

table[,1] <- first
table[,2] <- second

# colnames(table)[1] <- " " 
# colnames(table)[2] <- "  " 

colnames(table) <- c(" ", "  ", "Finland", "    ", "UK", "     ", "Spain", "      ", "Poland", "       ")

# add %
rows <- c(3,5,7:nrow(table))

for (i in rows) {
  
  for (j in 3:10) {
    
    # i = row
    # j = col
    
    value <- (as.numeric(table[i,j]) / as.numeric(table[1,j])) * 100
    
    pasted <- paste0(table[i,j], " (", round(value), ")")
    
    table[i,j] <- pasted
    
    
    
  }
  
}

# headers
table_flex <- flextable(table)
# table_flex <- delete_part(table_flex, part = "header")
table_flex <- add_header_row(table_flex, values = c("", "", rep(c("Discovery", "Replication"), 4)), top = F)

table_flex <- align(table_flex, part = "header", align = "left")

# merge vertical duplicated names
table_flex <- merge_v(table_flex, j = c(" ", "  "))
table_flex <- valign(table_flex, j = 1, valign = "top")

# add footnote
table_flex <- add_footer_lines(table_flex, "GvHD, graft-versus-host disease; aGvHD, acute GvHD; cGvHD, chronic GvHD")
table_flex <- footnote(table_flex, i = c(5,7), j = 2,
                       value = as_paragraph(
                         c("Missing ages were imputed, see Materials and methods")
                       ),
                       ref_symbols = c("a"),
                       part = "body")
table_flex <- footnote(table_flex, i = 38, j = 1,
                       value = as_paragraph(
                         c("Five most frequent diagnoses are presented here")
                       ),
                       ref_symbols = c("b"),
                       part = "body")

table_flex_onePage <- table_flex
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)

set_flextable_defaults(background.color = "white")
save_as_image(table_flex, path = "./results/for_paper/flextables/table_1_train_test.png", bg = "white")

# paginate(table_flex, init = TRUE, hdr_ftr = TRUE, group = " ")
save_as_docx(table_flex, path = "./results/for_paper/flextables/table_1_train_test_docx.docx", align = "left")
# and changing thd document to be in landscape format manually after saving the table
# -> wide enough to see all columns

save_as_docx(table_flex_onePage, path = "./results/for_paper/flextables/table_1_train_test_docx_onePage.docx", align = "left")

#--------------------------------------------------------------------------------------------------------------------------------------

# supple table 1: all NK receptor SNPs

table <- fread('./results/for_paper/supple_table_all_snps.txt', data.table=F) 

# remove one row with no rs code 
table <- table[!(is.na(table$rs_code)),]

# remove the name in our datasets
table <- table[,-3]

# rename the columns
colnames(table) <- c("Gene", "Rsid", "Chr", "Position", "Finland", "UK", "Spain", "Poland")

# make one masisve table from these
table_flex <- flextable(table)
table_flex <- align(table_flex, part = "header", align = "left")
table_flex <- footnote( table_flex, i = 1, j = 4,
                        value = as_paragraph(
                          c("Position in hg38")
                        ),
                        ref_symbols = c("a"),
                        part = "header")
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
set_flextable_defaults(background.color = "white")

save_as_docx(table_flex, path = "./results/for_paper/flextables/supple_table_all.docx", align = "left")


# divide these into separate tables to be on each page
tables <- list()
all <- seq(1, nrow(table), by=60)
for (i in 1:length(seq(1, nrow(table), by=60))) {
  
  beginning <- all[i]
  
  if (i == length(all)) {
    end <- nrow(table)
  } else {
    end <- all[i+1]-1
  }
  
  tables[[i]] <-table[beginning:end,]
  
}


# add formatting

for (i in 1:length(tables)) {
  
  table <- tables[[i]]
  
  table_flex <- flextable(table)
  table_flex <- align(table_flex, part = "header", align = "left")
  table_flex <- footnote( table_flex, i = 1, j = 4,
                          value = as_paragraph(
                            c("Position in hg38")
                          ),
                          ref_symbols = c("a"),
                          part = "header")
  table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
  
  save_as_image(table_flex, path = paste0("./results/for_paper/flextables/supple_table_1_", i, ".png"), bg = "white")
  
  save_as_docx(table_flex, path = paste0("./results/for_paper/flextables/supple_table_1_", i, ".docx"), align = "left")
  
  
}







