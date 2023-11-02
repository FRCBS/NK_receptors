library(data.table)
library(tidyverse)
library(readxl)

## Leenas's NK assay results data parsing
## 28 Oct 2022


## ---------------------------------------------------
## read excel data sheets
## Note: the sheet 'NK3' was manually edited to remove
## the first empty column
## ---------------------------------------------------

# get sheet names = sample names
sam <- excel_sheets('data/NK_assay_original/NK_Luc_Assays_cytotoxicity.xlsx')[3:16]
# read data from all donor samples
nk.dat <- map(sam, function(x) {
  read_xlsx('data/NK_assay_original/NK_Luc_Assays_cytotoxicity.xlsx', x) %>% data.frame
})
names(nk.dat) <- sam

# read CMV, KIR and HLA-C covariates
info <- read_xlsx('data/NK_assay_original/NK_Luc_Assays_cytotoxicity.xlsx', 1, range = "V3:Y18") %>% 
  na.omit %>% data.frame


## ---------------------------------------------------
## parse sheet contents
## ---------------------------------------------------

# omit NA values based on first column only
na.omit.tmp <- function(x) x[!is.na(x[, 1]), ]

# extract asay value data from sheets
nk.dat.parsed <- map(1:length(nk.dat), function(i) {
  print(i)
  # row where data begins & PRIESS cell line
  start.row <- which(grepl('PRIESS', nk.dat[[i]][, 1]))
  # find columns for the other two cell lines
  ssto <- which(grepl('SSTO', nk.dat[[i]][start.row, ]))
  k562 <- which(grepl('K562', nk.dat[[i]][start.row, ]))
  # extract dilution and measurement for PRIESS
  priess.val <- nk.dat[[i]][(start.row+1):(start.row+36), c(1, 6)] %>% na.omit.tmp
  # extract for the other two cell lines
  ssto.val <- nk.dat[[i]][(start.row+1):(start.row+36), c(ssto, ssto+5)] %>% na.omit.tmp
  k562.val <- nk.dat[[i]][(start.row+1):(start.row+36), c(k562, k562+5)] %>% na.omit.tmp
  # bind into data frame
  data.frame(DIL=priess.val[, 1],
             rbind(
               data.frame(VAL=priess.val[, 2], CELL='PRIESS'),
               data.frame(VAL=ssto.val[, 2], CELL='SSTO'),
               data.frame(VAL=k562.val[, 2], CELL='K562')
             ),
             SAMPLE=sam[i]
  ) %>% return
}) %>% do.call(rbind, .)

# join with covariates
nk.dat.parsed <- left_join(nk.dat.parsed, info, by=c('SAMPLE'='Donor.'))

# dilution to numeric
nk.dat.parsed$DIL_NUM <- nk.dat.parsed$DIL %>% str_split_fixed(., ':', 2) %>% .[, 1] %>% as.numeric

# edit col names
colnames(nk.dat.parsed)[c(5:6)] <- c('HLAC', 'KIR')

# write
fwrite(nk.dat.parsed, 'data/NK_assay_processed/NK_data_parsed.tsv', sep='\t')





