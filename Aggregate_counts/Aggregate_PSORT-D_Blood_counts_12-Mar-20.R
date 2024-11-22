# Load packages
library(tidyverse)
library(tximport)
library(mgsub)

# File paths to all abundance tables
file_path <- "F:/2019_09_Hs_PsO_PSORT_06-Jan-20/"
sample_dirs <- list.files(file_path)
files <- paste0(file_path,sample_dirs,"/abundance.tsv")

# Rename files so txi$counts has right column names
names(files) <- files
names(files) <- gsub("F:/2019_09_Hs_PsO_PSORT_06-Jan-20/", "", names(files))
names(files) <- gsub("/abundance.tsv", "", names(files))

# Load transcript to gene ID file
t2g <- readRDS('H:/my_documents/PhD/data/Data_12-Mar-20/other/Hs91.t2g.rds')

# Aggregate to gene level counts
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'lengthScaledTPM', ignoreTxVersion = T)

# Save txi object
saveRDS(txi, 
        paste0("H:/my_documents/PhD/data/Data_12-Mar-20/txi_objects/PSORT-D_Blood_txi_", 
               mgsub(as.character(format(Sys.time(), "%d %b %Y %X")), c(" ",":"), c("-", "-")),".rds"))

# Save cnts
saveRDS(txi$counts, 
        paste0("H:/my_documents/PhD/data/Data_12-Mar-20/gene_level_counts/PSORT-D_Blood_counts_", 
               mgsub(as.character(format(Sys.time(), "%d %b %Y %X")), c(" ",":"), c("-", "-")),".rds"))
