# Load packages
library(tidyverse)
library(tximport)
library(mgsub)

# Load PSORT-R fastq filenames and remove last 2 ("filenames.txt", "psort_index")
fastq_filenames <- read.delim("H:/my_documents/PhD/data/Data_12-Mar-20/PSORT-R_abundances_from_rocket_12-Mar-20/filenames.txt", header = F)
fastq_filenames <- fastq_filenames %>% dplyr::filter(V1 != "filenames.txt", V1 != "psort_index")
colnames(fastq_filenames) <- "File_name"

# Make column of number prefixes and sample IDs
fastq_filenames$Number_prefix <- gsub("_.*", "", fastq_filenames$File_name)
fastq_filenames$Sample_id <- gsub("^\\d{1,3}_{1}", "", fastq_filenames$File_name)
fastq_filenames$Sample_id <- gsub("_.*$", "", fastq_filenames$Sample_id)
fastq_filenames$Sample_id <- sub("(....)$", "-\\1", fastq_filenames$Sample_id)
fastq_filenames$Sample_id <- gsub("^[0]+", "", fastq_filenames$Sample_id)

# There are two read files oer sample - remove duplication
fastq_filenames <- fastq_filenames[grep("R1", fastq_filenames$File_name),]

# File paths to all abundance tables
files <- paste0("H:/my_documents/PhD/data/Data_12-Mar-20/PSORT-R_abundances_from_rocket_12-Mar-20/abundance_all/",
                list.files("H:/my_documents/PhD/data/Data_12-Mar-20/PSORT-R_abundances_from_rocket_12-Mar-20/abundance_all/"))

# Rename files so txi$counts has right column names
names(files) <- files
names(files) <- gsub("H:/my_documents/PhD/data/Data_12-Mar-20/PSORT-R_abundances_from_rocket_12-Mar-20/abundance_all/kallisto_out_", "", names(files))
names(files) <- gsub("_abundance.tsv", "", names(files))
for(name in names(files)){
  name_2 <- fastq_filenames$Sample_id[fastq_filenames$Number_prefix == name]
  names(files)[names(files) == name] <- name_2
}

# Load transcript to gene ID file
t2g <- readRDS('H:/my_documents/PhD/data/Data_12-Mar-20/other/Hs91.t2g.rds')

# Aggregate to gene level counts
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, countsFromAbundance = 'lengthScaledTPM', ignoreTxVersion = T)

# Save txi object
saveRDS(txi, 
        paste0("H:/my_documents/PhD/data/Data_12-Mar-20/txi_objects/PSORT-R_Skin_txi_", 
               mgsub(as.character(format(Sys.time(), "%d %b %Y %X")), c(" ",":"), c("-", "-")),".rds"))

# Save cnts
saveRDS(txi$counts, 
        paste0("H:/my_documents/PhD/data/Data_12-Mar-20/gene_level_counts/PSORT-R_Skin_counts_", 
               mgsub(as.character(format(Sys.time(), "%d %b %Y %X")), c(" ",":"), c("-", "-")),".rds"))

