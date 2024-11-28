library(ica)
library(MineICA)

for (r_file in list.files("tools/ica_latent_model/R/")) {
    source(paste0("tools/ica_latent_model/R/", r_file))
}

#--- load data
blood_counts <- readRDS("data/20230822_osf_repo/Blood normalised counts discovery//counts.rds")

blood_meta <- read.csv("data/20230822_osf_repo/PSORT-D_Blood_Clinical_Data_01-Apr-20.txt", sep = "\t")
rownames(blood_meta) <- blood_meta$Sample_id

# check meta lines up
blood_meta <- blood_meta[which(blood_meta$Sample_id %in% rownames(blood_counts)) ,]
stopifnot(all(rownames(blood_counts) == blood_meta$Sample_id))

#--- ica

# slow!
nc <- run_mstd(blood_counts, plot_path = "scripts/20230907_run_ica_both_tissues/mstd_blood.png", n_runs = 50, min_mean_stability = 0.85, max_components = 50)
nc <- 21  # Based on mstd plot

# run ICA in discovery
ica_blood <- run_ica(blood_counts, package = "ica", method = "imax", nc = nc)

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
saveRDS(mart, "mart.rds")
mart <- readRDS("mart.rds")

ica_blood <- create_ica_list(ica_blood, blood_meta, mart = mart)

saveRDS(ica_blood, "scripts/20230907_run_ica_both_tissues/ica_blood.rds")
