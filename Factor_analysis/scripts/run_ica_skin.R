library(ica)
library(MineICA)

for (r_file in list.files("tools/ica_latent_model/R/")) {
    source(paste0("tools/ica_latent_model/R/", r_file))
}

#--- load data
skin_d_counts <- readRDS("data/20230822_osf_repo/Skin normalised counts discovery/counts.rds")
skin_r_counts <- readRDS("data/20230822_osf_repo/Skin normalised counts replication/counts_r.rds")

skin_d_counts <- skin_d_counts[, which(colnames(skin_d_counts) %in% colnames(skin_r_counts))]
skin_r_counts <- skin_r_counts[, which(colnames(skin_r_counts) %in% colnames(skin_d_counts))]
stopifnot(all(colnames(skin_d_counts) == colnames(skin_r_counts)))

skin_d_meta <- read.csv("data/20230822_osf_repo/PSORT-D_Skin_Clinical_Data_01-Apr-20.txt", sep = "\t")
skin_r_meta <- read.csv("data/20230822_osf_repo/PSORT-R_Skin_Clinical_Data_01-Apr-20.txt", sep = "\t")
skin_r_meta <- skin_r_meta[which(skin_r_meta$Sample_id %in% rownames(skin_r_counts)) ,]

rownames(skin_d_meta) <- skin_d_meta$Sample_id
rownames(skin_r_meta) <- skin_r_meta$Sample_id

# check meta lines up
stopifnot(all(rownames(skin_d_counts) == skin_d_meta$Sample_id))
stopifnot(all(rownames(skin_r_counts) == skin_r_meta$Sample_id))

# check meta lines up
stopifnot(all(rownames(skin_d_counts) == skin_d_meta$Sample_id))
stopifnot(all(rownames(skin_r_counts) == skin_r_meta$Sample_id))
#--- ica

# slow!
nc <- run_mstd(skin_d_counts, plot_path = "scripts/20230907_run_ica_both_tissues/mstd_skin.png", n_runs = 50, min_mean_stability = 0.90, max_components = 50)
nc <- 24  # Based on mstd plot

# run ICA in discovery then project in replication
ica_skin_d <- run_ica(skin_d_counts, package = "ica", method = "imax", nc = nc)
ica_skin_r <- project_ica(skin_r_counts, ica_skin_d)

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
saveRDS(mart, "mart.rds")
mart <- readRDS("mart.rds")

ica_skin_d <- create_ica_list(ica_skin_d, skin_d_meta, mart = mart)
ica_skin_r <- create_ica_list(ica_skin_r, skin_r_meta, mart = mart)

saveRDS(ica_skin_d, "scripts/20230907_run_ica_both_tissues/ica_skin_d.rds")
saveRDS(ica_skin_r, "scripts/20230907_run_ica_both_tissues/ica_skin_r.rds")
