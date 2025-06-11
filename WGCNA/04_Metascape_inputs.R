#' ---
#' title: "Format inputs for Metascape"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 12, fig.height = 9)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' Here, we format the inputs for [Metascape](https://metascape.org/gp/index.html#/main/step1)
#' ([Zhou et al, 2019](https://www.nature.com/articles/s41467-019-09234-6)),
#' which will be used to carry out functional enrichment analysis of the WGCNA module and latent factor genes.
#' 
#' # Load packages

library(tidyverse)

#' # Create output directory

output_directory <- "results/WGCNA/04_Metascape_inputs"
dir.create(output_directory)

#' # Skin
#' 
#' ## Load data

# Module genes
modules <- read.delim("results/WGCNA/01_Module_identification/Skin/modules.txt")

# Factor genes
factors <- read.delim("data/latent_factors/factor_genes_skin.csv", sep = ",")
colnames(factors) <- c("X", "Factor", "EnsemblID", "Coef", "GeneSymbol")

# Enrichment background for factor genes
factor_background <- read.delim("data/latent_factors/enrichment_background.csv", sep = ",")
factor_background <- factor_background[[1]]

# Module-trait correlations
module_traits <- read.delim(
  "results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/Module-trait_correlations.txt"
)

#' ## Create Metascape input
#' 
#' ### Modules

genes <- vector(mode = "list", length = length(unique(modules$Module)))
names(genes) <- sort(unique(modules$Module))
for(i in 1:length(genes)){
  genes[[i]] <- modules %>%
    filter(Module == names(genes)[i]) %>%
    pull(EnsemblID) %>%
    paste(collapse = ",")
}
genes$BACKGROUND <- paste(modules$EnsemblID, collapse = ",")
genes <- stack(genes)
genes <- genes[,c(2,1)]
colnames(genes) <- c("Name", "Genes")
genes$Name <- as.character(genes$Name)
genes$Name[which(genes$Name == "BACKGROUND")] <- "_BACKGROUND"

module_traits <- module_traits %>% filter(padj_all_d <= 0.05, Module != "grey")

genes <- genes %>% dplyr::filter(Name %in% c(unique(module_traits$Module),"_BACKGROUND"))

write.table(genes, paste0(output_directory,"/Skin_modules_padj_all_0.05_in_disc.txt"), 
            sep = "\t", row.names = F, quote = F)

#' ### Factors

genes <- vector(mode = "list", length = length(unique(factors$Factor)))
names(genes) <- sort(unique(factors$Factor))
for(i in 1:length(genes)){
  genes[[i]] <- factors %>%
    filter(Factor == names(genes)[i]) %>%
    pull(EnsemblID) %>%
    paste(collapse = ",")
}
genes$BACKGROUND <- paste(factor_background, collapse = ",")
genes <- stack(genes)
genes <- genes[,c(2,1)]
colnames(genes) <- c("Name", "Genes")
genes$Name <- as.character(genes$Name)
genes$Name[which(genes$Name == "BACKGROUND")] <- "_BACKGROUND"

write.table(genes, paste0(output_directory,"/Skin_factors.txt"), sep = "\t", row.names = F, quote = F)

#' # Blood
#' 
#' ## Load data

# Module genes
modules <- read.delim("results/WGCNA/01_Module_identification/Blood/modules.txt")

# Factor genes
factors <- read.delim("data/latent_factors/factor_genes_blood.csv", sep = ",")
colnames(factors) <- c("X", "Factor", "EnsemblID", "Coef", "GeneSymbol")

#' ## Create Metascape input
#' 
#' ### Modules

genes <- vector(mode = "list", length = length(unique(modules$Module)))
names(genes) <- sort(unique(modules$Module))
for(i in 1:length(genes)){
  genes[[i]] <- modules %>%
    filter(Module == names(genes)[i]) %>%
    pull(EnsemblID) %>%
    paste(collapse = ",")
}
genes$BACKGROUND <- paste(modules$EnsemblID, collapse = ",")
genes <- stack(genes)
genes <- genes[,c(2,1)]
colnames(genes) <- c("Name", "Genes")
genes$Name <- as.character(genes$Name)
genes$Name[which(genes$Name == "BACKGROUND")] <- "_BACKGROUND"

genes <- genes %>% dplyr::filter(!Name %in% "grey")

write.table(genes, paste0(output_directory,"/Blood_modules.txt"), sep = "\t", row.names = F, quote = F)

#' ### Factors

genes <- vector(mode = "list", length = length(unique(factors$Factor)))
names(genes) <- sort(unique(factors$Factor))
for(i in 1:length(genes)){
  genes[[i]] <- factors %>%
    filter(Factor == names(genes)[i]) %>%
    pull(EnsemblID) %>%
    paste(collapse = ",")
}
genes$BACKGROUND <- paste(factor_background, collapse = ",")
genes <- stack(genes)
genes <- genes[,c(2,1)]
colnames(genes) <- c("Name", "Genes")
genes$Name <- as.character(genes$Name)
genes$Name[which(genes$Name == "BACKGROUND")] <- "_BACKGROUND"

write.table(genes, paste0(output_directory,"/Blood_factors.txt"), sep = "\t", row.names = F, quote = F)

#' # Session information

sessionInfo()