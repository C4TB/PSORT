#' ---
#' title: "Gene-level trait correlation plots"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 23, fig.height = 16)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' # Preliminaries
#' 
#' ## Load packages

library(tidyverse)
library(gridExtra)
library(ggh4x)
library(scales)
library(splines)
library(reshape2)
library(limma)
library(edgeR)

#' ## Create output directory

output_directory <- "results/WGCNA/08_Gene-level_trait_correlation_plots"
dir.create(output_directory)

#' # Skin
#' 
#' ## Load data
#' 
#' ### Clinical data

clin <- list(
  Discovery_skin = read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/clin.txt"),
  Replication_skin = read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/clin_r.txt")
)

#' ### Counts

cnts <- list(
  Discovery_skin = readRDS("data/gene_level_counts/PSORT-D_Skin_counts_01-Apr-2020-13-00-07.rds"),
  Replication_skin = readRDS("data/gene_level_counts/PSORT-R_Skin_counts_13-Mar-2020-15-35-24.rds")
)

intersect_genes <- intersect(rownames(cnts$Discovery_skin), rownames(cnts$Replication_skin))

cnts$Replication_skin <- cnts$Replication_skin[intersect_genes,]

#' ### DEGs

pasi_degs <- read.delim("data/PASI_DE_analysis/LS_vp_sig_fc_genes_module_anno.splinedf3.tsv", sep = "\t")
bmi_degs <- read.delim("data/BMI_DE_analysis/combined_res.csv", sep = ",")

#' ### Gene annotation data

anno <- read.delim("data/gene_annotation_data/Hs.anno.csv", sep = ",") %>%
  # Drop Description column
  select(EnsemblID, GeneSymbol) %>%
  # Replace gene symbols that are "" with corresponding Ensembl ID
  mutate(GeneSymbol = if_else(GeneSymbol == "", EnsemblID, GeneSymbol))

#' ## Filter and normalise counts

normCounts <- function(cnts_dat, clin_dat, anno_dat){
  # Add Drug-Tissue-Time interaction variable to clinical data
  clin_dat <- clin_dat %>%
    mutate(Drug.Tissue.Time = paste(Drug, Tissue, Time, sep = "."))
  # Filter counts
  gene_ids <- tibble(EnsemblID = rownames(cnts_dat)) %>%
    inner_join(anno_dat, by = 'EnsemblID') %>%
    dplyr::select(EnsemblID, GeneSymbol)
  cnts_dat <- cnts_dat[gene_ids$EnsemblID, clin_dat$Sample_id]
  y <- DGEList(cnts_dat, genes = gene_ids)
  keep <- rowSums(cpm(y) >= 1) >= nrow(clin_dat) / length(unique(clin_dat$Drug.Tissue.Time))
  y <- DGEList(y[keep,])
  # TMM normalisation
  y <- calcNormFactors(y)
  return(y)
}

cnts$Discovery_skin <- normCounts(
  cnts_dat = cnts$Discovery_skin, 
  clin_dat = clin$Discovery_skin,
  anno_dat = anno
)

cnts$Replication_skin <- normCounts(
  cnts_dat = cnts$Replication_skin, 
  clin_dat = clin$Replication_skin,
  anno_dat = anno
)

#' ## Combine data

clin[[1]]$Cohort <- "Discovery"
clin[[2]]$Cohort <- "Replication"
clin <- bind_rows(clin)

intersect_genes <- intersect(rownames(cnts$Discovery_skin), rownames(cnts$Replication_skin))

cnts$Discovery_skin <- cnts$Discovery_skin[intersect_genes,]
cnts$Replication_skin <- cnts$Replication_skin[intersect_genes,]

cnts <- cbind(cnts$Discovery_skin, cnts$Replication_skin)

cnts <- voom(cnts)$E

#' ## Plots

plotGeneTrait <- function(gene, gene_symbol, module, trait, cohort, tissue, time, drug, clin, cnts, x_title, cohort_facet, 
                          text_size = 20){
  # Assemble data
  clin <- clin[,c("Sample_id", "Cohort", "Tissue", "Time", "Drug", trait)] %>%
    filter(Cohort %in% cohort, Tissue %in% tissue, Time %in% time, Drug %in% drug)
  cnts <- melt(cnts[gene, clin$Sample_id])
  cnts <- cnts %>%
    rename(Expression = value) %>%
    rownames_to_column("Sample_id")
  dat <- merge(clin, cnts) 
  if(trait != "Time"){
    dat <- dat %>% rename(all_of(c(Trait = trait)))
  }else{
    dat <- dat %>% rename(Trait = Time.1)
  }
  # Gene facet variable
  dat$Gene_facet <- gene_symbol
  # Drug facet variable
  dat$Drug_facet <- dat$Drug
  if(length(drug) == 2){
    dat$Drug_facet <- "Both drugs"
  }
  # Cohort facet variable
  dat$Cohort_facet <- NA
  dat$Cohort_facet[which(dat$Cohort == "Discovery")] <- cohort_facet["Discovery"]
  dat$Cohort_facet[which(dat$Cohort == "Replication")] <- cohort_facet["Replication"]
  # Facet strip colours
  if(unique(dat$Tissue) == "Lesional"){
    tissue_strip_col <- "purple"
    tissue_text_col <- "white"
  }else if(unique(dat$Tissue) == "Nonlesional"){
    tissue_strip_col <- "pink"
    tissue_text_col <- "black"
  }else if(unique(dat$Tissue) == "Blood"){
    tissue_strip_col <- "red"
    tissue_text_col <- "white"
  }
  if(unique(dat$Drug_facet) == "Adalimumab"){
    drug_strip_col <- "darkslategray4"
    drug_text_col <- "white"
  }else if(unique(dat$Drug_facet) == "Ustekinumab"){
    drug_strip_col <- "darkorange3"
    drug_text_col <- "white"
  }else if(unique(dat$Drug_facet) == "Both drugs"){
    drug_strip_col <- "cyan"
    drug_text_col <- "black"
  }
  if(module == "white"){
    module_fill_col <- "white"
    module_col <- "black"
    module_text_col <- "black"
  }else if(module %in% c("black", "blue", "darkred", "navy", "darkslategrey", "firebrick", "darkorchid", "purple")){
    module_fill_col <- module
    module_col <- module
    module_text_col <- "white"
  }else{
    module_fill_col <- module
    module_col <- module
    module_text_col <- "black"
  }
  # Assemble facet strips
  strip_fun <- strip_nested(
    # Horizontal strips
    background_x = elem_list_rect(fill = c(tissue_strip_col, drug_strip_col, "white"), color = c("white", "white", "white")),
    text_x = elem_list_text(color = c(tissue_text_col, drug_text_col, "black")),
    by_layer_x = TRUE,
    # Vertical strips
    background_y = elem_list_rect(fill = module_fill_col, color = module_col),
    text_y = elem_list_text(color = module_text_col),
    by_layer_y = FALSE
  )
  # Plot
  p <- ggplot(data = dat, aes(x = Trait, y = Expression)) +
    scale_fill_manual(values = c("wk00" = "#1F77B4FF", "wk01" = "#2CA02CFF", "wk04" = "#D62728FF", "wk12" = "#FF7F0EFF")) +
    scale_color_manual(values = c("wk00" = "#1F77B4FF", "wk01" = "#2CA02CFF", "wk04" = "#D62728FF", "wk12" = "#FF7F0EFF")) +
    labs(x = x_title) +
    theme_minimal() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size),
          strip.placement = "outside",
          axis.title.y = element_blank(),
          legend.position = "none") +
    facet_nested(Gene_facet ~ Tissue + Drug_facet + Cohort_facet, strip = strip_fun, switch = "y")
  if(is.numeric(dat$Trait)){
    if(grepl("PASI", trait)){
      p <- p + stat_smooth(method="lm", formula = "y ~ ns(x, df=3)", color = "black")
    }else{
      p <- p + geom_smooth(method = "lm", formula = y ~ x, color = "black")
    }
    p <- p + geom_point(aes(color = Time), size = 2.5, alpha = 0.75)
  }else{
    p <- p + 
      geom_violin(aes(fill = Time, color = Time), width = 0.7, alpha = 0.5, draw_quantiles = 0.5, lwd = 0.5) +
      geom_point(aes(fill = Time, color = Time), position = position_jitter(0.1), size = 2.5, alpha = 0.75)
  }
  return(p)
}

genes = c("DNER", "PDE9A", "GNG7", "SCGB1D2", "MMP7")

plot_list <- list()
for(i in 1:length(genes)){
  # Gene
  gene_symbol <- genes[i]
  gene <- anno %>% filter(GeneSymbol == gene_symbol) %>% pull(EnsemblID)
  # BMI plot
  facet_d <- paste0(
    "Discovery\n",
    "log2FC=",bmi_degs %>% filter(ensembl_id == gene) %>% pull(logFC) %>% signif(2),
    "\nFDR=",bmi_degs %>% filter(ensembl_id == gene) %>% pull(adj.P.Val) %>% signif(2)
  )
  facet_r <- paste0(
    "Replication\n",
    "log2FC=",bmi_degs %>% filter(ensembl_id == gene) %>% pull(logFC_r) %>% signif(2),
    "\np-value=",bmi_degs %>% filter(ensembl_id == gene) %>% pull(P.Value_r) %>% signif(2)
  )
  p <- plotGeneTrait(gene = gene, 
                     gene_symbol = gene_symbol,
                     module = "lightyellow",
                     trait = "BMI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Nonlesional",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin,
                     cnts = cnts,
                     x_title = "BMI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
  plot_list[[length(plot_list) + 1]] <- p
  # ADA PASI plot
  facet_both <- paste0(
    "signed fit range=",
    pasi_degs %>% 
      filter(EnsemblID == gene, Drug == "ADA") %>% 
      pull(signed_fit_range) %>% 
      signif(2),
    "\nq-value=",
    pasi_degs %>% 
      filter(EnsemblID == gene, Drug == "ADA") %>% 
      pull(q.value) %>% 
      signif(2)
  )
  p <- plotGeneTrait(gene = gene, 
                     gene_symbol = gene_symbol,
                     module = "lightyellow",
                     trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     cnts = cnts,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_both, Replication = facet_both))
  plot_list[[length(plot_list) + 1]] <- p
  # UST PASI plot
  facet_both <- paste0(
    "signed fit range=",
    pasi_degs %>% 
      filter(EnsemblID == gene, Drug == "UST") %>% 
      pull(signed_fit_range) %>% 
      signif(2),
    "\nq-value=",
    pasi_degs %>% 
      filter(EnsemblID == gene, Drug == "UST") %>% 
      pull(q.value) %>% 
      signif(2)
  )
  p <- plotGeneTrait(gene = gene, 
                     gene_symbol = gene_symbol,
                     module = "lightyellow",
                     trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     cnts = cnts,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_both, Replication = facet_both))
  plot_list[[length(plot_list) + 1]] <- p
}

# Save
png(paste0(output_directory,"/BMI-PASI_trait_plots.png"), width = 32, height = 16, units = "in", res = 300)
grid.arrange(grobs = plot_list, ncol = 6)
dev.off()
pdf(paste0(output_directory,"/BMI-PASI_trait_plots.pdf"), width = 32, height = 16)
grid.arrange(grobs = plot_list, ncol = 6)
dev.off()

#' # Session information

sessionInfo()