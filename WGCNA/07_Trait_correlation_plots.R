#' ---
#' title: "Trait correlation plots"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 19, fig.height = 14)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' Here, we'll plot some exemplar module and factor-tait correlations.
#' 
#' # Preliminaries
#' 
#' ## Load packages

library(tidyverse)
library(gridExtra)
library(ggh4x)
library(scales)
library(splines)

#' ## Create output directory

output_directory <- "results/WGCNA/07_Trait_correlation_plots"
dir.create(output_directory)

#' # Skin
#' 
#' We'll start by creating some plots for significant associations with the skin modules and factors. First, we need to load
#' the data.
#' 
#' ## Load data
#' 
#' ### Clinical data
#' 
#' We'll load the clinical data for the PSORT-D and PSORT-R cohorts and bind it together.

clin <- list(
  read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/clin.txt"),
  read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/clin_r.txt")
)
clin[[1]]$Cohort <- "Discovery"
clin[[2]]$Cohort <- "Replication"
clin <- bind_rows(clin)

#' ### Eigengenes
#' 
#' We'll do the same for the module eigengenes...

eigen <- rbind(
  read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes.txt"),
  read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes_r.txt")
)

#' ### Latent factors
#' 
#' And the latent factors...

factors <- rbind(
  read.delim("data/latent_factors/skin_d_m.csv", sep = ","),
  read.delim("data/latent_factors/skin_r_m.csv", sep = ",")
)
colnames(factors)[1] <- "Sample_id"

#' ### Module and factor-trait correlations
#' 
#' We'll also load the module/factor-trait correlation data.

cor_dat <- rbind(
 read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/Module-trait_correlations.txt"),
 read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Skin/Factor-trait_correlations.txt")
)

#' ## Plots
#' 
#' Below we define a function that assembles all the above data into a correlation plot and use it to 
#' plot some exemplars in skin.

plotModuleTrait <- function(module, trait, cohort, tissue, time, drug, clin, eigen, x_title, cohort_facet, text_size = 12){
  # Assemble data
  clin <- clin[,c("Sample_id", "Cohort", "Tissue", "Time", "Drug", trait)] %>%
    filter(Cohort %in% cohort, Tissue %in% tissue, Time %in% time, Drug %in% drug)
  eigen <- eigen[,c("Sample_id", module)]
  dat <- merge(clin, eigen) 
  if(trait != "Time"){
    dat <- dat %>% rename(all_of(c(Trait = trait, Eigengene = module)))
  }else{
    dat <- dat %>% rename(Trait = Time.1) %>% rename(all_of(c(Eigengene = module)))
  }
  # Module facet variable
  if(grepl("factor",module)){
    if(unique(dat$Tissue) %in% c("Lesional", "Nonlesional")){
      dat$Module_facet <- gsub("_", " S", module)
    }else if(unique(dat$Tissue) == "Blood"){
      dat$Module_facet <- gsub("_", " B", module)
    }
  }else{
    dat$Module_facet <- module
  }
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
  if(grepl("factor", module)){
    module_fill_col <- "grey"
    module_col <- "grey"
    module_text_col <- "black"
  }else if(module == "white"){
    module_fill_col <- "white"
    module_col <- "black"
    module_text_col <- "black"
  }else if(module %in% c("black", "blue", "darkred", "navy", "darkslategrey", "firebrick", "darkorchid", "purple",
                         "darkmagenta")){
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
  p <- ggplot(data = dat, aes(x = Trait, y = Eigengene)) +
    scale_fill_manual(values = c("wk00" = "#1F77B4FF", "wk01" = "#2CA02CFF", "wk04" = "#D62728FF", "wk12" = "#FF7F0EFF")) +
    scale_color_manual(values = c("wk00" = "#1F77B4FF", "wk01" = "#2CA02CFF", "wk04" = "#D62728FF", "wk12" = "#FF7F0EFF")) +
    labs(x = x_title) +
    theme_minimal() +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size),
          strip.placement = "outside",
          axis.title.y = element_blank(),
          legend.position = "none") +
    facet_nested(Module_facet ~ Tissue + Drug_facet + Cohort_facet, strip = strip_fun, switch = "y")
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

plot_list <- list()

# Module and trait
module <- "lightyellow"
trait <- "NL_wk00_BMI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "BMI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Nonlesional",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin,
                     eigen = eigen,
                     x_title = "BMI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_9"
trait <- "NL_wk00_BMI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "BMI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Nonlesional",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin,
                     eigen = factors,
                     x_title = "BMI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "salmon"
trait <- "NL_wk00_Psoriatic_Arthritis"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "Psoriatic_Arthritis",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Nonlesional",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin,
                     eigen = eigen,
                     x_title = "PsA",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "turquoise"
trait <- "LS_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "turquoise"
trait <- "LS_UST_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "turquoise"
trait <- "NL_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Nonlesional",
                     time = c("wk00", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "blue"
trait <- "LS_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "blue"
trait <- "LS_UST_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "lightyellow"
trait <- "LS_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "lightyellow"
trait <- "LS_UST_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_1"
trait <- "LS_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = factors,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_9"
trait <- "LS_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
cor_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_r) %>% signif(2)
p_r <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(P.Value_r) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
facet_r <- paste0("Rep cor=",cor_r," (p=",p_r,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = c("Discovery", "Replication"),
                     tissue = "Lesional",
                     time = c("wk00", "wk01", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = factors,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d, Replication = facet_r))
plot_list[[length(plot_list) + 1]] <- p

#' # Blood
#' 
#' Now we'll do the same as above for blood.
#' 
#' ## Load data
#' 
#' ### Clinical data

clin <- read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Blood/clin.txt")
clin$Cohort <- "Discovery"

#' ### Eigengenes

eigen <- read.delim("results/WGCNA/01_Module_identification/Blood/eigengenes.txt")

#' ### Latent factors

factors <- read.delim("data/latent_factors/blood_m.csv", sep = ",")
colnames(factors)[1] <- "Sample_id"

#' ### Module and factor-trait correlations

cor_dat <- rbind(
  read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Blood/Module-trait_correlations.txt"),
  read.delim("results/WGCNA/03_Get_disease_and_disease_severity_correlations/Blood/Factor-trait_correlations.txt")
)

#' ## Plots

# Module and trait
module <- "khaki"
trait <- "BL_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = c("wk00", "wk01", "wk04", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "khaki"
trait <- "BL_UST_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = c("wk00", "wk01", "wk04", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     eigen = eigen,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_19"
trait <- "BL_wk00_Cw6_PosNeg"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "Cw6_PosNeg",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin[!is.na(clin$Cw6_PosNeg),],
                     eigen = factors,
                     x_title = "Cw6",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_8"
trait <- "BL_wk00_Cw6_PosNeg"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "Cw6_PosNeg",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = "wk00",
                     drug = c("Adalimumab", "Ustekinumab"),
                     clin = clin[!is.na(clin$Cw6_PosNeg),],
                     eigen = factors,
                     x_title = "Cw6",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_1"
trait <- "BL_ADA_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = c("wk00", "wk01", "wk04", "wk12"),
                     drug = "Adalimumab",
                     clin = clin,
                     eigen = factors,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

# Module and trait
module <- "factor_1"
trait <- "BL_UST_PASI"
# Cohort facet strips
cor_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(Cor_d) %>% signif(2)
padj_d <- cor_dat %>% filter(Module == module, Trait == trait) %>% pull(padj_all_d) %>% signif(2)
facet_d <- paste0("Disc cor=",cor_d," (padj=",padj_d,")")
# Plot
p <- plotModuleTrait(module = module, trait = "PASI",
                     cohort = "Discovery",
                     tissue = "Blood",
                     time = c("wk00", "wk01", "wk04", "wk12"),
                     drug = "Ustekinumab",
                     clin = clin,
                     eigen = factors,
                     x_title = "PASI",
                     cohort_facet = c(Discovery = facet_d))
plot_list[[length(plot_list) + 1]] <- p

#' # Assemble plot panel
#' 
#' Finally, we'll assemble the skin and blood plots into one panel and save.

grid.arrange(grobs = plot_list, ncol = 4)

# Save
png(paste0(output_directory,"/Skin_Blood_trait_plots.png"), width = 19, height = 14, units = "in", res = 300)
grid.arrange(grobs = plot_list, ncol = 4)
dev.off()
pdf(paste0(output_directory,"/Skin_Blood_trait_plots.pdf"), width = 19, height = 14)
grid.arrange(grobs = plot_list, ncol = 4)
dev.off()

#' # Session information

sessionInfo()