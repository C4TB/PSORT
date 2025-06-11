#' ---
#' title: "Module preservation analysis"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' # Preliminaries
#' 
#' ## Load packages

library(tidyverse)
library(edgeR)
library(limma)
library(WGCNA)
library(reshape2)

#' ## Load data

# Module genes
modules <- list(
  Skin = read.delim("results/WGCNA/01_Module_identification/Skin/modules.txt"),
  Blood = read.delim("results/WGCNA/01_Module_identification/Blood/modules.txt")
)

# Normalised counts
cnts <- list(
  Skin = readRDS("results/WGCNA/01_Module_identification/Skin/cnts.rds")$Discovery_skin,
  Blood = readRDS("results/WGCNA/01_Module_identification/Blood/cnts.rds")
)

# Create module colors
module_colors <- list(
  Skin = setNames(modules$Skin$Module, modules$Skin$EnsemblID),
  Blood = setNames(modules$Blood$Module, modules$Blood$EnsemblID)
)

# Make sure counts are in correct order
cnts$Skin <- cnts$Skin[,names(module_colors$Skin)]
cnts$Blood <- cnts$Blood[,names(module_colors$Blood)]

#' # Skin-blood

# Output directory
output_path <- "results/WGCNA/10_Module_preservation_analysis/Skin-blood"
dir.create(output_path, recursive = T)

# Module preservation
cor <- WGCNA::cor
mp <- WGCNA::modulePreservation(
  multiData = list(Set1 = list(data = cnts$Skin), Set2 = list(data = cnts$Blood)), 
  multiColor = list(Set1 = module_colors$Skin, Set2 = module_colors$Blood),
  dataIsExpr = T,
  networkType = "signed",
  referenceNetworks = 1,
  nPermutations = 200,
  randomSeed = 1,
  quickCor = 0,
  verbose = 3,
  goldName = "random_sample",
  maxModuleSize = max(table(module_colors$Skin))
)

# Save
saveRDS(mp, paste0(output_path,"/mp.rds"))

#' # Blood-skin

# Output directory
output_path <- "results/WGCNA/10_Module_preservation_analysis/Blood-skin"
dir.create(output_path, recursive = T)

# Module preservation
cor <- WGCNA::cor
mp <- WGCNA::modulePreservation(
  multiData = list(Set1 = list(data = cnts$Blood), Set2 = list(data = cnts$Skin)), 
  multiColor = list(Set1 = module_colors$Blood, Set2 = module_colors$Skin),
  dataIsExpr = T,
  networkType = "signed",
  referenceNetworks = 1,
  nPermutations = 200,
  randomSeed = 1,
  quickCor = 0,
  verbose = 3,
  goldName = "random_sample",
  maxModuleSize = max(table(module_colors$Blood))
)

# Save
saveRDS(mp, paste0(output_path,"/mp.rds"))

#' # Session information

sessionInfo()