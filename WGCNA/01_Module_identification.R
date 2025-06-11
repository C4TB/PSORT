#' ---
#' title: "Identification of co-expressed gene modules in skin and blood"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 12, fig.height = 9)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' Here, we use WGCNA, or **W**eighted **G**ene **C**o-expression **N**etwork **A**nalysis 
#' ([Langfelder & Horvath, 2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)),
#' to identify co-expressed gene modules in skin and blood.
#' 
#' # Preliminaries
#' 
#' ## Load packages

library(tidyverse)
library(edgeR)
library(limma)
library(WGCNA)
library(reshape2)

#' ## Load gene annotation data

anno <- read.delim("data/gene_annotation_data/Hs.anno.csv", sep = ",") %>%
  # Drop Description column
  select(EnsemblID, GeneSymbol) %>%
  # Replace gene symbols that are "" with corresponding Ensembl ID
  mutate(GeneSymbol = if_else(GeneSymbol == "", EnsemblID, GeneSymbol))

head(anno)

#' ## Create output directory

output_directory <- "results/WGCNA/01_Module_identification"
dir.create(output_directory)

#' # Skin
#'
#' We'll start by analysing the skin data.
#' 
#' ## Create output sub-directory

output_directory2 <- paste0(output_directory,"/Skin")
dir.create(output_directory2)

#' ## Load data
#' 
#' ### Clinical data
#' 
#' We load the clinical data for each cohort.

# Load clinical data for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
clin <- list()
clin$Discovery_skin <- read.delim("data/clinical_data/PSORT-D_Skin_Clinical_Data_01-Apr-20.txt") %>%
  select(Patient_id, Sample_id, Drug, Tissue, Time)
clin$Replication_skin <- read.delim("data/clinical_data/PSORT-R_Skin_Clinical_Data_01-Apr-20.txt") %>%
  select(Patient_id, Sample_id, Drug, Tissue, Time) %>% 
  dplyr::filter(!Sample_id %in% c("3A", "3B", "4A", "4B"))

dim(clin$Discovery_skin)

head(clin$Discovery_skin)

dim(clin$Replication_skin)

head(clin$Replication_skin)

#' Preliminary analysis revealed mislabelling of some lesional and non-lesional samples, which were swapped
#' for several patients. We correct these swaps here.

# Sample swaps based on consensus of S100A7 and S100A9 expression
clin$Discovery_skin <- clin$Discovery_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0005", "23012-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0006", "23012-1206"), "Lesional", Tissue))
clin$Replication_skin <- clin$Replication_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0005", "5031-1205", "6049-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0006", "5031-1206", "6049-1206"), "Lesional", Tissue))

#' ### Gene-level counts
#' 
#' We also load gene-level counts for each cohort.

# Load counts for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
cnts <- list()
cnts$Discovery_skin <- readRDS("data/gene_level_counts/PSORT-D_Skin_counts_01-Apr-2020-13-00-07.rds")
cnts$Replication_skin <- readRDS("data/gene_level_counts/PSORT-R_Skin_counts_13-Mar-2020-15-35-24.rds")

dim(cnts$Discovery_skin)

cnts$Discovery_skin[1:5, 1:5]

dim(cnts$Replication_skin)

cnts$Replication_skin[1:5, 1:5]

#' The PSORT-R counts contain some extra genes that are not present in the PSORT-D counts. We'll proceed with
#' just the genes that are present in both datasets.

intersect_genes <- intersect(rownames(cnts$Discovery_skin), rownames(cnts$Replication_skin))

cnts$Replication_skin <- cnts$Replication_skin[intersect_genes,]

dim(cnts$Discovery_skin)

dim(cnts$Replication_skin)

#' ## Filter and normalise counts
#' 
#' We now filter and normalise the counts. We use a filtering threshold that requires at least 1 CPM in at least
#' *n/k* samples, where *n* equals the number of samples and *k* equals the number of unique combinations of
#' drug, tissue type (i.e. lesional and non-lesional) and time point. We then normalise the counts using the
#' trimmed mean of m-values (TMM) method from 
#' [Robinson & Oshlack, 2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25). We wrap
#' this workflow in a function and apply it to the PSORT-D and PSORT-R data separately.

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

dim(cnts$Discovery_skin)

dim(cnts$Replication_skin)

#' Analysis will be conducted using log2-CPM counts, which we derive using the voom function from the limma package.

cnts$Discovery_skin <- voom(cnts$Discovery_skin)$E

cnts$Replication_skin <- voom(cnts$Replication_skin)$E

#' The WGCNA R package requires that the count data have samples as rows and genes as columns; therefore, we
#' transpose the count matrices here.

cnts$Discovery_skin <- t(cnts$Discovery_skin)

cnts$Replication_skin <- t(cnts$Replication_skin)

cnts$Discovery_skin[1:5, 1:5]

cnts$Replication_skin[1:5, 1:5]

#' We'll save the normalised counts for use in downstream analyses.

saveRDS(cnts, paste0(output_directory2,"/cnts.rds"))

#' ## Choose soft-thresholding power

#' Module identification will be carried out using the PSORT-D data. Prior to this we need to identify an
#' appropriate soft-thresholding power. We'll do this by plotting the candidate values 1-20 against R2, 
#' a measure of scale-free topology, and mean connectivity. We'll choose the lowest value that reaches 
#' an R2 threshold of 0.8.

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(
  cnts$Discovery_skin, 
  powerVector = powers, 
  verbose = 5, 
  networkType = "signed", 
  corFnc = WGCNA::cor
)

# Plot the results
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",
  type="n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
  labels = powers, 
  cex = cex1, 
  col="red"
)
# This line corresponds to using an R^2 cut-off of h
abline(h = 0.80, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1], 
  sft$fitIndices[,5],
  xlab="Soft Threshold (power)",
  ylab="Mean Connectivity", 
  type="n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[,1], 
  sft$fitIndices[,5], 
  labels = powers, 
  cex = cex1, 
  col = "red"
)

#' Based on these plots, we'll choose a soft-thresholding power of 12.

#' ## One-step network construction and module detection
#' 
#' Now we can construct the co-expression network and identify modules using the blockwiseModules function
#' from the WGCNA R package.
#' In brief, here we calculate the pearson correlations between each pair of genes
#' and raise these estimates to the selected soft-thresholding power of 12; this amplifies the differences
#' between high and low correlations. 
#' These correlations are then used to generate a topological overlap matrix (TOM) and hierarchical clustering 
#' of this matrix is used to group genes with similar expression profiles into modules. 
#' We specify a minimum module size of 30, and a dendrogram cut height (for merging of similar modules) of 0.1. 
#' We also use a signed network so that the correlations between genes are scaled to lie between 0 and 1.

net <- blockwiseModules(
  cnts$Discovery_skin, 
  power = 12, 
  networkType = "signed",
  corType = "pearson", 
  maxPOutliers = 0.1,
  TOMType = "signed", 
  minModuleSize = 30,
  reassignThreshold = 0, 
  mergeCutHeight = 0.1,
  numericLabels = TRUE, 
  pamRespectsDendro = FALSE,
  stabilityCriterion = "Individual fraction",
  saveTOMs = FALSE,
  verbose = 3, 
  maxBlockSize = ncol(cnts$Discovery_skin)
)

# Save
saveRDS(net, paste0(output_directory2,"/net.rds"))

#' Now we can visualise the module dendrogram.

# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  net$dendrograms[[1]], 
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05
)

#' ## Module assignments
#' 
#' Next, we record the assignments of genes to modules and save this to file.

modules <- bind_rows(lapply(
  X = unique(moduleColors),
  FUN = function(x) anno %>% 
    filter(EnsemblID %in% colnames(cnts$Discovery_skin)[moduleColors == x]) %>%
    mutate(Module = x)
))

head(modules)

# Save
write.table(modules, paste0(output_directory2,"/modules.txt"), sep = "\t", row.names = F, quote = F)

#' We can also examine the size of each module.

table(modules$Module)

#' ## Module eigengenes
#' 
#' Now we can use the module assignments to calculate module eigengenes for each module.

eigen <- moduleEigengenes(cnts$Discovery_skin, moduleColors)$eigengenes %>%
  rownames_to_column(var = "Sample_id") %>%
  rename_with(~ gsub("ME", "", .x, fixed = TRUE))

eigen[1:5, 1:5]

# Save
write.table(eigen, paste0(output_directory2,"/eigengenes.txt"), sep = "\t", row.names = F, quote = F)

#' We also do this for the PSORT-R cohort (using module assignments defined in the PSORT-D cohort). 
#' The number of genes retained in the filtered counts is different between PSORT-D and PSORT-R; therefore, we
#' use the overlapping genes to calculate eigengenes for the PSORT-R cohort. 

# Module labels for PSORT-R data
moduleColors_r <- net$colors[intersect(colnames(cnts$Replication_skin), names(net$colors))]

# Subset counts
cnts$Replication_skin <- cnts$Replication_skin[,names(moduleColors_r)]

# Convert labels to colors
moduleColors_r <- labels2colors(moduleColors_r)

# Calculate eigengenes
eigen_r <- moduleEigengenes(cnts$Replication_skin, moduleColors_r)$eigengenes %>%
  rownames_to_column(var = "Sample_id") %>%
  rename_with(~ gsub("ME", "", .x, fixed = TRUE))

eigen_r[1:5, 1:5]

# Save
write.table(eigen_r, paste0(output_directory2,"/eigengenes_r.txt"), sep = "\t", row.names = F, quote = F)

#' ## Module membership
#' 
#' Here we will calculate module membership for each gene with every module; this is defined as the correlation
#' of a gene's expression values with a module eigengene.

# Performs gene-module correlation
geneModuleCor <- function(cnts_dat, eigen_dat){
  eigen_dat <- eigen_dat %>% 
    remove_rownames %>% 
    column_to_rownames("Sample_id")
  eigen_dat <- eigen_dat[rownames(cnts_dat),]
  # Calculate gene eigengene correlations and associated p-values
  cp = WGCNA::corAndPvalue(cnts_dat, eigen_dat)
  # Merge correlation coefficients and p-values for each gene-eigen_datgene pair into one data frame
  cor_dat <- melt(cp$cor)
  p_dat <- melt(cp$p)
  dat <- merge(cor_dat, p_dat, by = c("Var1", "Var2"))
  colnames(dat) <- c("EnsemblID", "Module", "Cor", "P.Value")
  # Convert factor columns to character columns
  dat$EnsemblID <- as.character(dat$EnsemblID)
  dat$Module <- as.character(dat$Module)
  return(dat)
}

mm <- geneModuleCor(cnts_dat = cnts$Discovery_skin, eigen_dat = eigen)

head(mm)

# Save
write.table(mm, paste0(output_directory2,"/mm.txt"), sep = "\t", row.names = F, quote = F)

#' We'll also do this for the PSORT-R data.

mm_r <- geneModuleCor(cnts_dat = cnts$Replication_skin, eigen_dat = eigen_r)

head(mm_r)

# Save
write.table(mm_r, paste0(output_directory2,"/mm_r.txt"), sep = "\t", row.names = F, quote = F)

#' # Blood
#'
#' Now we will analyse the blood data.
#' 
#' ## Create output sub-directory

output_directory2 <- paste0(output_directory,"/Blood")
dir.create(output_directory2)

#' ## Load data
#' 
#' ### Clinical data
#' 
#' We load the clinical data for each cohort.

# Load clinical data for PSORT-D blood samples
clin <- read.delim("data/clinical_data/PSORT-D_Blood_Clinical_Data_01-Apr-20.txt") %>%
  select(Patient_id, Sample_id, Drug, Tissue, Time)

dim(clin)

head(clin)

#' Preliminary analysis revealed some outlier samples in the blood data that we don't want to include in the
#' analysis. There are also some smaples from other cohorts that we don't want to include.
#' Here we will read in a file containing the IDs of samples that we want to analyse; we'll use this to
#' subset the clinical data.

samples <- read.delim("data/clinical_data/PSORT-D_Blood_analysis_samples.txt")

clin <- clin %>% filter(Sample_id %in% samples$Sample_id)

dim(clin)

#' ### Gene-level counts
#' 
#' We also load gene-level counts for each cohort.

# Load counts for PSORT-D blood samples
cnts <- readRDS("data/gene_level_counts/PSORT-D_Blood_counts_01-Apr-2020-13-00-10.rds")

# Subset
cnts <- cnts[,clin$Sample_id]

dim(cnts)

cnts[1:5, 1:5]

#' ## Filter and normalise counts
#' 
#' We now filter and normalise the counts using the function defined in the skin section above.

cnts <- normCounts(cnts_dat = cnts, clin_dat = clin,anno_dat = anno)

#' And we derive log2-CPM counts using the coom function.

cnts <- voom(cnts)$E

#' We also need to transpose the counts.

cnts <- t(cnts)

cnts[1:5, 1:5]

#' Again, we'll save the normalised counts for use in downstream analyses.

saveRDS(cnts, paste0(output_directory2,"/cnts.rds"))

#' ## Choose soft-thresholding power

#' Again, we need to define an appropriate soft-thresholding power.

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(
  cnts, 
  powerVector = powers, 
  verbose = 5, 
  networkType = "signed", 
  corFnc = WGCNA::cor
)

# Plot the results
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",
  type="n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
  labels = powers, 
  cex = cex1, 
  col="red"
)
# This line corresponds to using an R^2 cut-off of h
abline(h = 0.80, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[,1], 
  sft$fitIndices[,5],
  xlab="Soft Threshold (power)",
  ylab="Mean Connectivity", 
  type="n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[,1], 
  sft$fitIndices[,5], 
  labels = powers, 
  cex = cex1, 
  col = "red"
)

#' Based on these plots, we'll choose a soft-thresholding power of 12.

#' ## One-step network construction and module detection

net <- blockwiseModules(
  cnts, 
  power = 5, 
  networkType = "signed",
  corType = "pearson", 
  maxPOutliers = 0.1,
  TOMType = "signed", 
  minModuleSize = 30,
  reassignThreshold = 0, 
  mergeCutHeight = 0.1,
  numericLabels = TRUE, 
  pamRespectsDendro = FALSE,
  stabilityCriterion = "Individual fraction",
  saveTOMs = FALSE,
  verbose = 3, 
  maxBlockSize = ncol(cnts)
)

# Save
saveRDS(net, paste0(output_directory2,"/net.rds"))

#' Now we can visualise the module dendrogram. We will also rename the blood modules to make them distinct from
#' the skin modules.

# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)

table(moduleColors)

# Create data frame that maps between old names and new names
new_names <- c(
  "darkorchid" = "black",
  "dodgerblue" = "blue",
  "gold" = "brown",
  "aquamarine" = "cyan",
  "sienna" = "darkgreen",
  "darkslategrey" = "darkgrey",
  "coral" = "darkorange",
  "firebrick" = "darkred",
  "khaki" = "darkturquoise",
  "limegreen" = "green",
  "palegreen" = "greenyellow",
  "grey" = "grey",
  "slategrey" = "grey60",
  "thistle" = "lightcyan",
  "rosybrown" = "lightgreen",
  "sandybrown" = "lightyellow",
  "maroon" = "magenta",
  "lightcoral" = "midnightblue",
  "orangered" = "orange",
  "plum" = "pink",
  "orchid" = "purple",
  "chocolate" = "red",
  "cadetblue" = "royalblue",
  "burlywood" = "salmon",
  "darkviolet" = "tan",
  "navy" = "turquoise",
  "olivedrab" = "yellow"
)
new_names <- data.frame(Old = new_names, New = names(new_names), row.names = NULL)

# Rename modules
for(i in 1:length(moduleColors)){
  moduleColors[i] <- new_names$New[which(new_names$Old == moduleColors[i])]
}

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  net$dendrograms[[1]], 
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05
)

#' ## Module assignments
#' 
#' Next, we record the assignments of genes to modules and save this to file.

modules <- bind_rows(lapply(
  X = unique(moduleColors),
  FUN = function(x) anno %>% 
    filter(EnsemblID %in% colnames(cnts)[moduleColors == x]) %>%
    mutate(Module = x)
))

head(modules)

# Save
write.table(modules, paste0(output_directory2,"/modules.txt"), sep = "\t", row.names = F, quote = F)

#' We can also examine the size of each module.

table(modules$Module)

#' ## Module eigengenes
#' 
#' Now we can use the module assignments to calculate module eigengenes for each module.

eigen <- moduleEigengenes(cnts, moduleColors)$eigengenes %>%
  rownames_to_column(var = "Sample_id") %>%
  rename_with(~ gsub("ME", "", .x, fixed = TRUE))

eigen[1:5, 1:5]

# Save
write.table(eigen, paste0(output_directory2,"/eigengenes.txt"), sep = "\t", row.names = F, quote = F)

#' ## Module membership
#' 
#' Here we will calculate module membership for each gene with every module.

mm <- geneModuleCor(cnts_dat = cnts, eigen_dat = eigen)

head(mm)

# Save
write.table(mm, paste0(output_directory2,"/mm.txt"), sep = "\t", row.names = F, quote = F)

#' # Session information

sessionInfo()