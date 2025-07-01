#' ---
#' title: "Cross-referencing with Broad single cell data"
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
library(WGCNA)
library(gridExtra)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(mgsub)
library(readxl)
library(data.table)

#' ## Load data

options(stringsAsFactors = F)

# Path to Broad SC data
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"

# Load expression file
dat <- fread(paste0(file_path,"/ExpressionFile.txt"))

# Load metadata
meta <- fread(paste0(file_path,"/MetadataUploadFinal.txt"))
meta <- meta[2:nrow(meta),]

# Skin modules
modules <- read.delim("results/WGCNA/01_Module_identification/Skin/modules.txt")

#' # Heatmaps

# Output directory
output_path <- "results/WGCNA/09_Broad_single_cell_heatmaps"
dir.create(output_path)

# Function to draw heatmaps
drawHeatmap <- function(expr_dat, meta_dat, pcnt_dat, circle_factor = 0.025){
  # Melt expression data
  expr_dat <- reshape2::melt(expr_dat)
  colnames(expr_dat) <- c("GENE", "NAME", "Expression")
  # Merge with metadata
  expr_dat <- merge(expr_dat, meta_dat, by = "NAME")
  # Calculate average expression for each gene in each cluster
  expr_dat <- expr_dat %>%
    group_by(GENE, Specific) %>%
    summarise(Mean = mean(Expression), Median = median(Expression))
  # Plot mean
  expr_dat <- tidyr::spread(expr_dat[,c("GENE", "Specific", "Mean")], Specific, Mean)
  # Coerce to matrix
  expr_dat <- as.data.frame(expr_dat)
  rownames(expr_dat) <- expr_dat$GENE
  expr_dat$GENE <- NULL
  expr_dat <- data.matrix(expr_dat)
  
  pcnt_dat <- tidyr::spread(pcnt_dat[,c("GENE", "Specific", "TRUE")], Specific, `TRUE`)
  pcnt_dat <- as.data.frame(pcnt_dat)
  rownames(pcnt_dat) <- pcnt_dat$GENE
  pcnt_dat$GENE <- NULL
  pcnt_dat <- data.matrix(pcnt_dat)
  
  # Module membership of markers
  modules <- modules %>% 
    select(GeneSymbol, Module) %>%
    filter(GeneSymbol %in% rownames(expr_dat))
  rownames(modules) <- modules$GeneSymbol
  # Order/subset rows
  expr_dat <- expr_dat[rownames(modules),]
  pcnt_dat <- pcnt_dat[rownames(modules),]
  # Row annotation to indicate module membership
  anno_col <- unique(modules$Module)
  names(anno_col) <- unique(modules$Module)
  side_anno = rowAnnotation(
    ModuleName = anno_text(modules$Module, gp = gpar(fontsize = 6)),
    ModuleColour = modules$Module,
    col = list(ModuleColour = anno_col),
    show_legend = F,
    show_annotation_name = F
  )
  # Change column names if looking at KT clusters
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-1")] <- "KC-1 (spinous)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-2")] <- "KC-2 (spinous)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-4")] <- "KC-4 (spinous)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-6")] <- "KC-6 (suprabasal)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-3")] <- "KC-3 (suprabasal)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-7")] <- "KC-7 (basal)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-8")] <- "KC-8 (basal)"
  colnames(expr_dat)[which(colnames(expr_dat) == "Keratinocyte-5")] <- "KC-5 (supra-spinous)"
  # Change column names if looking at KT clusters
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-1")] <- "KC-1 (spinous)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-2")] <- "KC-2 (spinous)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-4")] <- "KC-4 (spinous)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-6")] <- "KC-6 (suprabasal)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-3")] <- "KC-3 (suprabasal)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-7")] <- "KC-7 (basal)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-8")] <- "KC-8 (basal)"
  colnames(pcnt_dat)[which(colnames(pcnt_dat) == "Keratinocyte-5")] <- "KC-5 (supra-spinous)"
  # Draw heatmap
  col_ramp <- colorRamp2(breaks = c(0, max(expr_dat)), colors = c("white", "red"))
  hm <- Heatmap(
    expr_dat,
    col = col_ramp,
    rect_gp = gpar(type = "none"),
    width = ncol(expr_dat)*unit(15,"mm"),
    right_annotation = side_anno,
    row_names_gp = gpar(fontsize = 6),
    row_dend_side = "left",
    heatmap_legend_param = list(title = "Mean expression", legend_direction = "vertical"),
    cell_fun = function(j, i, x, y, w, h, col){
      grid.circle(x = x, y = y, r = pcnt_dat[i, j] * circle_factor, 
                  gp = gpar(fill = col_ramp(expr_dat[i, j]), col = "black",
                            alpha = 0.75))
    }
  )
  return(hm)
}

#' ## KC

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "KC") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
markers <- c(markers, "LTF")
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "KC", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- reshape2::melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Broad_KC_markers.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Fibro

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Fibro") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "Fibro", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Fibro.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## T-cells

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "T") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "T", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/T-cells.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Myeloid

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Myeloid") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "Myeloid", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Myeloid.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Endo

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Endo") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "Endo", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Endo.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Mast

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Mast") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "Mast", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Mast.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## VSMC

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "VSMC") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/VSMC.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## B-cells

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "B") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/B-cells.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Hair follicle

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "HairFollicle") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/HairFollicle.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Langerhans

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Langerhans") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Langerhans.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Lymphatic

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Lymphatic") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Lymphatic.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Melanocyte

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Melanocyte") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Melanocyte.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Plasma

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Plasma") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Plasma.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Schwann

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Schwann") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Schwann.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Sebocyte

# Choose marker genes to plot
file_path <- "C:/Users/b6054775A/Documents/PhD/data/Data_12-Mar-20/Broad_single_cell"
markers <- read_xlsx(paste0(file_path,"/mmc4.xlsx"), sheet = "GenericSignatures") %>%
  filter(cluster == "Sebocyte") %>%
  slice_max(avg_logFC, n = 80) %>%
  pull(gene)
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 %in% c("VSMC", "B", "HairFollicle", "Langerhans", "Lymphatic", 
                                                "Melanocyte", "Plasma", "Schwann", "Sebocyte"), 
                            Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Sebocyte.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' ## Fibro (new differentiating markers)

# Choose marker genes to plot
markers <- unique(c(
  "DCN", "COL3A1", "SPARC", "APOE", "CXCL12", "COL18A1", "LUM",
  "PDGFRA", "APCDD1", "CCL19", "COL8A1", "DPEP1", "LEF1", "MKX",
  "TNMD", "FGF18", "HHIP", "S100B", "SOX10", "MYL4", "COL5A1",
  "VCAM1", "POSTN", "TNC", "SFRP1", "PI16", "C7", "PLA2G2A",
  "CH25H", "TNFSF13B", "MEF2C", "IL15", "WIF1", "RSPO4", "CRABP1",
  "ITGA6", "PPARG", "COCH", "RGS5", "SLPI", "NKD2", "COL23A1",
  "IL33", "CD34", "PCOLCE2", "RSPO3", "CORIN", "COL13A1", "DPP4",
  "SCN7A", "CCL19", "TNC", "LRRC15", "PDLIM3", "LRRC17", "CTHRC1",
  "ACTA2", "SFRP4", "MMP1", "PIEZO2", "THBS4", "CXCL6", "CXCL5",
  "CXCL13", "ADAM12", "MMP3", "CDH2", "IL24", "ITGA10", "IL11",
  "SFRP1", "SFRP2", "SFRP3", "SFRP4"
))
# Subset metadata
meta_dat <- meta %>% filter(CellType1.50 == "Fibro", Condition == "Psoriasis")
# Subset expression data
cells <- c("GENE", meta_dat$NAME)
expr_dat <- as.data.frame(dat[GENE %in% markers, ..cells])
# Percentage expressed
pcnt_dat <- melt(expr_dat)
pcnt_dat <- pcnt_dat %>%
  rename(all_of(c(NAME = "variable", Expression = "value"))) %>%
  left_join(meta_dat, by = "NAME") %>%
  mutate(Expressed = if_else(Expression > 0, T, F))
pcnt_dat <- pcnt_dat %>%
  group_by(Specific, GENE) %>%
  summarise("TRUE" = mean(Expressed),
            "FALSE" = mean(!Expressed)) %>%
  as.data.frame()
# Draw heatmap
hm <- drawHeatmap(expr_dat = expr_dat, meta_dat = meta_dat, pcnt_dat = pcnt_dat)
#+fig.width=9,fig.height=8.5
draw(hm)
# Save
png(paste0(output_path,"/Fibro_new_diff_markers.png"), 
    width = 9, height = 8.5, units = "in", res = 300)
ComplexHeatmap::draw(hm)
dev.off()

#' # Session information

sessionInfo()