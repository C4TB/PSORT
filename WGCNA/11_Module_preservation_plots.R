#' ---
#' title: "Module preservation plots"
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

#' ## Load data

# Module genes
modules <- list(
  Skin = read.delim("results/WGCNA/01_Module_identification/Skin/modules.txt"),
  Blood = read.delim("results/WGCNA/01_Module_identification/Blood/modules.txt")
)

#' ## Output directory

output_path <- "results/WGCNA/11_Module_preservation_plots"
dir.create(output_path)

#' # Functions

overlapHeatmap <- function(mods1, mods2, title1, title2, show_grey = F){
  n_mat <- matrix(data = NA, nrow = length(mods1), ncol = length(mods2), dimnames = list(names(mods1), names(mods2)))
  p_mat <- n_mat
  df <- data.frame(EnsemblID = intersect(bind_rows(mods1)$EnsemblID, bind_rows(mods2)$EnsemblID))
  for(i in 1:length(mods1)){
    for(j in 1:length(mods2)){
      # N overlap matrix
      genes1 <- mods1[[rownames(n_mat)[i]]]$EnsemblID
      genes2 <- mods2[[colnames(n_mat)[j]]]$EnsemblID
      n_mat[i, j] <- signif((length(intersect(genes1, genes2))/length(unique(c(genes1, genes2))))*100,2)
      # P-value matrix
      df <- df %>% 
        mutate(Module1 = ifelse(EnsemblID %in% genes1, T, F)) %>%
        mutate(Module2 = ifelse(EnsemblID %in% genes2, T, F))
      p_mat[i, j] <- -log10(fisher.test(table(df$Module1, df$Module2), alternative = "greater")$p.value)
    }
  }
  
  if(show_grey == F){
    p_mat <- p_mat[rownames(p_mat)[rownames(p_mat) != "grey"], colnames(p_mat)[colnames(p_mat) != "grey"]]
    n_mat <- n_mat[rownames(n_mat)[rownames(n_mat) != "grey"], colnames(n_mat)[colnames(n_mat) != "grey"]]
  }
  
  n_mat[which(n_mat == 0)] <- ""
  
  
  # Row annotation
  anno_col <- rownames(p_mat)
  names(anno_col) <- rownames(p_mat)
  row_anno = rowAnnotation(
    Module = rownames(p_mat),
    col = list(Module = anno_col),
    show_legend = F,
    show_annotation_name = F
  )
  # Column annotation
  anno_col <- colnames(p_mat)
  names(anno_col) <- colnames(p_mat)
  column_anno = HeatmapAnnotation(
    Module = colnames(p_mat),
    col = list(Module = anno_col),
    show_legend = F,
    show_annotation_name = F
  )
  
  for(i in 1:nrow(p_mat)){
    rownames(p_mat)[i] <- paste0(
      rownames(p_mat)[i]," (",length(mods1[[rownames(p_mat)[i]]]$EnsemblID),")"
    )
    rownames(n_mat)[i] <- paste0(
      rownames(n_mat)[i]," (",length(mods1[[rownames(n_mat)[i]]]$EnsemblID),")"
    )
  }
  
  for(i in 1:ncol(p_mat)){
    colnames(p_mat)[i] <- paste0(
      colnames(p_mat)[i]," (",length(mods2[[colnames(p_mat)[i]]]$EnsemblID),")"
    )
    colnames(n_mat)[i] <- paste0(
      colnames(n_mat)[i]," (",length(mods2[[colnames(n_mat)[i]]]$EnsemblID),")"
    )
  }
  
  
  # Heatmap
  hm <- Heatmap(
    p_mat,
    col = colorRamp2(seq(0, 100, 25), brewer.pal(name = "YlOrRd", n = 5)),
    #col = colorRamp2(breaks = c(0, 25, 50, 75, 100), colors = c("white", "yellow", "orange", "red", "darkred")),
    heatmap_legend_param = list(title = "-log10 p-value", legend_direction = "vertical"),
    bottom_annotation = column_anno,
    left_annotation = row_anno,
    show_column_names = T,
    show_row_names = T,
    row_names_side = "left",
    row_dend_side = "right",
    row_title = title1,
    row_title_side = "left",
    column_title = title2,
    column_title_side = "bottom",
    cell_fun = function(j, i, x, y, w, h, col){
      grid.text(n_mat[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
    }
  )
  return(hm)
}

preservationPlots <- function(dat){
  dat <- cbind(dat$preservation$observed$ref.Set1$inColumnsAlsoPresentIn.Set2,
               dat$preservation$Z$ref.Set1$inColumnsAlsoPresentIn.Set2)
  dat <- dat[,!duplicated(colnames(dat))]
  dat <- rownames_to_column(dat, var = "Module")
  dat <- dat %>% filter(!Module %in% c("grey", "random_sample"))
  color_values <- dat$Module
  names(color_values) <- color_values
  p1 <- ggplot(data = dat, aes(x = moduleSize, y = Zsummary.pres, fill = Module, label = Module)) +
    geom_point(shape = 21, size = 5, color = "black") +
    geom_text(check_overlap = T, vjust = 0, nudge_y = 0.75) +
    scale_fill_manual(values = color_values) +
    geom_hline(yintercept = 2, color = "blue", linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = 10, color = "darkgreen", linetype = "dashed", linewidth = 0.7) +
    labs(x = "Module size", y = "Preservation Z summary") +
    theme_bw() +
    theme(legend.position = "none", 
          text = element_text(size = 20))
  p2 <- ggplot(data = dat, aes(x = moduleSize, y = medianRank.pres, fill = Module, label = Module)) +
    geom_point(shape = 21, size = 5, color = "black") +
    geom_text(check_overlap = T, vjust = 0, nudge_y = 0.75) +
    scale_y_reverse() +
    scale_fill_manual(values = color_values) +
    labs(x = "Module size", y = "Preservation median rank") +
    theme_bw() +
    theme(legend.position = "none", 
          text = element_text(size = 20))
  grid.arrange(p1, p2, ncol = 2)
}

#' # Heatmap

# Reformat module genes for drawing heatmap
modules_list <- list()
modules_list$Skin <- lapply(X = unique(modules$Skin$Module), FUN = function(x){
  modules$Skin %>% filter(Module == x)
})
names(modules_list$Skin) <- unique(modules$Skin$Module)
modules_list$Blood <- lapply(X = unique(modules$Blood$Module), FUN = function(x){
  modules$Blood %>% filter(Module == x)
})
names(modules_list$Blood) <- unique(modules$Blood$Module)

# Draw heatmap
hm <- overlapHeatmap(
  mods1 = modules_list$Skin, 
  mods2 = modules_list$Blood, 
  title1 = "Skin modules", 
  title2 = "Blood modules"
)

#+ fig.width = 14, fig.height = 9
draw(hm)

# Save
png(paste0(output_path,"/heatmap.png"), width = 14, height = 9, units = "in", res = 300)
draw(hm)
dev.off()

#' # Preservation plots

# Skin-blood
dat <- readRDS("results/WGCNA/10_Module_preservation_analysis/Skin-blood/mp.rds")

#+ fig.width = 12, fig.height = 7
preservationPlots(dat = dat)

# Save
png(paste0(output_path,"/Skin-blood.png"), width = 12, height = 7, units = "in", res = 300)
preservationPlots(dat = dat)
dev.off()

# Blood-skin
dat <- readRDS("results/WGCNA/10_Module_preservation_analysis/Blood-skin/mp.rds")

#+ fig.width = 12, fig.height = 7
preservationPlots(dat = dat)

# Save
png(paste0(output_path,"/Blood-skin.png"), width = 12, height = 7, units = "in", res = 300)
preservationPlots(dat = dat)
dev.off()

#' # Session information

sessionInfo()