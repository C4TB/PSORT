library(tidyverse)

output_directory <- "results/WGCNA/05_Module_factor_names"
dir.create(output_directory)

skin_names <- c(
  black = "Translation module",
  blue = "Muscle contraction/WNT signalling module",
  brown = "Regulatory module",
  darkgreen = "Cornification module 1",
  darkgrey = "Histone modification module",
  darkmagenta = "Ion transport module",
  darkolivegreen = "RORA-aligned module",
  darkred = "PI3K-AKT-mTOR signalling module",
  darkturquoise = "Chromatin modification module",
  green = "mRNA processing module",
  greenyellow = "Cell cycle module 1 (mitosis)",
  grey60 = "ECM module 1",
  lightcyan = "Notch signalling module",
  lightgreen = "Cell-cell interaction module",
  lightyellow = "Obesity-associated module",
  orange = "Lipid metabolism module",
  paleturquoise = "Sphingolipid metabolism module",
  purple = "DNA binding and repair module",
  red = "Cell proliferation module",
  salmon = "Cell cycle module 2 (checkpoints)",
  skyblue = "Cornification module 2",
  steelblue = "Innate/adaptive immunity module",
  tan = "Oxidative phosphorylation module",
  turquoise = "Cytokine signalling module",
  violet = "Complement module",
  white = "Antigen presentation module",
  yellow = "ECM module 2"
)
skin_names <- skin_names %>% as.data.frame() %>% rownames_to_column("Module")
colnames(skin_names)[2] <- "Name"

blood_names <- c(
  cadetblue = "Chromatin modification module",
  chocolate = "TLR signalling module",
  coral = "VEGF signalling module",
  darkorchid = "B-cell module",
  darkviolet = "Intracellular biosynthesis module",
  dodgerblue = "RNA processing module",
  firebrick = "Sex-linked module",
  gold = "Translation module",
  khaki = "Innate immune cell module",
  lightcoral = "Mitochondrial assembly module",
  navy = "Innate/adaptive immunity module 1",
  olivedrab = "Oxidative phosphorylation module",
  orangered = "IL17/TNF signalling module",
  plum = "Innate/adaptive immunity module 2",
  sandybrown = "Cell growth and metabolism module",
  sienna = "Autoimmunity module",
  slategrey = "IL17 signalling/T-cell activation module",
  thistle = "Inflammatory cell death module"
)
blood_names <- blood_names %>% as.data.frame() %>% rownames_to_column("Module")
colnames(blood_names)[2] <- "Name"

skin_factor_names <- c(
  factor_1 = "Antimicrobial protein factor",
  factor_2 = "Mixed inflammatory factor",
  factor_4 = "PPAR signalling factor",
  factor_5 = "HLA_aligned ECM factor",
  factor_6 = "ECM-enriched factor",
  factor_8 = "Sex-linked factor",
  factor_9 = "Obesity-associated factor",
  factor_10 = "Fatty acid metabolism factor",
  factor_12 = "Interferon signalling factor",
  factor_13 = "Cornified envelope factor",
  factor_14 = "HLA-enriched factor A",
  factor_20 = "HLA-enriched factor B",
  factor_21 = "HLA-enriched factor C",
  factor_22 = "HLA-enriched factor D"
)
skin_factor_names <- skin_factor_names %>% as.data.frame() %>% rownames_to_column("Module")
colnames(skin_factor_names)[2] <- "Name"

blood_factor_names <- c(
  factor_1 = "Inflammatory and HLA-aligned factor A",
  factor_8 = "Inflammatory and HLA-aligned factor B",
  factor_10 = "Sex-linked factor",
  factor_19 = "Inflammatory and HLA-aligned factor C"
)
blood_factor_names <- blood_factor_names %>% as.data.frame() %>% rownames_to_column("Module")
colnames(blood_factor_names)[2] <- "Name"

write.table(skin_names, paste0(output_directory,"/Skin_module_names.txt"), 
            sep = "\t", row.names = F, quote = F)
write.table(blood_names, paste0(output_directory,"/Blood_module_names.txt"), 
            sep = "\t", row.names = F, quote = F)
write.table(skin_factor_names, paste0(output_directory,"/Skin_factor_names.txt"), 
            sep = "\t", row.names = F, quote = F)
write.table(blood_factor_names, paste0(output_directory,"/Blood_factor_names.txt"), 
            sep = "\t", row.names = F, quote = F)