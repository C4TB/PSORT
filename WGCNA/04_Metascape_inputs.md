Format inputs for Metascape
================
Ashley Rider
2024-10-04

- <a href="#load-packages" id="toc-load-packages">Load packages</a>
- <a href="#create-output-directory"
  id="toc-create-output-directory">Create output directory</a>
- <a href="#skin" id="toc-skin">Skin</a>
  - <a href="#load-data" id="toc-load-data">Load data</a>
  - <a href="#create-metascape-input" id="toc-create-metascape-input">Create
    Metascape input</a>
    - <a href="#modules" id="toc-modules">Modules</a>
    - <a href="#factors" id="toc-factors">Factors</a>
- <a href="#blood" id="toc-blood">Blood</a>
  - <a href="#load-data-1" id="toc-load-data-1">Load data</a>
  - <a href="#create-metascape-input-1"
    id="toc-create-metascape-input-1">Create Metascape input</a>
    - <a href="#modules-1" id="toc-modules-1">Modules</a>
    - <a href="#factors-1" id="toc-factors-1">Factors</a>
- <a href="#session-information" id="toc-session-information">Session
  information</a>

Here, we format the inputs for
[Metascape](https://metascape.org/gp/index.html#/main/step1) ([Zhou et
al, 2019](https://www.nature.com/articles/s41467-019-09234-6)), which
will be used to carry out functional enrichment analysis of the WGCNA
module and latent factor genes.

# Load packages

``` r
library(tidyverse)
```

# Create output directory

``` r
output_directory <- "results/WGCNA/04_Metascape_inputs"
dir.create(output_directory)
```

# Skin

## Load data

``` r
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
```

## Create Metascape input

### Modules

``` r
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
```

### Factors

``` r
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
```

# Blood

## Load data

``` r
# Module genes
modules <- read.delim("results/WGCNA/01_Module_identification/Blood/modules.txt")

# Factor genes
factors <- read.delim("data/latent_factors/factor_genes_blood.csv", sep = ",")
colnames(factors) <- c("X", "Factor", "EnsemblID", "Coef", "GeneSymbol")
```

## Create Metascape input

### Modules

``` r
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
```

### Factors

``` r
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
```

# Session information

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.0    
    ##  [5] purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.0   
    ##  [9] ggplot2_3.4.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] pillar_1.8.1     compiler_4.2.3   tools_4.2.3      digest_0.6.31   
    ##  [5] timechange_0.2.0 evaluate_0.20    lifecycle_1.0.3  gtable_0.3.2    
    ##  [9] pkgconfig_2.0.3  rlang_1.1.0      cli_3.6.0        rstudioapi_0.14 
    ## [13] yaml_2.3.7       xfun_0.39        fastmap_1.1.1    withr_2.5.0     
    ## [17] knitr_1.42       generics_0.1.3   vctrs_0.6.0      hms_1.1.2       
    ## [21] rprojroot_2.0.3  grid_4.2.3       tidyselect_1.2.0 glue_1.6.2      
    ## [25] R6_2.5.1         fansi_1.0.4      rmarkdown_2.20   tzdb_0.3.0      
    ## [29] magrittr_2.0.3   scales_1.2.1     htmltools_0.5.4  ellipsis_0.3.2  
    ## [33] colorspace_2.1-0 utf8_1.2.3       stringi_1.7.12   munsell_0.5.0
