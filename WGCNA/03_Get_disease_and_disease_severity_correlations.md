Correlation of disease endotypes and disease severity with module
eigengenes and latent factors
================
Ashley Rider
2024-10-04

- <a href="#preliminaries" id="toc-preliminaries">Preliminaries</a>
  - <a href="#load-packages" id="toc-load-packages">Load packages</a>
  - <a href="#create-output-directory"
    id="toc-create-output-directory">Create output directory</a>
- <a href="#skin" id="toc-skin">Skin</a>
  - <a href="#create-output-sub-directory"
    id="toc-create-output-sub-directory">Create output sub-directory</a>
  - <a href="#load-data" id="toc-load-data">Load data</a>
    - <a href="#clinical-data" id="toc-clinical-data">Clinical data</a>
    - <a href="#module-eigengenes" id="toc-module-eigengenes">Module
      eigengenes</a>
    - <a href="#latent-factors" id="toc-latent-factors">Latent factors</a>
  - <a href="#create-traits" id="toc-create-traits">Create traits</a>
  - <a href="#trait-correlations" id="toc-trait-correlations">Trait
    correlations</a>
- <a href="#blood" id="toc-blood">Blood</a>
  - <a href="#create-output-sub-directory-1"
    id="toc-create-output-sub-directory-1">Create output sub-directory</a>
  - <a href="#load-data-1" id="toc-load-data-1">Load data</a>
    - <a href="#clinical-data-1" id="toc-clinical-data-1">Clinical data</a>
    - <a href="#module-eigengenes-1" id="toc-module-eigengenes-1">Module
      eigengenes</a>
    - <a href="#latent-factors-1" id="toc-latent-factors-1">Latent factors</a>
  - <a href="#create-traits-1" id="toc-create-traits-1">Create traits</a>
  - <a href="#trait-correlations-1" id="toc-trait-correlations-1">Trait
    correlations</a>
- <a href="#session-information" id="toc-session-information">Session
  information</a>

Here, we identify associations between clinical traits and gene
expression profiles by correlating these traits with module eigengenes
and latent factors.

# Preliminaries

## Load packages

``` r
library(tidyverse)
library(WGCNA)
library(reshape2)
```

## Create output directory

``` r
output_directory <- "results/WGCNA/03_Get_disease_and_disease_severity_correlations"
dir.create(output_directory)
```

# Skin

## Create output sub-directory

``` r
output_directory2 <- paste0(output_directory,"/Skin")
dir.create(output_directory2)
```

## Load data

### Clinical data

We load the clinical data for each cohort.

``` r
# Load clinical data for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
clin <- list()
clin$Discovery_skin <- read.delim("data/clinical_data/PSORT-D_Skin_Clinical_Data_01-Apr-20.txt")
clin$Replication_skin <- read.delim("data/clinical_data/PSORT-R_Skin_Clinical_Data_01-Apr-20.txt") %>%
  dplyr::filter(!Sample_id %in% c("3A", "3B", "4A", "4B"))
```

Preliminary analysis revealed mislabelling of some lesional and
non-lesional samples, which were swapped for several patients. We
correct these swaps here.

``` r
# Sample swaps based on consensus of S100A7 and S100A9 expression
clin$Discovery_skin <- clin$Discovery_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0005", "23012-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0006", "23012-1206"), "Lesional", Tissue))
clin$Replication_skin <- clin$Replication_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0005", "5031-1205", "6049-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0006", "5031-1206", "6049-1206"), "Lesional", Tissue))
```

Preliminary analysis also revealed that some PASI (**P**soriasis
**A**rea and **S**everity **I**ndex) scores were misrecorded in the
clinical data. We correct these here.

``` r
clin$Replication_skin <- clin$Replication_skin %>%
  mutate(PASI = if_else(Patient_id == "P.6044" & Time == "wk01", 5.8, PASI)) %>%
  mutate(PASI = if_else(Patient_id == "P.60006" & Time == "wk01", 7, PASI))
```

We’ll also define traits that indicate whether or not the patients are
biologic naive or, more specifically, anti-TNF naive.

``` r
for(i in 1:length(clin)){
  clin[[i]] <- clin[[i]] %>%
    mutate(Biologic_Naive = if_else(Adalimumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Biologic_Naive = if_else(Ustekinumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Biologic_Naive = if_else(Efalizumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Biologic_Naive = if_else(Etanercept == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Biologic_Naive = if_else(Infliximab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Biologic_Naive = if_else(Secukinumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
    mutate(Anti_TNF_Naive = "Yes") %>%
    mutate(Anti_TNF_Naive = if_else(Adalimumab == "Stopped_before_wk00", "No", Anti_TNF_Naive)) %>%
    mutate(Anti_TNF_Naive = if_else(Etanercept == "Stopped_before_wk00", "No", Anti_TNF_Naive))
}
```

Finally we select the columns we’ll need for the analysis.

``` r
for(i in 1:length(clin)){
  clin[[i]] <- clin[[i]] %>%
    select(Sample_id, Patient_id, Tissue, Drug, Time, PASI,
           Onset_type, Age_of_onset, Biologic_Naive, Anti_TNF_Naive,
           Psoriatic_Arthritis, Gender, Age, BMI, HLA_C_0602, Cw6_PosNeg,
           wk12_PASI)
}
```

We’ll save the processed clinical data to file for downstream use.

``` r
write.table(clin$Discovery_skin, paste0(output_directory2,"/clin.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(clin$Replication_skin, paste0(output_directory2,"/clin_r.txt"),
            sep = "\t", quote = F, row.names = F)
```

### Module eigengenes

Next, we’ll load the module eigengenes.

``` r
eigen <- list(
  Discovery_skin = read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes.txt"),
  Replication_skin = read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes_r.txt")
)
```

### Latent factors

Finally, we’ll load the latent factors.

``` r
factors <- list(
  Discovery_skin = read.delim("data/latent_factors/skin_d_m.csv", sep = ","),
  Replication_skin = read.delim("data/latent_factors/skin_r_m.csv", sep = ",")
)
colnames(factors$Discovery_skin)[1] <- "Sample_id"
colnames(factors$Replication_skin)[1] <- "Sample_id"
```

## Create traits

We now need to format the clinical data so that it can be correlated
with the module eigengenes and latent factors. We wrap the workflow to
do this in a function and apply it to the PSORT-D and PSORT-R data.

``` r
createTraitsSkin <- function(clin_dat){
  trait_anno <- list()
  # Disease traits
  trait_anno$Disease <- expand.grid(
    Trait = c("Age_of_onset", "Onset_type", "Anti_TNF_Naive", "Psoriatic_Arthritis", 
              "Gender", "Age", "BMI", "Cw6_PosNeg"),
    Type = "Disease",
    Tissue = c("Lesional", "Nonlesional"),
    Time = "wk00",
    Drug = "Both drugs",
    stringsAsFactors = F
  )
  trait_anno$Disease <- trait_anno$Disease %>%
    mutate(Name = paste0(Tissue,"_wk00_",Trait)) %>%
    mutate(Name = gsub("Lesional", "LS", Name)) %>%
    mutate(Name = gsub("Nonlesional", "NL", Name)) %>%
    mutate(Short_name = case_when(Trait == "Age_of_onset" ~ "Age of onset",
                                  Trait == "Onset_type" ~ "Onset type",
                                  Trait == "Anti_TNF_Naive" ~ "TNFi-naive",
                                  Trait == "Psoriatic_Arthritis" ~ "PsA",
                                  Trait == "Gender" ~ "Sex",
                                  Trait == "Age" ~ "Age",
                                  Trait == "BMI" ~ "BMI",
                                  Trait == "Cw6_PosNeg" ~ "Cw6"))
  for(i in 1:nrow(trait_anno$Disease)){
    name <- trait_anno$Disease[i,"Name"]
    trait <- trait_anno$Disease[i,"Trait"]
    tissue <- trait_anno$Disease[i,"Tissue"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Time == "wk00") %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Disease severity traits
  trait_anno$Disease_severity <- expand.grid(
    Trait = "PASI",
    Type = "Disease severity",
    Tissue = c("Lesional", "Nonlesional"),
    Time = "All times",
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Disease_severity <- trait_anno$Disease_severity %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_PASI")) %>%
    mutate(Name = gsub("Lesional", "LS", Name)) %>%
    mutate(Name = gsub("Nonlesional", "NL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = Trait)
  for(i in 1:nrow(trait_anno$Disease_severity)){
    name <- trait_anno$Disease_severity[i,"Name"]
    trait <- trait_anno$Disease_severity[i,"Trait"]
    tissue <- trait_anno$Disease_severity[i,"Tissue"]
    drug <- trait_anno$Disease_severity[i,"Drug"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Drug == drug) %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Bind together trait annotation data frames
  trait_anno <- bind_rows(trait_anno)
  # Select trait columns needed for analysis
  clin_dat <- clin_dat %>%
    select(all_of(c("Sample_id", trait_anno$Name)))
  # Re-level some field and make trait columns numeric
  clin_dat[clin_dat == "Late-onset"] <- 0
  clin_dat[clin_dat == "Early-onset"] <- 1
  clin_dat[clin_dat == "Yes"] <- 0
  clin_dat[clin_dat == "No"] <- 1
  clin_dat[clin_dat == "M"] <- 1
  clin_dat[clin_dat == "F"] <- 0
  clin_dat[clin_dat == "Positive"] <- 1
  clin_dat[clin_dat == "Negative"] <- 0
  for(i in 2:ncol(clin_dat)){
    if(is.character(clin_dat[,i])){
      clin_dat[,i] <- as.numeric(clin_dat[,i])
    }
  }
  return(list(traits = clin_dat, trait_anno = trait_anno))
}

traits <- list(
  Discovery_skin = createTraitsSkin(clin_dat = clin$Discovery_skin),
  Replication_skin = createTraitsSkin(clin_dat = clin$Replication_skin)
)
```

As part of this workflow we created an annotation table which gives
details about each trait, including the tissues, time points and drug
cohorts that comprise them. Let’s examine this table.

``` r
traits$Discovery_skin$trait_anno
```

    ##                  Trait             Type      Tissue      Time        Drug
    ## 1         Age_of_onset          Disease    Lesional      wk00  Both drugs
    ## 2           Onset_type          Disease    Lesional      wk00  Both drugs
    ## 3       Anti_TNF_Naive          Disease    Lesional      wk00  Both drugs
    ## 4  Psoriatic_Arthritis          Disease    Lesional      wk00  Both drugs
    ## 5               Gender          Disease    Lesional      wk00  Both drugs
    ## 6                  Age          Disease    Lesional      wk00  Both drugs
    ## 7                  BMI          Disease    Lesional      wk00  Both drugs
    ## 8           Cw6_PosNeg          Disease    Lesional      wk00  Both drugs
    ## 9         Age_of_onset          Disease Nonlesional      wk00  Both drugs
    ## 10          Onset_type          Disease Nonlesional      wk00  Both drugs
    ## 11      Anti_TNF_Naive          Disease Nonlesional      wk00  Both drugs
    ## 12 Psoriatic_Arthritis          Disease Nonlesional      wk00  Both drugs
    ## 13              Gender          Disease Nonlesional      wk00  Both drugs
    ## 14                 Age          Disease Nonlesional      wk00  Both drugs
    ## 15                 BMI          Disease Nonlesional      wk00  Both drugs
    ## 16          Cw6_PosNeg          Disease Nonlesional      wk00  Both drugs
    ## 17                PASI Disease severity    Lesional All times  Adalimumab
    ## 18                PASI Disease severity Nonlesional All times  Adalimumab
    ## 19                PASI Disease severity    Lesional All times Ustekinumab
    ## 20                PASI Disease severity Nonlesional All times Ustekinumab
    ##                           Name   Short_name
    ## 1         LS_wk00_Age_of_onset Age of onset
    ## 2           LS_wk00_Onset_type   Onset type
    ## 3       LS_wk00_Anti_TNF_Naive   TNFi-naive
    ## 4  LS_wk00_Psoriatic_Arthritis          PsA
    ## 5               LS_wk00_Gender          Sex
    ## 6                  LS_wk00_Age          Age
    ## 7                  LS_wk00_BMI          BMI
    ## 8           LS_wk00_Cw6_PosNeg          Cw6
    ## 9         NL_wk00_Age_of_onset Age of onset
    ## 10          NL_wk00_Onset_type   Onset type
    ## 11      NL_wk00_Anti_TNF_Naive   TNFi-naive
    ## 12 NL_wk00_Psoriatic_Arthritis          PsA
    ## 13              NL_wk00_Gender          Sex
    ## 14                 NL_wk00_Age          Age
    ## 15                 NL_wk00_BMI          BMI
    ## 16          NL_wk00_Cw6_PosNeg          Cw6
    ## 17                 LS_ADA_PASI         PASI
    ## 18                 NL_ADA_PASI         PASI
    ## 19                 LS_UST_PASI         PASI
    ## 20                 NL_UST_PASI         PASI

We’ll save this data to file for later use.

``` r
# Save
write.table(traits$Discovery_skin$traits, paste0(output_directory2,"/Traits.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$Replication_skin$traits, paste0(output_directory2,"/Traits_r.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$Discovery_skin$trait_anno, paste0(output_directory2,"/Traits_anno.txt"),
            sep = "\t", quote = F, row.names = F)
```

## Trait correlations

Now we can correlate the traits defined above with the module eigengenes
and latent factors. Again, we define a function to do this and apply it
to the modules and factors for PSORT-D and PSORT-R separately.
Specifically, this calculates pearson correlation.

``` r
# Performs module-trait correlation
modTraitCor <- function(eigen, traits){
  # Use sample IDs as row names
  eigen <- column_to_rownames(eigen, "Sample_id")
  traits <- column_to_rownames(traits, "Sample_id")
  traits <- traits[rownames(eigen),]
  # Calculate module trait correlations and associated p-values
  cp = WGCNA::corAndPvalue(eigen, traits)
  # Merge correlation coefficients and p-values for each module-trait pair into one data frame
  cor_dat <- melt(cp$cor)
  p_dat <- melt(cp$p)
  dat <- merge(cor_dat, p_dat, by = c("Var1", "Var2"))
  colnames(dat) <- c("Module", "Trait", "Cor", "P.Value")
  # Convert factor columns to character columns
  dat$Module <- as.character(dat$Module)
  dat$Trait <- as.character(dat$Trait)
  return(dat)
}

dat <- list(
  modTraitCor(eigen = eigen$Discovery_skin, traits = traits$Discovery_skin$traits),
  modTraitCor(eigen = eigen$Replication_skin, traits = traits$Replication_skin$traits),
  modTraitCor(eigen = factors$Discovery_skin, traits = traits$Discovery_skin$traits),
  modTraitCor(eigen = factors$Replication_skin, traits = traits$Replication_skin$traits)
)

# Trait types to merge with correlation data
trait_types <- traits$Discovery_skin$trait_anno %>% 
  select(Name, Type) %>% 
  rename(Trait = Name)

for(i in 1:length(dat)){
  # Merge with trait types
  dat[[i]] <- merge(dat[[i]], trait_types, by = "Trait")
  # Adjust all p-values
  dat[[i]]$padj_all <- p.adjust(dat[[i]]$P.Value, method = "fdr")
  # Adjust p-values for trait types separately
  dat[[i]] <- dat[[i]] %>% 
    group_by(Type) %>% 
    mutate(padj_type = p.adjust(P.Value, method = "fdr")) %>% 
    as.data.frame()
}

# Merge module-trait correlations for Discovery and Replication
module_dat <- merge(dat[[1]], dat[[2]], by = c("Module", "Trait", "Type"), suffixes = c("_d", "_r"))

# Merge module-trait correlations for Discovery and Replication
factor_dat <- merge(dat[[3]], dat[[4]], by = c("Module", "Trait", "Type"), suffixes = c("_d", "_r"))

# Save
write.table(module_dat, paste0(output_directory2,"/Module-trait_correlations.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(factor_dat, paste0(output_directory2,"/Factor-trait_correlations.txt"),
            sep = "\t", quote = F, row.names = F)
```

# Blood

## Create output sub-directory

``` r
output_directory2 <- paste0(output_directory,"/Blood")
dir.create(output_directory2)
```

## Load data

### Clinical data

We load the clinical data for each cohort.

``` r
# Load clinical data for PSORT-D blood samples
clin <- read.delim("data/clinical_data/PSORT-D_Blood_Clinical_Data_01-Apr-20.txt")
```

Preliminary analysis revealed some outlier samples in the blood data
that we don’t want to include in the analysis. There are also some
smaples from other cohorts that we don’t want to include. Here we will
read in a file containing the IDs of samples that we want to analyse;
we’ll use this to subset the clinical data.

``` r
samples <- read.delim("data/clinical_data/PSORT-D_Blood_analysis_samples.txt")

clin <- clin %>% filter(Sample_id %in% samples$Sample_id)
```

Preliminary analysis also revealed that some PASI scores were
misrecorded in the clinical data. We correct these here.

``` r
clin <- clin %>%
  mutate(PASI = if_else(Patient_id == "P.6040" & Time == "wk04", 6.3, PASI)) %>%
  mutate(PASI = if_else(Patient_id == "P.5003" & Time == "wk04", 17.4, PASI))
```

We’ll also define traits that indicate whether or not the patients are
biologic naive or, more specifically, anti-TNF naive.

``` r
clin <- clin %>%
  mutate(Biologic_Naive = if_else(Adalimumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Biologic_Naive = if_else(Ustekinumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Biologic_Naive = if_else(Efalizumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Biologic_Naive = if_else(Etanercept == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Biologic_Naive = if_else(Infliximab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Biologic_Naive = if_else(Secukinumab == "Stopped_before_wk00", "No", Biologic_Naive)) %>%
  mutate(Anti_TNF_Naive = "Yes") %>%
  mutate(Anti_TNF_Naive = if_else(Adalimumab == "Stopped_before_wk00", "No", Anti_TNF_Naive)) %>%
  mutate(Anti_TNF_Naive = if_else(Etanercept == "Stopped_before_wk00", "No", Anti_TNF_Naive))
```

Finally we select the columns we’ll need for the analysis.

``` r
clin <- clin %>%
  select(Sample_id, Patient_id, Tissue, Drug, Time, PASI,
         Onset_type, Age_of_onset, Biologic_Naive, Anti_TNF_Naive,
         Psoriatic_Arthritis, Gender, Age, BMI, HLA_C_0602, Cw6_PosNeg,
         wk12_PASI)
```

We’ll save the processed clinical data to file for downstream use.

``` r
write.table(clin, paste0(output_directory2,"/clin.txt"),
            sep = "\t", quote = F, row.names = F)
```

### Module eigengenes

Next, we’ll load the module eigengenes.

``` r
eigen <- read.delim("results/WGCNA/01_Module_identification/Blood/eigengenes.txt")
```

### Latent factors

Finally, we’ll load the latent factors.

``` r
factors <- read.delim("data/latent_factors/blood_m.csv", sep = ",")
colnames(factors)[1] <- "Sample_id"
```

## Create traits

We now need to format the clinical data so that it can be correlated
with the module eigengenes and latent factors.

``` r
createTraitsBlood <- function(clin_dat){
  trait_anno <- list()
  # Disease traits
  trait_anno$Disease <- expand.grid(
    Trait = c("Age_of_onset", "Anti_TNF_Naive", "Psoriatic_Arthritis", 
              "Gender", "Age", "BMI", "Cw6_PosNeg"),
    Type = "Disease",
    Tissue = "Blood",
    Time = "wk00",
    Drug = "Both drugs",
    stringsAsFactors = F
  )
  trait_anno$Disease <- trait_anno$Disease %>%
    mutate(Name = paste0(Tissue,"_wk00_",Trait)) %>%
    mutate(Name = gsub("Blood", "BL", Name)) %>%
    mutate(Short_name = case_when(Trait == "Age_of_onset" ~ "Age of onset",
                                  Trait == "Onset_type" ~ "Onset type",
                                  Trait == "Anti_TNF_Naive" ~ "TNFi-naive",
                                  Trait == "Psoriatic_Arthritis" ~ "PsA",
                                  Trait == "Gender" ~ "Sex",
                                  Trait == "Age" ~ "Age",
                                  Trait == "BMI" ~ "BMI",
                                  Trait == "Cw6_PosNeg" ~ "Cw6"))
  for(i in 1:nrow(trait_anno$Disease)){
    name <- trait_anno$Disease[i,"Name"]
    trait <- trait_anno$Disease[i,"Trait"]
    tissue <- trait_anno$Disease[i,"Tissue"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Time == "wk00") %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Disease severity traits
  trait_anno$Disease_severity <- expand.grid(
    Trait = "PASI",
    Type = "Disease severity",
    Tissue = "Blood",
    Time = "All times",
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Disease_severity <- trait_anno$Disease_severity %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_PASI")) %>%
    mutate(Name = gsub("Blood", "BL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = Trait)
  for(i in 1:nrow(trait_anno$Disease_severity)){
    name <- trait_anno$Disease_severity[i,"Name"]
    trait <- trait_anno$Disease_severity[i,"Trait"]
    tissue <- trait_anno$Disease_severity[i,"Tissue"]
    drug <- trait_anno$Disease_severity[i,"Drug"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Drug == drug) %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Bind together trait annotation data frames
  trait_anno <- bind_rows(trait_anno)
  # Select trait columns needed for analysis
  clin_dat <- clin_dat %>%
    select(all_of(c("Sample_id", trait_anno$Name)))
  # Re-level some field and make trait columns numeric
  clin_dat[clin_dat == "Late-onset"] <- 0
  clin_dat[clin_dat == "Early-onset"] <- 1
  clin_dat[clin_dat == "Yes"] <- 0
  clin_dat[clin_dat == "No"] <- 1
  clin_dat[clin_dat == "M"] <- 1
  clin_dat[clin_dat == "F"] <- 0
  clin_dat[clin_dat == "Positive"] <- 1
  clin_dat[clin_dat == "Negative"] <- 0
  for(i in 2:ncol(clin_dat)){
    if(is.character(clin_dat[,i])){
      clin_dat[,i] <- as.numeric(clin_dat[,i])
    }
  }
  return(list(traits = clin_dat, trait_anno = trait_anno))
}

traits <- createTraitsBlood(clin_dat = clin)
```

As part of this workflow we created an annotation table which gives
details about each trait, including the tissues, time points and drug
cohorts that comprise them. Let’s examine this table.

``` r
traits$trait_anno
```

    ##                 Trait             Type Tissue      Time        Drug
    ## 1        Age_of_onset          Disease  Blood      wk00  Both drugs
    ## 2      Anti_TNF_Naive          Disease  Blood      wk00  Both drugs
    ## 3 Psoriatic_Arthritis          Disease  Blood      wk00  Both drugs
    ## 4              Gender          Disease  Blood      wk00  Both drugs
    ## 5                 Age          Disease  Blood      wk00  Both drugs
    ## 6                 BMI          Disease  Blood      wk00  Both drugs
    ## 7          Cw6_PosNeg          Disease  Blood      wk00  Both drugs
    ## 8                PASI Disease severity  Blood All times  Adalimumab
    ## 9                PASI Disease severity  Blood All times Ustekinumab
    ##                          Name   Short_name
    ## 1        BL_wk00_Age_of_onset Age of onset
    ## 2      BL_wk00_Anti_TNF_Naive   TNFi-naive
    ## 3 BL_wk00_Psoriatic_Arthritis          PsA
    ## 4              BL_wk00_Gender          Sex
    ## 5                 BL_wk00_Age          Age
    ## 6                 BL_wk00_BMI          BMI
    ## 7          BL_wk00_Cw6_PosNeg          Cw6
    ## 8                 BL_ADA_PASI         PASI
    ## 9                 BL_UST_PASI         PASI

We’ll save this data to file for later use.

``` r
# Save
write.table(traits$traits, paste0(output_directory2,"/Traits.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$trait_anno, paste0(output_directory2,"/Traits_anno.txt"),
            sep = "\t", quote = F, row.names = F)
```

## Trait correlations

Now we can correlate the traits defined above with the module eigengenes
and latent factors.

``` r
dat <- list(
  modTraitCor(eigen = eigen, traits = traits$traits),
  modTraitCor(eigen = factors, traits = traits$traits)
)

# Trait types to merge with correlation data
trait_types <- traits$trait_anno %>% 
  select(Name, Type) %>% 
  rename(Trait = Name)

for(i in 1:length(dat)){
  # Merge with trait types
  dat[[i]] <- merge(dat[[i]], trait_types, by = "Trait")
  # Adjust all p-values
  dat[[i]]$padj_all <- p.adjust(dat[[i]]$P.Value, method = "fdr")
  # Adjust p-values for trait types separately
  dat[[i]] <- dat[[i]] %>% 
    group_by(Type) %>% 
    mutate(padj_type = p.adjust(P.Value, method = "fdr")) %>% 
    as.data.frame()
  # Change column names
  dat[[i]] <- dat[[i]] %>% select(Trait, Module, Type, Cor, P.Value, padj_all, padj_type)
  colnames(dat[[i]]) <- c("Trait", "Module", "Type", "Cor_d", "P.Value_d", "padj_all_d", "padj_type_d")
}

# Save
write.table(dat[[1]], paste0(output_directory2,"/Module-trait_correlations.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(dat[[2]], paste0(output_directory2,"/Factor-trait_correlations.txt"),
            sep = "\t", quote = F, row.names = F)
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
    ##  [1] reshape2_1.4.4        WGCNA_1.72-1          fastcluster_1.2.3    
    ##  [4] dynamicTreeCut_1.63-1 lubridate_1.9.2       forcats_1.0.0        
    ##  [7] stringr_1.5.0         dplyr_1.1.0           purrr_1.0.1          
    ## [10] readr_2.1.4           tidyr_1.3.0           tibble_3.2.0         
    ## [13] ggplot2_3.4.2         tidyverse_2.0.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7           matrixStats_0.63.0     bit64_4.0.5           
    ##  [4] doParallel_1.0.17      httr_1.4.5             rprojroot_2.0.3       
    ##  [7] GenomeInfoDb_1.34.9    tools_4.2.3            backports_1.4.1       
    ## [10] utf8_1.2.3             R6_2.5.1               rpart_4.1.19          
    ## [13] Hmisc_5.0-1            DBI_1.1.3              BiocGenerics_0.44.0   
    ## [16] colorspace_2.1-0       nnet_7.3-18            withr_2.5.0           
    ## [19] tidyselect_1.2.0       gridExtra_2.3          preprocessCore_1.60.2 
    ## [22] bit_4.0.5              compiler_4.2.3         cli_3.6.0             
    ## [25] Biobase_2.58.0         htmlTable_2.4.1        scales_1.2.1          
    ## [28] checkmate_2.1.0        digest_0.6.31          foreign_0.8-84        
    ## [31] rmarkdown_2.20         XVector_0.38.0         base64enc_0.1-3       
    ## [34] pkgconfig_2.0.3        htmltools_0.5.4        fastmap_1.1.1         
    ## [37] htmlwidgets_1.6.2      rlang_1.1.0            impute_1.72.3         
    ## [40] rstudioapi_0.14        RSQLite_2.3.0          generics_0.1.3        
    ## [43] RCurl_1.98-1.10        magrittr_2.0.3         GO.db_3.16.0          
    ## [46] GenomeInfoDbData_1.2.9 Formula_1.2-5          Matrix_1.6-1.1        
    ## [49] Rcpp_1.0.10            munsell_0.5.0          S4Vectors_0.36.2      
    ## [52] fansi_1.0.4            lifecycle_1.0.3        stringi_1.7.12        
    ## [55] yaml_2.3.7             zlibbioc_1.44.0        plyr_1.8.8            
    ## [58] grid_4.2.3             blob_1.2.4             parallel_4.2.3        
    ## [61] crayon_1.5.2           lattice_0.20-45        Biostrings_2.66.0     
    ## [64] splines_4.2.3          hms_1.1.2              KEGGREST_1.38.0       
    ## [67] knitr_1.42             pillar_1.8.1           codetools_0.2-19      
    ## [70] stats4_4.2.3           glue_1.6.2             evaluate_0.20         
    ## [73] data.table_1.14.8      png_0.1-8              vctrs_0.6.0           
    ## [76] tzdb_0.3.0             foreach_1.5.2          gtable_0.3.2          
    ## [79] cachem_1.0.7           xfun_0.39              survival_3.5-3        
    ## [82] iterators_1.0.14       AnnotationDbi_1.60.2   memoise_2.0.1         
    ## [85] IRanges_2.32.0         cluster_2.1.4          timechange_0.2.0      
    ## [88] ellipsis_0.3.2
