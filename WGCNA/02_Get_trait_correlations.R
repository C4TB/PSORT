#' ---
#' title: "Correlation of traits with module eigengenes and latent factors"
#' author: Ashley Rider
#' output:
#'    github_document:
#'      toc: TRUE
#' ---

#+ echo = FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 12, fig.height = 9)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

#' Here, we identify associations between clinical traits and gene expression profiles by correlating these traits
#' with module eigengenes and latent factors.
#' 
#' # Preliminaries
#' 
#' ## Load packages

library(tidyverse)
library(WGCNA)
library(reshape2)

#' ## Create output directory

output_directory <- "results/WGCNA/02_Get_trait_correlations"
dir.create(output_directory)

#' # Skin
#' 
#' ## Create output sub-directory

output_directory2 <- paste0(output_directory,"/Skin")
dir.create(output_directory2)

#' ## Load data
#' 
#' ### Clinical data

#' We load the clinical data for each cohort.

# Load clinical data for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
clin <- list()
clin$Discovery_skin <- read.delim("data/clinical_data/PSORT-D_Skin_Clinical_Data_01-Apr-20.txt")
clin$Replication_skin <- read.delim("data/clinical_data/PSORT-R_Skin_Clinical_Data_01-Apr-20.txt") %>%
  dplyr::filter(!Sample_id %in% c("3A", "3B", "4A", "4B"))

#' Preliminary analysis revealed mislabelling of some lesional and non-lesional samples, which were swapped
#' for several patients. We correct these swaps here.

# Sample swaps based on consensus of S100A7 and S100A9 expression
clin$Discovery_skin <- clin$Discovery_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0005", "23012-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0006", "23012-1206"), "Lesional", Tissue))
clin$Replication_skin <- clin$Replication_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0005", "5031-1205", "6049-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0006", "5031-1206", "6049-1206"), "Lesional", Tissue))

#' Preliminary analysis also revealed that some PASI (**P**soriasis **A**rea and **S**everity **I**ndex) scores
#' were misrecorded in the clinical data. We correct these here.

clin$Replication_skin <- clin$Replication_skin %>%
  mutate(PASI = if_else(Patient_id == "P.6044" & Time == "wk01", 5.8, PASI)) %>%
  mutate(PASI = if_else(Patient_id == "P.60006" & Time == "wk01", 7, PASI))

#' We'll also define traits that indicate whether or not the patients are biologic naive or, more specifically,
#' anti-TNF naive.

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

#' Prelinary analysis also revealed that the distribution of delta PASI is very skewed; we correct this here
#' with a double exponential transformation.
  
clin$Discovery_skin$dexpDeltaPASI <- exp(exp(clin$Discovery_skin$DeltaPASI))
clin$Replication_skin$dexpDeltaPASI <- exp(exp(clin$Replication_skin$DeltaPASI))

#' Finally we select the columns we'll need for the analysis.

for(i in 1:length(clin)){
  clin[[i]] <- clin[[i]] %>%
    select(Sample_id, Patient_id, Tissue, Drug, Time, PASI, DeltaPASI, dexpDeltaPASI,
           Onset_type, Age_of_onset, Biologic_Naive, Anti_TNF_Naive,
           Psoriatic_Arthritis, Gender, Age, BMI, HLA_C_0602, Cw6_PosNeg,
           wk12_PASI)
}

#' We'll save the processed clinical data to file for downstream use.

write.table(clin$Discovery_skin, paste0(output_directory2,"/clin.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(clin$Replication_skin, paste0(output_directory2,"/clin_r.txt"),
            sep = "\t", quote = F, row.names = F)

#' ### Module eigengenes
#' 
#' Next, we'll load the module eigengenes.

eigen <- list(
  Discovery_skin = read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes.txt"),
  Replication_skin = read.delim("results/WGCNA/01_Module_identification/Skin/eigengenes_r.txt")
)

#' ### Latent factors
#' 
#' Finally, we'll load the latent factors.

factors <- list(
  Discovery_skin = read.delim("data/latent_factors/skin_d_m.csv", sep = ","),
  Replication_skin = read.delim("data/latent_factors/skin_r_m.csv", sep = ",")
)
colnames(factors$Discovery_skin)[1] <- "Sample_id"
colnames(factors$Replication_skin)[1] <- "Sample_id"

#' ## Create traits
#' 
#' We now need to format the clinical data so that it can be correlated with the module eigengenes and 
#' latent factors. We wrap the workflow to do this in a function and apply it to the PSORT-D and PSORT-R data.

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
  # Response traits
  trait_anno$Response <- expand.grid(
    Trait = "dexpDeltaPASI",
    Type = "Response",
    Tissue = c("Lesional", "Nonlesional"),
    Time = c("wk00", "wk01", "wk12"),
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Response <- trait_anno$Response %>%
    filter(!(Tissue == "Nonlesional" & Time == "wk01")) %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_",Time,"_Resp")) %>%
    mutate(Name = gsub("Lesional", "LS", Name)) %>%
    mutate(Name = gsub("Nonlesional", "NL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = Time)
  for(i in 1:nrow(trait_anno$Response)){
    name <- trait_anno$Response[i,"Name"]
    trait <- trait_anno$Response[i,"Trait"]
    tissue <- trait_anno$Response[i,"Tissue"]
    drug <- trait_anno$Response[i,"Drug"]
    time <- trait_anno$Response[i,"Time"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Drug == drug, Time == time) %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Treatment traits
  trait_anno$Treatment <- expand.grid(
    Trait = NA,
    Type = "Treatment",
    Tissue = c("Lesional", "Nonlesional"),
    Time = c("wk01", "wk12"),
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Treatment <- trait_anno$Treatment %>%
    filter(!(Tissue == "Nonlesional" & Time == "wk01")) %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_",Time,"_vs_wk00")) %>%
    mutate(Name = gsub("Lesional", "LS", Name)) %>%
    mutate(Name = gsub("Nonlesional", "NL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = paste0(Time," vs wk00"))
  for(i in 1:nrow(trait_anno$Treatment)){
    name <- trait_anno$Treatment[i,"Name"]
    tissue <- trait_anno$Treatment[i,"Tissue"]
    drug <- trait_anno$Treatment[i,"Drug"]
    time <- trait_anno$Treatment[i,"Time"]
    clin_dat[[name]] <- NA
    wk00_samples <- clin_dat %>% filter(Tissue == tissue, Drug == drug, Time == "wk00") %>% pull(Sample_id)
    treated_samples <- clin_dat %>% filter(Tissue == tissue, Drug == drug, Time == time) %>% pull(Sample_id)
    clin_dat[[name]][which(clin_dat$Sample_id %in% wk00_samples)] <- 0
    clin_dat[[name]][which(clin_dat$Sample_id %in% treated_samples)] <- 1
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

#' As part of this workflow we created an annotation table which gives details about each trait, including the
#' tissues, time points and drug cohorts that comprise them. Let's examine this table.

traits$Discovery_skin$trait_anno

#' We'll save this data to file for later use.

# Save
write.table(traits$Discovery_skin$traits, paste0(output_directory2,"/Traits.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$Replication_skin$traits, paste0(output_directory2,"/Traits_r.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$Discovery_skin$trait_anno, paste0(output_directory2,"/Traits_anno.txt"),
            sep = "\t", quote = F, row.names = F)

#' ## Trait correlations
#' 
#' Now we can correlate the traits defined above with the module eigengenes and latent factors. Again, we define
#' a function to do this and apply it to the modules and factors for PSORT-D and PSORT-R separately. Specifically,
#' this calculates pearson correlation.

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

#' # Blood
#' 
#' ## Create output sub-directory

output_directory2 <- paste0(output_directory,"/Blood")
dir.create(output_directory2)

#' ## Load data
#' 
#' ### Clinical data

#' We load the clinical data for each cohort.

# Load clinical data for PSORT-D blood samples
clin <- read.delim("data/clinical_data/PSORT-D_Blood_Clinical_Data_01-Apr-20.txt")

#' Preliminary analysis revealed some outlier samples in the blood data that we don't want to include in the
#' analysis. There are also some smaples from other cohorts that we don't want to include.
#' Here we will read in a file containing the IDs of samples that we want to analyse; we'll use this to
#' subset the clinical data.

samples <- read.delim("data/clinical_data/PSORT-D_Blood_analysis_samples.txt")

clin <- clin %>% filter(Sample_id %in% samples$Sample_id)

#' Preliminary analysis also revealed that some PASI scores were misrecorded in the clinical data. 
#' We correct these here.

clin <- clin %>%
  mutate(PASI = if_else(Patient_id == "P.6040" & Time == "wk04", 6.3, PASI)) %>%
  mutate(PASI = if_else(Patient_id == "P.5003" & Time == "wk04", 17.4, PASI))

#' We'll also define traits that indicate whether or not the patients are biologic naive or, more specifically,
#' anti-TNF naive.

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

#' Prelinary analysis also revealed that the distribution of delta PASI is very skewed; we correct this here
#' with a double exponential transformation.

clin$dexpDeltaPASI <- exp(exp(clin$DeltaPASI))

#' Finally we select the columns we'll need for the analysis.

clin <- clin %>%
  select(Sample_id, Patient_id, Tissue, Drug, Time, PASI, DeltaPASI, dexpDeltaPASI,
         Onset_type, Age_of_onset, Biologic_Naive, Anti_TNF_Naive,
         Psoriatic_Arthritis, Gender, Age, BMI, HLA_C_0602, Cw6_PosNeg,
         wk12_PASI)

#' We'll save the processed clinical data to file for downstream use.

write.table(clin, paste0(output_directory2,"/clin.txt"),
            sep = "\t", quote = F, row.names = F)

#' ### Module eigengenes
#' 
#' Next, we'll load the module eigengenes.

eigen <- read.delim("results/WGCNA/01_Module_identification/Blood/eigengenes.txt")

#' ### Latent factors
#' 
#' Finally, we'll load the latent factors.

factors <- read.delim("data/latent_factors/blood_m.csv", sep = ",")
colnames(factors)[1] <- "Sample_id"

#' ## Create traits
#' 
#' We now need to format the clinical data so that it can be correlated with the module eigengenes and 
#' latent factors.

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
  # Response traits
  trait_anno$Response <- expand.grid(
    Trait = "dexpDeltaPASI",
    Type = "Response",
    Tissue = "Blood",
    Time = c("wk00", "wk01", "wk04", "wk12"),
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Response <- trait_anno$Response %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_",Time,"_Resp")) %>%
    mutate(Name = gsub("Blood", "BL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = Time)
  for(i in 1:nrow(trait_anno$Response)){
    name <- trait_anno$Response[i,"Name"]
    trait <- trait_anno$Response[i,"Trait"]
    tissue <- trait_anno$Response[i,"Tissue"]
    drug <- trait_anno$Response[i,"Drug"]
    time <- trait_anno$Response[i,"Time"]
    samples <- clin_dat %>% 
      filter(Tissue == tissue, Drug == drug, Time == time) %>% 
      pull(Sample_id)
    clin_dat[[name]] <- clin_dat[[trait]]
    clin_dat[[name]][which(!clin_dat$Sample_id %in% samples)] <- NA
  }
  # Treatment traits
  trait_anno$Treatment <- expand.grid(
    Trait = NA,
    Type = "Treatment",
    Tissue = "Blood",
    Time = c("wk01", "wk04", "wk12"),
    Drug = c("Adalimumab", "Ustekinumab"),
    stringsAsFactors = F
  )
  trait_anno$Treatment <- trait_anno$Treatment %>%
    mutate(Name = paste0(Tissue,"_",Drug,"_",Time,"_vs_wk00")) %>%
    mutate(Name = gsub("Blood", "BL", Name)) %>%
    mutate(Name = gsub("Adalimumab", "ADA", Name)) %>%
    mutate(Name = gsub("Ustekinumab", "UST", Name)) %>%
    mutate(Short_name = paste0(Time," vs wk00"))
  for(i in 1:nrow(trait_anno$Treatment)){
    name <- trait_anno$Treatment[i,"Name"]
    tissue <- trait_anno$Treatment[i,"Tissue"]
    drug <- trait_anno$Treatment[i,"Drug"]
    time <- trait_anno$Treatment[i,"Time"]
    clin_dat[[name]] <- NA
    wk00_samples <- clin_dat %>% filter(Tissue == tissue, Drug == drug, Time == "wk00") %>% pull(Sample_id)
    treated_samples <- clin_dat %>% filter(Tissue == tissue, Drug == drug, Time == time) %>% pull(Sample_id)
    clin_dat[[name]][which(clin_dat$Sample_id %in% wk00_samples)] <- 0
    clin_dat[[name]][which(clin_dat$Sample_id %in% treated_samples)] <- 1
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

#' As part of this workflow we created an annotation table which gives details about each trait, including the
#' tissues, time points and drug cohorts that comprise them. Let's examine this table.

traits$trait_anno

#' We'll save this data to file for later use.

# Save
write.table(traits$traits, paste0(output_directory2,"/Traits.txt"),
            sep = "\t", quote = F, row.names = F)
write.table(traits$trait_anno, paste0(output_directory2,"/Traits_anno.txt"),
            sep = "\t", quote = F, row.names = F)

#' ## Trait correlations
#' 
#' Now we can correlate the traits defined above with the module eigengenes and latent factors.

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

#' # Session information

sessionInfo()