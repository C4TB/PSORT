Identification of co-expressed gene modules in skin and blood
================
Ashley Rider
2024-10-01

- <a href="#preliminaries" id="toc-preliminaries">Preliminaries</a>
- <a href="#skin" id="toc-skin">Skin</a>
  - <a href="#load-data" id="toc-load-data">Load data</a>
    - <a href="#clinical-data" id="toc-clinical-data">Clinical data</a>
    - <a href="#gene-level-counts" id="toc-gene-level-counts">Gene-level
      counts</a>
  - <a href="#filter-and-normalise-counts"
    id="toc-filter-and-normalise-counts">Filter and normalise counts</a>

Here, we use WGCNA, or **W**eighted **G**ene **C**o-expression
**N**etwork **A**nalysis ([Langfelder & Horvath,
2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)),
to identify co-expressed gene modules in skin and blood.

# Preliminaries

``` r
# Load packages
library(tidyverse)
library(edgeR)
library(limma)
library(WGCNA)

# Load gene annotation data
anno <- read.delim("data/gene_annotation_data/Hs.anno.csv", sep = ",") %>%
  # Drop Description column
  select(EnsemblID, GeneSymbol) %>%
  # Replace gene symbols that are "" with corresponding Ensembl ID
  mutate(GeneSymbol = if_else(GeneSymbol == "", EnsemblID, GeneSymbol))

head(anno)
```

    ##         EnsemblID      GeneSymbol
    ## 1 ENSG00000252303       RNU6-280P
    ## 2 ENSG00000281771 ENSG00000281771
    ## 3 ENSG00000281256 ENSG00000281256
    ## 4 ENSG00000280864 ENSG00000280864
    ## 5 ENSG00000280792 ENSG00000280792
    ## 6 ENSG00000281822        RNU1-62P

# Skin

We’ll start by analysing the skin data.

## Load data

### Clinical data

We load the clinical data for each cohort.

``` r
# Load clinical data for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
clin <- list()
clin$Discovery_skin <- read.delim("data/clinical_data/PSORT-D_Skin_Clinical_Data_01-Apr-20.txt") %>%
  select(Patient_id, Sample_id, Drug, Tissue, Time)
clin$Replication_skin <- read.delim("data/clinical_data/PSORT-R_Skin_Clinical_Data_01-Apr-20.txt") %>%
  select(Patient_id, Sample_id, Drug, Tissue, Time)

dim(clin$Discovery_skin)
```

    ## [1] 400   5

``` r
head(clin$Discovery_skin)
```

    ##   Patient_id Sample_id       Drug      Tissue Time
    ## 1     P.1001 1001-0006 Adalimumab Nonlesional wk00
    ## 2     P.1001 1001-0105 Adalimumab    Lesional wk01
    ## 3     P.1001 1001-1206 Adalimumab Nonlesional wk12
    ## 4     P.1001 1001-1205 Adalimumab    Lesional wk12
    ## 5     P.1001 1001-0005 Adalimumab    Lesional wk00
    ## 6     P.1002 1002-0005 Adalimumab    Lesional wk00

``` r
dim(clin$Replication_skin)
```

    ## [1] 280   5

``` r
head(clin$Replication_skin)
```

    ##   Patient_id Sample_id        Drug      Tissue Time
    ## 1     P.1015 1015-0006  Adalimumab Nonlesional wk00
    ## 2     P.1015 1015-1205  Adalimumab    Lesional wk12
    ## 3     P.1015 1015-0105  Adalimumab    Lesional wk01
    ## 4     P.1015 1015-1206  Adalimumab Nonlesional wk12
    ## 5     P.1015 1015-0005  Adalimumab    Lesional wk00
    ## 6     P.1016 1016-0006 Ustekinumab Nonlesional wk00

Preliminary analysis revealed mislabelling of some lesional and
non-lesional samples, which were swapped for several patients. We
correct this here.

``` r
# Sample swaps based on consensus of S100A7 and S100A9 expression
clin$Discovery_skin <- clin$Discovery_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0005", "23012-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("6041-0006", "23012-1206"), "Lesional", Tissue))
clin$Replication_skin <- clin$Replication_skin %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0005", "5031-1205", "6049-1205"), "Nonlesional", Tissue)) %>%
  mutate(Tissue = if_else(Sample_id %in% c("5033-0006", "5031-1206", "6049-1206"), "Lesional", Tissue))
```

### Gene-level counts

We also load gene-level counts for each cohort.

``` r
# Load counts for PSORT-D (Discovery) and PSORT-R (Replication) and add to list
cnts <- list()
cnts$Discovery_skin <- readRDS("data/gene_level_counts/PSORT-D_Skin_counts_01-Apr-2020-13-00-07.rds")
cnts$Replication_skin <- readRDS("data/gene_level_counts/PSORT-R_Skin_counts_13-Mar-2020-15-35-24.rds")

dim(cnts$Discovery_skin)
```

    ## [1] 34947   400

``` r
cnts$Discovery_skin[1:5, 1:5]
```

    ##                  1001-0005 1002-0005 1003-0005 1005-0005 1013-0005
    ## ENSG00000000003 1175.44411 1562.9835 1252.3025 1361.3136  821.4603
    ## ENSG00000000005   24.27478  291.7607   14.9610  352.9931  198.2883
    ## ENSG00000000419 2669.42467 2493.0064 2185.6786 2649.1577 1656.7727
    ## ENSG00000000457  720.88201  950.2264 1023.5560  632.5259  777.8005
    ## ENSG00000000460  340.11979  375.3907  358.9917  335.8198  179.4848

``` r
dim(cnts$Replication_skin)
```

    ## [1] 39297   280

``` r
cnts$Replication_skin[1:5, 1:5]
```

    ##                  1015-0005 1022-1206 5045-0006 5045-0105 5045-1205
    ## ENSG00000000003 1324.26564 2870.1355 1547.0102 1260.4806 1137.9581
    ## ENSG00000000005   37.50147  146.5253  161.0882  231.0615  412.8736
    ## ENSG00000000419 3216.75517 2261.0516 1820.7281 1678.0153 1980.3527
    ## ENSG00000000457 1258.47998 1952.3299 1420.5515 1421.3044 1495.6756
    ## ENSG00000000460  477.45813  547.8123  403.2949  399.3427  473.1268

The PSORT-R counts contain some extra genes that are not present in the
PSORT-D counts. We’ll proceed with just the genes that are present in
both datasets.

``` r
intersect_genes <- intersect(rownames(cnts$Discovery_skin), rownames(cnts$Replication_skin))

cnts$Discovery_skin <- cnts$Discovery_skin[intersect_genes,]

cnts$Replication_skin <- cnts$Replication_skin[intersect_genes,]

dim(cnts$Discovery_skin)
```

    ## [1] 34947   400

``` r
dim(cnts$Replication_skin)
```

    ## [1] 34947   280

## Filter and normalise counts

We now filter and normalise the counts. We use a filtering threshold
that requires at least 1 CPM in at least *n/k* samples, where *n* equals
the number of samples and *k* equals the number of unique combinations
of drug, tissue type (i.e. lesional and non-lesional) and time point. We
then normalise the counts using the trimmed mean of m-values (TMM)
method from [Robinson & Oshlack,
2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25).
We wrap this filtering and normalisation workflow in a function and
apply it to the PSORT-D and PSORT-R data separately.

``` r
normCounts <- function(cnts_dat, clin_dat){
  # Add Drug-Tissue-Time interaction variable to clinical data
  clin_dat <- clin_dat %>%
    mutate(Drug.Tissue.Time = paste(Drug, Tissue, Time, sep = "."))
  # Filter counts
  gene_ids <- tibble(EnsemblID = rownames(cnts_dat))
  cnts_dat <- cnts_dat[gene_ids$EnsemblID, clin_dat$Sample_id]
  y <- DGEList(cnts_dat, genes = gene_ids)
  keep <- rowSums(cpm(y) >= 1) >= nrow(clin_dat) / length(unique(clin_dat$Drug.Tissue.Time))
  y <- DGEList(y[keep,])
  # TMM normalisation
  y <- calcNormFactors(y)
  return(y)
}

cnts$Discovery_skin <- normCounts(cnts_dat = cnts$Discovery_skin, clin_dat = clin$Discovery_skin)

cnts$Replication_skin <- normCounts(cnts_dat = cnts$Replication_skin, clin_dat = clin$Replication_skin)

dim(cnts$Discovery_skin)
```

    ## [1] 15797   400

``` r
dim(cnts$Replication_skin)
```

    ## [1] 15731   280
