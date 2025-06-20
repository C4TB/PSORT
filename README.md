# Code repository for: *Transcriptomic profiling and machine learning uncover gene signatures of psoriasis endotypes and disease severity*
Rider, Grantham HJ, Smith GR, Watson DS et al. on behalf of the [PSORT consortium](http://www.psort.org.uk/)
## ArrayExpress data link
The raw and adjusted gene count data from our RNA-seq analysis, along with associated clinical data are available at Array Express under accession number [E-MTAB-14509](https://www.ebi.ac.uk/biostudies/arrayexpress/studies?query=E-MTAB-14509).
## Shiny portal link
The RNA-Seq data may be visualised and further explored through a [R Shiny web interface](https://shiny-whri-c4tb.hpc.qmul.ac.uk/psort/).  
Username: psort  
Password: tower squid ramp
##  Graphical summary of study design 
<img src='./Images/PSORT schematic 11-06-25.png'>  

## Repository contents
In order to reproduce the analyses in this paper you will need to download the RNA-Seq data and clinical data from ArrayExpress (see link above) You will also need the supplementary data from the paper, including: the ICA factor loadings, metagenes and most aligned features; the WGCNA module genes and eigengenes. Some of this data is also available in this repository in the following directory: Cell_Type_Correlations/paper_data

This repository contains several R scripts and workbooks which cover the core analyses in this paper. These are oraganised into several directories and below we provide details about which scripts and workbooks correspond to which figures in the paper.

| Figure | Link to script or workbook | Details |
|------------------|---------------------------|---------------------------|
| 1A, 1B | [link](/WGCNA/03_Get_disease_and_disease_severity_correlations.md) | Calculates correlations between modules/factors and traits |
| 1A, 1B | [link](/WGCNA/04_Metascape_inputs.md) | Creates imput files for Metascape; functional annotations from Metascape were used to create descriptive module and factor names |
| 1A, 1B | [link](/WGCNA/05_Module_factor_names.R) | Creates descriptive module and factor names which are used as annotations in the module/factor-trait correlation heatmaps |
| 1A, 1B | [link](/WGCNA/06_Trait_correlation_heatmaps.md) | Creates module/factor-trait correlation heatmaps |
| 1C | [link](/WGCNA/07_Trait_correlation_plots.md) | Plots exemplar module/factor-trait correlations |
| 2A, 2B | [link](/Factor_analysis/factor_analysis.md) | Carries out BMI differential expression analysis and creates BMI and PASI association heatmap |
| 2C | [link](/WGCNA/08_Gene-level_trait_correlation_plots.R) | Creates exemplar gene-level BMI and PASI association plots |
