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

This repository contains several R scripts and associated workbooks which cover the core analyses in this paper. These are oraganised into several directories and below we provide details about which scripts and workbooks correspond to which figures in the paper.

### Figure 1

| Figure | Path to workbook                                                                    | Details                                                                      |
|------------------|---------------------------|---------------------------|
| 1A, 1B        | [link](/WGCNA/03_Get_disease_and_disease_severity_correlations.md)               | Calculates correlations between modules/factors and traits                                |
