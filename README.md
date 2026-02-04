# Cyanobacteria-functional-ecology
The code aims at reproducing the analyses and figures presented in the paper titled "Individual-level trait responses in cyanobacterial populations and communities". This code directly opens the dataset located at: doi.org/10.34894/9X9YMO.
The code uses R version 4.3.1.

## Content
  1. License
  2. README
  3. R scripts:
  - Material and Methods.R
  - Results_Flow cytometry traits.R
  - Results_Functional analyses from individual level.R
  - Results_Functional assessment of natural communities.R
  - Supplementary material.R

## Reproduction of analyses and figures 


#### Material and Methods

`Material and Methods.R`

The script details methodological steps necessar to prepare the data. The script produce figure 1 panels B to F.

Tab: IndTraits

#### Results

`Results_Flow cytometry traits.R`

The script produces the figure 2.

Tab: IndTraits

`Results_Functional analyses from individual level.R`

The script produces the figures 3 and 4. Initially, the code code ran only for one treatment starting from Line 175. An update has generalised the code for all treatments. 

Tab: IndTraits

`Results_Functional assessment of natural communities.R`

The script produces the figure 5.

Tabs: IndTraits, CommTraits, envData

#### Supplementary informations

`Supplementary material.R`

The script produces the figures S1 to S9.

Tabs: IndTraits, CommTraits, envData, Supp. S1, Supp. S6


## List of packages and versions

```r
# ==== Required Packages and Versions ====

# dplyr             1.1.3
# ggplot2           3.5.1
# ggrepel           0.9.4
# devtools          2.4.5
# dataverse         0.3.16
# PCAtest           0.0.1
# ggsci             3.0.0
# BAT               2.9.6
# reshape2          1.4.4
# pheatmap          1.0.12
# hypervolume       3.1.3
# ade4              1.7-22
# ggExtra           0.10.1
# tidyverse         2.0.0
# ggnewscale        0.4.10
# ggpmisc           0.5.4-1
# cowplot           1.1.3
# foreach           1.5.2
# doSNOW            1.0.20
# stringr           1.6.0
# factoextra        1.0.7
# rstatix           0.7.3
# ggpattern         1.2.1
# parallel          4.5.2
# plyr              1.8.9
# conover.test      1.1.5
# multcompView      0.1-9
# Matrix            1.6-1.1
# readxl            1.4.3
```
