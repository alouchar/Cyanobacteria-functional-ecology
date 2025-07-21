# Cyanobacteria-functional-ecology
The code aims at reproducing the analyses and figures presented in the paper titled Functional response of cyanobacteria to environmental changes. These code use the dataset located at: doi.org/10.34894/9X9YMO.
The code uses R version 4.3.1.

## Content
  1. License
  2. README
  3. R scripts:
  - Packages loading.R
  - Material and Methods.R
  - Results_Flow cytometry traits.R
  - Results_Functional analyses from individual level.R
  - Results_Functional assessment of natural communities.R

## Reproduction of analyses and figures 

#### Packages loading

`Packages loading.R`

#### Material and Methods

`Material and Methods.R`

The script details methodological steps necessar to prepare the data. The script produce figure 1 panels B to F.

Dataset: IndTraits

#### Results

`Results_Flow cytometry traits.R`

The script produces the figure 2.

Dataset: IndTraits

`Results_Functional analyses from individual level.R`

The script produces the figures 3 and 4. From Line 175 the code runs only for one treatment. 

Dataset: IndTraits

`Results_Functional assessment of natural communities.R`

The script produces the figure 5. Ensure that you read the data set "DATA" to produce the functional fingerprints.

Datasets: IndTraits, CommTraits

## List of packages and versions

<pre lang="markdown"> ## Required R Packages and Versions ```r # ==== Required Packages and Versions ==== # ggConvexHull 0.1.0 # svMisc 1.2.3 # dplyr 1.1.3 # ggplot2 3.5.1 # RColorBrewer 1.1-3 # ggrepel 0.9.4 # devtools 2.4.5 # PCAtest 0.0.1 # tidyr 1.3.0 # ggsci 3.0.0 # tictoc 1.2 # BAT 2.9.6 # reshape2 1.4.4 # ggpubr 0.6.0 # pheatmap 1.0.12 # hypervolume 3.1.3 # alphahull 2.5 # ade4 1.7-22 # ggExtra 0.10.1 # corrplot 0.92 # tidyverse 2.0.0 # ggnewscale 0.4.10 # ggpmisc 0.5.4-1 # scales 1.3.0 # cowplot 1.1.3 ``` </pre>
