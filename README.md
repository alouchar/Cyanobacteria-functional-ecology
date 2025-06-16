# Cyanobacteria-functional-ecology
The code aims at reproducing the analyses and figures presented in the paper titled Functional response of cyanobacteria to environmental changes. These code use the dataset located at: DOI.
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
