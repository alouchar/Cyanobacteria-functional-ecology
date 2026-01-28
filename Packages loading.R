#`---
#`Authors: Arnaud P. Louchart, Annemieke M. Drost, Chaohong Lin, Suzanne M.H. Wiezer, Zhipeng Duan, Elena Litchman, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 2025
#`Journal: Ecology Letters
#`Status: First submission
#`---

##############################################################################
############################## PACKAGES LOADING ##############################
##############################################################################

# Note: ggConvexHull and PCAtest must be installed from GitHub manually:
# install.packages("remotes") 
# install.packages("devtools")
# remotes::install_github("cmartin/ggConvexHull")
# devtools::install_github("arleyc/PCAtest")

install.packages("remotes")


packages = c("ggConvexHull","svMisc","dplyr","dataverse","ggplot2","RColorBrewer","ggrepel","devtools","PCAtest",
             "tidyr","ggsci","tictoc","BAT","reshape2","ggpubr","pheatmap","hypervolume","alphahull",
             "ade4","ggExtra","corrplot","tidyverse","ggnewscale","ggpmisc","scales","cowplot","foreach","doSNOW",
             "progress","plyr","conover.test","stringr","multcompView","Matrix","ggbreak", "readxl")

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
  require(p, character.only = T)
}

rm(p,packages)

options(scipen = 999)
Sys.setlocale("LC_ALL", "English")
