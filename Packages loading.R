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

install.packages("remotes")
remotes::install_github("cmartin/ggConvexHull")

packages = c("ggConvexHull","svMisc","dplyr","ggplot2","RColorBrewer","ggrepel","devtools","PCAtest",
             "tidyr","ggsci","tictoc","BAT","reshape2","ggpubr","pheatmap","hypervolume","alphahull",
             "ade4","ggExtra","corrplot","tidyverse","ggnewscale","ggpmisc","scales","cowplot")

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
  require(p, character.only = T)
}

rm(p,packages)

options(scipen = 999)
Sys.setlocale("LC_ALL", "English")
