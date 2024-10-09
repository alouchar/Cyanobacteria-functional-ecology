#`---
#`Authors: Arnaud Louchart, Annemieke Drost, Harley Lin, Suzanne Wiezer, Zhipeng Duan, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 202X
#`Journal: Functional Ecology
#`Status: Writing
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