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

# Note: ggConvexHull must be installed from GitHub manually:
# remotes::install_github("cmartin/ggConvexHull")

packages = c("ggConvexHull","svMisc","dplyr","dataverse","ggplot2","RColorBrewer","ggrepel","devtools","PCAtest",
             "tidyr","ggsci","tictoc","BAT","reshape2","ggpubr","pheatmap","hypervolume","alphahull",
             "ade4","ggExtra","corrplot","tidyverse","ggnewscale","ggpmisc","scales","cowplot","foreach","doSNOW",
            "progress","plyr","conover.test","stringr","multcompView","Matrix","ggbreak", "readxl")

missing_pkgs <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "The following packages are required but not installed:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

invisible(lapply(packages, library, character.only = TRUE))

rm(packages)

options(scipen = 999)
try(Sys.setlocale("LC_ALL", "English"), silent = TRUE)
