#`---
#`Authors: Arnaud P. Louchart, Annemieke M. Drost, Chaohong Lin, Suzanne M.H. Wiezer, Zhipeng Duan, Elena Litchman, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 2025
#`Journal: Ecology Letters
#`Status: First submission
#`---

##########################################################################################
############################## MATERIAL AND METHODS SECTION ##############################
##########################################################################################

#### 1. DATA LOADING
# IndTraits contains flow cytometry informations at the individual level for each treatment obtained through the PhytoCytoTraits GUI.
# Treatment contains informations on the different environmental conditions

dat <- readxl::read_excel("D:/Manuscript/Cyanobacterial traits - Louchart - 2025/Manuscript/IndCyano_Louchart_Limitation_experiement_June2025.xlsx", sheet = "IndTraits")

#### 2 DATA MANIPULATION
## 2.1. Log10 transformation
dat[,c(3:10)] <- log10(dat[,c(3:10)])

## 2.2. Z-score normalisation
dat[,c(3:10)] <- scale(dat[,c(3:10)], center = TRUE, scale = TRUE)
