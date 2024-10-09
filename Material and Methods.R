#`---
#`Authors: Arnaud Louchart, Annemieke Drost, Harley Lin, Suzanne Wiezer, Zhipeng Duan, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 202X
#`Journal: Functional Ecology
#`Status: Writing
#`---

##########################################################################################
############################## MATERIAL AND METHODS SECTION ##############################
##########################################################################################

#### 1. DATA IMPORTATION

# dat contains flow cytometry informations at the individual level for each treatment obtained through the PhytoCytoTraits GUI.
# Treatment contains informations on the different environmental conditions

dat <- readxl::read_excel("C:/Users/ArnaudL/Desktop/Zhipeng data/FCM data/Combined_data_Zhipeng.xlsx")
treatment <- readxl::read_excel("C:/Users/ArnaudL/Desktop/Zhipeng data/Treatment.xlsx", sheet = "Sheet2")

#### 2. CLEANING
## 2.1. Identify and rename cyanobacterial cluster. Quick and dirty step to improve.
{
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_001.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_002.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_003.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_004.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_007.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_008.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_014.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_015.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_017.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_019.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_020.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_021.fcs" & dat$db_FCM == 2] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_022.fcs" & dat$db_FCM == 2] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_023.fcs" & dat$db_FCM == 3] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_024.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_025.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_024.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_025.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_026.fcs" & dat$db_FCM == 5] <- "Cyanobacteria"
  
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_005.fcs" & dat$db_FCM == 2] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_006.fcs" & dat$db_FCM == 2] <- "Cyanobacteria"
  
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_009.fcs" & dat$db_FCM == 3] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_010.fcs" & dat$db_FCM == 3] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_013.fcs" & dat$db_FCM == 3] <- "Cyanobacteria"
  
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_011.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_012.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_016.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  dat$db_FCM[dat$filename == "Zhipeng DW cyanoMQ_018.fcs" & dat$db_FCM == 4] <- "Cyanobacteria"
  
}

## 2.2. Subset of data: working only with cyanobacteria cluster
dat <- 
  dat %>%
  filter(db_FCM == "Cyanobacteria")

# merge treatment and dat by filename
dat <- merge(dat, treatment, by.x = "filename",by.y = "Filename", all.x = TRUE)

## 2.3. Data are already log10 transformed. Back transformation to original values
dat[,c(3:24)] <- 10^(dat[,c(3:24)]-1)

#### 3. CREATION OF NEW TRAITS (see table 1)
dat$Fill_SSC <- dat$SSC/(4/3*(pi*(dat$`FSC [Par]`/2)^3))
dat$PC_size <- dat$`670/30[640]`/dat$`FSC [Par]`
dat$Chla_size <- dat$`692/40[488]`/dat$`FSC [Par]`
dat$PC_Chla <- dat$`670/30[640]`/dat$`692/40[488]`

#### 4. DATA CLEANING
## 4.1. Homogeneisation of labels

dat <-
  dat %>%
  dplyr::filter(Label == "Control" | Label == "Nitrogen" | Label == "Phosphorus" | Label == "CO2" | Label == "Light")

dat$Label <- gsub("CO2", "+CO2", dat$Label)
dat$Label <- stringr::str_replace_all(dat$Label, "Light", "-Light")
dat$Label <- gsub("Nitrogen", "-Nitrogen", dat$Label)
dat$Label <- gsub("Phosphorus", "-Phosphorus", dat$Label)

## 4.2. Removing the first nitrogen replicate (as being flagged)
dat <-
  dat %>%
  filter(!filename == "Zhipeng DW cyanoMQ_005.fcs")

## 4.3. Create a new column to determine replicate
dat["filename"][dat["filename"] == "Zhipeng DW cyanoMQ_001.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_006.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_010.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_014.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_018.fcs" ] <- "Replicate 1"

dat["filename"][dat["filename"] == "Zhipeng DW cyanoMQ_002.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_007.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_011.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_015.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_019.fcs" ] <- "Replicate 2"

dat["filename"][dat["filename"] == "Zhipeng DW cyanoMQ_003.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_008.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_012.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_016.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_020.fcs" ] <- "Replicate 3"

dat["filename"][dat["filename"] == "Zhipeng DW cyanoMQ_004.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_009.fcs" | 
                  dat["filename"] == "Zhipeng DW cyanoMQ_013.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_017.fcs" |
                  dat["filename"] == "Zhipeng DW cyanoMQ_021.fcs" ] <- "Replicate 4"

dat <-
  dat %>%
  mutate(Replicat = paste(Label,filename))

#### 5. DATA MANIPULATION
## 5.1. Subset 
DATA <- 
  dat[,c(3,5,9,17,20:22,27,26,1)]

## 5.2. Log10 transformation
DATA[,c(1:8)] <- log10(DATA[,c(1:8)])

## 5.3. Z-score normalisation
DATA[,c(1:8)] <- scale(DATA[,c(1:8)], center = TRUE, scale = TRUE)
