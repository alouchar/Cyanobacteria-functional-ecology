#`---
#`Authors: Arnaud P. Louchart, Annemieke M. Drost, Chaohong Lin, Suzanne M.H. Wiezer, Zhipeng Duan, Elena Litchman, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 2025
#`Journal: Ecology Letters
#`Status: First submission
#`---

#############################################################################
############################## RESULTS SECTION ##############################
#############################################################################

## 1. Load the datasets
df <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "CommTraits"),
  )
df$Treatment <- NA
df$Replicate <- NA

dat <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "IndTraits"),
  )
dat$Sample <- NA
dat <- dat[,c(1:10)]

env <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "EnvData"),
  )
env$Date <- as.Date(env$Date, format = "%d-%m-%Y")
env <- env[,c(1:8,10,9)]


## 2. Data preparation
# Log10 transformation
df[,c(2:9)] <- log10(df[,c(2:9)])

# Z-score normalisation
df[,c(2:9)] <- scale(df[,c(2:9)], center = TRUE, scale = TRUE)

# Log10 transformation
dat[,c(3:10)] <- log10(dat[,c(3:10)])

# Z-score normalisation
dat[,c(3:10)] <- scale(dat[,c(3:10)], center = TRUE, scale = TRUE)

# Merge biotic and abiotic field data
df <- merge(df, env, by = "Sample", all.x = TRUE, all.y = TRUE)

# Convert character to numeric
df[14:19] <- lapply(df[14:19], function(x) {
  if(is.factor(x)) x <- as.character(x)
  
  x[x == "NA"] <- NA
  
  as.numeric(x)
})

df <- df[complete.cases(df[, 2:9]), ]

df <- rbind.fill(dat, df)

## 3. Alignement trait-space between experiment-field (individual-level: Fig. S9)
# 3.1. Dimension reduction through a Principal Component Analysis
PCA <- dudi.pca(df[,c(3:10)], scannf = FALSE, nf = 7)

# 3.2. Quick visualisation of eigenvalues
screeplot(PCA, main = "Screeplot - Eigenvalues")

# 3.3. Eigenvalues contribution
(PCA$eig*100)/sum(PCA$eig)

# 3.4. Testing the number of axis to keep
PCAtest(PCA$tab, 999, 999, 0.001, varcorr=TRUE, counter=FALSE, plot=TRUE)

# 3.5. Contribution of each variable
rowSums(100*(get_pca_var(PCA)$cos2)[,c(1,2)])

df$Axis1 <- -PCA$li$Axis1
df$Axis2 <- PCA$li$Axis2

# Define the 5-point palette
desat_palette <- c(
  "#f7f7f7",  # very light gray
  "#cccccc",  # light gray
  "#969696",  # medium gray
  "#636363",  # dark gray
  "#252525"   # almost black
)

# Separation between field and experiment data
df <- df %>%
  mutate(Source = case_when(
    !is.na(Sample) ~ "field",
    !is.na(Treatment) ~ "Experiment",
    TRUE ~ NA_character_
  ))

PCA$co$Label <- rownames(PCA$co)


## 4. Data preparation - Fig. 5
# 4.1. Calculate barycenters (centroids) for each group
centroids <- df %>%
  group_by(Sample, Date) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

centroids_treat <- 
  df %>%
  group_by(Treatment) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

centroids_treat <- na.omit(centroids_treat)
centroids <- merge(centroids, env[,-3], by = "Sample", all.x = TRUE)

# Convert character to numeric
centroids[6:11] <- lapply(centroids[6:11], function(x) {
  if(is.factor(x)) x <- as.character(x)
  
  x[x == "NA"] <- NA
  
  as.numeric(x)
})

# Identification Nitrogen limitation in field data
N_lim <-
  centroids %>%
  filter(Date == as.Date("2023-06-07") |
           Date == as.Date("2023-06-21") |
           Date == as.Date("2023-07-26")|
           Date == as.Date("2023-08-09")|
           Date == as.Date("2023-08-30")|
           Date == as.Date("2023-09-21"))

# Plot
p <- ggplot() +
  geom_point(data = centroids, mapping = aes(x = Axis1, y = Axis2, color = `TN (mmol/L)`/I_m), size = 4) +
  scale_color_gradientn(colors = desat_palette, name = "TN:Light ratio") +
  new_scale_color() +
  geom_point(data = centroids_treat, mapping = aes(x = Axis1, y = Axis2, color = Treatment), size = 6, shape = 18) +
  scale_color_d3(name = "Treatment") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(data = N_lim, mapping = aes(x = Axis1, y = Axis2), col = "orange", size = 8, shape = 21,  stroke = 1.2) +
  geom_segment(
    data = PCA$co,
    aes(x = 0, y = 0, xend = -Comp1*3, yend = Comp2*3),
    col = "grey20", arrow = arrow(length = unit(0.02, "npc")),
    inherit.aes = FALSE
  ) +
  geom_label_repel(
    data = PCA$co,
    aes(x = -Comp1*4, y = Comp2*4, label = Label),
    col = "grey20", alpha = 0.7, size = 4, inherit.aes = FALSE
  ) +
  xlab("PC1 (49.0%)") +
  ylab("PC2 (28.7%)") +
  xlim(c(-4,4)) +
  ylim(c(-4,4)) +
  theme_linedraw() +
  theme(legend.position = "left",
        legend.text=element_text(size=12),
        axis.text=element_text(size=18),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p
ggsave(file="figure 5.svg", plot=p, width=20, height=14, units = "cm")
