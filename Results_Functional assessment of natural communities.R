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
require(ggnewscale)
require(ggpmisc)

## 1. Load the datasets
dat <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "CommTraits"),
  )
dat$Treatment <- NA
dat$Label <- NA

df <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "IndTraits"),
  )
df$sample <- NA

env <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "EnvData"),
  )
env$Date <- as.Date(env$Date, format = "%d-%m-%Y")
env <- env[,c(1:8,10,9)]


## 2. Data preparation - Fig. S9
# Log10 transformation
dat[,c(1:8)] <- log10(dat[,c(1:8)])

# Z-score normalisation
dat[,c(1:8)] <- scale(dat[,c(1:8)], center = TRUE, scale = TRUE)

# Log10 transformation
df[,c(3:10)] <- log10(df[,c(3:10)])

# Z-score normalisation
df[,c(3:10)] <- scale(df[,c(3:10)], center = TRUE, scale = TRUE)

# Merge biotic and abiotic data
dat <- merge(dat, env, by.x = "sample", by.y = "Sample", all.x = TRUE, all.y = TRUE)

# Convert character to numeric
dat[12:17] <- lapply(dat[12:17], function(x) {
  if(is.factor(x)) x <- as.character(x)
  
  x[x == "NA"] <- NA
  
  as.numeric(x)
})

dat <- dat[complete.cases(dat[, 2:9]), ]

dat <- plyr::rbind.fill(df, dat)

## 3. Alignement trait-space between experiment-field (individual-level)
# 3.1. Dimension reduction through a Principal Component Analysis
PCA <- dudi.pca(dat[,c(3:10)], scannf = FALSE, nf = 7)

# 3.2. Quick visualisation of eigenvalues
screeplot(PCA, main = "Screeplot - Eigenvalues")

# 3.3. Eigenvalues contribution
(PCA$eig*100)/sum(PCA$eig)

# 3.4. Quick visualisation of the PCA
s.corcircle(PCA$co)

# 3.5. Testing the number of axis to keep
PCAtest(PCA$tab, 999, 999, 0.001, varcorr=TRUE, counter=FALSE, plot=TRUE)

# 3.6. Contribution of each variable
rowSums(100*(factoextra::get_pca_var(PCA)$cos2)[,c(1,2)])

dat$Axis1 <- PCA$li$Axis1
dat$Axis2 <- PCA$li$Axis2

# Define the 5-point palette
desat_palette <- c(
  "#f7f7f7",  # very light gray
  "#cccccc",  # light gray
  "#969696",  # medium gray
  "#636363",  # dark gray
  "#252525"   # almost black
)

# Separation between field and experiment data
dat <- dat %>%
  mutate(Source = case_when(
    !is.na(sample) ~ "field",
    !is.na(Treatment) ~ "Experiment",
    TRUE ~ NA_character_
  ))

PCA$co$Label <- rownames(PCA$co)

# 3.7. Plot
p <- ggplot() +
  geom_point(data = subset(dat, Source == "field"), mapping = aes(x = Axis1, y = Axis2, color = `TN (mmol/L)`/I_m), size = 2) +
  scale_color_gradientn(colors = desat_palette, name = "TN:Light ratio") +
  new_scale_color() +
  geom_point(data = subset(dat, Source == "Experiment"), mapping = aes(x = Axis1, y = Axis2, color = Treatment), size = 2, shape = 18) +
  scale_color_d3(name = "Treatment") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
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
  theme_linedraw() +
  theme(legend.position = "left",
        legend.text=element_text(size=12),
        axis.text=element_text(size=18),
        axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Source)
p
ggsave(file="supp figure 5.svg", plot=p, width=26, height=18, units = "cm")


## 4. Data preparation - Fig. 5
# 4.1. Calculate barycenters (centroids) for each group
library(dplyr)

centroids <- dat %>%
  group_by(sample, Date) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

centroids_treat <- dat %>%
  group_by(Treatment) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

centroids_treat <- na.omit(centroids_treat)
centroids <- merge(centroids, env[,-3], by.x = "sample", by.y = "Sample", all.x = TRUE)

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


## Supp. Material S7
require(tidyr)

common_theme <- theme(
  axis.text  = element_text(size = 10),
  axis.title = element_text(size = 10)
) +
  theme_bw()

p1 <- 
  centroids %>%
  drop_na(`TN (mmol/L)`) %>%
  ggplot(mapping = aes(x = Date, y = `TN (mmol/L)`)) +
  geom_point() +
  geom_point(data = N_lim, mapping = aes(x = Date, `TN (mmol/L)`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_line() +
  xlab("Date") +
  ylab(expression("Total Nitrogen (mmol L"^{-1}*")" )) +
  common_theme

p2 <- 
  centroids %>%
  drop_na(`DIN (mmol/L)`) %>%
  ggplot(mapping = aes(x = Date, y = `DIN (mmol/L)`)) +
  geom_point() +
  geom_point(data = N_lim, mapping = aes(x = Date, y = `DIN (mmol/L)`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_line() +
  xlab("Date") +
  ylab(expression("Dissolved Inorganic Nitrogen (mmol L"^{-1}*")" )) +
  common_theme

p3 <- 
  centroids %>%
  drop_na(`TP (mmol/L)`) %>%
  ggplot(mapping = aes(x = Date, y = `TP (mmol/L)`)) +
  geom_point() +
  geom_point(data = N_lim, mapping = aes(x = Date, `TP (mmol/L)`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_line() +
  xlab("Date") +
  ylab(expression("Total Phosphorus (mmol L"^{-1}*")" )) +
  common_theme

p4 <- 
  centroids %>%
  drop_na(I_m) %>%
  ggplot(mapping = aes(x = Date, y = I_m)) +
  geom_point() +
  geom_point(data = N_lim, mapping = aes(x = Date, y = I_m), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_line() +
  xlab("Date") +
  ylab(expression("Irradiance of the mixed layer (" ~ mu * "mol photons" ~ m^2 ~ s^-1*")")) +
  common_theme

p5 <- 
  centroids %>%
  ggplot(mapping = aes(y = `TN (mmol/L)`/`TP (mmol/L)`, x = Date))+
  geom_point(data = N_lim, aes(y = `TN (mmol/L)`/`TP (mmol/L)`, x = Date) , col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_point() +
  ylab("TN:TP") +
  xlab("Date") +
  stat_smooth(se = TRUE, method = "lm", col = "grey20", linetype = "dashed") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4, label.x = 0.95) +
  common_theme

p6 <- 
  centroids %>%
  ggplot(mapping = aes(y = `Global irradiance`, x = Date))+
  geom_point(data = N_lim, aes(y = `Global irradiance`, x = Date) , col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  geom_point() +
  ylab(expression("Global irradiance (" ~ mu * "mol photons" ~ m^2 ~ s^-1*")")) + 
  xlab("Date") +
  stat_smooth(se = FALSE, method = "gam", col = "grey20", linetype = "dashed") +
  common_theme

# 7.9) Produce the Time series plot
top <- plot_grid(p1, p2, labels = c('A','B'))
mid <- plot_grid(p3, p4, labels = c('C','D'))
bot <- plot_grid(p5, p6, labels = c('E','F'))

p <- plot_grid(top, mid, bot, ncol =1)

ggsave(file="Supplementary mat TS env param.svg", plot=p, width=18, height=24, units = "cm")


## Supp. Material S8
traits_table <- reshape2::melt(
  dat[,c(1,3:10)],
  id.vars = "Treatment",
  variable.name = "Trait",
  value.name = "value"
)

traits_table <- na.omit(traits_table)

p <- traits_table %>%
  ggplot() +
  geom_density(aes(x = value, color = Treatment)) +
  facet_wrap(.~Trait, nrow = 2, scales = "free") +
  scale_color_d3(name = traits_table$Treatment) +
  xlab("Trait value") +
  ylab("Density distribution") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

ggsave(file="Supplementary mat density distribution.svg", plot=p, width=16, height=10, units = "cm")
