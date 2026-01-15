#`---
#`Authors: Arnaud Louchart, Annemieke Drost, Harley Lin, Suzanne Wiezer, Zhipeng Duan, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 202X
#`Journal: Functional Ecology
#`Status: Writing
#`---

#####################################################################################
############################## SUPPLEMENTARY MATERIALS ##############################
#####################################################################################
rm(list = ls(all=TRUE))

require(dplyr)
require(ggplot2)
require(tidyr)
require(cowplot)
require(patchwork)
require(rstatix)
require(ggpattern)
library(tidyr)
library(multcompView)
library(ggrepel)


## Supp. Material fig. S1
dat <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "Supp. S1"),
  )

p1 <-
  dat %>%
  ggplot(aes(color = as.factor(Cluster))) +
  geom_point(mapping = aes(x = FSC, y = `670/30 (640nm)`), size = 0.7) +
  facet_grid(.~Method) +
  xlab("") +
  scale_color_manual(values = c("#56B4E9","#CC79A7"), name = "Cluster") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))

legend <- get_plot_component(p1, 'guide-box-bottom', return_all = TRUE)

p2 <-
  dat %>%
  ggplot(aes(color = as.factor(Cluster))) +
  geom_point(mapping = aes(x = FSC, y = `692/40 (488nm)`), size = 0.7) +
  facet_grid(.~Method) +
  scale_color_manual(values = c("#56B4E9","#CC79A7"), name = "Cluster") +
  theme_bw() 


p <- plot_grid(p1 + theme(legend.position = "none"), 
          p2 + theme(legend.position = "none"),
          legend,
          nrow = 3, rel_heights = c(1,1,0.2))

ggsave(file="Algorithm_selection.svg", plot=p, width= 17.3, height=12, units = "cm")


## Supp. Material fig. S2
dat <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "IndTraits"),
  )

p1 <- dat %>%
  filter(Treatment == "Control") %>%
  filter(Replicate == "Replicate 1") %>%
  ggplot(mapping = aes(x = `Cell size`)) +
  xlab("Cell size (FSC)") +
  ylab("Number of events") +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(tag = "Distribution of cell size (original values)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.5, 1.05),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10))

p2 <- dat %>%
  filter(Treatment == "Control") %>%
  filter(Replicate == "Replicate 1") %>%
  ggplot(mapping = aes(x = log10(`Cell size`+1))) +
  xlab("log10 Cell size (FSC)") +
  ylab("Number of events") +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(tag = "Distribution of cell size (log10 values)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.5, 1.05),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10))

p3 <- dat %>%
  filter(Treatment == "-Nitrogen") %>%
  filter(Replicate == "Replicate 3") %>%
  ggplot(mapping = aes(x = Phycocyanin)) +
  xlab("Phycocyanin (670/30 [640 nm])") +
  ylab("Number of events") +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(tag = "Distribution of phycocyanin (original values)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.5, 1.05),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10))

p4 <- dat %>%
  filter(Treatment == "-Nitrogen") %>%
  filter(Replicate == "Replicate 3") %>%
  ggplot(mapping = aes(x = log10(Phycocyanin+1))) +
  xlab("log10 Phycocyanin (670/30 [640 nm])") +
  ylab("Number of events") +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(tag = "Distribution of phycocyanin (log10 values)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.5, 1.05),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10))


bot <- plot_grid(p2, p4)
top <- plot_grid(p1, p3)
p <- plot_grid(top, bot, ncol =1)
p
ggsave(file="Supplementary mat distri.svg", plot=p, width=36, height=36, units = "cm")


## Supp. Material fig. S3
p <- dat[,-c(1,2,11)] %>% cor_mat()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

get_upper_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

p <- as.matrix(p[,-1])
row.names(p) <- colnames(p)

p <- reorder_cormat(p)
p <- get_upper_tri(p)

melt_mat <- reshape2::melt(p, na.rm = TRUE)
melt_mat$value <- round(melt_mat$value,2)


img <- ggplot(data = melt_mat, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(high = "#B2182B", mid = "grey90", low = "#053061",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson \nCorrelation") +
  geom_tile_pattern(
    data = subset(melt_mat, value > 0.8), 
    aes(pattern = "point"),
    fill = NA,
    colour = "black", 
    pattern_density = 0.01,
    show.legend = FALSE) +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 3) +
  theme_bw() +
  xlab("") +
  ylab("") +
  coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(size = 12, family = "Arial", hjust = 1, angle = 35),
        axis.text.y = element_text(size = 12, family = "Arial"))
img
ggsave(file = "collinerity traits.svg", plot = img, width = 22, height = 22,  units = "cm")


## Supp. Material fig. S4 
# Prior this analysis, ensure to run "Material and Methods.R", and "Results_Flow cytometry traits.R" until L17
wss <- data.frame(
  k = 1:8,
  tot_withinss = sapply(1:8, function(k){
    kmeans(dat[,c(3:10)], centers = k, nstart = 2)$tot.withinss
  })
)

# Plot elbow
p <- ggplot(wss, aes(x = k, y = tot_withinss)) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Elbow Method for K-means",
    x = "Number of clusters (k)",
    y = "Total within-cluster sum of squares (WSS)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 8)
  )

ggsave(filename = "nb_kmeans_cluster.svg",plot = p, width=7, height=7,  units = "cm")


## Supp. Material fig. S5
# Run script "Packages Loading.R", "Material and Methods.R" and until line 26 of "Results_Functional analysis from individual level.R"
p1 <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis1, y = PCA$li$Axis3, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp1*3, yend = PCA$co$Comp3*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp1*4, y = PCA$co$Comp3*4, label = rownames(PCA$co)), col = "grey20", alpha = 0.7) +
  xlab("PC1 (51.8%)") +
  ylab("PC3 (11.4%)") +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p1 <- ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)

p2 <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis1, y = PCA$li$Axis4, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp1*3, yend = PCA$co$Comp4*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp1*4, y = PCA$co$Comp4*4, label = rownames(PCA$co)), col = "grey20", alpha = 0.7) +
  xlab("PC1 (51.8%)") +
  ylab("PC4 (9.3%)") +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2 <- ggMarginal(p2, groupColour = TRUE, groupFill = TRUE)

p3 <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = -PCA$li$Axis2, y = PCA$li$Axis3, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = -PCA$co$Comp2*3, yend = PCA$co$Comp3*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = -PCA$co$Comp2*4, y = PCA$co$Comp3*4, label = rownames(PCA$co)), col = "grey20", alpha = 0.7) +
  xlab("PC2 (27.4%)") +
  ylab("PC3 (11.4%)") +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3 <- ggMarginal(p3, groupColour = TRUE, groupFill = TRUE)

p4 <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = -PCA$li$Axis2, y = PCA$li$Axis4, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = -PCA$co$Comp2*3, yend = PCA$co$Comp4*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = -PCA$co$Comp2*4, y = PCA$co$Comp4*4, label = rownames(PCA$co)), col = "grey20", alpha = 0.7) +
  xlab("PC2 (27.4%)") +
  ylab("PC4 (9.3%)") +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p4 <- ggMarginal(p4, groupColour = TRUE, groupFill = TRUE)

p5 <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis3, y = PCA$li$Axis4, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp3*3, yend = PCA$co$Comp4*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp3*4, y = PCA$co$Comp4*4, label = rownames(PCA$co)), col = "grey20", alpha = 0.7, size = 1) +
  xlab("PC3 (11.4%)") +
  ylab("PC4 (9.3%)") +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
  theme(legend.position = "right",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

legend <- get_plot_component(p5, 'guide-box-right', return_all = TRUE)

p5 <- p5 + guides(color = "none")
p5 <- ggMarginal(p5, groupColour = TRUE, groupFill = TRUE)

p <- plot_grid(p1, p2, p3, p4, p5, legend, nrow = 2)
p
ggsave(file="PCA_axes.png", plot=p, width= 17.3, height=12, units = "cm")


## Supp. Material fig. S6
df <-
  get_dataframe_by_name(
    filename  = "IndCyano_Louchart_Limitation_experiment_June2025.xlsx",
    dataset   = "10.34894/9X9YMO",
    server = "dataverse.nl",
    .f = function(file) read_excel(file, sheet = "Supp. S6"),
  )

df_long <- df %>%
  pivot_longer(
    cols = c(Chla_start, Chla_end),
    names_to = "Phase",
    values_to = "Chla"
  ) %>%
  mutate(
    Phase = recode(Phase,
               Chla_start = "Start",
               Chla_end   = "End"),
    Phase = factor(Phase, levels = c("Start", "End"))
  )

df_sum <- df_long %>%
  group_by(`Start date`, Treatment, Phase) %>%
  summarise(
    mean = mean(Chla, na.rm = TRUE),
    SD   = sd(Chla, na.rm = TRUE),
    .groups = "drop"
  )


df_end <- df %>%
  select(`Start date`, Treatment, Chla_end) %>%
  filter(!is.na(Chla_end))

start_dates <- sort(unique(df_end$`Start date`))


letters_list <- list()

for (i in seq_along(start_dates)) {
  
  this_date <- start_dates[i]
  df_sub <- df_end %>% 
    filter(`Start date` == this_date)
  
  # Sécurité : au moins 2 traitements
  if (n_distinct(df_sub$Treatment) < 2) next
  
  fit <- aov(Chla_end ~ Treatment, data = df_sub)
  tuk <- TukeyHSD(fit)
  
  # Extraction des lettres
  letters <- multcompLetters4(fit, tuk)$Treatment
  
  letters_df <- data.frame(
    Treatment  = names(letters$Letters),
    Letters    = letters$Letters,
    `Start date` = this_date,
    stringsAsFactors = FALSE
  )
  
  letters_list[[i]] <- letters_df
}

letters_df <- bind_rows(letters_list)

letters_df <- letters_df %>%
  rename(`Start date` = Start.date)


df_sum <- df_sum %>%
  mutate(`Start date` = as.character(`Start date`)) %>%
  left_join(
    letters_df %>%
      mutate(`Start date` = as.character(`Start date`)),
    by = c("Start date", "Treatment")
  ) %>%
  mutate(
    Letters = ifelse(Phase == "End", Letters, NA_character_)
  )


p <- ggplot(df_sum,
       aes(x = Phase, y = mean,
           color = Treatment,
           group = Treatment)) +
  
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  geom_errorbar(
    aes(ymin = mean - SD, ymax = mean + SD),
    width = 0.05,
    size = 0.5
  ) +
  
  geom_text_repel(
    data = df_sum %>% filter(Phase == "End"),
    aes(
      label = Letters,
      color = Treatment,
    ),
    size = 5,
    direction = "both",
    force = 2,
    force_pull = 0.5,
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  
  facet_wrap(~ `Start date`, ncol = 3, scales = "free") +
  
  scale_color_manual(
    values = c(
      "Control" = "#8E6BBE",
      "N"       = "#FF7F0E",
      "NP"      = "#8C564B",
      "P"       = "#2CA02C"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Eq. chlorophyll-a concentration",
    color = "Treatment"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 11),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(filename = "figure_chla.svg", plot = p, width = 24, height = 18, units = "cm", dpi = 400)


## Supp. Material fig. S7
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


## Supp. Material fig. S8
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


## Supp. Material fig. S9 
# Prior the plot, load the datasets in  "Results_Functional assessment of natural communities.R" and run the script
p <- ggplot() +
  geom_point(data = subset(df, Source == "field"), mapping = aes(x = Axis1, y = Axis2, color = `TN (mmol/L)`/I_m), size = 2) +
  scale_color_gradientn(colors = desat_palette, name = "TN:Light ratio") +
  new_scale_color() +
  geom_point(data = subset(df, Source == "Experiment"), mapping = aes(x = Axis1, y = Axis2, color = Treatment), size = 2, shape = 18) +
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

