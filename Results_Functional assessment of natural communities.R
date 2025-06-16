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

## 6. Formating functional fingerprints
# 6.1. Reload the functional responses of Microcystis obtained at the individual level
dat <- readxl::read_excel("IndCyano_Louchart_Limitation_experiment_June2025.xlsx", sheet = "IndTraits")

# 6.2. Compute the average functional responses of Microcystis
Ref_transform_average <- 
  dat %>%
  group_by(Treatment) %>%
  dplyr::summarise(across(c('Cell size', Granularity,'Phycocyanin', 'Chlorophyll-a', 'Gas vesicle','PC per size unit', 'Chla per size unit', 'PC:Chla ratio'), mean))

# 6.3. Log10 transformation and normalisation
Ref_transform_average[,c(2:9)] <- log10(Ref_transform_average[,c(2:9)])
Ref_transform_average[,c(2:9)] <- scale(Ref_transform_average[,c(2:9)], center = TRUE, scale = TRUE)

## 7. Functional assessment of natural communities
CommTraits <- readxl::read_excel("IndCyano_Louchart_Limitation_experiment_June2025.xlsx", sheet = "CommTraits")
CommTraits$Date <- as.Date(CommTraits$Date, format = "%d-%m-%Y")

# 7.1. Log10 transformation and normalisation
CommTraits[,c(4:11)] <- log10(CommTraits[,c(4:11)])
CommTraits <- CommTraits[complete.cases(CommTraits),]
CommTraits[,c(4:11)] <- scale(CommTraits[,c(4:11)], center = TRUE, scale = TRUE)

# 7.2. Plots to produce figure 5
p1 <- 
  ggplot(data = CommTraits,aes(x = Date, y = Phycocyanin), shape = 16) +
  geom_rect(aes(xmin = min(CommTraits$Date), xmax = max(CommTraits$Date), ymin = -2.5, ymax = Ref_transform_average$Phycocyanin[5]), 
            fill =  "#FF7F0EFF", alpha = 0.01) + # Rectangle for "likely N limitation"
  geom_rect(aes(xmin = min(CommTraits$Date), xmax = max(CommTraits$Date), ymin = Ref_transform_average$Phycocyanin[5], ymax = 2.5), 
            fill = "#1F77B4FF", alpha = 0.01) + 
  annotate("text", x = as.Date("2023-05-15"), y = -2.25, 
           label = "Likely\nnitrogen limited", color = "black", size = 3) +
  annotate("text", x = as.Date("2023-05-15"), y = 2.25, 
           label = "Likely\nlight limited", color = "black", size = 3) +
  geom_point( ) +
  geom_line() +
  geom_point(data = subset(CommTraits, `N Limitation` == "Y"), aes(x = Date, y = Phycocyanin), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  geom_hline(yintercept = Ref_transform_average$Phycocyanin[5], linetype = "dotted") +
  stat_smooth(mapping = aes(y = Phycocyanin, x = Date), method = "lm", se = FALSE, col = "grey20", linetype = "dashed") +
  # stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4, label.x = 0.95) +
  ylim(-2.5,2.5) +
  theme_bw() +
  xlab("Date") +
  ylab("Phycocyanin fluorescence (a.u.)") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )
p1

p2 <-  
  ggplot(data = CommTraits, aes(x = Date, y = `Chlorophyll-a`), shape = 16) +
  geom_rect(aes(xmin = min(CommTraits$Date), xmax = max(CommTraits$Date), ymin = -2.5, ymax = Ref_transform_average$`Chlorophyll-a`[5]), 
            fill =  "#FF7F0EFF", alpha = 0.01) +  # Rectangle for "likely N limitation"
  geom_rect(aes(xmin = min(CommTraits$Date), xmax = max(CommTraits$Date), ymin = Ref_transform_average$`Chlorophyll-a`[5], ymax = 2.5), 
            fill = "#1F77B4FF", alpha = 0.01) + 
  annotate("text", x = as.Date("2023-05-15"), y = -2.25, 
           label = "Likely\nnitrogen limited", color = "black", size = 3) +
  annotate("text", x = as.Date("2023-05-15"), y = 2.25, 
           label = "Likely\nlight limited", color = "black", size = 3) +
  geom_hline(yintercept = Ref_transform_average$`Chlorophyll-a`[5], linetype = "dotted") +
  geom_point() +
  geom_line() +
  geom_point(data = subset(CommTraits, `N Limitation` == "Y"), aes(x = Date, y = `Chlorophyll-a`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  stat_smooth(mapping = aes(y = `Chlorophyll-a`, x = Date), method = "lm", se = FALSE, col = "grey20", linetype = "dashed") +
  # stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4, label.x = 0.95) +
  ylim(-2.5,2.5) +
  theme_bw() +
  xlab("Date") +
  ylab("Chlorophyll a fluorescence (a.u.)") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )
p2

p3 <- ggplot() +
  geom_point(data = CommTraits, mapping = aes(`Chlorophyll-a`, Phycocyanin , col = Date), shape = 16) +
  geom_point(data = subset(CommTraits, `N Limitation` == "Y"), mapping = aes(`Chlorophyll-a`, Phycocyanin), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_color_gradient(low = "grey90", high = "grey10", name = "Date", breaks = as.numeric(CommTraits$Date[seq(1, nrow(CommTraits), length.out = 4)]),
                       labels = c("12-Apr.","mid June","mid Aug.","25-Oct.")
  ) +
  new_scale_color() +
  geom_point(Ref_transform_average, mapping = aes(x = `Chlorophyll-a`, y = Phycocyanin, color = Treatment), size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(data = CommTraits, aes(x = `Chlorophyll-a`, y = Phycocyanin), se = FALSE, method = "lm", col = "grey20", linetype = "dashed") +
  scale_color_d3(name = "Treatment") +
  theme_bw() +
  xlim(-2.5,2.5) +
  ylim(-2.5,2.5) +
  xlab("Chlorophyll a fluorescence (a.u.)") +
  ylab("Phycocyanin fluorescence (a.u.)") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )
p3

bottom <- plot_grid(p3, labels = c('C'))
top_row <- plot_grid(p1, p2, labels = c('A','B'))
p <- plot_grid(top_row, bottom, ncol =1)
p

ggsave(file="PC_Chla over time.svg", plot=p, width=20, height=20, units = "cm")
