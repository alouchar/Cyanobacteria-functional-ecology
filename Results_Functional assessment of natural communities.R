#`---
#`Authors: Arnaud Louchart, Annemieke Drost, Harley Lin, Suzanne Wiezer, Zhipeng Duan, Dedmer B. Van de Waal
#`Title: Individual-level trait responses in natural cyanobacterial communities.
#`Year: 202X
#`Journal: Functional Ecology
#`Status: Writing
#`---

#############################################################################
############################## RESULTS SECTION ##############################
#############################################################################
Ref_transform_average <- 
  DATA %>%
  group_by(Label) %>%
  dplyr::summarise(across(c(`FSC [Par]`,SSC, `692/40[488]`, `670/30[640]`, Fill_SSC, PC_size, Chla_size, PC_Chla), mean))

#### FUNCTIONAL ASSESSMENT OF NATURAL COMMUNITIES
dat <- readxl::read_xlsx("C:/Users/ArnaudL/Desktop/BloomTox/Temporal FCM/Data/Natural community Netherlands 2023.xlsx", sheet = "Abundance")
Bloomtox <- readxl::read_xlsx("C:/Users/ArnaudL/Desktop/BloomTox/Temporal FCM/Data/Natural community Netherlands 2023.xlsx", sheet = "Trait")


## Pivot wide to long traits
Trait_table <- Bloomtox %>%
  pivot_longer(!c(Sample, Date, Lake, FCM_trait), names_to = "cluster", values_to = "values") %>%
  filter(!cluster == "All events" & !cluster == "CHLA AND NOT PHY" & !cluster == "CHLA AND PHY") %>%
  separate(cluster, into = c("size","pigment"), sep = "AND", remove = FALSE, extra = "merge") %>%
  rowwise() %>% 
  mutate(size = gsub("over", '>', size)) %>%
  mutate(size = gsub("below", '<', size)) %>%
  ungroup()

## Pivot wide to long abundance
Abundance_table <- dat %>%
  pivot_longer(!c(Sample, Date, Lake), names_to = "cluster", values_to = "Abundance") %>%
  filter(!cluster == "All events" & !cluster == "CHLA AND NOT PHY" & !cluster == "CHLA AND PHY") %>%
  separate(cluster, into = c("size","pigment"), sep = "AND", remove = FALSE, extra = "merge") %>%
  rowwise() %>% 
  mutate(size = gsub("over", '>', size)) %>%
  mutate(size = gsub("below", '<', size)) %>%
  ungroup()


merged_table <- merge(Trait_table, Abundance_table, by.x = c("Sample","size","pigment"), by.y = c("Sample","size","pigment"), all.x = TRUE)
merged_table <- merged_table[,-c(7,9,10,11)]
merged_table$Total_trait <- merged_table$values*merged_table$Abundance
names(merged_table)[4] <- "Date"
names(merged_table)[5] <- "Lake"


##
merged_table <-
  merged_table %>%
  filter(str_detect(size, "< 05|05-1|1-2|2-10")) %>%
  filter(str_detect(Lake, "Grote Plas")) %>%
  filter(str_detect(pigment, "CHLA AND PHY"))

##
Trait_table <-
  merged_table %>%
  dplyr::group_by(Sample, Date, Lake, FCM_trait) %>%
  dplyr::summarise(Abundance = sum(Abundance),
                   Total_trait = sum(values))

Trait_table <- Trait_table[is.finite(Trait_table$Total_trait),]
# 
# 
Trait_table$individual_trait <-
  Trait_table$Total_trait/Trait_table$Abundance


Trait_table <-
  Trait_table[,-6] %>% 
  pivot_wider(names_from = FCM_trait, values_from = individual_trait)

Trait_table$Fill_SSC <- Trait_table$SSC/(4/3*(pi*(Trait_table$FSC/2)^3))
Trait_table$PC_size <- Trait_table$`670/30`/Trait_table$FSC
Trait_table$Chla_size <- Trait_table$`692/40`/Trait_table$FSC
Trait_table$PC_Chla <- Trait_table$`670/30`/Trait_table$`692/40`

Trait_table[,c(5:12)] <- log10(Trait_table[,c(5:12)])

Trait_table <- Trait_table[complete.cases(Trait_table),]

Trait_table[,c(5:12)] <- scale(Trait_table[,c(5:12)], center = TRUE, scale = TRUE)


Lake_data <- readxl::read_excel("C:/Users/ArnaudL/Desktop/BloomTox/Temporal FCM/Data/Lake_data.xlsx")
Lake_data$Date <- as.Date(Lake_data$Date, format = "%d-%m-%Y")

dat <- merge(Trait_table, Lake_data, by.x = c("Lake","Date"), by.y = c("Lake","Date"))

N_lim <-
  Trait_table %>%
  filter(Lake == "Grote Plas") %>%
  filter(Date == as.Date("2023-06-07") |
           Date == as.Date("2023-06-21") |
           Date == as.Date("2023-07-12") |
           Date == as.Date("2023-07-27")|
           Date == as.Date("2023-08-09")|
           Date == as.Date("2023-08-30")|
           Date == as.Date("2023-09-20"))

Trait_table$Date <- as.Date(Trait_table$Date, format = "%d-%m-%Y")
N_lim$Date <- as.Date(N_lim$Date, format = "%d-%m-%Y")


p1 <- 
  ggplot(data = Trait_table,aes(x = Date, y = `670/30`), shape = 16) +
  geom_rect(aes(xmin = min(Trait_table$Date), xmax = max(Trait_table$Date), ymin = -2.5, ymax = Ref_transform_average$`670/30[640]`[5]), 
            fill =  "#FF7F0EFF", alpha = 0.01) +  # Rectangle for "likely N limitation"
  geom_rect(aes(xmin = min(Trait_table$Date), xmax = max(Trait_table$Date), ymin = Ref_transform_average$`670/30[640]`[5], ymax = 2.5), 
            fill = "#1F77B4FF", alpha = 0.01) + 
  annotate("text", x = as.Date("2023-05-15"), y = -2.25, 
           label = "Likely\nnitrogen limited", color = "black", size = 3) +
  annotate("text", x = as.Date("2023-05-15"), y = 2.25, 
           label = "Likely\nlight limited", color = "black", size = 3) +
  geom_point( ) +
  geom_line() +
  geom_point(data = N_lim, aes(x = Date, y = `670/30`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  geom_hline(yintercept = Ref_transform_average$`670/30[640]`[5], linetype = "dotted") +
  stat_smooth(mapping = aes(y = `670/30`, x = Date), method = "lm", se = FALSE, col = "grey20", linetype = "dashed") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4, label.x = 0.95) +
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
  ggplot(data = Trait_table, aes(x = Date, y = `692/40`), shape = 16) +
  geom_rect(aes(xmin = min(Trait_table$Date), xmax = max(Trait_table$Date), ymin = -2.5, ymax = Ref_transform_average$`692/40[488]`[5]), 
            fill =  "#FF7F0EFF", alpha = 0.01) +  # Rectangle for "likely N limitation"
  geom_rect(aes(xmin = min(Trait_table$Date), xmax = max(Trait_table$Date), ymin = Ref_transform_average$`692/40[488]`[5], ymax = 2.5), 
            fill = "#1F77B4FF", alpha = 0.01) + 
  annotate("text", x = as.Date("2023-05-15"), y = -2.25, 
           label = "Likely\nnitrogen limited", color = "black", size = 3) +
  annotate("text", x = as.Date("2023-05-15"), y = 2.25, 
           label = "Likely\nlight limited", color = "black", size = 3) +
  geom_hline(yintercept = Ref_transform_average$`692/40[488]`[5], linetype = "dotted") +
  geom_point() +
  geom_line() +
  geom_point(data = N_lim, aes(x = Date, y = `692/40`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  stat_smooth(mapping = aes(y = `692/40`, x = Date), method = "lm", se = FALSE, col = "grey20", linetype = "dashed") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4, label.x = 0.95) +
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
  geom_point(data = Trait_table, mapping = aes(`692/40`, `670/30` , col = Date), shape = 16) +
  geom_point(data = N_lim, mapping = aes(`692/40`, `670/30`), col = "orange", size = 3, shape = 21,  stroke = 1.2) +
  scale_color_gradient(low = "grey90", high = "grey10", name = "Date", breaks = as.numeric(Trait_table$Date[seq(1, nrow(Trait_table), length.out = 4)]),
                       labels = c("12-Apr.","mid June","mid Aug.","25-Oct.")
  ) +
  new_scale_color() +
  geom_point(Ref_transform_average, mapping = aes(x = `692/40[488]`, y = `670/30[640]`, color = Label), size = 5) +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(data = Trait_table, aes(x = `692/40`, y = `670/30`), se = FALSE, method = "lm", col = "grey20", linetype = "dashed") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "P")), size = 4) +
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

bottom <- plot_grid(p3, labels = c('C'))
top_row <- plot_grid(p1, p2, labels = c('A','B'))
p <- plot_grid(top_row, bottom, ncol =1)
p

ggsave(file="PC_Chla over time.svg", plot=p, width=20, height=20, units = "cm")