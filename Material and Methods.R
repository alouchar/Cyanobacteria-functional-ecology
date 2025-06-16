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

#### 1. CONCEPTUAL FIGURE
# Load the necessary library
library(ggplot2)
library(dplyr)

# Create 5 dummy data
set.seed(132)
data1 <- data.frame(
  Trait_1 = rnorm(5, mean = 7, sd = 1),
  Trait_2 = rnorm(5, mean = 7, sd = 1),
  Trait_3 = rnorm(5, mean = 4, sd = 0.9),
  Trait_4 = rnorm(5, mean = 2, sd = 1.4)
)

#
data2 <- data.frame(
  Trait_1 = rnorm(6, mean = 10, sd = 1),
  Trait_2 = rnorm(6, mean = 9, sd = 1),
  Trait_3 = rnorm(6, mean = 2, sd = 0.9),
  Trait_4 = rnorm(6, mean = 1, sd = 0.3)
)

#
data3 <- data.frame(
  Trait_1 = rnorm(10, mean = 17, sd = 1.2),
  Trait_2 = rnorm(10, mean = 11, sd = 1),
  Trait_3 = rnorm(10, mean = 8, sd = 0.7),
  Trait_4 = rnorm(10, mean = 10, sd = 1.4)
)

data1$Treatment <- "N"
data2$Treatment <- "P"
data3$Treatment <- "C"

dat <- rbind(data1, data2, data3) 

# Create the density contour plot for both datasets
p <-
  ggplot() +
  geom_density_2d(data = subset(dat, Treatment == "N"), aes(x = Trait_1, y = Trait_2), color = "#ff7f0e", bins = 4) +  # Contour for N data
  geom_point(data = subset(dat, Treatment == "N"), aes(x = Trait_1, y = Trait_2), color = "#ff7f0e", alpha = 0.6, size = 0.7) +
  geom_density_2d(data = subset(dat, Treatment == "P"), aes(x = Trait_1, y = Trait_2), color = "#2ca02c", bins = 4) +  # Contour for P data
  geom_point(data = subset(dat, Treatment == "P"), aes(x = Trait_1, y = Trait_2), color = "#2ca02c", alpha = 0.6, size = 0.7) +  # Contour for P data
  geom_density_2d(data = subset(dat, Treatment == "C"), aes(x = Trait_1, y = Trait_2), color = "#9467bd", bins = 4) +  # Contour for Control data
  geom_point(data = subset(dat, Treatment == "C"), aes(x = Trait_1, y = Trait_2), color = "#9467bd", alpha = 0.6, size = 0.7) +  # Contour for Control data
  xlim(5,21) +
  ylim(4,13.5) +
  labs(title = "Functional space \n(2 traits)",
       x = "Trait 1",
       y = "Trait 2") +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),                 # Remove axis ticks
    axis.title = element_text(size = 8),         # Adjust axis labels size
    plot.title = element_text(size = 8),         # Adjust title size
    axis.text = element_blank()                   # Remove axis tick labels
  )
p
ggsave(file="Functional space 2 traits.svg", plot=p, width=5, height=5, units = "cm")
  
## Dummy pca on the data
dat <- rbind(data1,data2,data3)
  
PCA <- ade4::dudi.pca(dat[,c(1:4)],  center = TRUE, scale =  TRUE, scannf = FALSE, nf = 7)
  
screeplot(PCA, main = "Screeplot - Eigenvalues")
  
(PCA$eig*100)/sum(PCA$eig)
  
# 7.2) Plot simple PCA
ade4::s.corcircle(PCA$co)

dat <- cbind(dat,PCA$li)


p <- PCA$li %>%
  ggplot(mapping = aes(x = Axis1, y = Axis2, color = as.factor(dat$Treatment))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 0.7) +
  scale_colour_manual(values = c("#9467bd","#ff7f0e","#2ca02c")) +
  geom_density_2d(data = subset(dat, Treatment == "P"), aes(x = Axis1, y = Axis2), color = "#2ca02c", bins = 4) +  # Contour for N data
  geom_density_2d(data = subset(dat, Treatment == "C"), aes(x = Axis1, y = Axis2), color = "#9467bd", bins = 4) +  # Contour for P data
  geom_density_2d(data = subset(dat, Treatment == "N"), aes(x = Axis1, y = Axis2), color = "#ff7f0e", bins = 4) +  # Contour for Control data
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2) +
  labs(title = "Functional space \n(n traits)",
       x = "Dim. 1",
       y = "Dim. 2") +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),                 # Remove axis ticks
    axis.title = element_text(size = 8),         # Adjust axis labels size
    plot.title = element_text(size = 8),         # Adjust title size
    axis.text = element_blank(), # Remove axis tick labels
    legend.position = "none"
  )
p
ggsave(file="Functional space n traits.svg", plot=p, width=5, height=5, units = "cm")

rm(data1,data2,data3)

DATA <- dat %>%
  group_by(Treatment) %>%
  summarise(across(c(Axis1, Axis2, Axis3, Axis4), mean))

##Extract centroid of func. space ------------------ N
p_N <- ggplot(subset(dat,Treatment == 'N'), mapping = aes(x = Axis1, y = Axis2)) +
  geom_density_2d(alpha = 0.4, bins = 4) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2)
  
contours_N <- layer_data(p_N)
contours_N <- contours_N[contours_N$piece == 5,]

contours <- NULL

contours$x <- mean(contours_N$x)
contours$y <- mean(contours_N$y)
contours$Treatment <- "N"

contours <- as.data.frame(unlist(contours))
centroid <- contours

##Extract centroid of func. space ------------------ P
p_P <- ggplot(subset(dat,Treatment == 'P'), mapping = aes(x = Axis1, y = Axis2)) +
  geom_density_2d(alpha = 0.4, bins = 4) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2)

contours_P <- layer_data(p_P)
contours_P <- contours_P[contours_P$piece == 1,]

contours <- NULL

contours$x <- mean(contours_P$x)
contours$y <- mean(contours_P$y)
contours$Treatment <- "P"

contours <- as.data.frame(unlist(contours))
centroid <- cbind(centroid,contours)

##Extract centroid of func. space ------------------ C
p_C <- ggplot(subset(dat,Treatment == 'C'), mapping = aes(x = Axis1, y = Axis2)) +
  geom_density_2d(alpha = 0.4, bins = 4) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2)

contours_C <- layer_data(p_C)
contours_C <- contours_C[contours_C$piece == 3,]

contours <- NULL

contours$x <- mean(contours_C$x)
contours$y <- mean(contours_C$y)
contours$Treatment <- "C"

contours <- as.data.frame(unlist(contours))
centroid <- cbind(centroid,contours)

centroid <- as.data.frame(t(centroid))
centroid$x <- as.numeric(as.character(centroid$x))
centroid$y <- as.numeric(as.character(centroid$y))

require(FNN)
require(spdep)

dat$id <- rownames(dat)

find_neighbors <- function(df) {
  # Compute pairwise distances
  dist_matrix <- as.matrix(dist(df[, c("Axis1", "Axis2")]))
  
  # Initialize edge list
  edges <- data.frame(from = integer(), to = integer(), stringsAsFactors = FALSE)
  
  # Find nearest neighbor for each point
  for (i in 1:nrow(df)) {
    distances <- dist_matrix[i, ]
    nearest_index <- which.min(distances[distances > 0])  # Exclude self distance
    edges <- rbind(edges, data.frame(from = df$id[i], to = df$id[nearest_index]))
  }
  
  return(edges)
}


edge_data <- dat %>%
  group_by(Treatment) %>%
  do(find_neighbors(.)) %>%
  ungroup()

# Merge edge_data with node data to get coordinates for plotting
edge_data <- merge(dat,edge_data, by.x = c("Treatment","id"), by.y = c("Treatment","from"))

edge_data <- edge_data %>%
  left_join(dat, by = c("from" = "id")) %>%
  rename(x_from = Axis1, y_from = Axis2) %>%
  left_join(dat, by = c("to" = "id")) %>%
  rename(x_to = Axis1, y_to = Axis2)

#### Functional Richness

p <- ggplot(dat, mapping = aes(x = Axis1, y = Axis2, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("#9467bd","#ff7f0e","#2ca02c")) +
  stat_density_2d(data = subset(dat, Treatment == "P"), aes(x = Axis1, y = Axis2, fill = ..level..), alpha = 0.4, geom = "polygon", bins = 4) +
  scale_fill_gradient(low = "grey70", high = "grey70") +  # Contour for N data
  stat_density_2d(data = subset(dat, Treatment == "C"), aes(x = Axis1, y = Axis2, fill = ..level..), alpha = 0.4, geom = "polygon", bins = 4) +
  stat_density_2d(data = subset(dat, Treatment == "N"), aes(x = Axis1, y = Axis2, fill = ..level..), alpha = 0.4, geom = "polygon", bins = 4) +  # Contour for Control data
  geom_point(alpha = 0.6, size = 0.7) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2) +
  labs(title = "Size of functional space",
       x = "Dim. 1",
       y = "Dim. 2") +
  annotate("text", x=-1.9, y=1.9, label= "F.Ric = High", color = "#ff7f0e", size = 1.7 ) +
  annotate("text", x=-1.9, y=1.6, label= "F.Ric = Medium", color = "#9467bd", size = 1.7 ) +
  annotate("text", x=-1.9, y=1.3, label= "F.Ric = Low", color = "#2ca02c", size = 1.7 ) +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),                 # Remove axis ticks
    axis.title = element_text(size = 8),         # Adjust axis labels size
    plot.title = element_text(size = 8),         # Adjust title size
    axis.text = element_blank(), # Remove axis tick labels
    legend.position = "none"
  )
p
ggsave(file="Functional Richness.svg", plot=p, width=5, height=5, units = "cm")

#### Functional Evenness
p <-
  ggplot(dat, mapping = aes(x = Axis1, y = Axis2, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("#9467bd","#ff7f0e","#2ca02c")) +
  geom_segment(data = edge_data, aes(x = x_from, y = y_from, xend = x_to, yend = y_to), color = "gray") +
  geom_point(alpha = 0.6, size = 0.7) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2) +
  labs(title = "Regularity of trait distribution \n in functional space)",
       x = "Dim. 1",
       y = "Dim. 2") +
  annotate("text", x=-1.9, y=1.9, label= "F.Eve = Low", color = "#ff7f0e", size = 1.7 ) +
  annotate("text", x=-1.9, y=1.6, label= "F.Eve = Medium", color = "#9467bd", size = 1.7) +
  annotate("text", x=-1.9, y=1.3, label= "F.Eve = High", color = "#2ca02c", size = 1.7 ) +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),                 # Remove axis ticks
    axis.title = element_text(size = 8),         # Adjust axis labels size
    plot.title = element_text(size = 8),         # Adjust title size
    axis.text = element_blank(), # Remove axis tick labels
    legend.position = "none"
  )
p
ggsave(file="Functional Evenness.svg", plot=p, width=5, height=5, units = "cm")

DATA <- merge(dat, centroid, by = "Treatment", all.x = TRUE)

#### Functional Dispersion
p <- 
  ggplot(DATA, mapping = aes(x = Axis1, y = Axis2, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("#9467bd","#ff7f0e","#2ca02c")) +
  geom_point(data = centroid, aes(x = x, y = y), color = 'black', alpha = 0.6, size = 0.9, shape = 4) +
  geom_segment(mapping = aes(x = x, y = y, xend = Axis1, yend = Axis2), linewidth = 0.2, linetype = "dashed") +
  # geom_density_2d(data = subset(dat, Treatment == "P"), aes(x = Axis1, y = Axis2), alpha = 0.4, bins = 4) +
  # geom_density_2d(data = subset(dat, Treatment == "C"), aes(x = Axis1, y = Axis2), alpha = 0.4, bins = 4) +
  # geom_density_2d(data = subset(dat, Treatment == "N"), aes(x = Axis1, y = Axis2), alpha = 0.4, bins = 4) +  # Contour for Control data
  geom_point(alpha = 0.6, size = 0.7) +
  theme_linedraw() +
  xlim(-2.6,3) +
  ylim(-2.6,2) +
  labs(title = "Deviation from the \n functional centroid",
       x = "Dim. 1",
       y = "Dim. 2") +
  annotate("text", x=-1.9, y=1.9, label= "F.Disp = High", color = "#ff7f0e", size = 1.7 ) +
  annotate("text", x=-1.9, y=1.6, label= "F.Disp = Medium", color = "#9467bd", size = 1.7) +
  annotate("text", x=-1.9, y=1.3, label= "F.Disp = Low", color = "#2ca02c", size = 1.7 ) +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),                 # Remove axis ticks
    axis.title = element_text(size = 8),         # Adjust axis labels size
    plot.title = element_text(size = 8),         # Adjust title size
    axis.text = element_blank(), # Remove axis tick labels
    legend.position = "none"
  )
ggsave(file="Functional Dispersion.svg", plot=p, width=5, height=5, units = "cm")
p

rm(list=setdiff(ls(), c("dat","edge_data", "DATA")))

#### 2. DATA LOADING
# IndTraits contains flow cytometry informations at the individual level for each treatment obtained through the PhytoCytoTraits GUI.
# Treatment contains informations on the different environmental conditions

dat <- readxl::read_excel("D:/IndCyano_Louchart_Limitation_experiment_June2025.xlsx", sheet = "IndTraits")

#### 3. DATA MANIPULATION
## 3.1. Log10 transformation
dat[,c(3:10)] <- log10(dat[,c(3:10)])

## 3.2. Z-score normalisation
dat[,c(3:10)] <- scale(dat[,c(3:10)], center = TRUE, scale = TRUE)
