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

#### FUNCTIONAL ANALYSES FROM INDIVIDUAL LEVEL

## 2. Computation of functional space 
# install_github("arleyc/PCAtest")

## 2.1. Dimension reduction through a Principal Component Analysis
PCA <- dudi.pca(DATA[,c(1:8)], scannf = FALSE, nf = 7)

## 2.2. Quick visualisation of eigenvalues
screeplot(PCA, main = "Screeplot - Eigenvalues")

## 2.3. Eigenvalues contribution
(PCA$eig*100)/sum(PCA$eig)

## 2.4. Quick visualisation of the PCA
s.corcircle(PCA$co)

## 2.5. Testing the number of axis to keep
PCAtest(PCA$tab, 999, 999, 0.001, varcorr=TRUE, counter=FALSE, plot=TRUE)

## 2.6. Contribution of each variable
rowSums(100*(factoextra::get_pca_var(PCA)$cos2)[,c(1,2)])


## 2.7. Plot PCA Biplot with density curve per environmental conditions (superimposition check between environmental conditions; Figure 3)
p <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis1, y = PCA$li$Axis2, color = as.factor(DATA$Label)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp1*3, yend = PCA$co$Comp2*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp1*4, y = PCA$co$Comp2*4, label = rownames(PCA$co)), col = "grey20") +
  xlab("PC1 (51.8%)") +
  ylab("PC2 (27.4%)") +
  xlim(c(-9,9)) +
  ylim(c(-9,9)) +
  theme_linedraw() +
  scale_color_d3(name = as.factor(DATA$Label)) +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        axis.text=element_text(size=12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p <- p + guides(colour = guide_legend(override.aes = list(size=14)))
pca_marg <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
pca_marg

ggsave(file="PCA_density.svg", plot=pca_marg, width=16, height=16, units = "cm", dpi = 400)

rm(heatmap, NOM, num_mat, treatment)

# 2.8. PCA eigenvalues dataframe binded with environmental conditions and replicate
reduced_dim <- as.data.frame(cbind(DATA, PCA$li$Axis1, PCA$li$Axis2, PCA$li$Axis3, PCA$li$Axis4))

names(reduced_dim)[11] <- "PCA1"
names(reduced_dim)[12] <- "PCA2"
names(reduced_dim)[13] <- "PCA3"
names(reduced_dim)[14] <- "PCA4"

my_theme <- theme(axis.text=element_text(size=10),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position.inside = c(0.90,0.30),
                  legend.direction = "vertical",
                  legend.key = element_rect(fill = "transparent"),
                  legend.key.size = unit(.7,"line"))

rm(list = setdiff(ls(),"reduced_dim"))

saveRDS(reduced_dim, "reduced_dim.rds")
reduced_dim <- readRDS("N:/NIOO-Data/Dep.AqE/Staging/van de Waal_group/Zhipeng data/reduced_dim.rds")

## 2.9. Delineation of multidimensional functional space

# The methodology relies on the probabilistic hypervolume since it considers abundance thus being less sensitive to outliers
bw_estimate <- hypervolume::estimate_bandwidth(reduced_dim[, c("PCA1","PCA2")], method = "cross-validation")


# 3. Computation of functional diversity indices based on functional hypervolume
# 3.1. Code optimization running foreach function on 4 cores
require(foreach)
require(doSNOW)
require(progress)

# adjust the number of core to use
num_cores <- 2

# Create a parallel cluster
cl <- makeCluster(num_cores)

# Save the cluster to use in foreach loop
registerDoSNOW(cl)

# Subset by environmental condition. Here an example with phosphorus limitation
temp <- reduced_dim %>%
  filter(Label == "-Phosphate")

# Define the unique replicates
replicat <- unique(temp$filename)

# Set a progression bar
pb <- txtProgressBar(max = length(replicat), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Loop parallelisation
results <- foreach(l = 1:length(replicat), .combine = c,.options.snow = opts) %dopar% {
  DATA <- temp[temp$filename == replicat[l], ]
  
  # Producing hypervolume for each replicate within the treatment
  hypervolume::hypervolume_gaussian(data = DATA[, c(12:13)], kde.bandwidth = bw_estimate, quantile.requested = 0.95, quantile.requested.type = "probability")
}

# Stop the progression bar
close(pb)

# Merge the results of the replicates into a single environmental condition
phosphate_vol <- hypervolume_join(results)

# Local save of the hypervolume data 
save(phosphate_vol, file = "Phosphate hypervolume.RData")

# Stop the parallel cluster
stopCluster(cl)

# When arrived here, re-run from line 108 with change of environmental condition.

# Load the data previously saved. Example with control conditions
load("Control hypervolume.Rdata")

## 3.2. Computation of functional diversity indices from the hypervolume
Control <- data.frame(
  Treatment = "Control",
  rich = kernel.alpha(Control_vol),
  even = kernel.evenness(Control_vol),
  div = kernel.dispersion(Control_vol)
)

# Save the indices on a csv file
write.table(Control, "Control.csv", row.names = FALSE, sep= ";", dec = ',')
