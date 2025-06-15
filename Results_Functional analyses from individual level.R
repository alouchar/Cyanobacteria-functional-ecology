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

#### FUNCTIONAL ANALYSES FROM INDIVIDUAL LEVEL

## 2. Computation of functional space 
# install_github("arleyc/PCAtest")

## 2.1. Dimension reduction through a Principal Component Analysis
PCA <- dudi.pca(dat[,c(3:10)], scannf = FALSE, nf = 7)

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


## 2.7. Plot PCA Biplot with density curve per environmental conditions (superimposition check between environmental conditions)
# Figure 3
p <- ggplot() +
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis1, y = PCA$li$Axis2, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp1*3, yend = PCA$co$Comp2*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp1*4, y = PCA$co$Comp2*4, label = rownames(PCA$co)), col = "grey20") +
  xlab("PC1 (51.8%)") +
  ylab("PC2 (27.4%)") +
  xlim(c(-9,9)) +
  ylim(c(-9,9)) +
  theme_linedraw() +
  scale_color_d3(name = as.factor(dat$Treatment)) +
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

rm(heatmap, NOM, num_mat)

# 2.8. PCA eigenvalues dataframe binded with environmental conditions and replicate
reduced_dim <- as.data.frame(cbind(dat, PCA$li$Axis1, PCA$li$Axis2, PCA$li$Axis3, PCA$li$Axis4))

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


## 2.9. Delineation of multidimensional functional space

# The methodology relies on the probabilistic hypervolume since it considers abundance thus being less sensitive to outliers
bw_estimate <- hypervolume::estimate_bandwidth(reduced_dim[, c("PCA1","PCA2")], method = "cross-validation")

#### The steps 3.1 and 3.2 display the code for the phosphorus treatment. Replace phosphate by nitrogen, control, light and CO2 to produce the results for these treatments

## 3. Computation of functional diversity indices based on functional hypervolume (step 3.1 and 3.2 are quite long, may takes several days for one treatment)
## 3.1. Code optimization running foreach function on 4 cores
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
  filter(Treatment == "-Phosphate")

# Define the unique replicates
replicat <- unique(temp$Replicate)

# Set a progression bar
pb <- txtProgressBar(max = length(replicat), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Loop parallelisation
results <- foreach(l = 1:length(replicat), .combine = c,.options.snow = opts) %dopar% {
  DATA <- temp[temp$filename == replicat[l], ]
  
  # Producing hypervolume for each replicate within the treatment
  hypervolume::hypervolume_gaussian(data = DATA[, c(11:12)], kde.bandwidth = bw_estimate, quantile.requested = 0.95, quantile.requested.type = "probability")
}

# Stop the progression bar
close(pb)
phosphorus_vol <- hypervolume_join(results)

# Saving the hypervolume data
save(phosphorus_vol, file = "Phosphorus hypervolume.RData")

# Stop the parallel cluster
stopCluster(cl)

# Load the data previously saved. Example with control conditions
load("Phosphorus hypervolume.Rdata")

## 3.2. Computation of functional diversity indices from the hypervolume
Phosphorus <- data.frame(
  Treatment = "Phosphorus",
  rich = kernel.alpha(phosphorus_vol),
  even = kernel.evenness(phosphorus_vol),
  div = kernel.dispersion(phosphorus_vol)
)

# Combine all into one data frame
Functional_rep_space_microcystis <- rbind(Phosphorus, Control, CO2, Light, Nitrogen)

## 3.3. Read the excel sheet that resume the results 
dat <- readxl::read_xlsx("IndCyano_Louchart_Limitation_experiement_June2025.xlsx", sheet = "Functional_rep_space_microcystis")


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

alpha_func <- data_summary(data = Functional_rep_space_microcystis, varname = "div", groupnames = "Treatment")

alpha_func$Treatment <- gsub("CO2", "+CO2", alpha_func$Treatment)
alpha_func$Treatment <- stringr::str_replace_all(alpha_func$Treatment, "Light", "-L")
alpha_func$Treatment <- gsub("Nitrogen", "-N", alpha_func$Treatment)
alpha_func$Treatment <- gsub("Phosphate", "-P", alpha_func$Treatment)
alpha_func$Treatment <- gsub("Control", "Cont", alpha_func$Treatment)


dat <- alpha_func

## Run from this line (l175) to l279 3 times: 1) functional richness, 2) functional evenness and 3) functional divergence.
# Replace the column names accordingly
dat <- cbind(dat,alpha_func)

dat <- dat[,-c(4,7)]

names(dat)[3] <- "sd_rich"
names(dat)[5] <- "sd_even"
names(dat)[7] <- "sd_div"

rm(list = setdiff(ls(),c('dat',"Functional_rep_space_microcystis")))


test_results <- kruskal.test(div ~ Treatment, data = alpha) #replace by rich, even or div
test_results$statistic <- round(test_results$statistic, 2)
test_results$p.value <- round(test_results$p.value, 3)


## 3.4. Conover post-hoc test
require(conover.test)
post_hoc <- conover.test(alpha$div, alpha$Treatment, kw = FALSE, alpha = 0.05, label = TRUE, list = FALSE, method="bh") #replace by rich, even or div


## 3.5. Extract pairwise significance as letters
v <- as.data.frame(cbind(comparisons = post_hoc$comparisons, p_val = post_hoc$P.adjusted*2))
require(stringr)

v[c('Treatment 1','Treatment 2')] <- str_split_fixed(v$comparisons, " - ", 2)
v <- v[,-1]

nameVals <- sort(unique(unlist(v[2:3])))
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
myMat[as.matrix(v[c("Treatment 1", "Treatment 2")])] <- v[["p_val"]]
diag(myMat) <- 1

myMat <- apply(myMat, 2, as.numeric)
rownames(myMat) <- colnames(myMat)

myMat <- as.matrix(Matrix::forceSymmetric(myMat, "U"))

require(multcompView)

if (!require(multcompView)){install.packages("multcompView")} 

# Matrix creation
nameVals <- sort(unique(unlist(v[2:3])))
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
myMat[as.matrix(v[c("Treatment 1", "Treatment 2")])] <- v[["p_val"]]
diag(myMat) <- 1

myMat <- apply(myMat, 2, as.numeric)
rownames(myMat) <- colnames(myMat)

# Define the type of matrix
myMat <- as.matrix(Matrix::forceSymmetric(myMat, "U"))

# Set letters with 0.05 as p-value
LET <- multcompLetters(myMat,
                       compare="<",
                       threshold=0.025,
                       Letters=letters,
                       reversed = FALSE)

# Display small letters
LET <- as.data.frame(LET[["Letters"]])

LET$Treatment <- row.names(LET)

LET$Treatment <- gsub("CO2", "+CO2", LET$Treatment)
LET$Treatment <- stringr::str_replace_all(LET$Treatment, "Light", "-L")
LET$Treatment <- gsub("Nitrogen", "-N", LET$Treatment)
LET$Treatment <- gsub("Phosphate", "-P", LET$Treatment)
LET$Treatment <- gsub("Control", "Cont", LET$Treatment)


# Choose location of the letters on the graph
require(ggbreak)
coords_letter <- as.data.frame(cbind(values = (dat$div + dat$sd_div), Treatment = dat$Treatment)) #replace by rich, even or div


coords_letter <- merge(LET, coords_letter, by = "Treatment")
names(coords_letter)[2] <- "p_val"
coords_letter$values <- as.numeric(coords_letter$values)


## 3.6. Producing histograms of the different alpha diversity indices (factor for letter placement 2% for richness = 2. Keep the same for the rest)
# Figures 4.A to 4.C
img3 <- ggplot(dat, aes(x = Treatment, y = div)) + #replace by rich, even or div
  geom_bar(stat="identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = div-sd_div, ymax = div+sd_div), width=.2, #replace by rich, even or div
                position = position_dodge(.9)) +
  annotate(geom = "text", x = 4, y = 1.2*max(dat$div), label = paste("H =", test_results$statistic,"; p-val <", test_results$p.value), size = 14/.pt) + #replace by rich, even or div
  annotate(geom = "text", x = coords_letter$Treatment, y = coords_letter$values+0.3, label = coords_letter$p_val, parse = TRUE, size = 18/.pt) +
  # ylim(0,60) + # Only for richness
  # ylim(0,0.40) + # Only for Evenness
  ylim(0,4.5)+ # Only for Dispersion
  xlab("Treatment") +
  ylab("Functional Dispersion") + #replace by rich, even or div
  # ylim(c(0,100)) +
  scale_x_discrete(limits=rev) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 18)
  )

ggsave(file="Functional dispersion.svg", plot=img3, width=12, height=12, units = "cm") #replace by rich, even or div

## 3.7. Combine the 3 graphs into a single one (Figure 4.D)
alpha_diversity <- ggpubr::ggarrange(img1, img2, img3, ncol = 3, nrow = 1)
ggsave(file = "alpha_diversity.svg", plot = alpha_diversity, width = 26, height = 26, units = "cm")


## 4. Dissimilarity across functional space
require(hypervolume)

## 4.1. Subset main dataset by environmental condition; here example with nitrogen
temp <- reduced_dim %>%
  filter(Label == "-Nitrogen")

## 4.2. Perform the hypervolume on the PCA eigenvalues
Nitrogen <- hypervolume::hypervolume_gaussian(data = temp[, c(12:13)], kde.bandwidth = bw_estimate, quantile.requested = 0.95, quantile.requested.type = "probability")

save(resultats, file = "Whole Nitrogen hypervolume.RData")

rm(reduced_dim, resultats, temp, bw_estimate)

## 4.3. Perform the functional dissimilarities across environmental conditions
hv <- hypervolume_set(hv1 = Nitrogen, hv2 = CO2, verbose = TRUE, check.memory = FALSE) ### Run the comparison one by one, between each treatment. Step to automatise
hypervolume_overlap_statistics(hv)

rm(list = ls())

# -> Feed manually an excel sheet with the values

## 4.5. Read the excel sheet previously written 
dat <- readxl::read_xlsx("IndCyano_Louchart_Limitation_experiement_June2025.xlsx", sheet = "Jaccard similarity")

## 4.6. Conversion into a matrix
mat <- as.matrix(dat[,-1])
row.names(mat) <- colnames(mat)

## 4.7. Creation a correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

## 4.8. Get the upper matrix
get_upper_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

## 4.9. Reorder the matrix
mat <- reorder_cormat(mat)
upper_tri <- get_upper_tri(mat)

melt_mat <- melt(upper_tri, na.rm = TRUE)
melt_mat$value <- round(melt_mat$value, 2)

# Name the variable
melt_mat$Var1 <- gsub("CO2", "+CO2", melt_mat$Var1)
melt_mat$Var1 <- stringr::str_replace_all(melt_mat$Var1, "Light", "-Light")
melt_mat$Var1 <- gsub("Nitrogen", "-Nitrogen", melt_mat$Var1)
melt_mat$Var1 <- gsub("Phosphorus", "-Phosphorus", melt_mat$Var1)

melt_mat$Var2 <- gsub("CO2", "+CO2", melt_mat$Var2)
melt_mat$Var2 <- stringr::str_replace_all(melt_mat$Var2, "Light", "-Light")
melt_mat$Var2 <- gsub("Nitrogen", "-Nitrogen", melt_mat$Var2)
melt_mat$Var2 <- gsub("Phosphorus", "-Phosphorus", melt_mat$Var2)

## 4.10. Graphical output
# Figure 4.D
img <- 
  ggplot(data = melt_mat, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
  scale_fill_gradient2(low = "snow2", mid = "paleturquoise2", high = "paleturquoise3",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Jaccard \nSimilarity") +
  theme_bw() +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme(axis.text.x = element_text(size = 12, family = "Arial", hjust = 1, angle = 35),
        axis.text.y = element_text(size = 12, family = "Arial")) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10,
                               title.hjust = 0.5, title.vjust = 5))
img
ggsave(file="Jaccard similarity multitraits.svg", plot=img, width=15, height=14, units = "cm")
