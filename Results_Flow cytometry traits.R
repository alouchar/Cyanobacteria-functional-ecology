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

#### PRODUCTION OF GRAPHICAL AND NUMERICAL OUTPUTS AS WELL AS INTERMEDIATE DATASETS
## 1. Flow cytometry traits
## 1.1. Correlation matrix between traits to identify collinearity (supplementary material figure 2)
res <- cor(DATA[,c(1:8)])
p.mat <- cor.mtest(DATA[,c(1:8)])

svg(filename = "Trait correlation matrix.svg", width=6, height=6, pointsize = 12, family = "sans")
corrplot(res, method = "shade", shade.col = NA,
         addCoef.col = "grey50", type = "upper",
         cl.pos = "n", order = "hclust", 
         p.mat = p.mat, sig.level = 0.001,
         tl.col = "black", tl.srt = 45)
dev.off()

rm(res, p.mat)

## 1.2. Heatmap of cyanobacterial traits across environmental conditions (Figure 2)
# Compute average values by environmental condition and replicate
dat <- DATA %>%
  group_by(Label, filename) %>%
  summarise_all(mean)

# Convert the data as a matrix
num_mat <- as.matrix(sapply(dat[,c(3:10)], as.numeric))

# Recode Label and filename
dat$Label <- recode(dat$Label, "-Light" = "-L", "-Nitrogen" = "-N","-Phosphorus" = "-P", "Control"  = "Cont")
dat$filename <- recode(dat$filename, "Replicate 1" = "R1", "Replicate 2" = "R2","Replicate 3" = "R3", "Replicate 4"  = "R4")

# Aggregate columns "Environmental condition" + "Replicate" with a space in between
dat$name <- apply(dat[,c("Label","filename")], 1, paste, collapse = " ")
NOM <- as.vector(dat[,11])

# Display row names of the matrix as "environmental condition" + "Replicate"
rownames(num_mat) <- NOM[["name"]]

# Producing the graphical outputs
heatmap <- pheatmap(num_mat, cutree_rows = 5, cutree_cols = 4, fontsize = 16, angle_col = "45", color = colorRampPalette(c("darkblue", "white", "darkred"))(200))

svg(filename = "Heatmap_filters_treatments.svg",width=7, height=7, family = "sans", pointsize = 12)
heatmap
dev.off()
