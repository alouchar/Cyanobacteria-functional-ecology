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

#### PRODUCTION OF GRAPHICAL AND NUMERICAL OUTPUTS AS WELL AS INTERMEDIATE DATASETS
# Compute average values by environmental condition and replicate
DAT <- dat %>%
  group_by(Treatment, Replicate) %>%
  summarise_all(mean)

# Convert the data as a matrix
num_mat <- as.matrix(sapply(DAT[,c(3:10)], as.numeric))

# Recode Label and filename
DAT$Treatment <- recode(DAT$Treatment, "-Light" = "-L", "-Nitrogen" = "-N","-Phosphorus" = "-P", "Control"  = "Control")
DAT$Replicate <- recode(DAT$Replicate, "Replicate 1" = "R1", "Replicate 2" = "R2","Replicate 3" = "R3", "Replicate 4"  = "R4")

# Aggregate columns "Environmental condition" + "Replicate" with a space in between
DAT$name <- apply(DAT[,c("Treatment","Replicate")], 1, paste, collapse = " ")
NOM <- as.vector(DAT[,11])

# Display row names of the matrix as "environmental condition" + "Replicate"
rownames(num_mat) <- NOM[["name"]]

# Producing the graphical outputs (figure 2)
heatmap <- pheatmap(num_mat, cutree_rows = 5, cutree_cols = 4, fontsize = 16, angle_col = "45", color = colorRampPalette(c("darkblue", "white", "darkred"))(200))

svg(filename = "Heatmap_filters_treatments.svg",width=7, height=7, family = "sans", pointsize = 12)
heatmap
dev.off()
