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
  geom_point(PCA$li, mapping = aes(x = PCA$li$Axis1, y = -PCA$li$Axis2, color = as.factor(dat$Treatment)), alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(PCA$co, mapping = aes(x = 0, y = 0, xend = PCA$co$Comp1*3, yend = -PCA$co$Comp2*3), col = "grey20",arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(size = 4, PCA$co, mapping = aes(x = PCA$co$Comp1*4, y = -PCA$co$Comp2*4, label = rownames(PCA$co)), col = "grey20") +
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
p
ggsave(file="PCA_density.svg", plot=pca_marg, width=16, height=16, units = "cm", dpi = 400)

rm(heatmap, NOM, num_mat)

## 2.8. PCA eigenvalues dataframe binded with environmental conditions and replicate
reduced_dim <- as.data.frame(cbind(dat, PCA$li$Axis1, -PCA$li$Axis2, PCA$li$Axis3, PCA$li$Axis4))

names(reduced_dim)[11] <- "PCA1"
names(reduced_dim)[12] <- "PCA2"
names(reduced_dim)[13] <- "PCA3"
names(reduced_dim)[14] <- "PCA4"

rm(list = setdiff(ls(),"reduced_dim"))


## 2.9. Delineation of multidimensional functional space

# The methodology relies on the probabilistic hypervolume since it considers abundance thus being less sensitive to outliers
bw_estimate <- hypervolume::estimate_bandwidth(reduced_dim[, c("PCA1","PCA2")], method = "cross-validation")


## 3. Computation of functional diversity indices based on functional hypervolumes
## Steps 3.1 and 3.2 are computationally intensive and may require several days per treatment depending on data size.

## 3.1. Hypervolume computation parallelized using foreach on multiple CPU cores

# List of environmental treatments to be analysed
treatments <- unique(reduced_dim$Treatment)

# Adjust the number of CPU cores according to available resources
num_cores <- 2

# Create and register a parallel cluster
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

compute_hypervolume_treatment <- function(tr, data, bw) {
  
  message("Processing treatment: ", tr)
  
  # Output file used to store intermediate results
  outfile <- paste0("hypervolume_", tr, ".RData")
  
  # Hypervolumes are computed only if they do not already exist
  if (!file.exists(outfile)) {
    
    # Subset data for the selected environmental treatment
    temp <- data %>% filter(Treatment == tr)
    
    # Identify unique biological replicates
    replicates <- unique(temp$Replicate)
    
    # Initialize a progress bar
    pb <- txtProgressBar(max = length(replicates), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Parallel loop: one hypervolume per replicate
    hv_list <- foreach(
      i = seq_along(replicates),
      .combine = c,
      .options.snow = opts,
      .packages = "hypervolume"
    ) %dopar% {
      
      DATA <- temp[temp$Replicate == replicates[i], ]
      
      # Construction of probabilistic Gaussian hypervolumes
      hypervolume::hypervolume_gaussian(
        data = DATA[, c("PCA1", "PCA2")],
        kde.bandwidth = bw,
        quantile.requested = 0.95,
        quantile.requested.type = "probability"
      )
    }
    
    # Close the progress bar
    close(pb)
    
    # Merge replicate-level hypervolumes into a treatment-level hypervolume
    hv <- hypervolume_join(hv_list)
    
    # Save hypervolume object to avoid recomputation
    save(hv, file = outfile)
    
  } else {
    # Load previously computed hypervolume
    load(outfile)
  }
  
  ## 3.2. Computation of functional diversity indices from the treatment-level hypervolume
  data.frame(
    Treatment = tr,
    rich = kernel.alpha(hv),        # Functional richness
    even = kernel.evenness(hv),     # Functional evenness
    div  = kernel.dispersion(hv)    # Functional dispersion
  )
}

reduced_sample <- reduced_dim


## 3.3. Run the analysis for all environmental treatments
Functional_rep_space_microcystis <- do.call(
  rbind,
  lapply(
    treatments,
    compute_hypervolume_treatment,
    data = reduced_dim,
    bw   = bw_estimate
  )
)

stopCluster(cl)


metrics <- c("rich", "even", "div")  # metrics to compute

## 3.4. Kruskal–Wallis tests
kw_results <- lapply(metrics, function(met){
  
  tmp <- Functional_rep_space_microcystis %>%
    select(Treatment, value = !!sym(met)) %>%
    filter(!is.na(value)) %>%
    group_by(Treatment) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  if (length(unique(tmp$Treatment)) < 2) {
    return(data.frame(
      metric = met,
      H = NA,
      p_value = NA
    ))
  }
  
  test <- kruskal.test(x = tmp$value, g = tmp$Treatment)
  
  data.frame(
    metric = met,
    H = round(as.numeric(test$statistic), 2),
    p_value = round(test$p.value, 3)
  )
})

kw_results <- do.call(rbind, kw_results)


## 3.5. Conover-lman tests
conover_results <- lapply(metrics, function(met){
  
  tmp <- Functional_rep_space_microcystis %>%
    select(Treatment, value = !!sym(met)) %>%
    filter(!is.na(value)) %>%
    group_by(Treatment) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  if (length(unique(tmp$Treatment)) < 2) return(NULL)
  
  conover.test(
    x = tmp$value,
    g = tmp$Treatment,
    kw = FALSE,
    method = "bh"
  )
})

names(conover_results) <- metrics

## 3.6. Pairwise significance letters
get_letters <- function(conover_obj){
  
  v <- data.frame(
    comparisons = conover_obj$comparisons,
    p_val = conover_obj$P.adjusted * 2
  )
  
  v[c("T1", "T2")] <- str_split_fixed(v$comparisons, " - ", 2)
  
  nameVals <- sort(unique(c(v$T1, v$T2)))
  mat <- matrix(1, length(nameVals), length(nameVals),
                dimnames = list(nameVals, nameVals))
  
  mat[as.matrix(v[, c("T1","T2")])] <- v$p_val
  mat <- as.matrix(forceSymmetric(mat, "U"))
  
  let <- multcompLetters(
    mat,
    compare = "<",
    threshold = 0.025,
    Letters = letters
  )
  
  data.frame(
    Treatment = names(let$Letters),
    Letters   = let$Letters
  )
}

letters_list <- lapply(conover_results, function(x){
  if (is.null(x)) return(NULL)
  get_letters(x)
})


## 3.7. Function to compute summary stats
data_summary <- function(data, varnames, groupnames){
  
  summary_list <- lapply(varnames, function(varname){
    
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm = TRUE),
        sd   = sd(x[[col]], na.rm = TRUE))
    }
    
    data_sum <- plyr::ddply(data, groupnames, .fun = summary_func, varname)
    data_sum$metric <- varname   # store metric name
    return(data_sum)
  })
  
  # Combine all metrics
  summary_all <- do.call(rbind, summary_list)
  
  return(summary_all)
}

# Produce the table
alpha_func <- data_summary(
  data = Functional_rep_space_microcystis,
  varnames = metrics,
  groupnames = "Treatment"
)


## 3.8. Produce the plots
plot_alpha <- function(metric_name, y_label, letters_df) {
  
  # Subset summary for this metric
  dat <- alpha_func[alpha_func$metric == metric_name, ]
  dat$Treatment <- as.character(dat$Treatment)
  
  Treatment <- c("Control", "+CO2", "-Phosphorus", "-Nitrogen", "-Light")
  
  # Merge letters
  coords <- merge(dat, letters_df, by = "Treatment", all.x = TRUE)
  coords$y <- 1.05 * coords$mean + coords$sd   # adjust height for letters
  
  # Assign colors for axis labels using ggsci palette
  treatment_colors <- ggsci::pal_d3("category20")(length(Treatment))
  names(treatment_colors) <- levels(dat$Treatment)
  
  # Replace +CO2 with proper subscript for axis labels
  axis_labels <- Treatment
  axis_labels[axis_labels == "+CO2"] <- expression("+CO"[2])
  
  # Base bar plot
  p <- ggplot(dat, aes(x = Treatment, y = mean, fill = Treatment)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_fill_manual(values = treatment_colors) +
    scale_x_discrete(labels = axis_labels, limits=rev) +
    ylab(y_label) +
    xlab("") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 14),
      legend.position = "none"
    )
  
  # Add Conover letters
  p <- p + geom_text(
    data = coords,
    aes(x = Treatment, y = y, label = Letters),
    inherit.aes = FALSE,
    size = 5
  )
  
  # Add Kruskal-Wallis annotation if available
  if (exists("kw_results")) {
    kw <- kw_results[kw_results$metric == metric_name, ]
    if (!is.na(kw$H)) {
      p <- p +
        annotate("text",
                 x = 1,
                 y = max(dat$mean + dat$sd) * 1.25,
                 label = paste("H =", kw$H, "; p =", kw$p_value),
                 hjust = 0,
                 size = 4)
    }
  }
  
  return(p)
}

img_rich <- plot_alpha("rich", "Functional richness", letters_df = letters_list[["rich"]])
img_even <- plot_alpha("even", "Functional evenness", letters_df = letters_list[["even"]])
img_div  <- plot_alpha("div",  "Functional divergence", letters_df = letters_list[["div"]])

img_rich

## 3.9. Dissimilarities between treatments
treatments <- unique(reduced_dim$Treatment)

hv_list <- setNames(
  lapply(treatments, function(tr) {
    temp <- reduced_dim %>% filter(Treatment == tr)
    hypervolume::hypervolume_gaussian(temp[, c("PCA1", "PCA2")],
                                      kde.bandwidth = bw_estimate,
                                      quantile.requested = 0.95,
                                      quantile.requested.type = "probability")
  }),
  treatments
)

pair_combinations <- combn(names(hv_list), 2, simplify=FALSE)

# compute pairwise overlap
overlap_results <- lapply(pair_combinations, function(pair) {
  
  hv1 <- hv_list[[pair[1]]]
  hv2 <- hv_list[[pair[2]]]
  
  # Compute hypervolume set
  hv_set <- hypervolume_set(hv1, hv2, verbose = FALSE, check.memory = FALSE)
  
  # Compute overlap statistics
  stats <- hypervolume_overlap_statistics(hv_set)
  
  # Extraction as numeric
  data.frame(
    Treatment1 = pair[1],
    Treatment2 = pair[2],
    Jaccard_similarity = as.numeric(stats["jaccard"])
  )
})

# Combine all results
overlap_table <- bind_rows(overlap_results)


# Define treatments and order
treatments <- c("-Light", "-Nitrogen", "-Phosphorus", "+CO2", "Control")

# Convert overlap_table to square matrix
mat <- matrix(NA, nrow = length(treatments), ncol = length(treatments),
              dimnames = list(treatments, treatments))

for(i in 1:nrow(overlap_table)){
  t1 <- overlap_table$Treatment1[i]
  t2 <- overlap_table$Treatment2[i]
  mat[t1, t2] <- overlap_table$Jaccard_similarity[i]
  mat[t2, t1] <- overlap_table$Jaccard_similarity[i]
}

diag(mat) <- 1

# Keep upper triangle
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

upper_tri <- get_upper_tri(mat)

# Melt for ggplot
melt_mat <- melt(upper_tri, na.rm = TRUE)
melt_mat$value <- round(melt_mat$value, 2)

# Factor levels for plotting
melt_mat$Var1 <- factor(melt_mat$Var1, levels = treatments)
melt_mat$Var2 <- factor(melt_mat$Var2, levels = treatments)

# Assign colors for axis labels using ggsci palette
axis_colors <- ggsci::pal_d3("category20")(length(treatments))
names(axis_colors) <- treatments

# Replace +CO2 with proper subscript for axis labels
axis_labels <- treatments
axis_labels[axis_labels == "+CO2"] <- expression("+CO"[2])

# Plot
img_jaccard <- ggplot(melt_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), color = "black", size = 5) +   # tile values
  scale_fill_gradient(low = "grey90", high = "grey40", name = "Jaccard\nSimilarity") +
  scale_x_discrete(labels = axis_labels) +
  scale_y_discrete(labels = axis_labels) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(color = axis_colors[treatments], size = 12, family = "Arial", hjust = 1, angle = 35),
    axis.text.y = element_text(color = axis_colors[treatments], size = 12, family = "Arial"),
    legend.position = "none",
    plot.margin = unit(c(0.3,0.3,0.1,0.1), "cm")
  )

# Figure and save
bot <- plot_grid(img_div, img_jaccard, labels = c("C", "D"), label_size = 14)
top <- plot_grid(img_rich, img_even, labels = c("A", "B"), label_size = 14)
p <- plot_grid(top, bot, ncol =1)
p

ggsave(file="Functional diversity.svg", plot=p, width=20, height=20, units = "cm")

