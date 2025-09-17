.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages_updated",.libPaths()))

library(monocle3)
library(dplyr)
library(ggplot2)

key <- read.csv(file = "ECMkey.csv", row.names = NULL, check.names = FALSE)

#### Human Stricture Bowel ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/HumanStricPicro/Bowel")
quantified <- read.csv(file = "Quantified.csv", row.names = NULL, check.names = FALSE)
key <- read.csv(file = "ECMkey.csv", row.names = NULL, check.names = FALSE)
key <- as.data.frame(t(key))
key <- key[["V2"]]
key <- key[-1]
#colnames(quantified) <- paste0("Parameter.", seq_len(ncol(quantified)))
colnames(quantified) <- key
colnames(quantified) <- paste0(seq(1:294))
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
metadata <- metadata %>%
  mutate(Region = ifelse(Region == "ADJ IF", "ADJ Bowel", Region))
unique(metadata$Region)

# data <- t(as.matrix(quantified)) #set rownames as features ("genes") and colnames as images ("cells")
# umap <- metadata[,c("UMAP_1", "UMAP_2")]
# cds <- new_cell_data_set(data,
#                          cell_metadata = metadata)
# cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
# reducedDim(cds, type = "PCA") <- umap
# cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(quantified)
# #cds <- preprocess_cds(cds, num_dim = 100)
# #cds <- reduce_dimension(cds)
# #plot_pc_variance_explained(cds)

quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) + 
  scale_color_manual(values = c("#e99e4e", "blue", "#D7191C")) + 
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/BowelUMAPall.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  scale_color_manual(values = c("#e99e4e", "blue", "#D7191C")) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 1)
ggsave("Figures/BowelUMAP.svg")

cds <- cluster_cells(cds, resolution=1e-2)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

rownames(cds) <- paste0("Parameter.", seq_len(nrow(cds)))

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_humbowel.svg")
saveRDS(cds, "humanstricbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "Figures/humanstricture.csv")
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("#e99e4e", "blue", "#D7191C")) +
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12), legend.position = "none")
ggsave("Figures/picropseudohumbowel.svg")

#Matrix features 


quantified["pseudotime"] <- data.pseudo$monocle3_pseudotime
# Compute Pearson correlations between each parameter column and pseudotime
cor_matrix <- apply(quantified[, -which(names(quantified) == "pseudotime")], 2, 
                    function(x) cor(x, quantified$pseudotime, method = "pearson"))

# Convert to a dataframe for heatmap plotting
cor_df <- data.frame(Parameter = names(cor_matrix), Correlation = cor_matrix)
row.names(cor_df) <- cor_df$Parameter  # Set row names
cor_matrix <- as.matrix(cor_df[, -1])  # Convert to matrix for pheatmap
cor_matrix[!is.finite(cor_matrix)] <- 0  # Replace with 0 (or use `na.omit()` to remove rows)

# # Select the top 10 highest and lowest correlations
# top_10 <- cor_df[order(-cor_df$Correlation), ][1:10, ]   # Top 10 highest
# bottom_10 <- cor_df[order(cor_df$Correlation), ][1:10, ] # Top 10 lowest
# 
# 
#Barplot
# Select only numeric columns, excluding 'pseudotime' and 'pseudotime_group'
quantified_numeric <- quantified[, sapply(quantified, is.numeric)]
quantified_numeric <- quantified_numeric[, -which(names(quantified_numeric) == "pseudotime")]

# Categorize pseudotime into low or high based on the median (or you can use other criteria)
median_pseudo <- median(quantified$pseudotime)
quantified$pseudotime_group <- ifelse(quantified$pseudotime <= median_pseudo, "Low", "High")


# Calculate Pearson correlation for each parameter within "Low" and "High" pseudotime groups
cor_low <- apply(quantified_numeric[quantified$pseudotime_group == "Low", ], 2, 
                 function(x) cor(x, quantified$pseudotime[quantified$pseudotime_group == "Low"], method = "pearson"))

cor_high <- apply(quantified_numeric[quantified$pseudotime_group == "High", ], 2, 
                  function(x) cor(x, quantified$pseudotime[quantified$pseudotime_group == "High"], method = "pearson"))

cor_low <- replace(cor_low, is.na(cor_low), 0)
cor_high <- replace(cor_high, is.na(cor_high), 0)
# Combine the correlations for both low and high pseudotime
cor_df_low_high <- data.frame(
  Parameter = names(cor_low),
  Correlation_Low = cor_low,
  Correlation_High = cor_high
)

# Calculate the difference between high and low pseudotime correlations (High - Low)
cor_df_low_high$Correlation_Difference <- cor_df_low_high$Correlation_High - cor_df_low_high$Correlation_Low

# Optionally, remove parameters with NA correlations (i.e., missing data for certain groups)
cor_df_low_high <- na.omit(cor_df_low_high)

# Select the top 10 highest and lowest correlation differences
top_10 <- cor_df_low_high[order(-cor_df_low_high$Correlation_Difference), ][1:10, ]
bottom_10 <- cor_df_low_high[order(cor_df_low_high$Correlation_Difference), ][1:10, ]

# Combine top and bottom 10 for plotting
cor_subset <- rbind(top_10, bottom_10)

# # Plot the top 10 highest and lowest correlation differences
# ggplot(cor_subset, aes(x = reorder(Parameter, Correlation_Difference), y = Correlation_Difference)) +
#   geom_bar(stat = "identity", aes(fill = Correlation_Difference > 0)) +
#   scale_fill_manual(values = c("blue", "red")) +
#   theme_classic() +
#   labs(title = "Top 10 Highest and Lowest Correlation Differences (High - Low) with Pseudotime",
#        x = "Parameter",
#        y = "Difference in Correlation (High - Low)") +
#   theme(axis.text.x = element_text(angle = 90, size = 10))

library(ggplot2)
library(ggrepel)

# Subset for low and high pseudotime Pearson correlation values greater than 0.75
cor_df_filtered <- cor_df_low_high[cor_df_low_high$Correlation_Low > 0.5 & cor_df_low_high$Correlation_High > 0.5, ]

# Calculate correlation difference (High - Low)
cor_df_filtered$Correlation_Difference <- cor_df_filtered$Correlation_High - cor_df_filtered$Correlation_Low

# Calculate the 95% quantile for the correlation difference
quantile_95 <- quantile(cor_df_filtered$Correlation_Difference, 0.95)
quantile_5 <- quantile(cor_df_filtered$Correlation_Difference, 0.05)

# Subset the top 5 highest and lowest correlation differences
top_5_positive <- cor_df_filtered[order(-cor_df_filtered$Correlation_Difference), ][1:5, ]
top_5_negative <- cor_df_filtered[order(cor_df_filtered$Correlation_Difference), ][1:5, ]
top_bottom_5 <- rbind(top_5_positive, top_5_negative)

# Create the plot
ggplot(cor_df_filtered, aes(x = Correlation_Low, y = Correlation_High)) +
  geom_point(aes(size = abs(Correlation_High), color = Correlation_Difference)) +  # Color by Correlation_Difference
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color gradient for Correlation Difference
  scale_size_continuous(range = c(2, 8)) +  # Adjust size range for clarity
  
  # Add a dotted line at y = x
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", size = 1) +
  
  # Set axis limits from 0.5 to 0.8
  xlim(0.55, 0.8) +
  ylim(0.55, 0.8) +
  
  theme_classic() +
  labs(title = "Bubble Plot: Low vs High Pseudotime Correlations",
       x = "Correlation with Low Pseudotime",
       y = "Correlation with High Pseudotime",
       color = "Correlation Difference", size = "Correlation Strength") +
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 15)) +
  theme(legend.position = "right") +
  
  # Add text labels only for the top 5 highest and lowest correlation differences
  geom_text_repel(data = top_bottom_5, aes(label = Parameter), color = "black", size = 6, fontface = "bold", box.padding = 0.5, max.overlaps = 15)


##### Parameters in human ####
# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Add pseudotime to the quantified data
quantified["pseudotime"] <- data.pseudo$monocle3_pseudotime

# Compute Pearson correlations between each parameter column and pseudotime
cor_matrix <- apply(quantified[, -which(names(quantified) == "pseudotime")], 2, 
                    function(x) cor(x, quantified$pseudotime, method = "pearson"))

# Convert to a dataframe for heatmap plotting
cor_df <- data.frame(Parameter = names(cor_matrix), Correlation = cor_matrix)
row.names(cor_df) <- cor_df$Parameter  # Set row names
cor_matrix <- as.matrix(cor_df[, -1])  # Convert to matrix for pheatmap
cor_matrix[!is.finite(cor_matrix)] <- 0  # Replace with 0 (or use `na.omit()` to remove rows)

# Select the top 10 highest and lowest correlations
# top_10 <- cor_df[order(-cor_df$Correlation), ][1:10, ]   # Top 10 highest
# bottom_10 <- cor_df[order(cor_df$Correlation), ][1:10, ] # Top 10 lowest

# Barplot: Select only numeric columns, excluding 'pseudotime' and 'pseudotime_group'
quantified_numeric <- quantified[, sapply(quantified, is.numeric)]
quantified_numeric <- quantified_numeric[, -which(names(quantified_numeric) == "pseudotime")]

# Categorize pseudotime into low or high based on the lower and upper quartiles
lower_quartile <- quantile(quantified$pseudotime, 0.25)
upper_quartile <- quantile(quantified$pseudotime, 0.75)
quantified$pseudotime_group <- ifelse(quantified$pseudotime <= lower_quartile, "Low", 
                                      ifelse(quantified$pseudotime >= upper_quartile, "High", "Mid"))

# Calculate Pearson correlation for each parameter within "Low" and "High" pseudotime groups
cor_low <- apply(quantified_numeric[quantified$pseudotime_group == "Low", ], 2, 
                 function(x) cor(x, quantified$pseudotime[quantified$pseudotime_group == "Low"], method = "pearson"))

cor_high <- apply(quantified_numeric[quantified$pseudotime_group == "High", ], 2, 
                  function(x) cor(x, quantified$pseudotime[quantified$pseudotime_group == "High"], method = "pearson"))

# Replace NAs with 0 in correlation results
cor_low <- replace(cor_low, is.na(cor_low), 0)
cor_high <- replace(cor_high, is.na(cor_high), 0)

# Combine the correlations for both low and high pseudotime
cor_df_low_high <- data.frame(
  Parameter = names(cor_low),
  Correlation_Low = cor_low,
  Correlation_High = cor_high
)

# Calculate the difference between high and low pseudotime correlations (High - Low)
cor_df_low_high$Correlation_Difference <- cor_df_low_high$Correlation_High - cor_df_low_high$Correlation_Low

# Optionally, remove parameters with NA correlations (i.e., missing data for certain groups)
cor_df_low_high <- na.omit(cor_df_low_high)

# Subset for low and high pseudotime Pearson correlation values greater than 0.75
cor_df_filtered <- cor_df_low_high[cor_df_low_high$Correlation_Low > 0.5 & cor_df_low_high$Correlation_High > 0.5ggplot(cor_df_filtered, aes(x = Correlation_Low, y = Correlation_High)) +
  geom_point(aes(size = abs(Correlation_High), color = Correlation_Difference)) +  # Color by Correlation_Difference
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color gradient for Correlation Difference
  scale_size_continuous(range = c(2, 8)) +  # Adjust size range for clarity
  
  # Add a dotted line at y = x
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", size = 1) +
  
  # Set axis limits from 0.5 to 0.8
  xlim(0.5, 0.8) +
  ylim(0.5, 0.8) +
  
  theme_classic() +
  labs(title = "Bubble Plot: Low vs High Pseudotime Correlations",
       x = "Correlation with Low Pseudotime",
       y = "Correlation with High Pseudotime",
       color = "Correlation Difference", size = "Correlation Strength") +
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(legend.position = "right") +
  
  # Add text labels for all points
  geom_text_repel(aes(label = Parameter), color = "black", size = 6, fontface = "bold", box.padding = 0.5, max.overlaps = 10), ]

# Subset the top 5 highest and lowest correlation differences
top_5_positive <- cor_df_filtered[order(-cor_df_filtered$Correlation_Difference), ][1:5, ]
top_5_negative <- cor_df_filtered[order(cor_df_filtered$Correlation_Difference), ][1:5, ]
top_bottom_5 <- rbind(top_5_positive, top_5_negative)

# Add a new column to store the maximum correlation (either Low or High)
cor_df_filtered$Max_Correlation <- pmax(cor_df_filtered$Correlation_Low, cor_df_filtered$Correlation_High)

# Create the bubble plot with custom color palette and adjusted size
# Add a new column to store the maximum correlation (either Low or High)
cor_df_filtered$Max_Correlation <- pmax(cor_df_filtered$Correlation_Low, cor_df_filtered$Correlation_High)

# Define custom color palette
custom_colors <- c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#e2e2e2", "#e1a692", "#de6e56", "#e14b31", "#c23728")

# Split data into low and high pseudotime groups
cor_df_low <- cor_df_filtered[cor_df_filtered$Correlation_Low > 0, ]
cor_df_high <- cor_df_filtered[cor_df_filtered$Correlation_High > 0, ]

# Select the top 7 points based on Correlation_Low for low pseudotime and Correlation_High for high pseudotime
top_7_low <- cor_df_low[order(-cor_df_low$Correlation_Low), ][1:7, ]
top_7_high <- cor_df_high[order(-cor_df_high$Correlation_High), ][1:7, ]

# Combine the top 7 points from both groups
top_7_all <- rbind(top_7_low, top_7_high)

# Create the bubble plot with custom color palette and adjusted size
ggplot(cor_df_filtered, aes(x = Correlation_Low, y = Correlation_High)) +
  geom_point(aes(size = Max_Correlation, color = Correlation_Difference)) +  # Size based on max correlation
  scale_color_gradientn(
    colors = custom_colors,
    values = scales::rescale(c(min(cor_df_filtered$Correlation_Difference), 0, max(cor_df_filtered$Correlation_Difference)))
  ) +  # Custom color palette, rescaled for better distribution of colors
  scale_size_continuous(range = c(2, 8)) +  # Adjust bubble size based on max correlation
  
  # Add a dotted line at y = x
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", size = 1) +
  
  # Set axis limits from 0.5 to 0.8
  xlim(0.55, 0.8) +
  ylim(0.55, 0.8) +
  
  theme_classic() +
  labs(title = "Bubble Plot: Low vs High Pseudotime Correlations",
       x = "Correlation with Low Pseudotime",
       y = "Correlation with High Pseudotime",
       color = "Correlation Difference", size = "Correlation Strength") +
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15)) +
  theme(legend.position = "right") +
  
  # Add text labels for the top 7 points from both low and high pseudotime groups
  geom_text_repel(data = top_7_all, aes(label = Parameter), color = "black", size = 4,  box.padding = 0.5, max.overlaps = 20)

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

plot_parameters <- function(data, param_names) {
  # Convert to long format for faceting
  data_long <- data %>%
    select(pseudotime, all_of(param_names)) %>%
    pivot_longer(cols = -pseudotime, names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = param_names))  # Ensure correct order
  
  # Create scatter plots with a common pseudotime color scale
  ggplot(data_long, aes(x = pseudotime, y = Value, color = pseudotime)) +
    geom_point(alpha = 0.7, size = 2) +  # Points with transparency
    geom_smooth(method = "loess", color = "black", linetype = "dashed", se = FALSE) +  # Add smooth trend line
    scale_color_viridis(option = "magma") +  # Apply magma color palette
    theme_classic() +
    labs(
      x = "Pseudotime",
      y = "Parameter Value",
      color = "Pseudotime"
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.background = element_blank(),  # Remove boxes around facet titles
      strip.text = element_text(size = 14, face = "bold")  # Keep facet titles readable
    ) +
    facet_wrap(~ Parameter, scales = "free_y", ncol = 2)  # Ensures left-right, top-bottom order
}

# Example usage
param_vector <- c(names(quantified)[2], names(quantified)[128], 
                  names(quantified)[3], names(quantified)[19])  # Example: Parameter 2 and 128
plot_parameters(quantified, param_vector)


#### Timecourse Bowel ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/TimeCourse/Bowel/Subset")
quantified <- read.csv(file = "Quantified_updated.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP_updated.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

# #Change labels 
# metadata$Region[metadata$Region == "ADJ Bowel" ] <- "ADJ"
# metadata$Region[metadata$Region == "ADJ IF" ] <- "ADJ"
# metadata$Region[metadata$Region == "Sham Bowel" ] <- "Sham"
# metadata$Region[metadata$Region == "POD3 IF" ] <- "POD 3"
# metadata$Region[metadata$Region == "POD3 Bowel" ] <- "POD 3"
# metadata$Region[metadata$Region == "POD 7 IF" ] <- "POD 7"
# metadata$Region[metadata$Region == "POD 7 Bowel" ] <- "POD 7"
# metadata$Region[metadata$Region == "POD 14 IF" ] <- "POD 14"
# metadata$Region[metadata$Region == "POD 14 Bowel" ] <- "POD 14"
# metadata$Region[metadata$Region == "POD 30 IF" ] <- "POD 30"
# metadata$Region[metadata$Region == "POD 30 Bowel" ] <- "POD 30"
# metadata$Region[metadata$Region == "POD 90 IF" ] <- "POD 90"
# metadata$Region[metadata$Region == "POD 90 Bowel" ] <- "POD 90"

bowel <- grep("_B_", rownames(metadata))
metadata$subregion <- NA
metadata$subregion[bowel] <- "Bowel"
bowel <- grep("Sham ", rownames(metadata))
metadata$subregion[bowel] <- "Bowel"

IF <- grep("_IF_", rownames(metadata))
metadata$subregion[IF] <- "IF"

sham <- grep("ShamBowel", rownames(metadata))
metadata$subregion[sham] <- "Bowel"

unique(metadata$Region)
unique(metadata$subregion)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ ", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("Sham", "ADJ", "POD3_COLO", "POD7_COLO", "POD14_COLO", "POD30_COLO", "POD90_COLO"))
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10)) +   theme(axis.line = element_line(linewidth = 10)) +   
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

ggsave("Figures/BowelUMAPall.pdf")

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + theme(axis.line = element_line(linewidth = 10)) +   
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) + 
  facet_wrap(~Region, nrow = 2)
ggsave("Figures/BowelUMAPall_split.pdf")

cds <- cluster_cells(cds, resolution=1e-2)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)
ggsave("Figures/TCMATUMAPall.svg")

#Plot by subregion 
groups <- c("POD7_COLO", "POD14_COLO", "POD30_COLO", "POD90_COLO")
pattern <- paste(groups, collapse = "|")
cols <- grep(pattern, colData(cds)$Region)
cds_sub <- cds[, cols]

cds_sub <- cluster_cells(cds_sub, resolution=1e-2)
plot_cells(cds_sub, color_cells_by = "subregion", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10)) +   
  theme(axis.line = element_line(linewidth = 15), 
   axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
    axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/BowelIFonly.pdf") 

plot_cells(cds_sub, color_cells_by = "subregion", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~subregion, ncol = 2) + 
  theme(axis.line = element_line(linewidth = 10)) +   
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) + 
  scale_color_manual(values = c("#0b5394", "#cc0000"))
ggsave("Figures/BowelIFonlysplit.pdf") 

#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0) + 
  theme(axis.line = element_line(linewidth = 10)) +   
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/m3_TCmat.pdf")
saveRDS(cds, "TCbowelmonocle3_updated.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "TCbowel_updated.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")

ggsave("Figures/picropseudoTCbp.svg")

#Pseudotime Plots 
#cds_sub <- learn_graph(cds_sub)
#cds_sub <- order_cells(cds_sub)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_TCbowelsub.pdf")
saveRDS(cds_sub, "TCbowelmonocle3sub_updated.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds_sub))
data.pseudo$monocle3_pseudotime <- pseudotime(cds_sub)
write.csv(data.pseudo, "TCbowelsub_updated.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = subregion)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")

ggsave("Figures/picropseudoTCbp.svg")


#### Timecourse MAT ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/TimeCourse/MAT")
quantified <- read.csv(file = "MAT_TC_Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
rows_with_pod30 <- grep("POD30", rownames(metadata))
metadata$Region[rows_with_pod30] <- "POD 30"
rows_with_pod30 <- grep("POD30_ADJ", rownames(metadata))
metadata$Region[rows_with_pod30] <- "ADJ"
metadata$Region[metadata$Region == "ADJ POD3" ] <- "ADJ"
metadata$Region[metadata$Region == "ADJ POD7" ] <- "ADJ"
metadata$Region[metadata$Region == "ADJ POD14" ] <- "ADJ"
metadata$Region[metadata$Region == "ADJ POD30" ] <- "ADJ"
metadata$Region[metadata$Region == "ADJ POD90" ] <- "ADJ"
metadata$Region[metadata$Region == "COLO POD3" ] <- "POD 3"
metadata$Region[metadata$Region == "COLO POD7" ] <- "POD 7"
metadata$Region[metadata$Region == "COLO POD14" ] <- "POD 14"
metadata$Region[metadata$Region == "COLO POD90" ] <- "POD 90"

unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ", "POD 3", "POD 7", "POD 14", "POD 30", "POD 90"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/TCMATUMAPall.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2)
ggsave("Figures/MATUMAPsplit.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_TCmat.svg")
saveRDS(cds, "TCbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")

#### Stretch MAT ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/StretchCF")
quantified <- read.csv(file = "CF_Stretch_Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "CF_Stretch_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
rows_with_pod30 <- grep("US_COLO", rownames(metadata))
metadata$Region[rows_with_pod30] <- "US Colo"
metadata$Region[metadata$Region == "Stretched Adjacent" ] <- "S ADJ"
metadata$Region[metadata$Region == "Unstretched Adjacent" ] <- "US ADJ"
metadata$Region[metadata$Region == "Stretched Colotomy" ] <- "S Colo"
unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("US ADJ", "S ADJ", "US Colo", "S Colo"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/MATstretchedumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2)
ggsave("Figures/MATUMAPsplit.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_stretchmat.svg")
saveRDS(cds, "stretchmatmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "stretchmatm3.csv")

#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudomatbp.svg")


#### Stretch Bowel ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/VehStretch_POD14/Subset")
quantified <- read.csv(file = "Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
metadata$Region[metadata$Region == "Stretched Adjacent" ] <- "S ADJ"
metadata$Region[metadata$Region == "Unstretched Adjacent" ] <- "US ADJ"
metadata$Region[metadata$Region == "Unstretched" ] <- "US Colo"
metadata$Region[metadata$Region == "Stretched" ] <- "S Colo"
unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("US ADJ", "S ADJ", "US Colo", "S Colo"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/IFstretchedumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))

ggsave("Figures/IFstretchsplit.svg")
#7.33 x 5.91

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/m3_stretchbowel.svg")
saveRDS(cds, "stretchbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "stretchbowelm3.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudobowelbp.svg")


#### DSS Bowel ####
setwd("/oak/stanford/groups/longaker/KEBR/Picro/DSS/Subset")
quantified <- read.csv(file = "Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
metadata$Region[metadata$Region == "Adjacent DSS" ] <- "ADJ DSS"
metadata$Region[metadata$Region == "Sham DSS" ] <- "Sham DSS"
metadata$Region[metadata$Region == "Adjacent Vehicle" ] <- "ADJ Veh"
metadata$Region[metadata$Region == "Colotomy DSS" ] <- "Colo DSS"
metadata$Region[metadata$Region == "Colotomy Vehicle" ] <- "Colo Veh"
unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Veh", "ADJ DSS", "Sham DSS", "Colo Veh", "Colo DSS"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10)) +
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/DSSumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2) +
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/DSSUMAPsplit.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0) +
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave("Figures/m3_DSSbowel.svg")
saveRDS(cds, "DSSbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "dssbpseudo.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudobowelbp.svg")
#Run 2 way Anova: 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
data.pseudo$injury <- "NA"
data.pseudo$injury[data.pseudo$Region == "ADJ DSS" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Sham DSS" ] <- "Sham"
data.pseudo$injury[data.pseudo$Region == "ADJ Veh" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Colo DSS" ] <- "Colo"
data.pseudo$injury[data.pseudo$Region == "Colo Veh" ] <- "Colo"
unique(data.pseudo$injury)

data.pseudo$treatment <- "NA"
data.pseudo$treatment[data.pseudo$Region == "ADJ DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Sham DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "ADJ Veh" ] <- "Veh"
data.pseudo$treatment[data.pseudo$Region == "Colo DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Colo Veh" ] <- "Veh"
unique(data.pseudo$treatment)
anova_result <- aov(monocle3_pseudotime ~ injury * treatment, data = data.pseudo)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
significant_pairings <- as.data.frame(tukey_result$injury)
significant_pairings <- significant_pairings[significant_pairings$`p adj` < 0.05, ]
significant_pairings <- significant_pairings[, c("diff", "p adj")]

# Add the significant pairs to a list for ggsignif
pair_list <- rownames(significant_pairings)
pair_list <- strsplit(pair_list, "-")
pair_list <- pair_list[!sapply(pair_list_region, function(x) any(is.na(x)))]





ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none") +
  geom_signif(
    comparisons = pair_list,
    map_signif_level = TRUE
  )

#### DSS MAT ####
setwd("/oak/stanford/groups/longaker/KEBR/Picro/DSS_MAT")
quantified <- read.csv(file = "CF_DSS_Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "CF_DSS_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
metadata$Region[metadata$Region == "DSS Sham" ] <- "Sham DSS"
metadata$Region[metadata$Region == "DSS Colo" ] <- "Colo DSS"
metadata$Region[metadata$Region == "Vehicle Colo" ] <- "Colo Veh"
metadata$Region[metadata$Region == "Vehicle Adj" ] <- "ADJ Veh"
metadata$Region[metadata$Region == "DSS Adj" ] <- "ADJ DSS"
unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Veh", "ADJ DSS", "Sham DSS", "Colo Veh", "Colo DSS"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/DSSMATumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2)
ggsave("Figures/DSSMATUMAPsplit.svg")

cds <- cluster_cells(cds, resolution=1e-2)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_DSSmat.svg")
saveRDS(cds, "DSSbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudomatbp.svg")
#Run 2 way Anova: 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
data.pseudo$injury <- "NA"
data.pseudo$injury[data.pseudo$Region == "ADJ DSS" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Sham DSS" ] <- "Sham"
data.pseudo$injury[data.pseudo$Region == "ADJ Veh" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Colo DSS" ] <- "Colo"
data.pseudo$injury[data.pseudo$Region == "Colo Veh" ] <- "Colo"
unique(data.pseudo$injury)

data.pseudo$treatment <- "NA"
data.pseudo$treatment[data.pseudo$Region == "ADJ DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Sham DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "ADJ Veh" ] <- "Veh"
data.pseudo$treatment[data.pseudo$Region == "Colo DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Colo Veh" ] <- "Veh"
unique(data.pseudo$treatment)
anova_result <- aov(monocle3_pseudotime ~ injury * treatment, data = data.pseudo)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result


significant_pairings <- as.data.frame(tukey_result$injury)
significant_pairings <- significant_pairings[significant_pairings$`p adj` < 0.05, ]
significant_pairings <- significant_pairings[, c("diff", "p adj")]

# Add the significant pairs to a list for ggsignif
pair_list <- rownames(significant_pairings)
pair_list <- strsplit(pair_list, "-")
pair_list <- pair_list[!sapply(pair_list_region, function(x) any(is.na(x)))]





ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none") +
  geom_signif(
    comparisons = pair_list,
    map_signif_level = TRUE
  )




#### KO Bowel ####
setwd("/oak/stanford/groups/longaker/KEBR/Picro/YAP KO IF/Subset")
quantified <- read.csv(file = "Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels
metadata$Region[metadata$Region == "YAPHet_ADJ" ] <- "YAP Het ADJ"
metadata$Region[metadata$Region == "YAPHet" ] <- "YAP Het Colo"
metadata$Region[metadata$Region == "YAPKO_ADJ" ] <- "YAP KO ADJ"
metadata$Region[metadata$Region == "YAPKO" ] <- "YAP KO Colo"
unique(metadata$Region)

quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("YAP Het ADJ", "YAP KO ADJ", "YAP Het Colo", "YAP KO Colo"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10)) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/Hydrogelumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/Hydrogelsplit.svg")

cds <- cluster_cells(cds, resolution=1e-1)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE, 
           show_trajectory_graph = F) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/m3_YAPHydrogel.svg")
saveRDS(cds, "m3_YAPHydrogel.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "bowelKOpseudo.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudobowelYAPKObp.svg")
#Run 2 way Anova: 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
data.pseudo$injury <- "NA"
data.pseudo$injury[data.pseudo$Region == "ADJ DSS" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Sham DSS" ] <- "Sham"
data.pseudo$injury[data.pseudo$Region == "ADJ Veh" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Colo DSS" ] <- "Colo"
data.pseudo$injury[data.pseudo$Region == "Colo Veh" ] <- "Colo"
unique(data.pseudo$injury)

data.pseudo$treatment <- "NA"
data.pseudo$treatment[data.pseudo$Region == "ADJ DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Sham DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "ADJ Veh" ] <- "Veh"
data.pseudo$treatment[data.pseudo$Region == "Colo DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Colo Veh" ] <- "Veh"
unique(data.pseudo$treatment)
anova_result <- aov(monocle3_pseudotime ~ injury * treatment, data = data.pseudo)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result


significant_pairings <- as.data.frame(tukey_result$injury)
significant_pairings <- significant_pairings[significant_pairings$`p adj` < 0.05, ]
significant_pairings <- significant_pairings[, c("diff", "p adj")]

# Add the significant pairs to a list for ggsignif
pair_list <- rownames(significant_pairings)
pair_list <- strsplit(pair_list, "-")
pair_list <- pair_list[!sapply(pair_list_region, function(x) any(is.na(x)))]





ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none") +
  geom_signif(
    comparisons = pair_list,
    map_signif_level = TRUE
  )


#### KO MAT ####
setwd("/oak/stanford/groups/longaker/KEBR/Picro/YAPKO/MAT3/")
quantified <- read.csv(file = "Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
# metadata$Region[metadata$Region == "DSS Sham" ] <- "Sham DSS"
# metadata$Region[metadata$Region == "DSS Colo" ] <- "Colo DSS"
# metadata$Region[metadata$Region == "Vehicle Colo" ] <- "Colo Veh"
# metadata$Region[metadata$Region == "Vehicle Adj" ] <- "ADJ Veh"
# metadata$Region[metadata$Region == "DSS Adj" ] <- "ADJ DSS"
# unique(metadata$Region)
quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
#metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("HET ADJ MAT", "KO ADJ MAT", "HET COLO MAT", "KO COLO MAT"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/KOMATumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2)
ggsave("Figures/KOMATUMAPsplit.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_KOmat.svg")
saveRDS(cds, "KOMATmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudomatbp.svg")
#Run 2 way Anova: 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
data.pseudo$injury <- "NA"
data.pseudo$injury[data.pseudo$Region == "ADJ DSS" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Sham DSS" ] <- "Sham"
data.pseudo$injury[data.pseudo$Region == "ADJ Veh" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Colo DSS" ] <- "Colo"
data.pseudo$injury[data.pseudo$Region == "Colo Veh" ] <- "Colo"
unique(data.pseudo$injury)

data.pseudo$treatment <- "NA"
data.pseudo$treatment[data.pseudo$Region == "ADJ DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Sham DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "ADJ Veh" ] <- "Veh"
data.pseudo$treatment[data.pseudo$Region == "Colo DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Colo Veh" ] <- "Veh"
unique(data.pseudo$treatment)
anova_result <- aov(monocle3_pseudotime ~ injury * treatment, data = data.pseudo)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result


significant_pairings <- as.data.frame(tukey_result$injury)
significant_pairings <- significant_pairings[significant_pairings$`p adj` < 0.05, ]
significant_pairings <- significant_pairings[, c("diff", "p adj")]

# Add the significant pairs to a list for ggsignif
pair_list <- rownames(significant_pairings)
pair_list <- strsplit(pair_list, "-")
pair_list <- pair_list[!sapply(pair_list_region, function(x) any(is.na(x)))]





ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none") +
  geom_signif(
    comparisons = pair_list,
    map_signif_level = TRUE
  )


#### KO Hydrogel ####
setwd("/oak/stanford/groups/longaker/KEBR/Picro/Hydrogels/Subset")
quantified <- read.csv(file = "Bowel_TC_Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "Bowel_TC_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels
metadata$Region[metadata$Region == "Vert_ADJ" ] <- "Vert ADJ"
metadata$Region[metadata$Region == "Empty_ADJ" ] <- "Empty ADJ"
metadata$Region[metadata$Region == "Vert Colo" ] <- "Vert Colo"
metadata$Region[metadata$Region == "Empty Colo" ] <- "Empty Colo"
unique(metadata$Region)

quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap
#colData(cds)$Region <- factor(colData(cds)$Region, levels = c("ADJ Bowel", "POD3 Bowel", "POD3 IF", "POD 7 Bowel", "POD 7 IF", 
#                                                              "POD 14 Bowel", "POD 14 IF", "POD 30 Bowel", "POD 30 IF", "POD 90 Bowel", "POD 90 IF"))
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("Empty ADJ", "Vert ADJ", "Empty Colo", "Vert Colo"))

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) +
  theme(axis.line = element_line(linewidth = 10)) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/Hydrogelumap.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 2) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/Hydrogelsplit.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE, 
           show_trajectory_graph = F) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))

# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0) + 
  theme(axis.line = element_line(linewidth = 15), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), 
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
ggsave("Figures/m3_hydrogelbowel.svg")
saveRDS(cds, "m3_hydrogelbowel.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "hydrogelpseudo.csv")
#data.pseudo.sub <- data.pseudo[, data.pseudo$Region = c("ADJ Bowel", "ADJ IF")]
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none")
ggsave("Figures/pseudobowelhydrogelbp.svg")
#Run 2 way Anova: 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
data.pseudo$injury <- "NA"
data.pseudo$injury[data.pseudo$Region == "ADJ DSS" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Sham DSS" ] <- "Sham"
data.pseudo$injury[data.pseudo$Region == "ADJ Veh" ] <- "ADJ"
data.pseudo$injury[data.pseudo$Region == "Colo DSS" ] <- "Colo"
data.pseudo$injury[data.pseudo$Region == "Colo Veh" ] <- "Colo"
unique(data.pseudo$injury)

data.pseudo$treatment <- "NA"
data.pseudo$treatment[data.pseudo$Region == "ADJ DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Sham DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "ADJ Veh" ] <- "Veh"
data.pseudo$treatment[data.pseudo$Region == "Colo DSS" ] <- "DSS"
data.pseudo$treatment[data.pseudo$Region == "Colo Veh" ] <- "Veh"
unique(data.pseudo$treatment)
anova_result <- aov(monocle3_pseudotime ~ injury * treatment, data = data.pseudo)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
tukey_result


significant_pairings <- as.data.frame(tukey_result$injury)
significant_pairings <- significant_pairings[significant_pairings$`p adj` < 0.05, ]
significant_pairings <- significant_pairings[, c("diff", "p adj")]

# Add the significant pairs to a list for ggsignif
pair_list <- rownames(significant_pairings)
pair_list <- strsplit(pair_list, "-")
pair_list <- pair_list[!sapply(pair_list_region, function(x) any(is.na(x)))]





ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(
    axis.text.x = element_text(size = 10),  # Right align and center x-axis labels angle = 90, hjust = 1, vjust = 0.5
    axis.text.y = element_text(size = 10),                         # Adjust y-axis text
    axis.title.y = element_text(size = 12), 
    legend.position = "none") +
  geom_signif(
    comparisons = pair_list,
    map_signif_level = TRUE
  )


#### Human Visium H&E ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/Viz_H&E_Integration/Image_Files/369x277 Tiles_trimed/")
quantified <- read.csv(file = "Quantified.csv", row.names = NULL, check.names = FALSE)
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

quantified <- t(quantified)
cds <- new_cell_data_set(as(quantified,"dgCMatrix"), cell_metadata = metadata)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

#Add in UMAP
umap <- metadata[,c("UMAP_1", "UMAP_2")]
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap #add umap embedding to cds object
reducedDim(cds, type = "PCA") <- umap

plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8, show_trajectory_graph = FALSE) + 
  scale_color_manual(values = c("#e99e4e", "blue", "#D7191C")) + 
  theme(axis.line = element_line(linewidth = 10))
ggsave("Figures/BowelUMAPall.svg")
plot_cells(cds, color_cells_by = "Region", 
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           graph_label_size=1.5,
           label_cell_groups=FALSE,
           cell_size = 0.8) + 
  scale_color_manual(values = c("#e99e4e", "blue", "#D7191C")) + 
  theme(axis.line = element_line(linewidth = 10)) + 
  facet_wrap(~Region, nrow = 1)
ggsave("Figures/BowelUMAP.svg")

cds <- cluster_cells(cds, resolution=1e-2)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)


# #Find features in each cluster 
# marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group = 100,  speedglm.maxiter = 1*10^100)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
#Pseudotime Plots 
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
saveRDS(cds, "vizHEmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "Figures/humanstricture.csv")
ggplot(data.pseudo, aes(x = Region, y = monocle3_pseudotime, fill = Region)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("#e99e4e", "blue", "#D7191C")) +
  labs(y = "Pseudotime", x = NULL) +    # Change y-axis label
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12), legend.position = "none")
ggsave("Figures/picropseudohumbowel.svg")