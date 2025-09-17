setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales); library(readxl); library(RColorBrewer); library(cowplot)
library(ggpubr)
options(Seurat.object.assay.version = "v5")

#### Global cells ####
global <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Annotated_Global.rds")
Idents(global) <- "subpop"

## Course annotation plots 
levels <- c("Fibroblasts", "Pericytes", "SMCs", "Vascular ECs", "Lymphatic ECs", "Glia Cells", 
            "T Cells", "B Cells", "Plasma Cells", "Macrophages", "Mast Cells", "Proliferating Immune Cells", 
            "Epithelial Cells")
global$subpop <- factor(global$subpop, levels = newlevels)



##Fine annotation plots 

newlevels <- c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib", "GREM1+ ssFib",  "KCNN3+ ssFib" ,  "PI16+ ssFib", "CD74+ apFib","F3+ ssFib", "ADAMDEC1+ ssFib", 
               "LUM+ Fibroblast-Like SMCs", "PRKG1+ ssSMCs" , "JUNBhi Stress-Responsive SMCs"  , 
               "EBF1+ Steady-State Peri", "JUNhi Activated Peri" , 
               "NOTCH4hi Arterial ECs", "CA4hi Capillary ECs", "ACKR1+ Venous ECs", "CXCL12hi Microvascular ECs", 
               "STAB1hi Medullary Sinus LECs", "CLDN11+ Valve LECs", "CCL21hi Paracortical Sinus LECs"  , "NTShi Ceiling LECs", "CXCL3hi Floor LECs" ,   
               "NRXN1hi ssGlia", "SOCS3hi Submucosal Glia" , "HLA-DRAhi Activated Glia" , 
               "LYVE1+ Macrophages", "C1QChi Macrophages", "CXCL2hi Macrophages", "CD300E+ Monocytes", "CPA3+ Mast Cells", "CSF3R+ Neutrophils", "IL4I1 Activated DC", "GZMB+ pDC", "LAMP3+ mregDC", "CLEC9A+ cDC1" , "CD1C+ cDC2"  ,
               "CD8A+ CD8 T Cells", "CCR6+ TH17 T Cells", "CCR7+ TH1 T Cells"  , "IL7R+ ILC" , "FOX3P+ Tregs", "AFF3+ Naïve B Cells"  , "TNFRSF13B+ Memory B Cells" , "IGHD+ Follicular B Cells"  , "AICDA+ GC B Cells"  , "FCER1G+ NK Cells"  , "IGHG1+ Plasma Cells" ,"IGHA1+ Plasma Cells", "MKI67+ PI Cells" , 
               "OLFM4+ ISCs" , "MKI67+ TA Cells", "DEFA5+ Paneth Cells", "MUC2+ Goblet Cells"   , "CHGB+ Enteroendocrine",  "FABP2+ Enterocytes", "BEST4+ Enterocytes", "CEACAM7+ Colonocytes" ,    "SH2D6+ Tuft Cells"  
               )
   

global$subpop <- factor(global$subpop, levels = newlevels)

Idents(global) <- "subpop"


library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)

# Generate a high-contrast color palette with 59 distinct colors
set1_colors <- brewer.pal(9, "Set1")
set3_colors <- brewer.pal(12, "Set3")
colors <- c(colorRampPalette(set1_colors)(30), colorRampPalette(set3_colors)(29))

# Assuming `global` is your Seurat object containing the UMAP results
umap_data <- as.data.frame(Embeddings(global, "full.umap.harmony"))
umap_data$cluster <- global$subpop

# Rename UMAP columns to match your dataframe
colnames(umap_data)[1:2] <- c("fullumapharmony_1", "fullumapharmony_2")
# Calculate cluster centers using median
cluster_centers <- umap_data %>%
  group_by(cluster_num) %>%
  summarize(
    fullumapharmony_1 = median(fullumapharmony_1),
    fullumapharmony_2 = median(fullumapharmony_2),
    .groups = 'drop'
  )

# Number of clusters for reference
num_clusters <- nrow(cluster_centers)

umap_plot <- DimPlot(global, group.by = "subpop", cols = colors, raster = TRUE, pt.size = 1, alpha = 0.2) +
  geom_text_repel(
    data = cluster_centers,
    aes(x = fullumapharmony_1, y = fullumapharmony_2, label = cluster_num),
    size = 3.5,                                  # Reduced label size
    label.padding = 0.1,                         # Padding around label text
    box.padding = 0.5,                           # Reduced padding around labels
    point.padding = 0.3,                         # Reduced space between point and text
    min.segment.length = 0,                      # Ensure segments if necessary
    max.overlaps = num_clusters + 10,            # Allow labels to overlap within reason
    force = 30,                                  # Increase force for separation
    force_pull = 3,                              # Stronger pull force to keep labels close
    direction = "both",                          # Allow both X and Y adjustments
    seed = 42,                                   # Ensure reproducibility
    segment.color = "grey40",
    fill = alpha("white", 0.8)
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Display the plot
print(umap_plot)
ggsave("Figures/globalumapupdated.pdf")

#generate legend
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(gridExtra)

# Calculate cluster order and labels
num_clusters <- length(levels(umap_data$cluster_num))
clusters_ordered <- paste0("C", 0:(num_clusters-1))

# Create legend data frame
legend_data <- data.frame(
  cluster_num = factor(clusters_ordered, levels = clusters_ordered),
  label = as.character(levels(factor(umap_data$cluster))),
  color = colors
)

# Calculate column layout
num_cols <- 5
rows_per_col <- ceiling(num_clusters / num_cols)

# Create column groups with reversed row order
legend_data <- legend_data %>%
  mutate(
    column_group = (row_number() - 1) %/% rows_per_col,
    row_in_column = (rows_per_col - (row_number() - 1) %% rows_per_col - 1)
  ) %>%
  group_by(column_group) %>%
  mutate(y = row_number()) %>%
  ungroup()

# Calculate dynamic column positions based on label widths
column_groups <- sort(unique(legend_data$column_group))
x_positions <- numeric(length(column_groups))
current_x <- 0
char_width <- 0.13  # Width per character in plot units
column_padding <- 0.3  # Minimum space between columns
label_offset <- 0.2   # Space between dot and text

for(i in seq_along(column_groups)) {
  group <- column_groups[i]
  group_labels <- legend_data$label[legend_data$column_group == group]
  
  # Calculate max width for current column
  max_width <- max(nchar(group_labels)) * char_width
  
  # Set position for current column
  x_positions[i] <- current_x
  
  # Calculate starting position for next column (current position + content width + padding)
  if(i < length(column_groups)) {
    current_x <- current_x + max_width + label_offset + column_padding
  }
}

# Create position mapping and join to legend data
position_map <- data.frame(
  column_group = column_groups,
  x_position = x_positions
)

legend_data <- legend_data %>%
  left_join(position_map, by = "column_group") %>%
  mutate(
    x = x_position,
    y = y * 0.08  # Vertical row spacing
  )

# Calculate plot boundaries
max_x <- max(legend_data$x + nchar(legend_data$label) * char_width + label_offset)
max_y <- max(legend_data$y)

# Create legend plot
custom_legend <- ggplot(legend_data, aes(x = x, y = y)) +
  geom_point(
    shape = 21,
    size = 5,
    aes(fill = color)
  ) +
  geom_text(
    aes(label = cluster_num),
    size = 2.2,
    hjust = 0.5,
    vjust = 0.5
  ) +
  geom_text(
    aes(
      x = x + label_offset,
      label = label
    ),
    size = 3,
    hjust = 0,
    color = "black"
  ) +
  scale_fill_identity() +
  scale_y_reverse() +
  theme_void() +
  coord_cartesian(
    xlim = c(-0.2, max_x),
    ylim = c(max_y + 0.05, -0.05),
    clip = "off"
  ) +
  theme(plot.margin = margin(10, 10, 10, 10, "mm"))

print(custom_legend)
ggsave("custom_cluster_legend.pdf")

#Reassignment of global cluster names
global$manual.annotation[global$manual.annotation == "Proliferating Immune Cells"] <- "PI Cells"
global$manual.annotation[global$manual.annotation == "Macrophages"] <- "Myeloid Cells"

##UMAPS across conditions


##Heatmap of marker genes 
Idents(global) <- "manual.annotation"
levels <- c("Fibroblasts", "Pericytes", "SMCs", "Vascular ECs", 
            "Lymphatic ECs", "Glia Cells", "T Cells", "B Cells", 
            "Plasma Cells", "Myeloid Cells", "Mast Cells", "PI Cells", 
            "Epithelial Cells")
global$manual.annotation <- factor(global$manual.annotation, levels = levels)

Idents(global) <- "manual.annotation"

# Subset the data
set.seed(123)  # Ensure reproducibility
sub <- subset(global, downsample = 1000)
sub[["RNA"]]$data <- as(object = sub[["RNA"]]$data, Class = "dgCMatrix")

# Find all markers with specified parameters
markers <- FindAllMarkers(sub, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

# Select the top 100 differentially expressed genes based on average log FC for each cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()

# Subset the data for the top 100 genes
expr_data <- GetAssayData(sub, slot = "data")[top100_genes, ]

# Calculate mean expression for each gene in each cluster
cluster_means <- AverageExpression(sub, features = top100_genes, group.by = "manual.annotation")$RNA

# Reorder the rows of cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    cluster_expression <- matrix[level, , drop=FALSE]
    ordered_genes <- c(ordered_genes, colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)])
  }
  return(ordered_genes)
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Define the row splitting
row_labels <- rownames(cluster_means_z_t)
row_split <- factor(row_labels, levels = levels)

# Plot the heatmap
ht <- Heatmap(cluster_means_z_t, 
              name = "Z-score",
              col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
              show_row_names = TRUE,
              show_column_names = FALSE, 
              cluster_rows = FALSE,  # Maintain predefined row order
              cluster_columns = FALSE,  # Do not cluster columns 
              row_title = "Clusters",  # Label rows 
              column_title = "Top 100 Genes",  # Label columns 
              #row_split = row_split,  # Split rows by factor clusters
              border = TRUE
              #row_gap = unit(2, "mm")  # Fixed gap between split rows
)
pdf(file="Figures/globalfibrohm.pdf")
draw(ht)
dev.off()

# Reorder the rows of cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    cluster_expression <- matrix[level, , drop=FALSE]
    ordered_genes <- c(ordered_genes, colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)])
  }
  return(ordered_genes)
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Define the row splitting
row_labels <- rownames(cluster_means_z_t)
row_split <- factor(row_labels, levels = levels)

# Plot the heatmap
ht <- Heatmap(cluster_means_z_t, 
              name = "Z-score",
              col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
              show_row_names = TRUE,
              show_column_names = FALSE, 
              cluster_rows = TRUE,  # Cluster rows
              cluster_columns = FALSE,  # Do not cluster columns 
              row_title = "Clusters",  # Label rows 
              column_title = "Top 100 Genes",  # Label columns 
              row_split = row_split,  # Split rows by factor clusters
              row_gap = unit(0, "mm"),  # Remove gaps between rows
              border = FALSE)  # Remove outer border

# Draw the heatmap
draw(ht, heatmap_legend_side = "right")

pdf(file="Figures/globalfibrohmclustered.pdf", height = 2, width = 4.5)
draw(ht)
dev.off()

# Adjusting the Seurat FeaturePlot to fit the desired specifications
features = c("PDGFRA", "ACTA2",  "RGS5","PECAM1", "PROX1", "PLP1",  "CD68", "CPA3", "CD3E", "MS4A1", "IGHA1", "MKI67", "EPCAM")

FeaturePlot(global, features = features, order = TRUE, raster = TRUE, max.cutoff = "q95", ncol = 5, pt.size = 3, alpha = 0.2) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1) # Ensure square aspect ratio

ggsave("Figures/global_fplots.pdf")


#### Global Fibroblasts DimPlots ####
fibro <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Fibro4/union_fibro4_harmony_clusters_0.3.rds")
fibro$diseasestate[fibro$study == "LongakerBowel" & fibro$diseasestate == "Involved"] <- "Stricture"
fibro <- readRDS("globalfibroannotated_0.3.rds")
Idents(fibro) <- "seurat_clusters"
fibro$fibro_annotation <- "tbd"
fibro$fibro_annotation[fibro$seurat_clusters == 0] <- "ADAMDEC1+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 1] <- "THBS1+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 2] <- "PI16+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 3] <- "F3+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 4] <- "GREM1+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 5] <- "KCNN3+ ssFib"
fibro$fibro_annotation[fibro$seurat_clusters == 6] <- "CTHRC1+ mFib"
fibro$fibro_annotation[fibro$seurat_clusters == 7] <- "LRRC7+ mFib"
fibro$fibro_annotation[fibro$seurat_clusters == 8] <- "MMP3+ iFib"
fibro$fibro_annotation[fibro$seurat_clusters == 9] <- "CD74+ apFib"
fibro$fibro_annotation[fibro$seurat_clusters == 10] <- "ADAMDEC1+ ssFib"
Idents(fibro) <- "fibro_annotation"
levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
                    "GREM1+ ssFib","KCNN3+ ssFib",
                    "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib" , 
                    "THBS1+ ssFib", "PI16+ ssFib")
fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = levels)
Idents(fibro) <- "fibro_annotation"
palette = rev(c(
  "#7b3294",  # deeper purple
  "#9467bd",  # original purple
  "#1f77b4",  # original blue
  "#17becf",  # original cyan
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff9900",  # more vivid orange
  "#ff7f0e",  # original orange
  "#db4d6d",  # new warm rose/red (fills the gap)
  #  "#e377c2",  # original pink
  "#8b0000"   # original dark red
))
DimPlot(fibro, cols = palette, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )
ggsave("Figures/globalumapupdatedgeneralDS.pdf") #5.2 x 3.15

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, 
                                   margin = margin(t = 5, unit = "pt")), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


StackedVlnPlot(fibro, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, cols = palette) +
  theme( axis.text.x = element_text(
  size = 15,
  angle = 45,      # Rotate 45°
  hjust = 1,       # Right-align labels with hashes
  vjust = 1,       # Baseline-align labels with hashes
  margin = margin(t = 5, unit = "pt")  # Small buffer above labels
))
ggsave("Figures/QCmetrics.pdf")


library(Seurat)
library(RColorBrewer)
library(scales)

Idents(fibro) <- "patientID"
# Assume `seurat_object` is your Seurat object
number_of_clusters <- length(unique(Idents(fibro)))

# Generate a color palette based on the number of clusters
generate_colors <- function(n) {
  set1_colors <- brewer.pal(9, "Set1")
  set3_colors <- brewer.pal(12, "Set3")
  combined_colors <- c(colorRampPalette(set1_colors)(30), colorRampPalette(set3_colors)(29))
  
  if(n > length(combined_colors)) {
    extended_colors <- c(combined_colors, colorRampPalette(brewer.pal(8, "Dark2"))(30), colorRampPalette(brewer.pal(11, "Spectral"))(31))
  } else {
    extended_colors <- combined_colors
  }
  
  if(n > length(extended_colors)) {
    extended_colors <- c(extended_colors, colorRampPalette(extended_colors)(n - length(extended_colors)))
  }
  
  return(extended_colors[1:n])
}

# Preview the generated palette
colors <- generate_colors(number_of_clusters)

DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none",  # Removes the legend
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_bypatientID.png")

Idents(fibro) <- "orig.ident"
# Assume `seurat_object` is your Seurat object
number_of_clusters <- length(unique(Idents(fibro)))

# Generate a color palette based on the number of clusters
generate_colors <- function(n) {
  set1_colors <- brewer.pal(9, "Set1")
  set3_colors <- brewer.pal(12, "Set3")
  combined_colors <- c(colorRampPalette(set1_colors)(30), colorRampPalette(set3_colors)(29))
  
  if(n > length(combined_colors)) {
    extended_colors <- c(combined_colors, colorRampPalette(brewer.pal(8, "Dark2"))(30), colorRampPalette(brewer.pal(11, "Spectral"))(31))
  } else {
    extended_colors <- combined_colors
  }
  
  if(n > length(extended_colors)) {
    extended_colors <- c(extended_colors, colorRampPalette(extended_colors)(n - length(extended_colors)))
  }
  
  return(extended_colors[1:n])
}

# Preview the generated palette
colors <- generate_colors(number_of_clusters)

DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none",  # Removes the legend
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_byorigident.jpg")

Idents(fibro) <- "study"
# Assume `seurat_object` is your Seurat object
number_of_clusters <- length(unique(Idents(fibro)))

# Generate a color palette based on the number of clusters
generate_colors <- function(n) {
  set1_colors <- brewer.pal(9, "Set1")
  set3_colors <- brewer.pal(12, "Set3")
  combined_colors <- c(colorRampPalette(set1_colors)(30), colorRampPalette(set3_colors)(29))
  
  if(n > length(combined_colors)) {
    extended_colors <- c(combined_colors, colorRampPalette(brewer.pal(8, "Dark2"))(30), colorRampPalette(brewer.pal(11, "Spectral"))(31))
  } else {
    extended_colors <- combined_colors
  }
  
  if(n > length(extended_colors)) {
    extended_colors <- c(extended_colors, colorRampPalette(extended_colors)(n - length(extended_colors)))
  }
  
  return(extended_colors[1:n])
}

# Preview the generated palette
colors <- generate_colors(number_of_clusters)

DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_bystudy.jpg")

Idents(fibro) <- "platform"
# Assume `seurat_object` is your Seurat object
number_of_clusters <- length(unique(Idents(fibro)))

generate_dynamic_colors <- function(n) {
  # Basic color palettes from RColorBrewer
  set1_colors <- brewer.pal(9, "Set1")
  paired_colors <- brewer.pal(12, "Paired")
  dark2_colors <- brewer.pal(8, "Dark2")
  spectral_colors <- brewer.pal(11, "Spectral")
  
  # If few levels, use a diverse and high-contrast palette
  if (n <= 9) {
    colors <- set1_colors[1:n]
  } else if (n <= 12) {
    colors <- paired_colors[1:n]
  } else if (n <= 20) {
    base_colors <- c(set1_colors, dark2_colors)
    colors <- colorRampPalette(base_colors)(n)
  } else {
    # For even larger number of levels, combine and extend palettes for variety
    base_colors <- c(set1_colors, paired_colors, dark2_colors, spectral_colors)
    colors <- colorRampPalette(base_colors)(n)
    
    # Ensure we have enough distinct colors
    if (n > length(base_colors)) {
      extra_colors <- colorRampPalette(base_colors)(n - length(base_colors))
      colors <- c(base_colors, extra_colors)
    }
  }
  
  return(colors)
}

DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_byplatform.jpg")

Idents(fibro) <- "tissue"
number_of_clusters <- length(unique(Idents(fibro)))
colors <- generate_dynamic_colors(number_of_clusters)
show_col(colors)
DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_bytissue.jpg")

Idents(fibro) <- "specimen_type"
number_of_clusters <- length(unique(Idents(fibro)))
colors <- generate_dynamic_colors(number_of_clusters)
show_col(colors)
DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_byspetype.jpg")

Idents(fibro) <- "treated"
treated <- subset(fibro, treated == "NA", invert = T)
number_of_clusters <- length(unique(Idents(treated)))
colors <- generate_dynamic_colors(number_of_clusters)
show_col(colors)
DimPlot(treated, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_bytreatment.jpg")

Idents(fibro) <- "diseasestate"
number_of_clusters <- length(unique(Idents(fibro)))
colors <- generate_dynamic_colors(number_of_clusters)
show_col(colors)
DimPlot(fibro, cols = colors, raster = F, pt.size = 0, alpha = 0.2) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave("Figures/fibro_dstate.jpg")



#### Global Fibroblast heatmap ####
Idents(fibro) <- "fibro_annotation"

# Subset the data
set.seed(123)  # Ensure reproducibility
sub <- subset(fibro, downsample = 1000)
sub[["RNA"]]$data <- as(object = sub[["RNA"]]$data, Class = "dgCMatrix")

# Find all markers with specified parameters
markers <- FindAllMarkers(sub, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
                    "GREM1+ ssFib","KCNN3+ ssFib",
                    "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib" , 
                    "THBS1+ ssFib", "PI16+ ssFib")
# Select the top 100 differentially expressed genes based on average log FC for each cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()

# Subset the expression data for the top 100 genes
expr_data <- GetAssayData(sub, slot = "data")[top100_genes, ]

# Calculate mean expression for each gene in each cluster
cluster_means <- AverageExpression(sub, features = top100_genes, group.by = "fibro_annotation")$RNA

# Reorder the rows of the cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    if (level %in% rownames(matrix)) {
      cluster_expression <- matrix[level, , drop = FALSE]
      ordered_genes_level <- colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)]
      ordered_genes <- c(ordered_genes, ordered_genes_level)
    } else {
      warning(paste("Level", level, "not found in the matrix row names."))
    }
  }
  return(unique(ordered_genes))
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z_t, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Transpose the matrix to get genes as rows and clusters as columns
heatmap_matrix <- t(cluster_means_z_t)

# Define genes of interest
genes_of_interest <- c("CTHRC1", "POSTN", "FAP", 
                       "CD74", "CCL19", "HLA-DRB1", 
                       "ADAMDEC1", "HAPLN1", "EDIL3", 
                       "F3", "FRZB", "BMP4", 
                       "GREM1", "MGP", "SFRP2", 
                       "KCNN3", "SPARCL1", "THBS4", 
                       "LRRC7", "AUTS2", "SDK1", 
                       "MMP3", "CHI3L1", "INHBA", 
                       "PI16", "SEMA3C", "CD55", 
                       "SOD2", "THBS1", "FMO2"
                      )

ha <- rowAnnotation(
  gene_labels = anno_mark(
    at = which(rownames(heatmap_matrix) %in% genes_of_interest),
    labels = rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% genes_of_interest],
    labels_gp = gpar(fontsize = 8)  # Set fontsize within anno_mark
  )
)

# Define colour mapping
col_fun <- colorRamp2(c(-2, 0, 2), brewer.pal(11, "RdYlBu")[c(11, 6, 1)])

# Plot transposed heatmap
ht <- Heatmap(
  heatmap_matrix,
  name = "Z-score",
  col = col_fun,
  show_row_names = FALSE,        # Hide gene names (we'll use annotations)
  show_column_names = TRUE,      # Show cluster names
  column_names_side = "bottom",  # Cluster names below columns
  cluster_rows = TRUE,           # Preserve gene order
  cluster_columns = TRUE,        # Preserve cluster order
  row_title = "Genes",
  column_title = "Clusters",
  right_annotation = ha,         # Gene annotations on right
  border = TRUE, 
  show_column_dend = TRUE,
  show_row_dend = FALSE
)

# Draw the heatmap
pdf(file="Figures/globalfibrohm.pdf", height = 7, width = 4)
draw(ht)
dev.off()

####
#### FeaturePlots of fibro markers ####
require(pals)
require(reshape2)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Adjusting the Seurat FeaturePlot to fit the desired specifications
features = c("CTHRC1", "LRRC7", "MMP3", "GREM1", "KCNN3", "CD74", "F3", "ADAMDEC1", "THBS1", "PI16")

FeaturePlot(fibro, features = features, order = TRUE, raster = TRUE, max.cutoff = "q95", ncol = 5) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1) # Ensure square aspect ratio
    
ggsave("Figures/glob_fibro_fplots.pds")

#Pathogenic fibroblast identification 


# Printing the plot
print(p)
#### Main signature characteristics ####
Yapgenes <- c('ACAT1','ADAM9','CCN1','KLF5','ACTA2','COL4A1','LIF','LACTB',
              'AMOTL2','IL11','SORCS3','ANKRD1','CCN2','C11orf96','LBH','ANLN',
              'HBEGF','EDN1','MAP2K6','AURKB','HSPG2','ATF3','ARHGAP42','BIRC5',
              'LAMC1','F3','ADRA1B','CASP3','LPL','SERPINE1','CRIM1','CASP8','IER3',
              'C18orf63','CAVIN2','THBS1','TRIB1','CDH2','CCNF','TIMP2','FAM178A','CDC42',
              'VCAM1','KLF10','SLC38A2','CDK1','CALM2','SMAD7','ABI3BP','CDK6','COL1A1','FZD8',
              'TBX18','COL1A2','RUNX1','DYNLRB2','COL5A2','FGF2','NRD1','DIAPH3','COL6A2',
              'ZRANB1','DICER1','CXCL12','TNS3','E2F1','FBN1','MUSK','INHBA','USP44','EGFR',
              'LAMA2','PMS1','FGF1','LAMB1','PRELID2','FN1','TFPI','JDP2','FOXM1','VEGFA','AL035681.1',
              'GAB1','BDNF','BASP1','IGF1R','CCL2','QSOX1','LATS2','COL4A2','DDHD1','LDHA','CSF1','CITED2',
              'LIFR','SEMA3E','DUSP10','MCM6','ANXA1','OR51B6','MYC','BGN','EFHC1','MYL9','COL11A1','MAGI2',
              'NDUFB3','COL14A1','JARID2','NEGR1','COL5A1','CCDC132','NEK2','COL6A1','ARNT2','NF2','ANGPT1',
              'PELO','PLK1','ADAM12','CENPV','PMAIP1','CCL9','RPL22L1','POLA1','PTN','SYNDIG1','PPP1R3B','SEMA3A',
              'RAD18','PTEN','THBS2','TSPAN3','RAD51','HMGB1','TBC1D2','SLC2A3','IRS1','SNAI2','NUAK2','TAGLN','TRPC3',
              'TEAD4','GAS6','VIM','SNTB1','WWC1','MTMR10','ZEB1','BMP5','FRY','IGFBP3','AC112715.2','SUSD1','RUNX2',
              'ROBO2','AC023590.1','FBXO47','FAM107A','RBMS3','SH3RF1','SYT16','XPNPEP1','PIK3C2G',
              'MED13L','RND3','CDH4','TSPAN18','NEDD9','KPNA3','SH2D4A','TRIO','DPYD','LRRTM1','ZNF469',
              'RP11-553A10.1','GLIS3','SLITRK6','IQCJ-SCHIP1','FGF5','NTRK3','ADAMTS8','NAALADL2',
              'PTGER4','ANXA3','ERCC2','FAM172A','GNB1L','PSMG2','NMBR','SLC38A4','ROR1','PARD3B','BET1', 
              "FOXF2", "CCDC80", "MYOF","FJX1","PTPN14","ARHGEF17", "DOCK5","NT5E", "ASAP1", "GADD45A", "TGFB2", "AXL")
Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
revfibro <- AddModuleScore(revfibro, features = list(Yapsub), name = "Yapsig") #Get score for YAP targets
#2. YAP components signature
components <- list(c("YAP1", "CTNNB1", "MST1",
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))

#Here we calculate AddModuleScores to characterize global fibroblasts 
migration <- c("CTHRC1", "TNC", "CEMIP", "CDH2", "CD44", "ARID5A", "DDR2", "IQGAP1", "SLC8A1", "TNS1", "FGF2", "PLAU", "PRKCE", "PODXL", "ID1", 
               "CNN1", "SNAI2", "ACTC1", "MEOX2", "BMP2", "LUM", "FGF18", 
               "PIK3CA", "EDN1", "GRB2", "THBS1", "ITGB1", "ITGA3", 
               "ZEB1", "ZEB2", "TWIST1", "TWIST2", 
               "ARID5B", "SGPL1", "AKT1", "TNS1", "ZFAND5", "CCN3", "MACIR", "ILK1", 
               "TMEM201", "FUT8", "NHERF1", "PIP5K1A", "IQGAP1", "PDLIM1", "FAM114A1", "PLEC", "AQP1", "PML", "CD248", "LAMTOR2") #combo of "Fibroblast Migration GO + known migration genes

inflam <- c("ABCA1","ABI1","ACVR1B","ACVR2A","ADGRE1","ADM","ADORA2B","ADRM1","AHR","APLNR","AQP9",
           "ATP2A2","ATP2B1","ATP2C1","AXL","BDKRB1","BEST1","BST2","BTG2","C3AR1","C5AR1","CALCRL",
           "CCL17","CCL2","CCL20","CCL22","CCL24","CCL5","CCL7","CCR7","CCRL2","CD14","CD40","CD48",
           "CD55","CD69","CD70","CD82","CDKN1A","CHST2","CLEC5A","CMKLR1","CSF1","CSF3","CSF3R","CX3CL1",
           "CXCL10","CXCL11","CXCL6","CXCL8","CXCL9","CXCR6","CYBB","DCBLD2","EBI3","EDN1","EIF2AK2","EMP3",
           "EREG","F3","FFAR2","FPR1","FZD5","GABBR1","GCH1","GNA15","GNAI3","GP1BA","GPC3","GPR132","GPR183",
           "HAS2","HBEGF","HIF1A","HPN","HRH1","ICAM1","ICAM4","ICOSLG","IFITM1","IFNAR1","IFNGR2","IL10",
           "IL10RA","IL12B","IL15","IL15RA","IL18","IL18R1","IL18RAP","IL1A","IL1B","IL1R1","IL2RB","IL4R",
           "IL6","IL7R","INHBA","IRAK2","IRF1","IRF7","ITGA5","ITGB3","ITGB8","KCNA3","KCNJ2","KCNMB2","KIF1B",
           "KLF6","LAMP3","LCK","LCP2","LDLR","LIF","LPAR1","LTA","LY6E","LYN","MARCO","MEFV","MEP1A","MET",
           "MMP14","MSR1","MXD1","MYC","NAMPT","NDP","NFKB1","NFKBIA","NLRP3","NMI","NMUR1","NOD2","NPFFR2",
           "OLR1","OPRK1","OSM","OSMR","P2RX4","P2RX7","P2RY2","PCDH7","PDE4B","PDPN","PIK3R5","PLAUR","PROK2",
           "PSEN1","PTAFR","PTGER2","PTGER4","PTGIR","PTPRE","PVR","RAF1","RASGRP1","RELA","RGS1","RGS16","RHOG",
           "RIPK2","RNF144B","ROS1","RTP4","SCARF1","SCN1B","SELE","SELENOS","SELL","SEMA4D","SERPINE1","SGMS2",
           "SLAMF1","SLC11A2","SLC1A2","SLC28A2","SLC31A1","SLC31A2","SLC4A4","SLC7A1","SLC7A2","SPHK1","SRI",
           "STAB1","TACR1","TACR3","TAPBP","TIMP1","TLR1","TLR2","TLR3","TNFAIP6","TNFRSF1B","TNFRSF9","TNFSF10",
           "TNFSF15","TNFSF9","TPBG","VIP")

"TNFRSF1B", "CCL2", "IL7R", "CXCL9", 

# inflam <- c("PLAU", "CHI3L1", "MMP3", "IL1R1", "IL13RA2", "TNFSF11", "MMP10", "OSMR", "IL11", "STRA6", "FAP", 
#             "WNT2", "TWIST1", "IL24")

#"GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"
antigen <- c("ABCB9","ACE","AP3B1","AP3D1","ARL8B","ATG5","B2M","CALR","CCL19","CCL21","CCR7","CD1A","CD1B","CD1C",
             "CD1D","CD1E","CD209","CD68","CD74","CD8A","CLEC4A","CLEC4M","CTSD","CTSE","CTSF","CTSH","CTSL","CTSS",
             "CTSV","DNM2","ERAP1","ERAP2","EXT1","FCER1G","FCER2","FCGR1A","FCGR2B","FGL2","GBA1","HFE","HLA-A","HLA-B",
             "HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
             "HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F","HLA-G","HLA-H","ICAM1",
             "IDE","IFI30","IKBKB","KDM5D","LGMN","LILRB2","LNPEP","MARCHF1","MARCHF8","MFSD6","MICA","MICB","MPEG1",
             "MR1","NOD1","NOD2","PDIA3","PIKFYVE","PSMB8","PSME1","PYCARD","RAB10","RAB27A","RAB32","RAB33A","RAB34",
             "RAB35","RAB3B","RAB3C","RAB4A","RAB5B","RAB6A","RAB8B","RAET1E","RAET1G","RAET1L","RELB","RFTN1","SAR1B",
             "SLC11A1","TAP1","TAP2","TAPBP","TAPBPL","THBS1","TRAF6","TREM2","TREX1","ULBP1","ULBP2","ULBP3","UNC93B1",
             "WAS","WDFY4","YTHDF1")


Yapgenes <- c('ACAT1','ADAM9','CCN1','KLF5','ACTA2','COL4A1','LIF','LACTB',
              'AMOTL2','IL11','SORCS3','ANKRD1','CCN2','C11orf96','LBH','ANLN',
              'HBEGF','EDN1','MAP2K6','AURKB','HSPG2','ATF3','ARHGAP42','BIRC5',
              'LAMC1','F3','ADRA1B','CASP3','LPL','SERPINE1','CRIM1','CASP8','IER3',
              'C18orf63','CAVIN2','THBS1','TRIB1','CDH2','CCNF','TIMP2','FAM178A','CDC42',
              'VCAM1','KLF10','SLC38A2','CDK1','CALM2','SMAD7','ABI3BP','CDK6','COL1A1','FZD8',
              'TBX18','COL1A2','RUNX1','DYNLRB2','COL5A2','FGF2','NRD1','DIAPH3','COL6A2',
              'ZRANB1','DICER1','CXCL12','TNS3','E2F1','FBN1','MUSK','INHBA','USP44','EGFR',
              'LAMA2','PMS1','FGF1','LAMB1','PRELID2','FN1','TFPI','JDP2','FOXM1','VEGFA','AL035681.1',
              'GAB1','BDNF','BASP1','IGF1R','CCL2','QSOX1','LATS2','COL4A2','DDHD1','LDHA','CSF1','CITED2',
              'LIFR','SEMA3E','DUSP10','MCM6','ANXA1','OR51B6','MYC','BGN','EFHC1','MYL9','COL11A1','MAGI2',
              'NDUFB3','COL14A1','JARID2','NEGR1','COL5A1','CCDC132','NEK2','COL6A1','ARNT2','NF2','ANGPT1',
              'PELO','PLK1','ADAM12','CENPV','PMAIP1','CCL9','RPL22L1','POLA1','PTN','SYNDIG1','PPP1R3B','SEMA3A',
              'RAD18','PTEN','THBS2','TSPAN3','RAD51','HMGB1','TBC1D2','SLC2A3','IRS1','SNAI2','NUAK2','TAGLN','TRPC3',
              'TEAD4','GAS6','VIM','SNTB1','WWC1','MTMR10','ZEB1','BMP5','FRY','IGFBP3','AC112715.2','SUSD1','RUNX2',
              'ROBO2','AC023590.1','FBXO47','FAM107A','RBMS3','SH3RF1','SYT16','XPNPEP1','PIK3C2G',
              'MED13L','RND3','CDH4','TSPAN18','NEDD9','KPNA3','SH2D4A','TRIO','DPYD','LRRTM1','ZNF469',
              'RP11-553A10.1','GLIS3','SLITRK6','IQCJ-SCHIP1','FGF5','NTRK3','ADAMTS8','NAALADL2',
              'PTGER4','ANXA3','ERCC2','FAM172A','GNB1L','PSMG2','NMBR','SLC38A4','ROR1','PARD3B','BET1', 
              "FOXF2", "CCDC80", "MYOF","FJX1","PTPN14","ARHGEF17", "DOCK5","NT5E", "ASAP1", "GADD45A", "TGFB2", "AXL",
              "YAP1", "CTNNB1", "MST1",
              "PTK2", "ITGA5", "ITGAV", "ITGB1",  
              "SAV1", "STK3", "LATS1", "LATS2", 
              "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
              "VIM", "VCL", "VASP")

ecmgenes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL24A1", 
              "LUM", "DCN", "BGN", "HSPG2", "AGRN", 
              "FN1", "VTN", "TGFBI", "NID1", "NID2", "LAMC1", "LAMA1", "LAMA2", "SPARC", 
              "FBLN1", "FBLN2", "FBLN5", "LTBP1", "LTBP5", "EMILIN1", "MFAP4", "EFEMP1", 
              "FBN1", "FBN2", "TNC")

migration <- c("CTHRC1", "TNC", "CEMIP", "CDH2", "CD44", "ARID5A", "DDR2", "IQGAP1", "SLC8A1", "TNS1", "FGF2", "PLAU", "PRKCE", "PODXL", "ID1", 
               "CNN1", "SNAI2", "ACTC1", "MEOX2", "BMP2", "LUM", "FGF18", 
               "PIK3CA", "EDN1", "GRB2", "THBS1", "ITGB1", "ITGA3", 
               "ZEB1", "ZEB2", "TWIST1", "TWIST2")

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
fibro <- AddModuleScore(fibro, features = list(inflam), name = "Inflam")
fibro <- AddModuleScore(fibro, features = list(Yapsub), name = "Mechano")
fibro <- AddModuleScore(fibro, features = list(ecmgenes), name = "ECMgen")
fibro <- AddModuleScore(fibro, features = list(antigen), name = "antigen")
fibro <- AddModuleScore(fibro, features = list(migration), name = "Mig")

fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = rev(levels))
Idents(fibro) <- "fibro_annotation"
DotPlot(fibro, features = c("ECMgen1", "Mechano1", "Mig1" , "Inflam1", "antigen1"), group.by = "fibro_annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 15,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 15
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "Mechanical", "Migration", "Inflammatory", "Antigen-\nPresenting"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  )
ggsave("Figures/fibropopsignatures.pdf")

#Global New Resolution 
fibro <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/Fibro_Annotated_0.28.rds")
fibro <- FindClusters(fibro, resolution = 0.3 )
DimPlot(fibro, raster = TRUE, label = T) #Saving 5.24 x 3.27 in image


#### Subsignatures for YAP and Dotplots + global YAP sig plot  ####

levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
           "GREM1+ ssFib","KCNN3+ ssFib",
           "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib" , 
           "THBS1+ ssFib", "PI16+ ssFib")

revfibro <- fibro
Idents(revfibro) <- "fibro_annotation"
revfibro$fibro_annotation <- factor(revfibro$fibro_annotation, levels = rev(levels))
Idents(revfibro) <- "fibro_annotation"

Yapgenes <- c('ACAT1','ADAM9','CCN1','KLF5','ACTA2','COL4A1','LIF','LACTB',
              'AMOTL2','IL11','SORCS3','ANKRD1','CCN2','C11orf96','LBH','ANLN',
              'HBEGF','EDN1','MAP2K6','AURKB','HSPG2','ATF3','ARHGAP42','BIRC5',
              'LAMC1','F3','ADRA1B','CASP3','LPL','SERPINE1','CRIM1','CASP8','IER3',
              'C18orf63','CAVIN2','THBS1','TRIB1','CDH2','CCNF','TIMP2','FAM178A','CDC42',
              'VCAM1','KLF10','SLC38A2','CDK1','CALM2','SMAD7','ABI3BP','CDK6','COL1A1','FZD8',
              'TBX18','COL1A2','RUNX1','DYNLRB2','COL5A2','FGF2','NRD1','DIAPH3','COL6A2',
              'ZRANB1','DICER1','CXCL12','TNS3','E2F1','FBN1','MUSK','INHBA','USP44','EGFR',
              'LAMA2','PMS1','FGF1','LAMB1','PRELID2','FN1','TFPI','JDP2','FOXM1','VEGFA','AL035681.1',
              'GAB1','BDNF','BASP1','IGF1R','CCL2','QSOX1','LATS2','COL4A2','DDHD1','LDHA','CSF1','CITED2',
              'LIFR','SEMA3E','DUSP10','MCM6','ANXA1','OR51B6','MYC','BGN','EFHC1','MYL9','COL11A1','MAGI2',
              'NDUFB3','COL14A1','JARID2','NEGR1','COL5A1','CCDC132','NEK2','COL6A1','ARNT2','NF2','ANGPT1',
              'PELO','PLK1','ADAM12','CENPV','PMAIP1','CCL9','RPL22L1','POLA1','PTN','SYNDIG1','PPP1R3B','SEMA3A',
              'RAD18','PTEN','THBS2','TSPAN3','RAD51','HMGB1','TBC1D2','SLC2A3','IRS1','SNAI2','NUAK2','TAGLN','TRPC3',
              'TEAD4','GAS6','VIM','SNTB1','WWC1','MTMR10','ZEB1','BMP5','FRY','IGFBP3','AC112715.2','SUSD1','RUNX2',
              'ROBO2','AC023590.1','FBXO47','FAM107A','RBMS3','SH3RF1','SYT16','XPNPEP1','PIK3C2G',
              'MED13L','RND3','CDH4','TSPAN18','NEDD9','KPNA3','SH2D4A','TRIO','DPYD','LRRTM1','ZNF469',
              'RP11-553A10.1','GLIS3','SLITRK6','IQCJ-SCHIP1','FGF5','NTRK3','ADAMTS8','NAALADL2',
              'PTGER4','ANXA3','ERCC2','FAM172A','GNB1L','PSMG2','NMBR','SLC38A4','ROR1','PARD3B','BET1', 
              "FOXF2", "CCDC80", "MYOF","FJX1","PTPN14","ARHGEF17", "DOCK5","NT5E", "ASAP1", "GADD45A", "TGFB2", "AXL")

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
revfibro <- AddModuleScore(revfibro, features = list(Yapsub), name = "Yapsig") #Get score for YAP targets
#2. YAP components signature
components <- list(c("YAP1", "CTNNB1", "MST1",
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))
revfibro <- AddModuleScore(revfibro, features = components, name = "Components")
revfibro <- AddModuleScore(revfibro, features = list(ecmgenes), name = "ECMgen")

DotPlot(revfibro, features = c("ECMgen1", "Yapsig1", "Components1"), group.by = "fibro_annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 10,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 10
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "YAP\nPathway", "Target\nGenes"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  ) + 
coord_fixed(ratio = 0.5) +
  guides(
    color = guide_colorbar(title = "Average\nExpression"),
    size = guide_legend(title = "Percent\nExpressed")
  ) +
  theme(
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(hjust = 0.5)
  )+ 
coord_fixed(ratio = 0.5) +
  guides(
    color = guide_colorbar(title = "Average\nExpression"),
    size = guide_legend(title = "Percent\nExpressed")
  ) +
  theme(
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(hjust = 0.5)
  )


ggsave("Figures/fibropopsubsignatures.pdf")

#Global signature plots 
subglobal <- subset(global, downsample = 1000)

revsubglobal <- subglobal
revsubglobal$manual.annotation <- factor(subglobal$manual.annotation, levels = rev(levels))
Idents(revsubglobal) <- "manual.annotation"
levels(revsubglobal)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(revsubglobal)] #retain YAP genes only in rownames of revsubglobal 
genes.to.keep <- Matrix::rowSums(revsubglobal[["RNA"]]$counts > 0) >= floor(0.1 * ncol(revsubglobal[["RNA"]]$counts))
counts.sub <- revsubglobal[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of subglobalblasts
revsubglobal <- AddModuleScore(revsubglobal, features = list(Yapsub), name = "Yapsig") #Get score for YAP targets
revsubglobal <- AddModuleScore(revsubglobal, features = list(ecmgenes), name = "ECM") #Get score for YAP targets
revsubglobal <- AddModuleScore(revsubglobal, features = components, name = "Components") #Get score for YAP targets

DotPlot(revsubglobal, features = c("ECM1", "Components1", "Yapsig1"), group.by = "manual.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 10,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 10
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "YAP\nPathway", "Target\nGenes"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  ) +coord_fixed(ratio = 0.4) +
  guides(
    color = guide_colorbar(title = "Average\nExpression"),
    size = guide_legend(title = "Percent\nExpressed")
  ) +
  theme(
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(hjust = 0.5)
  )


ggsave("Figures/globalyapsig.pdf")


#### DotPlot of rep marker genes for signatures ####
genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", 
           "YAP1", "PTK2", "TEAD1", "LATS1", 
           "CCN1", "CCN2", "SERPINE1", "ADAM12", 
           "TNFRSF1B", "TNFAIP6", "IL7R", "INHBA", 
           "ZEB1", "ZEB2", "TWIST1", "CNN1", 
           "CD74", "HLA-DRA", "HLA-DRB5", "HLA-DRB1", 
           "CDH11", "FAP", "TIMP1", "SPARC")

levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
           "GREM1+ ssFib","KCNN3+ ssFib",
           "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib" , 
           "THBS1+ ssFib", "PI16+ ssFib")
fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = rev(levels) )

DotPlot(fibro, features = genes, group.by = "fibro_annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 10,
      angle = 90,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 0.5,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 10
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) +coord_fixed(ratio = 1) +
  guides(
    color = guide_colorbar(title = "Average\nExpression"),
    size = guide_legend(title = "Percent\nExpressed")
  ) +
  theme(
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(hjust = 0.5)
  )



ggsave("Figures/repdotplot.pdf")
#### Generate a controlled dataset ####
#To enable controlled comparisons across conditions
sub <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel" | diseasestate == "Stricture Inflamed", invert = T)
sub <- subset(sub, specimen_type == "Resection")
sub <- subset(sub, tissue == "Colon", invert = T)
#sub <- subset(sub, idents = c("F3+ ssFib", "ADAMDEC1+ ssFib", "CD74+ apFib", "MMP3+ iFib"), invert = T) #Just subserosal fibroblasts 
saveRDS(sub, "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/controlledfibrocomparison.updated.rds", compress = F)

Idents(sub) <- "fibro_annotation"

# #Create a barplot showing the percentage changes 
# sub$newlabels <- "NA"
# sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
# sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
# sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
# sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
# sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninflamed Bowel"
# #sub$newlabels[sub$diseasestate == "Stricture Inflamed"] <- "Inflamed Bowel"
# sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"
# 
# # sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel","Inflamed Bowel" , "Stricture", 
# #                                                       "Normal MAT", "Uninvolved MAT", "CF"))
# sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel", "Stricture", 
#                                                       "Normal MAT", "Uninvolved MAT", "CF"))
# Idents(sub) <- "newlabels"
# saveRDS(sub, "updatedcontrolleddataset_fibro.rds", compress = F)

####General comparison across fibrosis ###
#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved"
sub$newlabels[sub$diseasestate == "Involved"] <- "Stricture"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved"
#sub$newlabels[sub$diseasestate == "Stricture Inflamed"] <- "Inflamed Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"

# sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel","Inflamed Bowel" , "Stricture", 
#                                                       "Normal MAT", "Uninvolved MAT", "CF"))
# sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel", "Stricture", 
#                                                       "Normal MAT", "Uninvolved MAT", "CF"))
sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal", "Uninvolved", "Stricture"))
Idents(sub) <- "newlabels"

palette = c(
  "#8b0000",  # dark red
  "#d62728",  # red
  "#bcbd22",  # yellow-green
  "#ff7f0e",  # orange
  "#f5b800",  # bright yellow
  "#2ca02c",  # green
  "#1f77b4",  # blue
  "#e377c2",  # pink
  "#9467bd",  # purple
  "#17becf"   # cyan/teal
)
DimPlot(sub, raster = T, cols = palette, group.by = "fibro_annotation", split.by = "newlabels")

# Count cells per annotation per condition
fibro_counts <- sub@meta.data %>%
  group_by(newlabels, fibro_annotation) %>%
  summarise(n = n()) %>%
  ungroup()

# Optional: get % of each annotation within each condition
fibro_counts <- fibro_counts %>%
  group_by(newlabels) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

#Plot overall 
overall <- ggplot(fibro_counts, aes(x = newlabels, y = freq, fill = fibro_annotation)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = palette) +
  ylab("Fraction of Fibroblasts") +
  xlab("Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


results <- read.csv("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/anndata/results.csv")


#Plot split by organ 
fibro_counts <- sub@meta.data %>%
  group_by(diseasestate, fibro_annotation) %>%
  summarise(n = n()) %>%
  ungroup()
fibro_counts <- fibro_counts %>%
  group_by(diseasestate) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()
fibro_counts$diseasestate <- factor(fibro_counts$diseasestate, levels = c("Normal Bowel", "Stricture Uninflamed", "Stricture", 
                                                                          "Normal MAT", 'Uninvolved', "Involved"))


ggplot(fibro_counts, aes(x = diseasestate, y = freq, fill = fibro_annotation)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = palette) +
  ylab("Fraction of Fibroblasts") +
  xlab("Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 14, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 16, 
      face = "bold"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
      angle = 45, 
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal",
      "Stricture Uninflamed" = "Uninvolved", 
      "Stricture" = "Stricture", 
      "Normal MAT" = "Normal", 
      "Uninvolved" = "Uninvolved", 
      "Involved" = "CF"
    )
  )


# #Equal numbers of cells 
# disease_idents <- unique(sub$newlabels)
# cell.list <- WhichCells(sub, idents = disease_idents, downsample = min(table(sub$newlabels)))
# sub_small <- sub[, cell.list]
# table(sub_small$newlabels)
# Idents(sub_small) <- "newlabels"

# palette = rev(c("#9467bd", "#1f77b4", "#17becf", "#2ca02c", "#bcbd22",
#                 "#f5b800", "#ff7f0e","#e377c2", "#d62728", "#8b0000"))

ggsave("Figures/proptablefibroblastsall.pdf") #ggsave("Figures/proptablefibroblasts.svg", width = 8, height = 4, units = "in" )

#Remove inflamed condition for consistency through the paper 
# Remove rows where Var1 contains "Inflamed Bowel"
library(dplyr)
# ptab_filtered <- ptab %>% 
#   filter(Var1 != "Inflamed Bowel")
ggplot(ptab_filtered, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Percent (%)", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  facet_grid(~tissue, scales = "free_x", space = "free_x") + 
  theme(
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 15, margin = margin(b = 2)),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
      angle = 45, 
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal",
      "Uninflamed Bowel" = "Uninvolved",  # Keep this if still needed
      "Stricture Bowel" = "Stricture", 
      "Normal MAT" = "Normal", 
      "Uninvolved MAT" = "Uninvolved", 
      "CF MAT" = "CF"
    )
  )

ggsave("Figures/proptablefibroblasts3conditions.pdf")

DotPlot(sub, features = c("ECMgen1", "Mechano1", "Inflam1", "antigen1"), group.by = "newlabels", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15, angle = 45),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size))
  scale_x_discrete(labels=c("ECM", "Mechanical", "Inflammatory", "Antigen-Presenting")) + 
  scale_y_discrete(labels = c(
    "Normal Bowel" = "Normal",
    "Uninflamed Bowel" = "Uninvolved", 
    "Inflamed Bowel" = "Inflamed", 
    "Stricture Bowel" = "Stricture", 
    "Normal MAT" = "Normal", 
    "Uninvolved MAT" = "Uninvolved", 
    "CF MAT" = "CF"))
ggsave("newHuman/Figures/globalsignatures.pdf")

#### Plot CTHRC1 Activation only ####

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)


sub <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/controlledfibrocomparison.updated.rds")
Idents(sub) <- "fibro_annotation"
CTH <- subset(sub, idents = "CTHRC1+ mFib")

CTH <- AddModuleScore(CTH, features = list(ecmgenes), name = "ECMgen")
CTH <- AddModuleScore(CTH, features = components, name = "Components")
CTH <- AddModuleScore(CTH, features = list(Yapsub), name = "YAPgenes")

# Generate DotPlot
dotplot <- DotPlot(CTH, features = c("ECMgen1", "YAPgenes1", "Components1"), group.by = "fibro_annotation", split.by = "diseasestate", cols = "RdYlBu", scale.min = 0, scale.max = 100)

# Extract data from the DotPlot
dotplot_data <- dotplot$data

# Extract diseasestate from id
dotplot_data <- dotplot_data %>%
  separate(id, into = c("fibro_annotation", "diseasestate"), sep = "_", extra = "merge")

# Reorder disease states for the Y-axis
diseasestate_levels <- c("Normal Bowel", "Stricture Uninflamed", "Stricture", "Normal MAT", "Uninvolved", "Involved")
dotplot_data$diseasestate <- factor(dotplot_data$diseasestate, levels = diseasestate_levels)

# Add tissue type based on diseasestate mapping
diseasestate_to_tissue <- c("Normal Bowel" = "Bowel",
                            "Stricture Uninflamed" = "Bowel",
                            "Stricture" = "Bowel",
                            "Normal MAT" = "MAT",
                            "Uninvolved" = "MAT",
                            "Involved" = "MAT")

dotplot_data <- dotplot_data %>%
  mutate(tissue = recode(diseasestate, !!!diseasestate_to_tissue))

# Reorder features (x-axis) and rename "ECMgen1" to "ECM"
dotplot_data <- dotplot_data %>%
  mutate(features = factor(features.plot, levels = c("ECMgen1", "YAPgenes1", "Components1"), labels = c("ECM", "YAP\nTargets", "YAP\nPathway")))

# Overlay Bowel and MAT data
combined_data <- dotplot_data %>%
  filter(diseasestate %in% names(diseasestate_to_tissue))

# Define reversed RdYlBu palette
RdYlBu_palette <- rev(brewer.pal(11, "RdYlBu"))

# Define color gradient limits by scaled average expression
color_limits <- range(dotplot_data$avg.exp.scaled, na.rm = TRUE)

# Define custom y-axis labels
custom_y_labels <- c("Normal", "Uninvolved", "Stricture", "Normal", "Uninvolved", "CF")

# Plot combined graph
combined_plot <- ggplot(combined_data, aes(x = features, y = diseasestate, color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  scale_color_gradientn(colors = RdYlBu_palette, limits = color_limits, name = "Average<br>Expression") +
  scale_size(range = c(1, 6), limits = c(0, 100), name = "Percent<br>Expressed") +  # Scaled down dot size in dotplot
  scale_y_discrete(labels = custom_y_labels) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),  # Decreased tick length to half of previous
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.title = ggtext::element_markdown(size = 8, color = "black", hjust = 0.5),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(0.1, "lines"),  # Decreased size of the legend keys even more
    legend.margin = margin(grid::unit(0.5, "lines")),  # Decreased margin around legends
    plot.title = element_text(size = 16, color = "black", hjust = 0.5),
    panel.border = element_rect(color = 'black', size = 1, fill = NA),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    aspect.ratio = 1
  ) +
  guides(color = guide_colorbar(order = 1, barwidth = unit(0.2, "cm"), barheight = unit(1.5, "cm")),  # Decreased color bar size
         size = guide_legend(order = 2, keywidth = unit(0.15, "cm"), keyheight = unit(0.15, "cm"))) +  # Decreased size legend more
  labs(title = "Combined Bowel and MAT", x = "Signatures", y = "Condition")

# Print the combined plot
print(combined_plot)

ggsave("Figures/CTHRC1activation.pdf")



# DotPlot(sub, features = c("ECMgen1", "Mechano1", "Inflam1", "antigen1"), 
#         group.by = "newlabels", 
#         scale.max = 100, 
#         scale.min = 0, 
#         cols = "RdYlBu") + 
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(
#       size = 15,
#       angle = 45,
#       hjust = 1,   # Right-align angled text
#       vjust = 1,    # Anchor text to the top
#       margin = margin(t = 10, b = 5, unit = "pt")  # Add space above/below labels
#     ),
#     axis.text.y = element_text(
#       hjust = 1, 
#       vjust = 1, 
#       size = 15
#     ),
#     # Increase bottom margin to push labels away
#     plot.margin = margin(t = 10, r = 10, b = 60, l = 10, unit = "pt"),
#     legend.key.size = unit(0.8, 'cm'),
#     legend.key.height = unit(0.4, 'cm'),
#     legend.key.width = unit(0.4, 'cm'),
#     legend.title = element_text(size = 9),
#     legend.text = element_text(size = 9)
#   ) +
#   scale_x_discrete(
#     labels = c(
#       "ECM", 
#       "Mechanical", 
#       "Inflammatory", 
#       "Antigen-\nPresenting"  # Newline preserved
#     ),
#     guide = guide_axis(n.dodge = 1)  # Remove `offset` argument
#   ) + 
#   scale_y_discrete(
#     labels = c(
#       "Normal Bowel" = "Normal",
#       "Uninflamed Bowel" = "Uninvolved", 
#       "Stricture Bowel" = "Stricture", 
#       "Normal MAT" = "Normal", 
#       "Uninvolved MAT" = "Uninvolved", 
#       "CF MAT" = "CF"
#     )
#   )

ggsave("Figures/signaturesbyregion.pdf")










#### Bowel subanalysis ####
bowel <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/Fibro_Bowel/union_fibro_harmony_clusters_0.35.rds")
Idents(bowel) <- "seurat_clusters"
bowel$bowel.annotation <- "tbd"
bowel$bowel.annotation[bowel$seurat_clusters == 0] <- "ADAMDEC1+ lpFib"
bowel$bowel.annotation[bowel$seurat_clusters == 1] <- "PI16+ ssFRC"
bowel$bowel.annotation[bowel$seurat_clusters == 2] <- "F3+ Telocytes"
bowel$bowel.annotation[bowel$seurat_clusters == 3] <- "SOD2+ ssFRC"
bowel$bowel.annotation[bowel$seurat_clusters == 4] <- "ADAMDEC1+ lpFib"
bowel$bowel.annotation[bowel$seurat_clusters == 5] <- "KCNN3+ ssFRC"
bowel$bowel.annotation[bowel$seurat_clusters == 6] <- "GREM1+ ssFRC"
bowel$bowel.annotation[bowel$seurat_clusters == 7] <- "CTHRC1+ mFib"
bowel$bowel.annotation[bowel$seurat_clusters == 8] <- "MMP3+ iFib"
bowel$bowel.annotation[bowel$seurat_clusters == 9] <- "LRRC7+ mFib"
bowel$bowel.annotation[bowel$seurat_clusters == 10] <- "CD74+ apFib"
bowel$bowel.annotation[bowel$seurat_clusters == 11] <- "CD74+ apFib"
bowel$bowel.annotation[bowel$seurat_clusters == 12] <- "ADAMDEC1+ lpFib"
Idents(bowel) <- "bowel.annotation"
levels = c("CTHRC1+ mFib", "LRRC7+ mFib", 
  "MMP3+ iFib", "CD74+ apFib", "F3+ Telocytes", "ADAMDEC1+ lpFib",
  "KCNN3+ ssFRC", "GREM1+ ssFRC", "SOD2+ ssFRC", "PI16+ ssFRC")
bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = levels)
Idents(bowel) <- "bowel.annotation"
palette = rev(c(
    "#7b3294",  # deeper purple
    "#9467bd",  # original purple
    "#1f77b4",  # original blue
    "#17becf",  # original cyan
    "#2ca02c",  # original green
    "#bcbd22",  # original olive
    "#f5b800",  # original yellow
    "#ff9900",  # more vivid orange
    "#ff7f0e",  # original orange
    "#db4d6d",  # new warm rose/red (fills the gap)
    #  "#e377c2",  # original pink
    "#8b0000"   # original dark red
  ))
DimPlot(bowel, cols= palette, raster = F, pt.size = 0.01) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave("Figures/bowelumap.pdf")

saveRDS(bowel, "Fibro_Bowel_Annotated_0.35.rds", compress = F)

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(bowel)] #retain YAP genes only in rownames of bowel 
genes.to.keep <- Matrix::rowSums(bowel[["RNA"]]$counts > 0) >= floor(0.1 * ncol(bowel[["RNA"]]$counts))
counts.sub <- bowel[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
bowel <- AddModuleScore(bowel, features = list(inflam), name = "Inflam")
bowel <- AddModuleScore(bowel, features = list(Yapsub), name = "Mechano")
bowel <- AddModuleScore(bowel, features = list(ecmgenes), name = "ECMgen")
bowel <- AddModuleScore(bowel, features = list(antigen), name = "antigen")

bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = rev(levels))
DotPlot(bowel, features = c("ECMgen1", "Mechano1", "Inflam1", "antigen1"), group.by = "bowel.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 10,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 10
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "Mechanical", "Inflammatory", "Antigen-\nPresenting"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  )
ggsave("Figures/bowelsigplot.pdf")
#Generate proptable and barplots of fibroblasts in resection 


bowel <- subset(bowel, diseasestate == "Ileostomy Bowel" | diseasestate == "Stricture Inflamed", invert = T)
bowel <- subset(bowel, specimen_type == "Resection")
bowel <- subset(bowel, tissue == "SB")

bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = levels)

tab = table(bowel$diseasestate, bowel$bowel.annotation)
ptab = as.data.frame(prop.table(tab, 1))
ptab$Var1 <- factor(ptab$Var1, levels = c("Normal Bowel", "Stricture Uninflamed", "Stricture"))

ptab = as.data.frame(ptab)
ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Fibroblast Enrichment", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 12, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 15, 
      colour = "black"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
      angle = 45,
      hjust = 1, # Right-align labels
      vjust = 1 # Baseline alignment
      
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal", 
      "Stricture Uninflamed" = "Uninvolved", 
      "Stricture" = "Stricture"
    )
  )

ggsave("Figures/bowelproptable.pdf")

# Adjusting the Seurat FeaturePlot to fit the desired specifications
features = c("CTHRC1", "LRRC7", "MMP3", "CD74", "ADAMDEC1", "F3", "KCNN3", "GREM1", "THBS1", "PI16")

FeaturePlot(bowel, features = features, order = TRUE, raster = TRUE, max.cutoff = "q95", ncol = 5) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1) # Ensure square aspect ratio

ggsave("Figures/bowel_fibro_fplots.pdf")

#### MAT subanalysis ####
mat <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/MAT_Fibro2/union_mat_harmony_clusters_0.2.rds")

Idents(mat) <- "seurat_clusters"
mat$mat.annotation <- NA
mat$mat.annotation[mat$seurat_clusters == 0] <- "FMO2+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 1] <- "DPP4+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 2] <- "ICAM1+ iAPC"
mat$mat.annotation[mat$seurat_clusters == 3] <- "PPARG+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 4] <- "CTHRC1+ mFAP"
mat$mat.annotation[mat$seurat_clusters == 5] <- "VIT+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 6] <- "LRRC7+ mFAP"
mat$mat.annotation[mat$seurat_clusters == 7] <- "NOTCH3+ Pericytes"

mat$mat.annotation <- factor(mat$mat.annotation, levels = c("CTHRC1+ mFAP", "LRRC7+ mFAP", "ICAM1+ iAPC", "PPARG+ ssAPC", 
                                                            "VIT+ ssAPC", "FMO2+ ssAPC", "DPP4+ ssAPC", "NOTCH3+ Pericytes"))

Idents(mat) <- "mat.annotation"

palette = rev(c(
  "#7b3294",  # deeper purple
  "#1f77b4",  # original blue
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff7f0e",  # original orange
  "#e377c2",  # original pink
  "#8b0000"   # original dark red
))

DimPlot(mat, cols= palette, raster = F, pt.size = 0.01) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave("Figures/matumapupdated.pdf")
saveRDS(mat, "matfibroblast.annotated.rds", compress = F)

mat <- readRDS("matfibroblast.annotated.rds")
#Generate proptable and barplots of fibroblasts in resection 
mat <- subset(mat, diseasestate == "Ileostomy MAT", invert = T)
tab = table(mat$diseasestate, mat$mat.annotation)
ptab = as.data.frame(prop.table(tab, 1))
ptab$Var1 <- factor(ptab$Var1, levels = c("Normal MAT", "Uninvolved", "Involved"))

ptab = as.data.frame(ptab)
ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Fibroblast Enrichment", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 12, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 15, 
      colour = "black"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
        angle = 45,
        hjust = 1, # Right-align labels
        vjust = 1 # Baseline alignment
        
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal MAT" = "Normal", 
      "Uninvolved" = "Uninvolved", 
      "Involved" = "Creeping\nFat"
    )
  )

ggsave("Figures/MATproptable.pdf")

#Proportion per patient 
metadata <- mat@meta.data
metadata <- metadata %>%
  mutate(unique_patientID = paste(patientID, diseasestate, tissue, sep = "_"))

#Generate proptable and barplots of fibroblasts in resection 
tab = table(sub$newlabels, sub$fibro_annotation)
ptab = prop.table(tab, 1)*100
#tissue <- c("Bowel", "Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")

ptab = as.data.frame(ptab)
ptab$tissue <- tissue
ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Percent (%)", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  facet_grid(~tissue, scales = "free_x", space = "free_x") + 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 14, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 16, 
      face = "bold"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
      angle = 45, 
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal",
      "Uninflamed Bowel" = "Uninvolved", 
      "Inflamed Bowel" = "Inflamed", 
      "Stricture Bowel" = "Stricture", 
      "Normal MAT" = "Normal", 
      "Uninvolved MAT" = "Uninvolved", 
      "CF MAT" = "CF"
    )
  )


# Calculate the proportion of each fibroblast matpopulation per patientID
fibroblast_proportions <- metadata %>%
  group_by(unique_patientID, mat.annotation, diseasestate) %>%
  tally() %>%
  group_by(unique_patientID, diseasestate) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Average across duplicate patientIDs
average_fibroblast_proportions <- fibroblast_proportions %>%
  group_by(diseasestate, mat.annotation) %>%
  summarise(mean_proportion = mean(proportion)) %>%
  ungroup()

# Plotting the data
p <- ggplot(average_fibroblast_proportions, aes(x=factor(mat.annotation), y=mean_proportion)) +
  geom_bar(stat="identity") +
  facet_wrap(~ diseasestate, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Fibroblast matpopulations by Disease State",
       x = "Fibroblast mattype",
       y = "Mean Proportion")

print(p)

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(mat)] #retain YAP genes only in rownames of mat 
genes.to.keep <- Matrix::rowSums(mat[["RNA"]]$counts > 0) >= floor(0.1 * ncol(mat[["RNA"]]$counts))
counts.sub <- mat[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of matblasts
mat <- AddModuleScore(mat, features = list(inflam), name = "Inflam")
mat <- AddModuleScore(mat, features = list(Yapsub), name = "Mechano")
mat <- AddModuleScore(mat, features = list(ecmgenes), name = "ECMgen")
mat <- AddModuleScore(mat, features = list(antigen), name = "antigen")
mat <- AddModuleScore(mat, features = list(migration), name = "migration")

mat$mat.annotation <- factor(mat$mat.annotation, levels = rev(levels))
Idents(mat) <- "mat.annotation"
DotPlot(mat, features = c("ECMgen1", "Mechano1", "Inflam1", "antigen1"), group.by = "mat.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 10,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 10
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "Mechanical", "Inflammatory", "Antigen-\nPresenting"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  )

ggsave("Figures/MATsignatures.pdf")

require(pals)
require(reshape2)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Adjusting the Seurat FeaturePlot to fit the desired specifications
p <- FeaturePlot(
  fibro, 
  features = c("CTHRC1", "LRRC7", "ICAM1", "PPARG", "F3", "FMO2", "DPP4", "NOTCH3"),
  order = TRUE, 
  raster = TRUE, 
  max.cutoff = "q95"
) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1)),
    plot.background = element_rect(fill = NA, colour = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )


features = c("CTHRC1", "LRRC7", "ICAM1", "PPARG", "VIT", "FMO2", "DPP4", "NOTCH3")

FeaturePlot(mat, features = features, order = TRUE, raster = TRUE, max.cutoff = "q95", ncol = 4, pt.size = 2) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1) # Ensure square aspect ratio

ggsave("Figures/mat_fibro_fplots.pdf")
#### Generate a controlled dataset ####
sub <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
sub <- subset(sub, specimen_type == "Resection")
sub <- subset(sub, tissue == "Colon", invert = T)
Idents(sub) <- "fibro_annotation"

#### Global prop table (all cells) ####

#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninflamed Bowel"
sub$newlabels[sub$diseasestate == "Stricture Inflamed"] <- "Inflamed Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"

sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel","Inflamed Bowel" , "Stricture", 
                                                      "Normal MAT", "Uninvolved MAT", "CF"))
Idents(sub) <- "newlabels"


#### Individual fibroblast activation ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Original DotPlot
plot = DotPlot(
  sub, 
  features = c("Mechano1"),
  group.by = "diseasestate",
  split.by = "fibro_annotation",
  scale.max = 100,
  scale.min = 0,
  cols = "RdYlBu"
)

# Extract the plot data
plot_data <- plot$data

# Inspect the plot_data to understand its structure
head(plot_data)

# Separate `id` into `diseasestate` and `fibro_annotation`
plot_data <- plot_data %>%
  separate(id, into = c("diseasestate", "fibro_annotation"), sep = "_", extra = "merge")

# Ensure the new columns are factors with correct levels based on your provided lists
plot_data <- plot_data %>%
  mutate(
    diseasestate = factor(diseasestate, levels = c("Normal Bowel", "Stricture Uninflamed", "Stricture", "Normal MAT", "Uninvolved", "Involved")),
    fibro_annotation = factor(fibro_annotation, levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
                                                           "GREM1+ ssFib", "KCNN3+ ssFib", "THBS1+ ssFib", "PI16+ ssFib",
                                                           "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib"))
  )

# Debug by checking the unique values
unique(plot_data$diseasestate)
unique(plot_data$fibro_annotation)

# Now create a new plot with the adjusted data
adjusted_plot <- ggplot(plot_data, aes(x = diseasestate, y = fibro_annotation, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # Align the x-axis labels with hashes
    axis.text.x = element_text(
      size = 15,
      angle = 45,
      hjust = 1, # Right-align labels
      vjust = 1 # Baseline alignment
    ),
    # Y-axis labels
    axis.text.y = element_text(
      hjust = 1, # Right-align text
      vjust = 0.5, # Center vertically
      size = 15
    ),
    # Remove the gray background
    panel.background = element_blank(),
    # Add lines to x and y axes
    axis.line = element_line(color = "black")
  ) +
  # Make sure percent expressed goes from 0 to 100
  scale_size_continuous(name = "Percent Expressed", range = c(0, 10), limits = c(0, 100)) +
  # Reverse the palette so red goes with higher values and blue with smaller
  scale_color_gradientn(name = "Average Expression (scaled)", colors = rev(brewer.pal(n = 9, name = "RdYlBu"))) +
  
  # Rename x-axis labels (disease states)
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal",
      "Stricture Uninflamed" = "Uninvolved", 
      "Stricture" = "Stricture", 
      "Normal MAT" = "Normal",
      "Uninvolved" = "Uninvolved",
      "Involved" = "CF"
    )
  )

# Print the adjusted plot
print(adjusted_plot)

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
fibro <- AddModuleScore(fibro, features = list(inflam), name = "Inflam")
fibro <- AddModuleScore(fibro, features = list(Yapsub), name = "Mechano")
fibro <- AddModuleScore(fibro, features = list(ecmgenes), name = "ECMgen")
fibro <- AddModuleScore(fibro, features = list(antigen), name = "antigen")
fibro <- AddModuleScore(fibro, features = list(migration), name = "migration")

fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = rev(levels))
Idents(fibro) <- "mat.annotation"

levels(fibro) <- rev()
DotPlot(fibro, features = c("ECMgen1", "Mechano1", "Inflam1", "antigen1"), group.by = "mat.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # Align rotated labels with hashes (angle + hjust/vjust)
    axis.text.x = element_text(
      size = 15,
      angle = 45,      # Rotate 45°
      hjust = 1,       # Right-align labels with hashes
      vjust = 1,       # Baseline-align labels with hashes
      margin = margin(t = 5, unit = "pt")  # Small buffer above labels
    ),
    axis.text.y = element_text(
      hjust = 1,       # Right-align y-labels
      vjust = 0.5,     # Center vertically
      size = 15
    ), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) + 
  scale_x_discrete(
    labels = c("ECM", "Mechanical", "Inflammatory", "Antigen-\nPresenting"),
    guide = guide_axis(n.dodge = 1)  # Prevents label crowding
  )

ggsave("Figures/MATsignatures.pdf")

#### ComplexHeatMap of fibroblast DEGs ####
library(ComplexHeatmap)
library(circlize)

# Load necessary libraries
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

#Global fibroblasts 
# Set identifiers and specific order for fibro_annotation
Idents(fibro) <- "fibro_annotation"
levels <- c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
            "GREM1+ ssFib", "KCNN3+ ssFib",
            "THBS1+ ssFib", "PI16+ ssFib",
            "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib")
fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = levels)
Idents(fibro) <- "fibro_annotation"

# Subset the data
set.seed(123)  # Ensure reproducibility
sub <- subset(fibro, downsample = 1000)
sub[["RNA"]]$data <- as(object = sub[["RNA"]]$data, Class = "dgCMatrix")

# Find all markers with specified parameters
markers <- FindAllMarkers(sub, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

# Select the top 100 differentially expressed genes based on average log FC for each cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()

# Subset the data for the top 100 genes
expr_data <- GetAssayData(sub, slot = "data")[top100_genes, ]

# Calculate mean expression for each gene in each cluster
cluster_means <- AverageExpression(sub, features = top100_genes, group.by = "fibro_annotation")$RNA

# Reorder the rows of cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    cluster_expression <- matrix[level, , drop=FALSE]
    ordered_genes <- c(ordered_genes, colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)])
  }
  return(ordered_genes)
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Define the row splitting
row_labels <- rownames(cluster_means_z_t)
row_split <- factor(row_labels, levels = levels)

# Plot the heatmap
ht <- Heatmap(cluster_means_z_t, 
        name = "Z-score",
        col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
        show_row_names = TRUE,
        show_column_names = FALSE, 
        cluster_rows = FALSE,  # Maintain predefined row order
        cluster_columns = FALSE,  # Do not cluster columns 
        row_title = "Clusters",  # Label rows 
        column_title = "Top 100 Genes",  # Label columns 
        #row_split = row_split,  # Split rows by factor clusters
        border = TRUE
        #row_gap = unit(2, "mm")  # Fixed gap between split rows
)
pdf(file="Figures/globalfibrohm.pdf")
draw(ht)
dev.off()

# Reorder the rows of cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    cluster_expression <- matrix[level, , drop=FALSE]
    ordered_genes <- c(ordered_genes, colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)])
  }
  return(ordered_genes)
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Define the row splitting
row_labels <- rownames(cluster_means_z_t)
row_split <- factor(row_labels, levels = levels)

# Plot the heatmap
ht <- Heatmap(cluster_means_z_t, 
              name = "Z-score",
              col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
              show_row_names = TRUE,
              show_column_names = FALSE, 
              cluster_rows = TRUE,  # Cluster rows
              cluster_columns = FALSE,  # Do not cluster columns 
              row_title = "Clusters",  # Label rows 
              column_title = "Top 100 Genes",  # Label columns 
              row_split = row_split,  # Split rows by factor clusters
              row_gap = unit(0, "mm"),  # Remove gaps between rows
              border = FALSE)  # Remove outer border

# Draw the heatmap
draw(ht, heatmap_legend_side = "right")

pdf(file="Figures/globalfibrohmclustered.pdf", height = 2, width = 4.5)
draw(ht)
dev.off()


#MAT heatmap 
# Set identifiers and specific order for mat.annotation
levels = c("CTHRC1+ mFAP", "LRRC7+ mFAP", "ICAM1+ iAPC", "PPARG+ ssAPC", 
           "VIT+ ssAPC", "FMO2+ ssAPC", "DPP4+ ssAPC", "NOTCH3+ Pericytes")
mat <- readRDS("matfibroblast.annotated.rds")
Idents(mat) <- "mat.annotation"
mat$mat.annotation <- factor(mat$mat.annotation, levels = c("CTHRC1+ mFAP", "LRRC7+ mFAP", "ICAM1+ iAPC", "PPARG+ ssAPC", 
                                                            "VIT+ ssAPC", "FMO2+ ssAPC", "DPP4+ ssAPC", "NOTCH3+ Pericytes"))

Idents(mat) <- "mat.annotation"

# Subset the data
set.seed(123)  # Ensure reproducibility
sub <- subset(mat, downsample = 1000)
sub[["RNA"]]$data <- as(object = sub[["RNA"]]$data, Class = "dgCMatrix")

# Find all markers with specified parameters
markers <- FindAllMarkers(sub, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

# Select the top 100 differentially expressed genes based on average log FC for each cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()

# Subset the expression data for the top 100 genes
expr_data <- GetAssayData(sub, slot = "data")[top100_genes, ]

# Calculate mean expression for each gene in each cluster
cluster_means <- AverageExpression(sub, features = top100_genes, group.by = "mat.annotation")$RNA

# Reorder the rows of the cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    if (level %in% rownames(matrix)) {
      cluster_expression <- matrix[level, , drop = FALSE]
      ordered_genes_level <- colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)]
      ordered_genes <- c(ordered_genes, ordered_genes_level)
    } else {
      warning(paste("Level", level, "not found in the matrix row names."))
    }
  }
  return(unique(ordered_genes))
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z_t, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Transpose the matrix to get genes as rows and clusters as columns
heatmap_matrix <- t(cluster_means_z_t)

# Define genes of interest
genes_of_interest <- c("FMO2", "MGP", "PTN", "CTHRC1", "POSTN", "FAP", 
                       "DPP4", "PI16", "CD55", "ICAM1", "CCL2", "CXCL2", 
                       "LRRC7", "BNC2", "SDK1", "NOTCH3", "FRZB", "TAGLN", 
                       "PPARG", "RORA", "EBF1", "VIT", "MEOX2", "APOD")

# Create annotation for genes (now rows)
ha <- rowAnnotation(
  gene_labels = anno_mark(
    at = which(rownames(heatmap_matrix) %in% genes_of_interest),
    labels = rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% genes_of_interest]
  )
)

# Plot transposed heatmap
ht <- Heatmap(
  heatmap_matrix,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
  show_row_names = FALSE,        # Hide gene names (we'll use annotations)
  show_column_names = TRUE,      # Show cluster names
  column_names_side = "bottom",  # Cluster names below columns
  cluster_rows = TRUE,          # Preserve gene order
  cluster_columns = TRUE,       # Preserve cluster order
  row_title = "Genes",
  column_title = "Clusters",
  right_annotation = ha,         # Gene annotations on right
  border = TRUE, 
  show_column_dend = TRUE,
  show_row_dend = FALSE
)

pdf(file="Figures/matfibrohm.pdf")
draw(ht)
dev.off()


#Bowel heatmap 
# Set identifiers and specific order for mat.annotation
bowel <- readRDS("Fibro_Bowel_Annotated_0.35.rds")
levels = c("CTHRC1+ mFib", "LRRC7+ mFib", 
           "MMP3+ iFib", "CD74+ apFib", "F3+ Telocytes", "ADAMDEC1+ lpFib",
           "KCNN3+ ssFRC", "GREM1+ ssFRC", "SOD2+ ssFRC", "PI16+ ssFRC")
bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = levels)

Idents(bowel) <- "bowel.annotation"
bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = levels)

Idents(bowel) <- "bowel.annotation"

# Subset the data
set.seed(123)  # Ensure reproducibility
sub <- subset(bowel, downsample = 1000)
sub[["RNA"]]$data <- as(object = sub[["RNA"]]$data, Class = "dgCMatrix")

# Find all markers with specified parameters
markers <- FindAllMarkers(sub, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

# Select the top 100 differentially expressed genes based on average log FC for each cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()

# Subset the expression data for the top 100 genes
expr_data <- GetAssayData(sub, slot = "data")[top100_genes, ]

# Calculate mean expression for each gene in each cluster
cluster_means <- AverageExpression(sub, features = top100_genes, group.by = "bowel.annotation")$RNA

# Reorder the rows of the cluster_means matrix based on top100_genes order
cluster_means <- cluster_means[top100_genes, ]

# Normalize the mean expression using Z-scores
cluster_means_z <- t(scale(t(cluster_means)))

# Transpose the matrix so that clusters are on rows and genes are on columns
cluster_means_z_t <- t(cluster_means_z)

# Function to order genes by Z-scores within each cluster
order_genes_by_zscore <- function(matrix, levels) {
  ordered_genes <- c()
  for (level in levels) {
    if (level %in% rownames(matrix)) {
      cluster_expression <- matrix[level, , drop = FALSE]
      ordered_genes_level <- colnames(cluster_expression)[order(cluster_expression, decreasing = TRUE)]
      ordered_genes <- c(ordered_genes, ordered_genes_level)
    } else {
      warning(paste("Level", level, "not found in the matrix row names."))
    }
  }
  return(unique(ordered_genes))
}

# Order genes within each cluster
ordered_genes_by_cluster <- order_genes_by_zscore(cluster_means_z_t, levels)

# Reorder columns (genes) by the computed order
cluster_means_z_t <- cluster_means_z_t[, ordered_genes_by_cluster]

# Transpose the matrix to get genes as rows and clusters as columns
heatmap_matrix <- t(cluster_means_z_t)

# Define genes of interest
genes_of_interest <- c("FMO2", "MGP", "PTN", "CTHRC1", "POSTN", "FAP", 
                       "DPP4", "PI16", "CD55", "ICAM1", "CCL2", "CXCL2", 
                       "LRRC7", "BNC2", "SDK1", "NOTCH3", "FRZB", "TAGLN", 
                       "PPARG", "RORA", "EBF1", "VIT", "MEOX2", "APOD")

# Create annotation for genes (now rows)
ha <- rowAnnotation(
  gene_labels = anno_mark(
    at = which(rownames(heatmap_matrix) %in% genes_of_interest),
    labels = rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% genes_of_interest]
  )
)

# Plot transposed heatmap
ht <- Heatmap(
  heatmap_matrix,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
  show_row_names = FALSE,        # Hide gene names (we'll use annotations)
  show_column_names = TRUE,      # Show cluster names
  column_names_side = "bottom",  # Cluster names below columns
  cluster_rows = TRUE,          # Preserve gene order
  cluster_columns = TRUE,       # Preserve cluster order
  row_title = "Genes",
  column_title = "Clusters",
  right_annotation = ha,         # Gene annotations on right
  border = TRUE, 
  show_column_dend = TRUE,
  show_row_dend = FALSE
)

pdf(file="Figures/bowelfibrohm.pdf")
draw(ht)
dev.off()





#### Comparison UMAPs ####
bowel <- subset(fibro, tissue == "SB" | tissue == "Colon")
states <- subset(fibro, diseasestate == "Normal Bowel" | diseasestate == "Normal MAT" | 
                 diseasestate == "Ileostomy Bowel" |  diseasestate == "Ileostomy MAT")
library(Seurat)
library(ggplot2)

# Creating the DimPlot with specified options
DimPlot(bowel, split.by = "specimen_type", cols = palette, raster = TRUE) + theme(aspect.ratio = 1) + 
  theme(axis.title = element_blank(),  # Remove axis titles
        axis.text = element_blank(),   # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank())   # Remove axis lines

ggsave("Figures/specimen_type_comparison.pdf")


states$status <- "NA"
states$status[states$diseasestate == "Ileostomy Bowel"] <- "Ileostomy"
states$status[states$diseasestate == "Ileostomy MAT"] <- "Ileostomy"

states$status[states$diseasestate == "Normal Bowel"] <- "Normal"
states$status[states$diseasestate == "Normal MAT"] <- "Normal"

Idents(states) <- "fibro_annotation"


DimPlot(states, split.by = "status", cols = palette, raster = TRUE) + theme(aspect.ratio = 1) + 
  theme(axis.title = element_blank(),  # Remove axis titles
        axis.text = element_blank(),   # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank())   # Remove axis lines

ggsave("Figures/specimen_type_comparison_ileo.pdf")

#### Patient Giant Colormap ####
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggdendro)

# Modified color bar creator with horizontal legend formatting
# Define the function to create color bars

# Assume fibro is your Seurat object
metadata <- fibro@meta.data

# Recalculate proportions for fibroblast subtypes per patient
prop_table <- metadata %>%
  count(patientID, fibro_annotation) %>%
  group_by(patientID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(fibro_annotation = factor(fibro_annotation, 
                                   levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib", 
                                              "GREM1+ ssFib", "KCNN3+ ssFib", 
                                              "THBS1+ ssFib", "PI16+ ssFib", 
                                              "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib")))

# Prepare fibroblast colors with exact level matching
fibro_palette <- c(
  "#d62728", "#ff7f0e", "#e377c2", "#bcbd22", "#7f7f7f", 
  "#1f77b4", "#9467bd", "#8c564b", "#2ca02c", "#17becf"
)
fibro_colors <- setNames(fibro_palette, levels(prop_table$fibro_annotation))

# Average proportions for duplicate patient IDs
prop_table_avg <- prop_table %>%
  group_by(patientID, fibro_annotation) %>%
  summarise(prop = mean(prop), .groups = 'drop')

# Create fibroblast matrix for clustering ensuring unique patient IDs
fibro_matrix <- prop_table_avg %>%
  pivot_wider(names_from = fibro_annotation, values_from = prop, values_fill = list(prop = 0)) %>%
  column_to_rownames("patientID")

# Hierarchical clustering
dist_matrix <- dist(as.matrix(fibro_matrix), method = "euclidean")
hclust_res <- hclust(dist_matrix, method = "complete")

# Generate dendrogram data
dendro_data <- dendro_data(as.dendrogram(hclust_res), type = "rectangle")

# Reorder patient IDs based on clustering
ordered_patients <- hclust_res$labels[hclust_res$order]
prop_table_avg <- prop_table_avg %>% mutate(patientID = factor(patientID, levels = ordered_patients))

# Create dendrogram plot without labels or axes
p1 <- ggplot(ggdendro::segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lineend = "round", size = 0.4) +  # Correcting the axes
  theme_void() +  # Ensure we have void theme (no extra margins or axes)
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

# Create bar plot aligned with the dendrogram labels with added outlines
p2 <- ggplot(prop_table_avg, aes(x = patientID, y = prop, fill = fibro_annotation)) +
  geom_bar(stat = "identity", width = 1, position = position_stack(reverse = FALSE), color = "black", size = 0.2) +
  scale_fill_manual(values = fibro_colors) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0), legend.position = "bottom", legend.direction = "horizontal", 
        legend.title = element_text(hjust = 0.5, size = 8, margin = margin(b = 2)), 
        legend.text = element_text(size = 7, margin = margin(l = 6)), 
        legend.spacing.x = unit(0.3, "cm"), legend.key.size = unit(0.4, "cm"), 
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
  guides(fill = guide_legend(title.position = "top", label.position = "right", byrow = TRUE, nrow = 3, ncol = 4))

make_color_bar <- function(data, colname, col_colors) {
  data[[colname]] <- factor(data[[colname]], levels = names(col_colors))
  
  ggplot(data, aes(x = factor(patientID), y = 1, fill = !!sym(colname))) + 
    geom_tile(height = 0.2, width = 1, color = "black", size = 0.1) + 
    scale_fill_manual(values = col_colors, name = colname) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),  # No space around the plot
      panel.spacing = unit(0, "lines"),  # No space between panels
      legend.position = "none"           # Remove legends
    )
}


metadata_columns <- metadata # Assuming 'metadata_columns' is same as 'metadata' here

# Convert patientID in metadata_columns to a factor with the levels in the desired order
metadata_columns <- metadata_columns %>%
  mutate(patientID = factor(patientID, levels = ordered_patients))

# Order the metadata_columns dataframe by patientID
metadata_columns <- metadata_columns %>%
  arrange(patientID)

# Generate color bars with horizontal legends
specimen_type_barplot <- make_color_bar(metadata_columns, "specimen_type", specimen_type_colors)
tissue_barplot <- make_color_bar(metadata_columns, "tissue", tissue_colors)
treated_barplot <- make_color_bar(metadata_columns, "treated", treated_colors)
sex_barplot <- make_color_bar(metadata_columns, "sex", sex_colors)
age_barplot <- make_color_bar(metadata_columns, "age", age_colors)
diseasestate_barplot <- make_color_bar(metadata_columns, "diseasestate", diseasestate_colors)

# Adjust legends for specific plots as requested
diseasestate_barplot <- diseasestate_barplot + guides(fill = guide_legend(title.position = "top", label.position = "right", byrow = TRUE, nrow = 3, ncol = 5))
treated_barplot <- treated_barplot + guides(fill = guide_legend(title.position = "top", label.position = "right", byrow = TRUE, nrow = 3, ncol = 4))

# Combine dendrogram, bar plot, and metadata color bars vertically
combined_plots <- wrap_plots(
  p1 + theme(plot.margin = margin(0, 0, 0, 0)),
  p2,
  diseasestate_barplot,
  specimen_type_barplot,
  tissue_barplot,
  treated_barplot,
  sex_barplot,
  age_barplot,
  ncol = 1,
  heights = c(0.4, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = 5, b = 5),
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(t = 8, b = 8)
  ) &
  guides(fill = guide_legend(
    nrow = 3,
    byrow = TRUE,
    direction = "horizontal",
    title.position = "top",
    title.hjust = 0.5,
    label.position = "right",
    label.hjust = 0.5,
    keywidth = unit(0.25, "cm"),
    keyheight = unit(0.25, "cm"),
    override.aes = list(
      shape = 22,
      size = 1.5,
      color = NA
    ),
    title.theme = element_text(size = 8),
    label.theme = element_text(
      size = 7,
      margin = margin(l = 1, r = 1),
      family = "sans"
    )
  ))

print(combined_plots)

# To save with proper dimensions
ggsave("final_plot_colors.pdf", final_plot, 
       width = 15, height = 6, units = "in")


#### Patient Comparisons at Global and Organ Level ####
##Global 
metadata <- sub@meta.data
metadata <- metadata %>%
  mutate(unique_patientID = paste(patientID, diseasestate, tissue, sep = "_"))

# Calculate the proportion of each fibroblast subpopulation per patientID
fibroblast_proportions <- metadata %>%
  group_by(unique_patientID, fibro_annotation, diseasestate) %>%
  tally() %>%
  group_by(unique_patientID, diseasestate) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Average across duplicate patientIDs
average_fibroblast_proportions <- fibroblast_proportions %>%
  group_by(diseasestate, fibro_annotation) %>%
  summarise(mean_proportion = mean(proportion)) %>%
  ungroup()

levels <- c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib", 
              "GREM1+ ssFib", "KCNN3+ ssFib", 
              "THBS1+ ssFib", "PI16+ ssFib", 
              "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib")
# Plotting the data
p <- ggplot(average_fibroblast_proportions, aes(x=factor(fibro_annotation, levels=levels), y=mean_proportion)) +
  geom_bar(stat="identity") +
  facet_wrap(~ diseasestate, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Fibroblast Subpopulations by Disease State",
       x = "Fibroblast Subtype",
       y = "Mean Proportion")

print(p)

# Filter to include only CTHRC1+ mFib
cthrc1_mfib_data <- fibroblast_proportions %>%
  filter(fibro_annotation == "CTHRC1+ mFib")

write.csv(cthrc1_mfib_data, "cthrc1.csv", row.names = FALSE)

##MAT fibroblasts only 
mat <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/MAT_Fibro2/union_mat_harmony_clusters_0.2.rds")

Idents(mat) <- "seurat_clusters"
mat$mat.annotation <- NA
mat$mat.annotation[mat$seurat_clusters == 0] <- "FMO2+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 1] <- "DPP4+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 2] <- "ICAM1+ iAPC"
mat$mat.annotation[mat$seurat_clusters == 3] <- "PPARG+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 4] <- "CTHRC1+ mFAP"
mat$mat.annotation[mat$seurat_clusters == 5] <- "VIT+ ssAPC"
mat$mat.annotation[mat$seurat_clusters == 6] <- "LRRC7+ mFAP"
mat$mat.annotation[mat$seurat_clusters == 7] <- "NOTCH3+ Pericytes"

mat$mat.annotation <- factor(mat$mat.annotation, levels = c("CTHRC1+ mFAP", "LRRC7+ mFAP", "ICAM1+ iAPC", "PPARG+ ssAPC", 
                                                    "VIT+ ssAPC", "FMO2+ ssAPC", "DPP4+ ssAPC", "NOTCH3+ Pericytes"))

Idents(mat) <- "mat.annotation"

palette <- c(
  "#d62728",  # muted red
  "#ff7f0e",  # muted orange
  "#e377c2",  # muted pink
  "#bcbd22",  # muted yellow-green
  "#7f7f7f",  # muted grey
  "#1f77b4",  # muted blue
  "#9467bd",  # muted purple
  "#2ca02c",  # muted green
  "#17becf"   # muted cyan
)
DimPlot(mat, raster = F, cols = palette)+
  theme(
    aspect.ratio = 1,  # Make the aspect ratio square
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_blank(),  # Remove axis lines
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank()  # Remove axis labels
  )
ggsave("Figures/matumapupdated.pdf")
saveRDS(mat, "matfibroblast.annotated.rds", compress = F)

#Generate proptable and barplots of fibroblasts in resection 
mat <- subset(mat, diseasestate == "Ileostomy MAT", invert = T)
tab = table(mat$diseasestate, mat$mat.annotation)
ptab = prop.table(tab, 1)
ptab$Var1 <- factor(ptab$Var1, levels = c("Normal MAT", "Uninvolved", "Involved"))

ptab = as.data.frame(ptab)
ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Fibroblast Enrichment", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 12, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 15, 
      colour = "black"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black"
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal MAT" = "N", 
      "Uninvolved" = "UI", 
      "Involved" = "CF"
    )
  )

ggsave("Figures/MATproptable.pdf")

#Proportion per patient 
metadata <- mat@meta.data
metadata <- metadata %>%
  mutate(unique_patientID = paste(patientID, diseasestate, tissue, sep = "_"))

#Generate proptable and barplots of fibroblasts in resection 
tab = table(sub$newlabels, sub$fibro_annotation)
ptab = prop.table(tab, 1)*100
#tissue <- c("Bowel", "Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")

ptab = as.data.frame(ptab)
ptab$tissue <- tissue
ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(
    x = "Disease State", 
    y = "Percent (%)", 
    fill = "Cluster"
  ) + 
  scale_fill_manual(values = palette) + 
  facet_grid(~tissue, scales = "free_x", space = "free_x") + 
  theme(
    # Regular y-axis numbers (removed face = "bold")
    axis.text.y = element_text(
      size = 14, 
      colour = "black"  # No bold here
    ),
    # Bold y-axis label only
    axis.title.y = element_text(
      size = 16, 
      face = "bold"
    ),
    # Rest remains the same
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 15, 
      margin = margin(b = 2)
    ),
    axis.text.x = element_text(
      size = 12, 
      colour = "black", 
      angle = 45, 
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  ) + 
  scale_x_discrete(
    labels = c(
      "Normal Bowel" = "Normal",
      "Uninflamed Bowel" = "Uninvolved", 
      "Inflamed Bowel" = "Inflamed", 
      "Stricture Bowel" = "Stricture", 
      "Normal MAT" = "Normal", 
      "Uninvolved MAT" = "Uninvolved", 
      "CF MAT" = "CF"
    )
  )


# Calculate the proportion of each fibroblast matpopulation per patientID
fibroblast_proportions <- metadata %>%
  group_by(unique_patientID, mat.annotation, diseasestate) %>%
  tally() %>%
  group_by(unique_patientID, diseasestate) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Average across duplicate patientIDs
average_fibroblast_proportions <- fibroblast_proportions %>%
  group_by(diseasestate, mat.annotation) %>%
  summarise(mean_proportion = mean(proportion)) %>%
  ungroup()

# Plotting the data
p <- ggplot(average_fibroblast_proportions, aes(x=factor(mat.annotation), y=mean_proportion)) +
  geom_bar(stat="identity") +
  facet_wrap(~ diseasestate, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Fibroblast matpopulations by Disease State",
       x = "Fibroblast mattype",
       y = "Mean Proportion")

print(p)

# Filter to include only CTHRC1+ mFib
cthrc1_mfib_data <- fibroblast_proportions %>%
  filter(fibro.annotation == "mCDF CTHRC1+")

write.csv(cthrc1_mfib_data, "cthrc1matonly.csv", row.names = FALSE)

##Bowel only fibroblasts 
bowel <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Fibro_Bowel_Annotated_0.35.rds")
metadata <- bowel@meta.data
metadata <- metadata %>%
  mutate(unique_patientID = paste(patientID, diseasestate, tissue, sep = "_"))

# Calculate the proportion of each fibroblast bowelpopulation per patientID
fibroblast_proportions <- metadata %>%
  group_by(unique_patientID, fibro_annotation, diseasestate) %>%
  tally() %>%
  group_by(unique_patientID, diseasestate) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Average across duplicate patientIDs
average_fibroblast_proportions <- fibroblast_proportions %>%
  group_by(diseasestate, fibro_annotation) %>%
  summarise(mean_proportion = mean(proportion)) %>%
  ungroup()

# Plotting the data
p <- ggplot(average_fibroblast_proportions, aes(x=factor(fibro_annotation), y=mean_proportion)) +
  geom_bar(stat="identity") +
  facet_wrap(~ diseasestate, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Fibroblast bowelpopulations by Disease State",
       x = "Fibroblast boweltype",
       y = "Mean Proportion")

print(p)

# Filter to include only CTHRC1+ mFib
cthrc1_mfib_data <- fibroblast_proportions %>%
  filter(fibro_annotation == "CTHRC1+ mFib")

write.csv(cthrc1_mfib_data, "cthrc1bowelonly.csv", row.names = FALSE)

#### Covariation analysis for treatment ####
fibro_cthrc1 <- subset(fibro, fibro_annotation == "CTHRC1+ mFib")
fibro_cthrc1 <- subset(fibro_cthrc1, treated == "NA", invert = T)
table(fibro_cthrc1$fibro_annotation, fibro_cthrc1$treated)

# Create a data frame summarizing the proportions
data_summary_cthrc1 <- fibro_cthrc1@meta.data %>%
  group_by(tissue, fibro_annotation, diseasestate, treated) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(diseasestate, treated) %>%
  mutate(proportion = count / sum(count))

# Transform categories to factors if they are not already
data_summary_cthrc1$fibro_annotation <- as.factor(data_summary_cthrc1$fibro_annotation)
data_summary_cthrc1$tissue <- as.factor(data_summary_cthrc1$tissue)
data_summary_cthrc1$diseasestate <- as.factor(data_summary_cthrc1$diseasestate)
data_summary_cthrc1$treated <- as.factor(data_summary_cthrc1$treated)

# Display the structure of the data before analysis
str(data_summary_cthrc1)

ancova_results_cthrc1 <- aov(proportion ~ diseasestate * treated + diseasestate + treated + tissue, data = data_summary_cthrc1)

summary(ancova_results_cthrc1)



#### Ligand Analysis for Pericyte to CTHRC1 transition ####
setwd("/home/kbauerro/Oak/IBDAFs/MetaUpdate/nichenetr")
library(nichenetr)
library(tidyverse)

#Load in networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

#Get target cell expression
receiver = "CTHRC1+ mFAP"
expressed_genes_receiver <- get_expressed_genes(receiver, mat, pct = 0.10)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("NOTCH3+ Pericytes")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, mat, 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

#DEGs between CTHRC1+ fibroblasts and NOTCH3low Pericytes
cluster.markers <- FindMarkers(mat, ident.1 = "NOTCH3+ Pericytes", ident.2 = "CTHRC1+ mFAP", only.pos = TRUE, min.pct = 0.10) %>% rownames_to_column("gene")
#cluster.markers <- FindMarkers(mat, ident.1 = "NOTCH3+ Pericytes", ident.2 = "CTHRC1+ mFAP", only.pos = TRUE)
write.csv(cluster.markers, "clustermarkers_CTHRC1_Pericytes.csv") 

DE_table_receiver <- read.csv("clustermarkers_CTHRC1_Pericytes.csv")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes <- rownames(mat)

#Run LR analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)


(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected))))
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>% bind_rows()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
order_targets <- order_targets[1:50]

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Differentiation Ligands", "DEGs",
                                               color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target_network

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Differentiation Ligands", "Ligand \n Activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank())
p_ligand_aupr


p_ligand_target_network <- p_ligand_target_network +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(), 
    plot.margin = unit(c(1, 0, 0, 0), "cm")
  )

p_ligand_aupr <- p_ligand_aupr +
  theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))

figures_without_legend = plot_grid(
  ncol = 2,
  nrow = 1,
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_target_network + theme(legend.position = "none"), 
  align = "hv", axis = "tb",
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_target))+1,
  rel_heights = c(nrow(vis_ligand_aupr)+3)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_aupr)),
  as_ggplot(get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h") + theme(plot.margin = unit(c(0,0,2,0),'cm'))

plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,1), nrow = 2, align = "hv")
ggsave("Figures/CTHRC1_Pericyte_activation.pdf")


#### CTHRC1+ spatially-informed nichenet analysis ####
global <- readRDS("Reference_Global.rds")
Idents(global) <- "subpop"
mat <- subset(global, tissue == "MAT")
#Create general disease state categories
mat$newlabels <- "NA"
mat$newlabels[mat$diseasestate == "Normal MAT"] <- "Non-Fibrotic"
mat$newlabels[mat$diseasestate == "Uninvolved"] <- "Non-Fibrotic"
mat$newlabels[mat$diseasestate == "Involved"] <- "Fibrotic"
sub <- subset(mat, diseasestate == "Ileostomy MAT", invert = T)
niche <- c("EBF1+ Steady-State Peri", "CD74+ apFib", "CCR7+ TH1 T Cells", "AFF3+ Naïve B Cells", 
           "CA4hi Capillary ECs", "JUNhi Activated Peri", "CXCL12hi Microvascular ECs", "GZMB+ pDC", 
           "TNFRSF13B+ Memory B Cells", "LYVE1+ Macrophages", "PI16+ ssFib", "IGHG1+ Plasma")
mat_cells <- intersect(levels(mat), levels(global))
mat_niche <- intersect(levels(mat), niche)
receiver = "CTHRC1+ mFib"
expressed_genes_receiver <- get_expressed_genes(receiver, sub, pct = 0.10)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- mat_niche

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, mat, 0.10)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

#Set conditions of interest and reference
condition_oi <-  "Fibrotic"
condition_reference <- "Non-Fibrotic"

Idents(sub) <- "subpop"
seurat_obj_receiver <- subset(sub, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "newlabels",
                                  test.use = "MAST", 
                                  min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

order_targets <- order_targets[1:50]
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")