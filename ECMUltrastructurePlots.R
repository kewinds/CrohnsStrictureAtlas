library(monocle3)
library(dplyr)
library(ggplot2)

#### Human Stricture Bowel ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/HumanStricPicro/Bowel")
quantified <- read.csv(file = "Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "Stricture_bowel_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
for (i in 1:length(metadata)) {
  if (metadata$Region[i] == "DS Bowel") {
    metadata$Region[i] <-  "Stricture Bowel"
  }
  if (metadata$Region[i] == "DS CF") {
    metadata$Region[i] <- "Stricture CF"
  }
  if (metadata$Region[i] == "DS IF") {
    metadata$Region[i] <- "Stricture IF"
  }
  if (metadata$Region[i] == "ADJ IF") {
    metadata$Region[i] <- "ADJ Bowel"
  }
}
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
#### Timecourse Bowel ####
#Set directory and import files 
setwd("/oak/stanford/groups/longaker/KEBR/Picro/TimeCourse/Bowel/Subset")
quantified <- read.csv(file = "Bowel_TC_Quantified.csv", row.names = 1, check.names = FALSE) #Picro quantification 
#umap <- read.csv(file = "UMAP.csv") #UMAP coordinates
metadata <- read.csv(file = "Bowel_TC_UMAP.csv", row.names = 1)
rownames(quantified) <- rownames(metadata)
unique(metadata$Region)

#Change labels 
metadata$Region[metadata$Region == "ADJ Bowel" ] <- "ADJ"
metadata$Region[metadata$Region == "ADJ IF" ] <- "ADJ"
metadata$Region[metadata$Region == "Sham Bowel" ] <- "Sham"
metadata$Region[metadata$Region == "POD3 IF" ] <- "POD 3"
metadata$Region[metadata$Region == "POD3 Bowel" ] <- "POD 3"
metadata$Region[metadata$Region == "POD 7 IF" ] <- "POD 7"
metadata$Region[metadata$Region == "POD 7 Bowel" ] <- "POD 7"
metadata$Region[metadata$Region == "POD 14 IF" ] <- "POD 14"
metadata$Region[metadata$Region == "POD 14 Bowel" ] <- "POD 14"
metadata$Region[metadata$Region == "POD 30 IF" ] <- "POD 30"
metadata$Region[metadata$Region == "POD 30 Bowel" ] <- "POD 30"
metadata$Region[metadata$Region == "POD 90 IF" ] <- "POD 90"
metadata$Region[metadata$Region == "POD 90 Bowel" ] <- "POD 90"

bowel <- grep("_B_", rownames(metadata))
metadata$subregion <- NA
metadata$subregion[bowel] <- "Bowel"
bowel <- grep("Sham ", rownames(metadata))
metadata$subregion[bowel] <- "Bowel"

IF <- grep("_IF_", rownames(metadata))
metadata$subregion[IF] <- "IF"

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
colData(cds)$Region <- factor(colData(cds)$Region, levels = c("Sham", "ADJ", "POD 3", "POD 7", "POD 14", "POD 30", "POD 90"))
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

ggsave("Figures/BowelUMAPall.svg")

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
ggsave("Figures/BowelUMAPall_split.svg")

cds <- cluster_cells(cds, resolution=1e-3)
plot_cells(cds, graph_label_size=1.5,
           cell_size = 1.0, 
           group_cells_by="cluster", 
           label_groups_by_cluster=TRUE, 
           label_cell_groups=FALSE)
ggsave("Figures/TCMATUMAPall.svg")

#Plot by subregion 
groups <- c("POD 7", "POD 14", "POD 30", "POD 90")
pattern <- paste(groups, collapse = "|")
cols <- grep(pattern, colData(cds)$Region)
cds_sub <- cds[, cols]

cds_sub <- cluster_cells(cds_sub, resolution=1e-3)
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
ggsave("Figures/BowelIFonly.svg") 

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
ggsave("Figures/BowelIFonlysplit.svg") 

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
ggsave("Figures/m3_TCmat.svg")
saveRDS(cds, "TCbowelmonocle3.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds))
data.pseudo$monocle3_pseudotime <- pseudotime(cds)
write.csv(data.pseudo, "TCbowel.csv")
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
cds_sub <- learn_graph(cds_sub)
cds_sub <- order_cells(cds_sub)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE,
           graph_label_size=1.5,
           cell_size = 1.0)
ggsave("Figures/m3_TCbowelsub.svg")
saveRDS(cds_sub, "TCbowelmonocle3sub.rds")
#Boxplot of pseudotime values 
data.pseudo <- as.data.frame(colData(cds_sub))
data.pseudo$monocle3_pseudotime <- pseudotime(cds_sub)
write.csv(data.pseudo, "TCbowelsub.csv")
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
ggsave("Figures/KOBowelumap.svg")
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
ggsave("Figures/KObowelsplit.svg")

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
ggsave("Figures/m3_YAPKObowel.svg")
saveRDS(cds, "m3_YAPKObowel.rds")
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

