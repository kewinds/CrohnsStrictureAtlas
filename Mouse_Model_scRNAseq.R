setwd("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy/r_packages",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers)
library(scales); library(SeuratObject); library(ComplexHeatmap); library(circlize); library(tidyverse)
options(Seurat.object.assay.version = "v5")

#### Create Seurat objects ####
#Create Seurat objects from the two studies: 
obj1 = makeSeuratHash('KB-19228_Sham', useHash=F); obj1$orig.ident = 'Sham_1'
obj2 = makeSeuratHash('KB-19228_Colo', useHash=F); obj2$orig.ident = 'Colotomy_1'

saveRDS(obj1, "readyForSeurat/Sham_1_preprocess.rds")
saveRDS(obj2, "readyForSeurat/Colotomy_1_preprocess.rds")

#### QC and merge objects ####
#Perform QC individually on objects
#Preprocess the objects and QC 
obj_list <- list.files(paste0("readyForSeurat/"),"*.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-15)
write.csv(obj_list, "processed/obj_list_01_preQC.csv")
for (i in 1:length(obj_list)) {
  objName <- obj_list[i]
  x <- readRDS(paste0("readyForSeurat/",objName, "_preprocess.rds"))
  # QC
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-01-pre-filt.jpg"))
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(x$nFeature_RNA, 0.95) & percent.mt < 10)
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-02-post-filt.jpg"))
  saveRDS(x, paste0("processed/",objName, "_01_qc.rds"), compress = FALSE)
}

obj_list <- list.files(paste0("processed/"),"*qc.rds", recursive = TRUE)
spatial.list <- list()
for (obj in obj_list) {
  x <- readRDS(paste0("processed/",obj))
  spatial.list <- append(spatial.list, x)
}
merged <- merge(spatial.list[[1]],spatial.list[2:length(spatial.list)]) #Matrix package sensitive. 
merged <- JoinLayers(merged)
saveRDS(merged, "processed/merged.global.processed.rds")

#Add in additional metadata
merged$tissue <- NA
merged$tissue[merged$orig.ident == "Colotomy_1" & merged$HTO_maxID == "HTO-1"] <- "Bowel"
merged$tissue[merged$orig.ident == "Colotomy_1" & merged$HTO_maxID == "HTO-2"] <- "MAT"
merged$tissue[merged$orig.ident == "Sham_1" & merged$HTO_maxID == "HTO-1"] <- "Bowel"
merged$tissue[merged$orig.ident == "Sham_1" & merged$HTO_maxID == "HTO-2"] <- "MAT"

merged$condition <- NA
merged$condition[merged$orig.ident == "Colotomy_1"] <- "Colotomy"
merged$condition[merged$orig.ident == "Sham_1"] <- "Sham"

saveRDS(merged, "processed/merged.global.processed.rds")

#### SCTransform and CCA R1 ####
global <- readRDS("mouseMerged_global_harmony_umap_res0.1.ANNOTATEDTRUE.rds")
Idents(global) <- "global.annotation"
fibro <- subset(global, idents = "Fibroblasts")

# re-generate seurat fibroect and remove extraneous data
options(Seurat.object.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
#saveRDS( fibro, 'SCT_CCA/newFibro1/SCT_CCA_fibro1.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
#write_matrix_dir(mat = counts.mat, dir = "SCT_CCA/newFibro1/SCT_CCA_fibro1_counts", overwrite = T )
#counts.mat <- open_matrix_dir(     dir = "SCT_CCA/newFibro1/SCT_CCA_fibro1_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
#saveRDS( fibro, 'SCT_CCA/newFibro1/SCT_CCA_fibro1_counts.rds', compress = F )

rm(list = ls())
fibro <- readRDS('SCT_CCA/newFibro1/SCT_CCA_fibro1_counts.rds')
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- SCTransform(fibro, assay = "RNA") #Fails if matrixStats is not 1.1.0. Make sure you downgrade. 
#saveRDS(fibro, 'SCT_CCA/newFibro1/SCT_CCA_transformed_fibro1_counts.rds', compress = F )
fibro <- RunPCA(fibro)
#saveRDS(fibro, 'SCT_CCA/newFibro1/SCT_CCA_transformed_pca_fibro1_counts.rds', compress = F )
fibro <- IntegrateLayers(object = fibro, method = CCAIntegration, normalization.method = "SCT", verbose = F)
#saveRDS(fibro, 'SCT_CCA/newFibro1/SCT_CCA_transformed_pca_cca_fibro1_counts.rds', compress = F )
fibro <- RunUMAP(fibro, dims = 1:30, reduction = "integrated.dr")
#saveRDS(fibro, 'SCT_CCA/newFibro1/SCT_CCA_transformed_pca_cca_umap_fibro1_counts.rds', compress = F )
fibro <- FindNeighbors(fibro, reduction = "integrated.dr", dims = 1:30)
fibro <- JoinLayers(fibro)
saveRDS(fibro, 'SCT_CCA/newFibro1/SCT_CCA_transformed_pca_cca_umap_joined_fibro1_counts.rds', compress = F )

#for (i in seq(0.03, 0.1, 0.01)) { 
for (i in seq(0.1, 1.5, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, label = TRUE)
  ggsave(paste0("SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_",i,".rds"), compress = FALSE)
}
DimPlot(fibro, reduction = "umap")

rm(list =ls())
fibro <- readRDS("SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_0.3.rds")

fibro <- PrepSCTFindMarkers(fibro, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_0.3_topmarkers.csv")


#### SCTransform and CCA R2 ####
rm(list = ls())
fibro <- readRDS("SCT_CCA/newFibro1/fibro1_SCT_CCA_clusters_0.3.rds")
Idents(fibro) <- "seurat_clusters"
fibro <- subset(fibro, idents = c(2, 4, 5, 9), invert = T)

# re-generate seurat fibroect and remove extraneous data
options(Seurat.object.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( fibro, 'SCT_CCA/newFibro2/SCT_CCA_fibro2.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "SCT_CCA/newFibro2/SCT_CCA_fibro2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "SCT_CCA/newFibro2/SCT_CCA_fibro2_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
saveRDS( fibro, 'SCT_CCA/newFibro2/SCT_CCA_fibro2_counts.rds', compress = F )

rm(list = ls())
fibro <- readRDS('SCT_CCA/newFibro2/SCT_CCA_fibro2_counts.rds')
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- SCTransform(fibro, assay = "RNA") #Fails if matrixStats is not 1.1.0. Make sure you downgrade. 
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_fibro2_counts.rds', compress = F )
fibro <- RunPCA(fibro)
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_fibro2_counts.rds', compress = F )
fibro <- IntegrateLayers(object = fibro, method = CCAIntegration, normalization.method = "SCT", verbose = F)
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_cca_fibro2_counts.rds', compress = F )
fibro <- RunUMAP(fibro, dims = 1:30, reduction = "integrated.dr")
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_cca_umap_fibro2_counts.rds', compress = F )
fibro <- FindNeighbors(fibro, reduction = "integrated.dr", dims = 1:30)
fibro <- JoinLayers(fibro)
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_cca_umap_joined_fibro2_counts.rds', compress = F )


#Run downstream 
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
options(future.globals.maxSize = 3e+09)
fibro <- SCTransform(fibro)
fibro <- RunPCA(fibro, npcs = 30, verbose = F)
fibro <- IntegrateLayers(
  object  = fibro,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", 
  new.reduction = "cca.integration",
  verbose = T, 
  assay = "SCT"
)

fibro <- RunUMAP(fibro, dims = 1:30, reduction = "cca.integration")
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_cca_umap_fibro2_counts.rds', compress = F )
fibro <- FindNeighbors(fibro, dims = 1:30, reduction = "cca.integration")
fibro <- JoinLayers(fibro)
saveRDS(fibro, 'SCT_CCA/newFibro2/SCT_CCA_transformed_pca_cca_umap_joined_fibro2_counts.rds', compress = F )

#for (i in seq(0.03, 0.1, 0.01)) { 
for (i in seq(0.1, 1.5, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, label = TRUE)
  ggsave(paste0("SCT_CCA/newFibro2/fibro2_SCT_CCA_clusters_",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("SCT_CCA/newFibro2/fibro2_SCT_CCA_clusters_",i,".rds"), compress = FALSE)
}
DimPlot(fibro, reduction = "umap")

rm(list =ls())

#Round 3 
fibro <- readRDS("SCT_CCA/newFibro2/fibro2_SCT_CCA_clusters_0.6.rds")
Idents(fibro) <- "seurat_clusters"
fibro <- subset(fibro, idents = c(0, 3, 7, 9, 12), invert = T)

# re-generate seurat fibroect and remove extraneous data
options(Seurat.object.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( fibro, 'SCT_CCA/newFibro3/SCT_CCA_fibro3.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "SCT_CCA/newFibro3/SCT_CCA_fibro3_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "SCT_CCA/newFibro3/SCT_CCA_fibro3_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
saveRDS( fibro, 'SCT_CCA/newFibro3/SCT_CCA_fibro3_counts.rds', compress = F )

fibro <- readRDS('SCT_CCA/newFibro3/SCT_CCA_fibro3_counts.rds')
#Run downstream 
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- SCTransform(fibro)
fibro <- RunPCA(fibro, npcs = 30, verbose = F)
fibro <- IntegrateLayers(object  = fibro, method = CCAIntegration, orig.reduction = "pca", new.reduction = "cca.integration",
                         verbose = T, normalization.method = "SCT")
fibro <- RunUMAP(fibro, dims = 1:30, reduction = "cca.integration")
saveRDS(fibro, 'SCT_CCA/newFibro3/SCT_CCA_transformed_pca_cca_umap_fibro3_counts.rds', compress = F )
fibro <- FindNeighbors(fibro, dims = 1:30, reduction = "cca.integration")
saveRDS(fibro, 'SCT_CCA/newFibro3/SCT_CCA_transformed_pca_cca_umap_joined_fibro3_counts.rds', compress = F )

#for (i in seq(0.03, 0.1, 0.01)) { 
for (i in seq(0.1, 1.5, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, label = TRUE)
  ggsave(paste0("SCT_CCA/newFibro3/fibro3_SCT_CCA_clusters_",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("SCT_CCA/newFibro3/fibro3_SCT_CCA_clusters_",i,".rds"), compress = FALSE)
}
DimPlot(fibro, reduction = "umap")

rm(list =ls())

fibro <- readRDS("SCT_CCA/newFibro2/fibro3_SCT_CCA_clusters_0.7.rds")

fibro <- PrepSCTFindMarkers(fibro, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "SCT_CCA/newFibro3/fibro3_SCT_CCA_clusters_0.7_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "SCT_CCA/newFibro3/fibro3_SCT_CCA_clusters_0.7_topmarkers.csv")

#Annotations
fibro$fibro.annotation <- NA
fibro$fibro.annotation[fibro$seurat_clusters == 0] <- "Pi16+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 1] <- "Col15a1+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 2] <- "Abca8a+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 3] <- "Cthrc1+ mFib"
fibro$fibro.annotation[fibro$seurat_clusters == 4] <- "Grem1+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 5] <- "Mgp+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 6] <- "Cthrc1+ mFib"
fibro$fibro.annotation[fibro$seurat_clusters == 7] <- "Pi16+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 8] <- "Eln+ mFib"
fibro$fibro.annotation[fibro$seurat_clusters == 9] <- "Adamdec1+ ssFib"
fibro$fibro.annotation[fibro$seurat_clusters == 10] <- "Cd74+ apFib"
fibro$fibro.annotation[fibro$seurat_clusters == 11] <- "Timp1+ iFib"

Idents(fibro) <- "fibro.annotation"

fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = c("Cthrc1+ mFib", "Eln+ mFib", "Timp1+ iFib",
                                                                        "Cd74+ apFib", "Adamdec1+ ssFib", "Mgp+ ssFib", 
                                                                        "Grem1+ ssFib", "Abca8a+ ssFib", 
                                                                        "Col15a1+ ssFib", "Pi16+ ssFib"))
Idents(fibro) <- "fibro.annotation"
DimPlot(fibro, cols = palette)
DimPlot(fibro, cols = palette, split.by = "group")

Idents(fibro) <- "fibro.annotation"
colors = c("#B2182B", "#ef8a62", "grey92", "#d1e5f0", "#67a9cf", "#2166ac")
dimplot <- DimPlot(fibro, cols = colors ) + theme(legend.key.size = unit(0.5, "in"))
levels(fibro) <- rev(c("Mechanosensitive Postn+", "Steady-State Abca8a+", "Steady-State Mki67+",
                       "Inflammatory Adamdec1+", "Steady-State Pi16+"))

palette = rev(c("#9467bd", "#1f77b4", "#17becf", "#2ca02c", "#bcbd22",
                "#f5b800", "#ff7f0e","#e377c2", "#d62728", "#8b0000"))

saveRDS(fibro, "mousefibro.CCA_0.7.annotated.rds")

#### Phenotyping Signatures ####
fibro <- readRDS("mousefibro.CCA_0.7.annotated.rds")
palette = rev(c(
  "#7b3294",  # deeper purple
  "#1f77b4",  # original blue
  "#17becf",  # original cyan
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff7f0e",  # original orange
  "#db4d6d",  # new warm rose/red (fills the gap)
  "#e377c2",  # original pink
  "#8b0000"   # original dark red
))

DimPlot(fibro, cols = palette, raster = F, pt.size = 1) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box around the plot
  )

ggsave(fibro, "mousedimplotupdated.pdf")

# Find all markers with specified parameters
fibro <- PrepSCTFindMarkers(fibro, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = TRUE)
significant_genes <- markers[markers$p_val_adj < 0.05, ]

levels = c("Cthrc1+ mFib", "Eln+ mFib", "Timp1+ iFib",
           "Cd74+ apFib", "Adamdec1+ ssFib", "Mgp+ ssFib", 
           "Grem1+ ssFib", "Abca8a+ ssFib", 
           "Col15a1+ ssFib", "Pi16+ ssFib")

# Get top 100 genes per cluster
top100_genes <- significant_genes %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(factor(cluster, levels = levels)) %>%
  pull(gene) %>%
  unique()
write.csv(significant_genes, "mousesig_genes.csv")
# Ensure that SCT is the default assay
DefaultAssay(fibro) <- "SCT"

# Extract SCT features
sct_features <- rownames(fibro[["SCT"]])

# Ensure top100_genes contains genes available in the SCT assay
valid_top100_genes <- intersect(top100_genes, sct_features)
if (length(valid_top100_genes) < length(top100_genes)) {
  warning("Some genes from top100_genes are not found in the SCT assay. Using only valid genes.")
}

cvalid_top100_genes <- c(valid_top100_genes)
# Calculate mean expression matrix using the counts slot of the SCT assay
cluster_means <- AverageExpression(
  fibro, 
  features = valid_top100_genes, 
  group.by = "fibro.annotation", 
  slot = "counts"  # Use the counts slot for counts data
)$SCT %>% 
  .[valid_top100_genes, levels]  # Enforce cluster order
# Z-score normalization
cluster_means_z <- t(scale(t(cluster_means)))

# New ordering strategy: Sort genes by max Z-score in their primary cluster
gene_order <- lapply(levels, function(cluster) {
  cluster_genes <- significant_genes %>% 
    filter(cluster == !!cluster) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene) %>% 
    head(100)
  
  # Get Z-scores for this cluster
  z_scores <- cluster_means_z[cluster_genes, cluster]
  
  # Order genes by Z-score descending
  cluster_genes[order(z_scores, decreasing = TRUE)]
}) %>% unlist() %>% unique()


# Apply the new ordering
ordered_matrix <- cluster_means_z[gene_order, ]

# Define genes of interest
genes_of_interest <- c("Cthrc1", "Postn", "Sparc", 
                       "Timp1", "Mmp3", 
                       "Abca8a", "Junb", "Pdgfra", 
                       "Pi16", "Sema3c", "Cd55", 
                       "Grem1", "Fmo2", "Aspn",
                       "Piezo2", "Eln", "Igf1", 
                       "Mgp", "Apoe", "Mdk",
                       "Col15a1", "Gas1" ,"Aebp1", 
                       "Adamdec1", "Edil3", "Hapln1",
                       "Cd74", "Mki67", "H2-Eb1")

# Create annotation
ha <- rowAnnotation(
  gene_labels = anno_mark(
    at = which(rownames(ordered_matrix) %in% genes_of_interest),
    labels = rownames(ordered_matrix)[rownames(ordered_matrix) %in% genes_of_interest],
    labels_gp = gpar(fontsize = 8)
  )
)

# Create heatmap
ht <- Heatmap(
  ordered_matrix,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  cluster_rows = FALSE,  # Crucial: disable row clustering
  cluster_columns = TRUE,  # Keep predefined cluster order
  row_title = "Genes",
  column_title = "Clusters",
  right_annotation = ha,
  border = TRUE,
  show_column_dend = TRUE,
  show_row_dend = FALSE,
  column_names_rot = 90
)

draw(ht)

# Draw the heatmap
pdf(file="Figures/fibromouse_hm.pdf", height = 7, width = 7)
draw(ht)
dev.off()

#Feature plots
require(pals)
require(reshape2)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Adjusting the Seurat FeaturePlot to fit the desired specifications
genes_of_interest <- c( "Postn","Timp1",  "Abca8a", "Pi16", "Grem1", "Mgp",  "Eln", 
                        "Col15a1","Adamdec1", "Cd74")

FeaturePlot(fibro, features = genes_of_interest, order = FALSE, raster = FALSE, ncol = 5) &
  scale_color_gradientn(colours = brewer.pal(9, "Reds"), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1)+
  theme(plot.background = element_rect(fill = "white", color = NA)) # Ensure square aspect ratio

ggsave("Figures/fibro_fplots.pdf")


#### Mouse YAP signatures ####
#1. Yap Gene Signatures
Yapgenes <- c("Ccn1", "Ccn2", "Foxf2", "Igfbp3", "Ccdc80", "F3", "Myof", "Fjx1", "Rbms3",
              "Ptpn14", "Ankrd1", "Serpine1", "Lats2", "Arhgef17", "Nuak2", "Amotl2", 
              "Dock5", "Crim1", "Nt5e", "Asap1", "Gadd45a", "Tgfb2", "Axl", 
              "Yap1", "Ctnnb1", 
              "Ptk2", "Itga5", "Itgav", "Itgb1",  
              "Sav1", "Stk3", "Lats1", "Lats2", 
              "Tead1", "Tead2", "Tead3", "Tead4", "Runx1", 
              "Vim", "Vcl", "Vasp") #see AKH Loe et al. 

#2. YAP components signature
components <- list(c("Yap1", "Ctnnb1", 
                     "Ptk2", "Itga5", "Itgav", "Itgb1",  
                     "Sav1", "Stk3", "Lats1", "Lats2", 
                     "Tead1", "Tead2", "Tead3", "Tead4", "Runx1", 
                     "Vim", "Vcl", "Vasp"))

#3. YAP overall 
overallYAP <- c(Yapgenes, components)

#4. ECM protein signature 
ecmgenes <- list(c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col24a1", 
                   "Lum", "Dcn", "Bgn", "Hspg2", "Agrn", 
                   "Fn1", "Vtn", "Tgfbi", "Nid1", "Nid2", "Lamc1", "Lama1", "Lama2", "Sparc", 
                   "Fbln1", "Fbln2", "Fbln5", "Ltbp1", "Ltbp5", "Emilin1", "Mfap4", "Efemp1", "Fbn1", "Fbn2", "Tnc"))

#5. Inflammatory signature
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
inflam <- str_to_sentence(tolower(inflam))

#6. Migration signature
migration <- c("CTHRC1", "TNC", "CEMIP", "CDH2", "CD44", "ARID5A", "DDR2", "IQGAP1", "SLC8A1", "TNS1", "FGF2", "PLAU", "PRKCE", "PODXL", "ID1", 
               "CNN1", "SNAI2", "ACTC1", "MEOX2", "BMP2", "LUM", "FGF18", 
               "PIK3CA", "EDN1", "GRB2", "THBS1", "ITGB1", "ITGA3", 
               "ZEB1", "ZEB2", "TWIST1", "TWIST2", 
               "ARID5B", "SGPL1", "AKT1", "TNS1", "ZFAND5", "CCN3", "MACIR", "ILK1", 
               "TMEM201", "FUT8", "NHERF1", "PIP5K1A", "IQGAP1", "PDLIM1", "FAM114A1", "PLEC", "AQP1", "PML", "CD248", "LAMTOR2") #combo of "Fibroblast Migration GO + known migration genes
migration <- str_to_sentence(tolower(migration))

#7. Antigen-presentation signature GO:0019882
antigen <- c("Abcc1","Ap3b1", "Ap3d1","Arl8b","Atg5","Azgp1",
             "B2m","Bag6","Calr","Ccl19","Ccl21a","Ccr7",
             "Cd1d1","Cd1d2","Cd68","Cd74","Clec4a2","Clec4a3", "Clec4a4", "Clec4b2","Cst7", "Ctse",
             "Ctsl","Ctss",  "Ext1","Fam3d", "Fcer1g","Fcer2a", "Fcgr1", "Fcgr2b","Fcgr3",
             "Fcgr4","Fgl2", "Flt3","Gba1","H2-Aa","H2-Ab1","H2-D1", "H2-DMa","H2-DMb1",
             "H2-DMb2","H2-Ea","H2-Eb1","H2-Eb2","H2-K1","H2-L","H2-M1",
             "H2-M2","H2-M3","H2-M5","H2-M9","H2-M10.1","H2-M10.2","H2-M10.3", "H2-M10.4","H2-M10.5",
             "H2-M10.6","H2-M11","H2-Oa","H2-Ob","H2-Q1","H2-Q2", "H2-Q4","H2-Q6","H2-Q7", "H2-Q8",
             "H2-Q9","H2-Q10","H2-T3","H2-T5","H2-T9","H2-T13","H2-T15","H2-T22", "H2-T23",
             "H2-T24","H60b","H60c","Hfe", "Icam1","Ide","Ifi30","Ifng","Ighe","Ighg2a",
             "Ighm","Kdm5d","Lgmn","Marchf1","Marchf8","Mfsd6","Mpeg1",
             "Mr1","Nod1","Nod2","Pdia3","Pikfyve","Psap","Psmb8","Psmb9","Psme1",
             "Psme2","Ptpn22","Pycard","Rab3b","Rab3c","Rab4a","Rab5b","Rab6a",
             "Rab8b","Rab10","Rab27a", "Rab32","Rab33a","Rab34","Rab35","Raet1a",
             "Raet1b","Raet1c","Raet1d","Raet1e","Relb","Rftn1","Slc11a1","Tap1","Tap2",
             "Tapbp","Tapbpl","Thbs1","Traf6","Trem2","Treml4","Trex1","Ulbp1",
             "Unc93b1","Was", "Washc1","Wdfy4","Ythdf1")

#Overall signatures 
overallYAP <- overallYAP[overallYAP %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
antigen  <- antigen[antigen %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 

#genes.to.keep <- Matrix::rowSums(fibro[["SCT"]]$counts > 0) >= floor(0.05 * ncol(fibro[["SCT"]]$counts))
#counts.sub <- fibro[["SCT"]]$counts[genes.to.keep,]
#Yapsub <- overallYAP[overallYAP %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
fibro <- AddModuleScore(fibro, features = overallYAP, name = "Mechanical")
fibro <- AddModuleScore(fibro, features = ecmgenes, name = "ECMgenes")
fibro <- AddModuleScore(fibro, features = list(inflam), name = "Inflam")
fibro <- AddModuleScore(fibro, features = list(antigen), name = "antigen")
fibro <- AddModuleScore(fibro, features = list(migration), name = "migration")

#Specific YAP signatures
# Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
# genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
# counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
# Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
fibro <- AddModuleScore(fibro, features = list(Yapgenes), name = "Yapsig")
fibro <- AddModuleScore(fibro, features = components, name = "Components")


Idents(fibro) <- "fibro.annotation"

fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = rev(c(levels = c("Cthrc1+ mFib", "Eln+ mFib", "Timp1+ iFib",
                                                                                       "Cd74+ apFib", "Adamdec1+ ssFib", "Mgp+ ssFib", 
                                                                                       "Grem1+ ssFib", "Abca8a+ ssFib", 
                                                                                       "Col15a1+ ssFib", "Pi16+ ssFib"))))
Idents(fibro) <- "fibro.annotation"
DotPlot(fibro, features = c("ECMgenes1", "Components1", "Yapsig1"), scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),  # center justified
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.title = element_text(size=9),
        legend.text = element_text(size=9)) +
  scale_x_discrete(labels = c('ECM\n', 'YAP\nPathway', 'YAP\nTargets'))
ggsave("Figures/YAPsignatures.pdf")

DotPlot(fibro, features = c("ECMgenes1", "Components1", "Yapsig1"), group.by = "group", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),  # center justified
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.title = element_text(size=9),
        legend.text = element_text(size=9), 
        aspect.ratio = 1) +
  scale_x_discrete(labels = c('ECM\n', 'YAP\nPathway', 'YAP\nTargets'))
ggsave("Figures/YAPsignaturesbygroup.pdf")


fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = rev(levels(fibro)))
Idents(fibro) <- "fibro.annotation"
DotPlot(fibro, features = c("ECMgenes1", "Mechanical1", "migration1", "Inflam1", "antigen1"), scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c("ECM", "Mechanical", "Migration", 'Inflammation', "Antigen-\nPresentation"))

ggsave("Figures/ECMsig_by_fibroblast.pdf")

fibro$group <- factor(x = fibro$group, levels = c("Sham Bowel", "Colotomy Bowel", "Sham MAT", "Colotomy MAT"))
Idents(fibro) <- "group"

DotPlot(fibro, features = c("ECMgenes1", "Mechanical1", "migration1", 'Inflam1', "antigen1"), group.by = "group", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c("ECM", "Mechanics", 'Migration', 'Inflammation', "Antigen\n-Presentation"))

#Create metadata category for each specific experimental group

fibro$group <- NA
fibro$group[fibro$tissue == "MAT" & fibro$condition == "Sham"] <- "Sham MAT"
fibro$group[fibro$tissue == "MAT" & fibro$condition == "Colotomy"] <- "Colotomy MAT"
fibro$group[fibro$tissue == "Bowel" & fibro$condition == "Sham"] <- "Sham Bowel"
fibro$group[fibro$tissue == "Bowel" & fibro$condition == "Colotomy"] <- "Colotomy Bowel"
Idents(fibro) <- "group"
#Compare conditions
cell.list <- WhichCells(fibro, idents = c("Sham MAT", "Colotomy MAT", "Sham Bowel", "Colotomy Bowel"), downsample = min(table(fibro$group)))
fibro_small <- fibro[, cell.list]
table(fibro_small$group)
Idents(fibro_small) <- "fibro.annotation"
fibro_small$group <- factor(x = fibro_small$group, levels = c("Sham Bowel", "Colotomy Bowel", "Sham MAT", "Colotomy MAT"))
Idents(fibro_small) <- "group"
Idents(fibro_small) <- "fibro.annotation"
DimPlot(fibro_small, split.by = "group", cols = colors)
ggsave("Figures/tissuesplitdimplot.jpg", width = 10.1, height = 2.76)

#Group by tissue
fibro$group <- NA
fibro$group[fibro$tissue == "MAT" & fibro$condition == "Sham"] <- "Sham MAT"
fibro$group[fibro$tissue == "MAT" & fibro$condition == "Colotomy"] <- "Colotomy MAT"
fibro$group[fibro$tissue == "Bowel" & fibro$condition == "Sham"] <- "Sham Bowel"
fibro$group[fibro$tissue == "Bowel" & fibro$condition == "Colotomy"] <- "Colotomy Bowel"
Idents(fibro) <- "group"
fibro$group <- factor(x = fibro$group, levels = c("Sham Bowel", "Colotomy Bowel", "Sham MAT", "Colotomy MAT"))

DotPlot(fibro, features = c("Components1", "Yapsig1", "ECMgenes1"), group.by = "group", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c("YAP Pathway Genes", "YAP Target Genes", 'ECM Genes'))
ggsave("Figures/YAPsignaturesbyregion.svg")

DotPlot(fibro, features = c("ECMgenes1", "Yapsig1", "migration1", "Inflam1", "antigen1"), group.by = "group", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c("ECM", "Mechanical", "Migration",  'Inflammation', 'Antigen-\n Presentation'))+ 
  scale_y_discrete(labels=c("Sham", "Colotomy", 'Sham', 'Colotomy'))
ggsave("Figures/YAPsignaturesbyregion.pdf")

Idents(fibro) <- "fibro.annotation"
genes <- c("Col1a1", "Col1a2", "Col3a1", "Col4a1", 
           "Tead1", "Tead2", "Ptk2", "Yap1", 
           "Ccn1", "Ccn2", "Serpine1", "Tgfb2", 
           "Fap", "Acta2", "Cdh11", "Sparc")
DotPlot(fibro, features = genes, group.by = "fibro.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) #+  #change legend text font size)) 
ggsave("Figures/repfibrogenes.pdf")


#Generate proptable and barplots of fibroblasts in mouse
Idents(fibro) <- "fibro.annotation"

fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = c(levels = c("Cthrc1+ mFib", "Eln+ mFib", "Timp1+ iFib",
                                                                                   "Cd74+ apFib", "Adamdec1+ ssFib", "Mgp+ ssFib", 
                                                                                   "Grem1+ ssFib", "Abca8a+ ssFib", 
                                                                                   "Col15a1+ ssFib", "Pi16+ ssFib")))
Idents(fibro) <- "fibro.annotation"


tab = table(fibro$group, fibro$fibro.annotation)
ptab = prop.table(tab, 1)
ptab = as.data.frame(ptab)
ptab$condition <- NA
ptab$condition[ptab$Var1 == "Sham Bowel" | ptab$Var1 == "Sham MAT"] <- "Sham"
ptab$condition[ptab$Var1 == "Colotomy Bowel" | ptab$Var1 == "Colotomy MAT"] <- "Colotomy"
ptab$tissue[ptab$Var1 == "Sham Bowel" | ptab$Var1 == "Colotomy Bowel"] <- "Bowel"
ptab$tissue[ptab$Var1 == "Sham MAT" | ptab$Var1 == "Colotomy MAT"] <- "MAT"

ptab$Var1 <- factor(ptab$Var1, levels = c("Sham Bowel", "Colotomy Bowel", "Sham MAT", "Colotomy MAT"))


ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col() +
  theme_classic() + 
  labs(y = "Enrichment", 
       fill = "Cluster") + 
  scale_fill_manual(values = palette) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15)) + 
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_discrete(labels=c('Sham', 'Colotomy', 'Sham', 'Colotomy')) + 
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size = 15, margin = margin(0, 0, 5, 0)),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),  # Increasing legend text size
        legend.key.size = unit(1.5, "line")) +  
  theme(panel.spacing = unit(0, "lines"))
ggsave("Figures/proptablemousefibroblasts.pdf")

#Plot YAP target gene activation 
sub <- subset(fibro, idents = c("Cthrc1+ mFib"))
sub$group <- NA
sub$group[sub$tissue == "MAT" & sub$condition == "Sham"] <- "Sham MAT"
sub$group[sub$tissue == "MAT" & sub$condition == "Colotomy"] <- "Colotomy MAT"
sub$group[sub$tissue == "Bowel" & sub$condition == "Sham"] <- "Sham Bowel"
sub$group[sub$tissue == "Bowel" & sub$condition == "Colotomy"] <- "Colotomy Bowel"

Idents(sub) <- "group"
dotplot <- DotPlot(sub, features = c("Components1", "Yapsig1", "ECMgenes1"), 
                   cols = "RdYlBu", split.by = "fibro.annotation", scale.min = 0, scale.max = 100) 
#Rearrange the data 
dotplot$data$id <- factor(dotplot$data$id, levels = c("Sham Bowel_Cthrc1+ mFib", "Colotomy Bowel_Cthrc1+ mFib", 
                                                      "Sham MAT_Cthrc1+ mFib", "Colotomy MAT_Cthrc1+ mFib"))

dotplot +
  scale_y_discrete(labels = c('Sham', 'Colotomy', 'Sham', 'Colotomy')) +
  scale_x_discrete(labels = c('YAP\nPathway', 'YAP\nTargets', 'ECM')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggsave("Figures/cthrc1_activation.pdf")


#### Global YAP signature ####
global <- AddModuleScore(global, features = list(Yapgenes), name = "Yapsig")
global <- AddModuleScore(global, features = components, name = "Components")
global <- AddModuleScore(global, features = ecmgenes,name = "ECM")

revglobal <- global 
revglobal$global.annotation <- factor(x = revglobal$global.annotation, levels = rev(c("Fibroblasts", "Smooth Muscle Cells", "Vascular ECs", "Lymphatic ECs", "Mesothelial Cells", 
                                                                                      "B Cells", "Plasma Cells", "T Cells", "Macrophages", "Myeloid Cells", "Proliferating Cells", 
                                                                                      "Dendritic Cells", "Plasmacytoid DCs", 
                                                                                      "Epithelial Cells", "Tuft Cells", "Erythrocytes")))
Idents(revglobal) <- "global.annotation"


DotPlot(revglobal, features = c("ECM1", "Components1", "Yapsig1"), group.by = "global.annotation", scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
  scale_x_discrete(labels = c("ECM", "YAP\n Targets", "YAP\n Pathway", "YAP\n Targets")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_text(size = 10, hjust = 1, vjust = 1, angle = 45),
    axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
    legend.key.size = unit(0.8, 'cm'),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )
ggsave("Figures/Yapsig_global.pdf")






# sub <- subset(fibro, idents = c("Cthrc1+ mFib", "Pi16+ ssFib"))
# sub$group <- NA
# sub$group[sub$tissue == "MAT" & sub$condition == "Sham"] <- "Sham MAT"
# sub$group[sub$tissue == "MAT" & sub$condition == "Colotomy"] <- "Colotomy MAT"
# sub$group[sub$tissue == "Bowel" & sub$condition == "Sham"] <- "Sham Bowel"
# sub$group[sub$tissue == "Bowel" & sub$condition == "Colotomy"] <- "Colotomy Bowel"
# 
# Idents(sub) <- "group"
# dotplot <- DotPlot(sub, features = c("Components1", "Yapsig1", "ECMgenes1"), 
#                    cols = "RdYlBu", split.by = "fibro.annotation", scale.min = 0, scale.max = 100) 
# #Rearrange the data 
# dotplot$data$id <- factor(dotplot$data$id, levels = c("Sham Bowel_Pi16+ ssFib", "Colotomy Bowel_Pi16+ ssFib", 
#                                                       "Sham MAT_Pi16+ ssFib", "Colotomy MAT_Pi16+ ssFib", 
#                                                       "Sham Bowel_Cthrc1+ mFib", "Colotomy Bowel_Cthrc1+ mFib", 
#                                                       "Sham MAT_Cthrc1+ mFib", "Colotomy MAT_Cthrc1+ mFib"))
# 
# dotplot + scale_y_discrete(labels=c('Sham', "Colotomy", 
#                                     'Sham', "Colotomy", 
#                                     'Sham', "Colotomy",
#                                     'Sham', "Colotomy")) + 
#   scale_x_discrete(labels = c("YAP Component Genes", "YAP Target Genes", "ECM Genes")) + 
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave("SuppFigures/YAPactivationinPOSTN.svg")

#Human mouse integration
library(patchwork)
#Anchor label transfer 
rm(list = ls())
human = readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman/finalfibroblastannotations.rds")
mouse = readRDS("mousefibro.annotated.rds")

DefaultAssay(human) = 'RNA'
temp2 = convertHumanToMouse(human)
saveRDS(temp2, "MouseHumanInt/mousetohumangenes.rds")

temp2 <- readRDS("MouseHumanInt/mousetohumangenes.rds")
data = human[["RNA"]]$counts
ii_toUse = !is.na(temp2)  #remove NAs
data_toUse = data[ ii_toUse, ]
rownames( data_toUse ) = temp2[ ii_toUse ]
human[[ 'ortholog_mouse' ]] = CreateAssayObject(data_toUse)
DefaultAssay(human) <- "ortholog_mouse"

# #Rerun NormalizeData and FVSs using the new "ortholog_mouse" assay
# human[["ortholog_mouse"]] <- split(human[["ortholog_mouse"]], f = human$orig.ident) 
# human <- NormalizeData(human)
# human <- FindVariableFeatures(human) #Find variable features in the new assay
# saveRDS(human, "MouseHumanInt/humanVFs.rds", compress = F)

human <- SCTransform(human, verbose = T, assay='ortholog_mouse')
mouse <- SCTransform(mouse, verbose = T, assay='RNA') #originalexp
DefaultAssay(mouse) = "SCT"
DefaultAssay(human) = "SCT"
saveRDS(human, "MouseHumanInt/SCThuman.rds", compress = F)
saveRDS(mouse, "MouseHumanInt/SCTmouse.rds", compress = F)


#Set the reference as human to map mouse onto human
rm(list = ls())

reference = readRDS("MouseHumanInt/SCThuman.rds")
query = readRDS("MouseHumanInt/SCTmouse.rds")

reference = readRDS("MouseHumanInt/SCThuman.rds")
query = readRDS("MouseHumanInt/SCTmouse.rds")

anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT",
                               reference.reduction = "pca", dims = 1:30, k.anchor = 15)


saveRDS(anchors, "MouseHumanInt/humanmouseanchorsVFs_k30.rds", compress = F)
saveRDS(anchors, "MouseHumanInt/humanmouseanchorsSCT_k15.rds", compress = F)

predictions <- TransferData(anchorset = anchors, refdata = reference$fibro.annotation, dims = 1:30, prediction.assay = T)
saveRDS(predictions, "MouseHumanInt/predictionsSCT_k15.rds")

query <- AddMetaData(query, metadata = predictions)
query$prediction.match <- query$predicted.id == query$fibro.annotation
query[["predictions"]] <- predictions
DefaultAssay(query) <- "predictions"
FeaturePlot(query, features = c("mCDF CTHRC1+"))

p1 = DimPlot(query, reduction = "umap.harmony", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(query, reduction = "umap.harmony", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

FeaturePlot(query, features = c("mCDF CTHRC1+"),  reduction = "umap.harmony", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))

predictions <- TransferData(anchorset = anchors, refdata = mouse$fibro.annotation, dims = 1:30, prediction.assay = T)
human[["predictions"]] <- predictions
DefaultAssay(human) <- "predictions"
FeaturePlot(human, features = c("mCDF CTHRC1+"))


human <- AddMetaData(human, metadata = predictions)
human$predicted.id <- factor(x = human$predicted.id, levels = rev(c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
                                                                    "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
                                                                    "mCDF LRRC7+", "mCDF CTHRC1+")))

Idents(mouse) = "predicted.id"
num_colors <- length(levels(mouse))
colors <- brewer.pal(num_colors, "RdYlBu")
dimplot <- DimPlot(mouse, cols = colors) + ggtitle("Predicted Metadata Clusters")
featureplot1 <- FeaturePlot(mouse, features = c("mCDF CTHRC1+")) + ggtitle("mCDF CTHRC1+") + guides(color = guide_legend(title = "Probability")) + 
  theme(legend.title =element_text(size=10)) + 
  scale_colour_gradientn(name = "Score", colors = c("#0d0887", "#cc4778", "#f0f921"))
featureplot2 <-FeaturePlot(mouse, features = c("mCDF LRRC7+")) + ggtitle("mCDF LRRC7+") + guides(color = guide_legend(title = "Probability")) + 
  theme(legend.title =element_text(size=10)) + 
  scale_colour_gradientn(name = "Score", colors = c("#0d0887", "#cc4778", "#f0f921"))
dimplot + featureplot1 + featureplot2 + plot_layout(guides = "collect")







