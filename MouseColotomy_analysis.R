setwd("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy/r_packages",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers)
library(scales); library(SeuratObject)
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

#### Downstream processing of global object ####
rm(list = ls())
merged <- readRDS("processed/merged.global.processed.rds")
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$orig.ident)
# run standard anlaysis workflow
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged, features = rownames(merged))
merged <- RunPCA(merged, npcs = 30, verbose = T)
merged <- IntegrateLayers(merged, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
merged <- JoinLayers(merged)
merged <- FindNeighbors(merged, dims = 1:30, reduction = "pca")
for (i in seq(0.01, 0.5, 0.01)) { 
  merged <- FindClusters(merged, resolution = i)
  DimPlot(merged, raster = TRUE, label = TRUE)
  ggsave(paste0("sketch_v5_all/global_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(merged, paste0("sketch_v5_all/global_harmony_umap_res",i,".rds"), compress = FALSE)
}

merged <- JoinLayers(merged)
markers <- FindAllMarkers(merged, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/mouseMerged_global_harmony_umap_res0.21_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/mouseMerged_global_harmony_umap_res0.21_topmarkers.csv")

merged <- readRDS("sketch_v5_all/global_harmony_umap_res0.21.rds")
merged$global.annotation <- NA
merged$global.annotation[merged$seurat_clusters == 0] <- "B Cells"
merged$global.annotation[merged$seurat_clusters == 1] <- "Fibroblasts"
merged$global.annotation[merged$seurat_clusters == 2] <- "T Cells"
merged$global.annotation[merged$seurat_clusters == 3] <- "Macrophages"
merged$global.annotation[merged$seurat_clusters == 4] <- "Mki67+ Cells"
merged$global.annotation[merged$seurat_clusters == 5] <- "Epithelial Cells"
merged$global.annotation[merged$seurat_clusters == 6] <- "Plasma Cells"
merged$global.annotation[merged$seurat_clusters == 7] <- "Erythrocytes"
merged$global.annotation[merged$seurat_clusters == 8] <- "Vascular ECs"
merged$global.annotation[merged$seurat_clusters == 9] <- "Mesothelial Cells"
merged$global.annotation[merged$seurat_clusters == 10] <- "Dendritic Cells"
merged$global.annotation[merged$seurat_clusters == 11] <- "Macrophages"
merged$global.annotation[merged$seurat_clusters == 12] <- "Myeloid Cells"
merged$global.annotation[merged$seurat_clusters == 13] <- "Smooth Muscle Cells"
merged$global.annotation[merged$seurat_clusters == 14] <- "Lymphatic ECs"
merged$global.annotation[merged$seurat_clusters == 15] <- "Plasmacytoid DCs"
merged$global.annotation[merged$seurat_clusters == 16] <- "Tuft Cells"
Idents(merged) <- "global.annotation"
colors = c("#BE0032", "#875692", "#E68FAC", "#8DB600", "#0067A5", 
                    "#F38400", "#F7E1A0", "#A1CAF1", "#C2B280","#DCD300",
                    "#B5EFB5", "#7ED7D1",
                    "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", 
"#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", "#5A5156", "#E4E1E3", 
                    "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
                    "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
                    "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6",
                    "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79",
                    "#66B0FF", "#3B00FB")
                    
merged$global.annotation <- factor(x = merged$global.annotation, levels = c("Fibroblasts", "Smooth Muscle Cells", "Vascular ECs", "Lymphatic ECs", "Mesothelial Cells", 
                                                                            "B Cells", "Plasma Cells", "T Cells", "Macrophages", "Myeloid Cells", "Mki67+ Cells", 
                                                                            "Dendritic Cells", "Plasmacytoid DCs", 
                                                                            "Epithelial Cells", "Tuft Cells", "Erythrocytes"))
Idents(merged) <- "global.annotation"
DimPlot(merged, cols = colors)
ggsave("Figures/globalumap.jpg", width = 6.22, height = 3.93)
saveRDS(merged, "mouseMerged_global_harmony_umap_res0.1.ANNOTATEDTRUE.rds")

#Dotplot of major markers 
global <- merged 
revglobal <- global 
revglobal$global.annotation <- factor(x = revglobal$global.annotation, levels = rev(c("Fibroblasts", "Smooth Muscle Cells", "Vascular ECs", "Lymphatic ECs", "Mesothelial Cells", 
                                                                            "B Cells", "Plasma Cells", "T Cells", "Macrophages", "Myeloid Cells", "Mki67+ Cells", 
                                                                            "Dendritic Cells", "Plasmacytoid DCs", 
                                                                            "Epithelial Cells", "Tuft Cells", "Erythrocytes")))
Idents(revglobal) <- "global.annotation"
genes = c("Pdgfra", "Des", "Pecam1", "Prox1", "Msln", "Bank1", "Jchain", "Cd3e", 
          "Cd68", "S100a8", "Mki67", "Flt3", "Tcf4", "Epcam", "Dclk1", "Hba-a1")

DotPlot(revglobal, features = genes, cols = "RdYlBu") +   theme(axis.title.y = element_blank(), 
                                                                         axis.title.x = element_blank(),
                                                                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                         legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                         legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                         legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                         legend.title = element_text(size=6), #change legend title font size
                                                                         legend.text = element_text(size=6)) #change legend text font size)
ggsave("SuppFigures/globalmarkers.svg")


#### Subclustering fibroblasts Round 1 #### 
rm(list = ls())
global <- readRDS("mouseMerged_global_harmony_umap_res0.1.ANNOTATEDTRUE.rds")
Idents(global) <- "global.annotation"
fibro <- subset(global, idents = "Fibroblasts")

DefaultAssay(fibro) <- "RNA"
fibro <- subset(global, idents = "Fibroblasts")
# re-generate seurat fibroect and remove extraneous data
options(Seurat.fibroect.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( fibro, 'newFibro1/fibro1.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro1/fibro1_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro1/fibro1_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
saveRDS( fibro, 'newFibro1/fibro1_counts.rds', compress = F )

# run downstream analysis
fibro <- NormalizeData(fibro)
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
saveRDS( fibro, 'newFibro1/fibro1_split.rds', compress = F )
fibro <- FindVariableFeatures(fibro, verbose = T )
saveRDS( fibro, 'newFibro1/fibro1_var.rds'  , compress = F )
fibro <- ScaleData(fibro, features = rownames(fibro), verbose = T)
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
saveRDS( fibro, 'newFibro1/fibro1_pca.rds'  , compress = F )
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( fibro, 'newFibro1/fibro1_harmony.rds', compress = F ) 
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( fibro, 'newFibro1/fibro1_harmony_umap.rds', compress = F ) 
fibro <- JoinLayers(fibro)
saveRDS( fibro, 'newFibro1/fibro1_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(fibro, features = c("Pecam1", "Mki67", "Ptprc", "Plp1","Igf2", "Dach1", "Pdgfra", "Dpt", "Epcam", "Krt19", "Kcnj8", "Rgs5", "Myh11"), order = T)
ggsave("newFibro1/fibro1_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

fibro <- FindNeighbors(fibro, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 0.5, 0.1) ) {
  fibro <- FindClusters(fibro, resolution = res )
  DimPlot(fibro)
  ggsave(paste0("newFibro1/fibro1_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( fibro, paste0( 'newFibro1/fibro1_harmony_clusters_', res, '.rds' ), compress = F ) 
}

fibro <- JoinLayers(fibro)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro1/mousefibro_fibro_harmony_umap_res0.5_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro1/mousefibro_fibro_harmony_umap_res0.5_topmarkers.csv")

#### Round 2 reclustering of fibroblasts ####

#Remove low quality clusters/ contaminants and re-perform clustering analysis 
Idents(fibro) <- "seurat_clusters"
fibro <- subset(fibro, idents = c(2, 4, 5, 9, 10, 11), invert = T)

# re-generate seurat object and remove extraneous data
options(Seurat.fibroect.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( fibro, 'newFibro2/fibro2.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro2/fibro2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro2/fibro2_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
saveRDS( fibro, 'newFibro2/fibro2_counts.rds', compress = F )

# run downstream analysis
fibro <- NormalizeData(fibro)
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
saveRDS( fibro, 'newFibro2/fibro2_split.rds', compress = F )
fibro <- FindVariableFeatures(fibro, verbose = T )
saveRDS( fibro, 'newFibro2/fibro2_var.rds'  , compress = F )
fibro <- ScaleData(fibro, features = rownames(fibro), verbose = T)
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
saveRDS( fibro, 'newFibro2/fibro2_pca.rds'  , compress = F )
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( fibro, 'newFibro2/fibro2_harmony.rds', compress = F ) 
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( fibro, 'newFibro2/fibro2_harmony_umap.rds', compress = F ) 
fibro <- JoinLayers(fibro)
saveRDS( fibro, 'newFibro2/fibro2_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(fibro, features = c("Pecam1", "Mki67", "Ptprc", "Plp1","Igf2", "Dach1", "Pdgfra", "Dpt", "Epcam", "Krt19", "Kcnj8", "Rgs5", "Myh11"), order = T)
ggsave("newFibro2/fibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

fibro <- FindNeighbors(fibro, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1.0, 0.1) ) {
  fibro <- FindClusters(fibro, resolution = res )
  DimPlot(fibro)
  ggsave(paste0("newFibro2/fibro2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( fibro, paste0( 'newFibro2/fibro2_harmony_clusters_', res, '.rds' ), compress = F ) 
}

fibro <- JoinLayers(fibro)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro2/mousefibro2_fibro_harmony_umap_res0.9_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro2/mousefibro2_fibro_harmony_umap_res0.9_topmarkers.csv")

#### Fibroblast clustering round 3 ####
#Remove low quality clusters/ contaminants and re-perform clustering analysis 
Idents(fibro) <- "seurat_clusters"
fibro <- subset(fibro, idents = c(4,9), invert = T)

# re-generate seurat object and remove extraneous data
options(Seurat.fibroect.assay.version = "v5")
fibro <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( fibro, 'newFibro3/fibro3.rds', compress = F )

# re-perform BPcells
counts.mat = fibro[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro3/fibro3_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro3/fibro3_counts" )
fibro = CreateSeuratObject( counts.mat, meta.data = fibro@meta.data )
saveRDS( fibro, 'newFibro3/fibro3_counts.rds', compress = F )

# run downstream analysis
fibro <- NormalizeData(fibro)
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
saveRDS( fibro, 'newFibro3/fibro3_split.rds', compress = F )
fibro <- FindVariableFeatures(fibro, verbose = T )
saveRDS( fibro, 'newFibro3/fibro3_var.rds'  , compress = F )
fibro <- ScaleData(fibro, features = rownames(fibro), verbose = T)
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
saveRDS( fibro, 'newFibro3/fibro3_pca.rds'  , compress = F )
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( fibro, 'newFibro3/fibro3_harmony.rds', compress = F ) 
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( fibro, 'newFibro3/fibro3_harmony_umap.rds', compress = F ) 
fibro <- JoinLayers(fibro)
saveRDS( fibro, 'newFibro3/fibro3_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(fibro, features = c("Pecam1", "Mki67", "Ptprc", "Plp1","Igf2", "Dach1", "Pdgfra", "Dpt", "Epcam", "Krt19", "Kcnj8", "Rgs5", "Myh11"), order = T)
ggsave("newFibro3/fibro3_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

fibro <- FindNeighbors(fibro, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1.0, 0.1) ) {
  fibro <- FindClusters(fibro, resolution = res )
  DimPlot(fibro)
  ggsave(paste0("newFibro3/fibro3_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( fibro, paste0( 'newFibro3/fibro3_harmony_clusters_', res, '.rds' ), compress = F ) 
}

fibro <- JoinLayers(fibro)
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro3/mousefibro3_fibro_harmony_umap_res0.7_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro3/mousefibro3_fibro_harmony_umap_res0.7_topmarkers.csv")

#Annotations
fibro$fibro.annotation <- NA
fibro$fibro.annotation[fibro$seurat_clusters == 0] <- "Mechanosensitive Postn+"
fibro$fibro.annotation[fibro$seurat_clusters == 1] <- "Steady-State Pi16+"
fibro$fibro.annotation[fibro$seurat_clusters == 2] <- "Mechanosensitive Postn+"
fibro$fibro.annotation[fibro$seurat_clusters == 3] <- "Steady-State Pi16+"
fibro$fibro.annotation[fibro$seurat_clusters == 4] <- "Steady-State Pi16+"
fibro$fibro.annotation[fibro$seurat_clusters == 5] <- "Mechanosensitive Postn+"
fibro$fibro.annotation[fibro$seurat_clusters == 6] <- "Inflammatory Adamdec1+"
fibro$fibro.annotation[fibro$seurat_clusters == 7] <- "Steady-State Abca8a+"
fibro$fibro.annotation[fibro$seurat_clusters == 8] <- "Mechanosensitive Postn+"
fibro$fibro.annotation[fibro$seurat_clusters == 9] <- "Steady-State Mki67+"
fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = c("Mechanosensitive Postn+", "Steady-State Abca8a+", "Steady-State Mki67+",
                                                                        "Inflammatory Adamdec1+", "Steady-State Pi16+"))
Idents(fibro) <- "fibro.annotation"
colors = c("#B2182B", "#ef8a62", "grey92", "#d1e5f0", "#67a9cf", "#2166ac")
dimplot <- DimPlot(fibro, cols = colors ) + theme(legend.key.size = unit(0.5, "in"))
levels(fibro) <- rev(c("Mechanosensitive Postn+", "Steady-State Abca8a+", "Steady-State Mki67+",
                   "Inflammatory Adamdec1+", "Steady-State Pi16+"))
dotplot <- DotPlot(fibro, features = c("Postn", "Abca8a", "Mki67", "Adamdec1", "Pi16"), cols = "RdBu") + 
  theme(axis.title.y = element_blank(), 
      axis.title.x = element_blank(),
      axis.text.y = element_blank(), 
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      legend.key.size = unit(0.4, 'cm'), #change legend key size
      legend.key.height = unit(0.2, 'cm'), #change legend key height
      legend.key.width = unit(0.2, 'cm'), #change legend key width
      legend.title = element_text(size=6), #change legend title font size
      legend.text = element_text(size=6)) #change legend text font size)

combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -35, unit = "pt")))) + dotplot 
combined_plot
ggsave(paste0("Figures/mousefibroblastUMAP.jpg"), width = 10, height = 3.7, units = "in")
saveRDS(fibro, "mousefibro.annotated.rds")

#Cluster markers 
genes = c("Postn", "Cthrc1", "Sparc", "Wnt4", "Clu", 
          "Abca8a", "Cp", "Dpep1", "Lpl", "Cxcl12",
          "Mki67", "Pclaf", "Top2a", "Coro1a", "Stmn1",
          "Adamdec1", "Bmp4", "Bmp5", "Grem1", "Fbln1", 
          "Pi16", "Mfap5", "Dpp4", "Cd55", "Cd248")
DotPlot(fibro, features = genes , cols = "RdYlBu") + theme(axis.title.y = element_blank(), 
                                                           axis.title.x = element_blank(),
                                                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                           legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                           legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                           legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                           legend.title = element_text(size=6), #change legend title font size
                                                           legend.text = element_text(size=6)) #change legend text font size)


#Fibrotic Cluster markers 
genes = c("Timp1", "Cdh11", "Bgn", "Fap", "Sfrp1", "Runx1")
DotPlot(fibro, features = genes , cols = "RdYlBu", scale.min = 0, scale.max = 100) + theme(axis.title.y = element_blank(), 
                                                           axis.title.x = element_blank(),
                                                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                           legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                           legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                           legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                           legend.title = element_text(size=6), #change legend title font size
                                                           legend.text = element_text(size=6)) #change legend text font size)

#Representative signature genes 
genes = c("Ccn1", "Ccn2", "Serpine1", "Tgfb2",
          "Ptk2", "Tead1", "Tead2", 
          "Col1a1", "Col1a2", "Col3a1")
DotPlot(fibro, features = genes , cols = "RdYlBu", scale.min = 0, scale.max = 100) + theme(axis.title.y = element_blank(), 
                                                                                           axis.title.x = element_blank(),
                                                                                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                           legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                           legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                           legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                           legend.title = element_text(size=6), #change legend title font size
                                                                                           legend.text = element_text(size=6)) #change legend text font size)

ggsave("SuppFigures/Mousesiggenes.svg")


#### Mouse YAP signatures ####
#1. YAP gene signature 
Yapgenes <- c("Ccn1", "Ccn2", "Foxf2", "Igfbp3", "Ccdc80", "F3", "Myof", "Fjx1", "Rbms3",
                            "Ptpn14", "Ankrd1", "Serpine1", "Lats2", "Arhgef17", "Nuak2", "Amotl2", 
                            "Dock5", "Crim1", "Nt5e", "Asap1", "Gadd45a", "Tgfb2", "Axl") #see AKH Loe et al. 


Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
fibro <- AddModuleScore(fibro, features = Yapsub, name = "Yapsig")

#2. YAP components signature
components <- list(c("Yap1", "Ctnnb1", 
                     "Ptk2", "Itga5", "Itgav", "Itgb1",  
                     "Sav1", "Stk3", "Lats1", "Lats2", 
                     "Tead1", "Tead2", "Tead3", "Tead4", "Runx1", 
                     "Vim", "Vcl", "Vasp"))
fibro <- AddModuleScore(fibro, features = components, name = "Components")
#3. ECM protein signature 
ecmgenes <- list(c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col24a1", 
              "Lum", "Dcn", "Bgn", "Hspg2", "Agrn", 
              "Fn1", "Vtn", "Tgfbi", "Nid1", "Nid2", "Lamc1", "Lama1", "Lama2", "Sparc", 
              "Fbln1", "Fbln2", "Fbln5", "Ltbp1", "Ltbp5", "Emilin1", "Mfap4", "Efemp1", "Fbn1", "Fbn2", "Tnc"))
fibro <- AddModuleScore(fibro, features = ecmgenes, name = "ECMgenes")
Idents(fibro) <- "fibro.annotation"
fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = rev(c("Mechanosensitive Postn+", "Steady-State Abca8a+", "Steady-State Mki67+",
                       "Inflammatory Adamdec1+", "Steady-State Pi16+")))
Idents(fibro) <- "fibro.annotation"
DotPlot(fibro, features = c("Components1", "Yapsig1", "ECMgenes1"), scale.max = 100, scale.min = 0, cols = "RdYlBu") + 
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
ggsave("Figures/YAPsignatures.svg")

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



#Generate proptable and barplots of fibroblasts in mouse
tab = table(fibro_small$group, fibro_small$fibro.annotation)
ptab = prop.table(tab, 1)*100
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
  labs(y = "Percent (%)", 
       fill = "Cluster") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold")) + 
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_discrete(labels=c('Sham', 'Colotomy', 'Sham', 'Colotomy')) + 
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold", margin = margin(0, 0, 5, 0)),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),  # Increasing legend text size
        legend.key.size = unit(1.5, "line")) +  
  theme(panel.spacing = unit(0, "lines"))
ggsave("Figures/proptablemousefibroblasts.jpg", width = 7.5, height = 3.89 )

#Plot YAP target gene activation 
sub <- subset(fibro, idents = c("Mechanosensitive Postn+", "Steady-State Pi16+"))
sub$group <- NA
sub$group[sub$tissue == "MAT" & sub$condition == "Sham"] <- "Sham MAT"
sub$group[sub$tissue == "MAT" & sub$condition == "Colotomy"] <- "Colotomy MAT"
sub$group[sub$tissue == "Bowel" & sub$condition == "Sham"] <- "Sham Bowel"
sub$group[sub$tissue == "Bowel" & sub$condition == "Colotomy"] <- "Colotomy Bowel"

Idents(sub) <- "group"
dotplot <- DotPlot(sub, features = c("Components1", "Yapsig1", "ECMgenes1"), 
                   cols = "RdYlBu", split.by = "fibro.annotation", scale.min = 0, scale.max = 100) 
#Rearrange the data 
dotplot$data$id <- factor(dotplot$data$id, levels = c("Sham Bowel_Steady-State Pi16+", "Colotomy Bowel_Steady-State Pi16+", 
                                                      "Sham MAT_Steady-State Pi16+", "Colotomy MAT_Steady-State Pi16+", 
                                                      "Sham Bowel_Mechanosensitive Postn+", "Colotomy Bowel_Mechanosensitive Postn+", 
                                                      "Sham MAT_Mechanosensitive Postn+", "Colotomy MAT_Mechanosensitive Postn+"))

dotplot1 <- dotplot[dotplot$data$id == "Sham Bowel_Steady-State Pi16+" | 
                      dotplot$data$id == "Colotomy Bowel_Steady-State Pi16+" |
                      dotplot$data$id == "Sham MAT_Steady-State Pi16+" |
                      dotplot$data$id == "Colotomy MAT_Steady-State Pi16+"]

dotplot + scale_y_discrete(labels=c('Sham', "Colotomy", 
                                    'Sham', "Colotomy", 
                                    'Sham', "Colotomy",
                                    'Sham', "Colotomy")) + 
  scale_x_discrete(labels = c("YAP Component Genes", "YAP Target Genes", "ECM Genes")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  
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




















#### SCT Harmony ####
rm(list = ls())
global <- readRDS("Shruti_global_labeled_cells.RDS")
global <- UpdateSeuratObject(global) #make object compatible with Seurat V5
Idents(global) <- "final_clusters"
fibro <- subset(global, idents = "Fibroblasts")
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- SCTransform(fibro, assay = "RNA") #Fails if matrixStats is not 1.1.0. Make sure you downgrade. 
DefaultAssay(fibro) <- "SCT"
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
saveRDS( fibro, 'newFibro1/SCT_CCA/SCT_CCA_fibro1_pca.rds'  , compress = F )
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( fibro, 'newFibro1/SCT_CCA/SCT_CCA_fibro1_harmony.rds', compress = F ) 
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( fibro, 'newFibro1/SCT_CCA/SCT_CCA_fibro1_harmony_umap.rds', compress = F ) 
fibro <- JoinLayers(fibro)
saveRDS( fibro, 'newFibro1/SCT_CCA/SCT_CCA_fibro1_harmony_umap_joined.rds', compress = F ) 

#### Setup immune subsets ####
global <- readRDS("mousemyeloid_global_harmony_umap_res0.1.ANNOTATEDTRUE.rds")
Idents(global) <- "global.annotation"
lymphoid <- subset(global, idents = c("T Cells", "B Cells", "Plasma Cells", "Proliferating Cells"))
myeloid <- subset(global, idents = c("Macrophages", "Myeloid Cells", "Dendritic Cells", "Plasmacytoid DCs"))

#### Lymphoid ####
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
lymphoid <- CreateSeuratObject(counts = lymphoid[["RNA"]]$counts, meta.data = lymphoid@meta.data)
saveRDS( lymphoid, "Lymphoid/Lymphoid.rds", compress = F )

# re-perform BPcells
counts.mat = lymphoid[["RNA"]]$counts
#counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as(counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "Lymphoid/LymphoidNew_counts", overwrite = T )
counts.mat <- open_matrix_dir(dir = "Lymphoid/LymphoidNew_counts" )
lymphoid = CreateSeuratObject( counts.mat, meta.data = lymphoid@meta.data )
saveRDS( lymphoid, "Lymphoid/LymphoidNew_counts.rds", compress = F )

# run downstream analysis
lymphoid <- NormalizeData(lymphoid)
lymphoid[["RNA"]] <- split(lymphoid[["RNA"]], f = lymphoid$orig.ident)
saveRDS( lymphoid, "Lymphoid/LymphoidNew_split.rds", compress = F )
lymphoid <- FindVariableFeatures(lymphoid, verbose = T )
saveRDS( lymphoid, "Lymphoid/LymphoidNew_var.rds"  , compress = F )
lymphoid <- ScaleData(lymphoid, verbose = T)
lymphoid <- RunPCA(lymphoid, npcs = 30, verbose = T)
saveRDS( lymphoid, "Lymphoid/LymphoidNew_pca.rds"  , compress = F )
lymphoid <- IntegrateLayers(lymphoid, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( lymphoid, "Lymphoid/LymphoidNew_harmony.rds", compress = F ) 
lymphoid <- RunUMAP(lymphoid, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( lymphoid, "Lymphoid/LymphoidNew_harmony_umap.rds", compress = F ) 
lymphoid <- JoinLayers(lymphoid)
saveRDS( lymphoid, "Lymphoid/LymphoidNew_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(lymphoid, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("Lymphoid/LymphoidNew_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
lymphoid <- FindNeighbors(lymphoid, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  lymphoid <- FindClusters(lymphoid, resolution = res )
  DimPlot(lymphoid, raster = TRUE, label = T, repel = T)
  ggsave(paste0("Lymphoid/LymphoidNew_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( lymphoid, paste0( "Lymphoid/LymphoidNew_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  lymphoid[["RNA"]]$data <- as(object = lymphoid[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(lymphoid, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
  write.csv(markers, paste0("Lymphoid/LymphoidNew_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("Lymphoid/LymphoidNew_harmony_clusters_topmarkers_",res,".csv"))
}

DefaultAssay(lymphoid) <- "RNA"

##Lymphoid Round 2
lymphoid <- readRDS("Lymphoid/LymphoidNew_harmony_clusters_0.6.rds")
lymphoid <- myeloidset(lymphoid, idents = c(1, 8, 12), invert = T)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
lymphoid <- CreateSeuratObject(counts = lymphoid[["RNA"]]$counts, meta.data = lymphoid@meta.data)
saveRDS( lymphoid, "Lymphoid2/Lymphoid2.rds", compress = F )

# re-perform BPcells
counts.mat = lymphoid[["RNA"]]$counts
#counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as(counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "Lymphoid2/Lymphoid2New_counts", overwrite = T )
counts.mat <- open_matrix_dir(dir = "Lymphoid2/Lymphoid2New_counts" )
lymphoid = CreateSeuratObject( counts.mat, meta.data = lymphoid@meta.data )
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_counts.rds", compress = F )

# run downstream analysis
lymphoid <- NormalizeData(lymphoid)
lymphoid[["RNA"]] <- split(lymphoid[["RNA"]], f = lymphoid$orig.ident)
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_split.rds", compress = F )
lymphoid <- FindVariableFeatures(lymphoid, verbose = T )
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_var.rds"  , compress = F )
lymphoid <- ScaleData(lymphoid, verbose = T)
lymphoid <- RunPCA(lymphoid, npcs = 30, verbose = T)
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_pca.rds"  , compress = F )
lymphoid <- IntegrateLayers(lymphoid, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_harmony.rds", compress = F ) 
lymphoid <- RunUMAP(lymphoid, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_harmony_umap.rds", compress = F ) 
lymphoid <- JoinLayers(lymphoid)
saveRDS( lymphoid, "Lymphoid2/Lymphoid2New_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(lymphoid, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("Lymphoid2/Lymphoid2New_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
lymphoid <- FindNeighbors(lymphoid, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  lymphoid <- FindClusters(lymphoid, resolution = res )
  DimPlot(lymphoid, raster = TRUE, label = T, repel = T)
  ggsave(paste0("Lymphoid2/Lymphoid2New_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( lymphoid, paste0( "Lymphoid2/Lymphoid2New_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  lymphoid[["RNA"]]$data <- as(object = lymphoid[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(lymphoid, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
  write.csv(markers, paste0("Lymphoid2/Lymphoid2New_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("Lymphoid2/Lymphoid2New_harmony_clusters_topmarkers_",res,".csv"))
}

DefaultAssay(lymphoid) <- "RNA"

lymphoid <- readRDS("Lymphoid2/Lymphoid2_harmony_clusters_0.5.rds")

lymphoid$immune.annotation <- NA
lymphoid$immune.annotation[lymphoid$seurat_clusters == 0] <- "Resting B Cells Ighd1+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 1] <- "CD8 T Cells Cd8b1+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 2] <- "Cycling B Cells Aicda+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 3] <- "Resting B Cells Ighd1+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 4] <- "Plasma Cells Igha+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 5] <- "NKT Cells Ccl5+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 6] <- "Treg Foxp3+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 7] <- "Cycling Plasma Cells Jchain+"
lymphoid$immune.annotation[lymphoid$seurat_clusters == 8] <- "CD4 T Cells Gata3+"

Idents(lymphoid) <- "immune.annotation"
lymphoid$immune.annotation <- factor(lymphoid$immune.annotation, levels = c("Resting B Cells Ighd1+", "Cycling B Cells Aicda+", "Plasma Cells Igha+", "Cycling Plasma Cells Jchain+", 
                                                                            "CD8 T Cells Cd8b1+", "NKT Cells Ccl5+", "Treg Foxp3+", "CD4 T Cells Gata3+"))
Idents(lymphoid) <- "immune.annotation"
colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
                    "#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
                    
dimplot <- DimPlot(lymphoid, cols = colors)                             

revlymphoid <- lymphoid
revlymphoid$immune.annotation <- factor(revlymphoid$immune.annotation, levels = rev(c("Resting B Cells Ighd1+", "Cycling B Cells Aicda+", "Plasma Cells Igha+", "Cycling Plasma Cells Jchain+", 
                                                                                      "CD8 T Cells Cd8b1+", "NKT Cells Ccl5+", "Treg Foxp3+", "CD4 T Cells Gata3+")))
Idents(revlymphoid) <- "immune.annotation"

genes <- c("Ighd",  "Aicda", "Igha", "Jchain", "Cd8b1", "Ccl5", "Foxp3", "Gata3", "Cd19", "Mki67", "Cd3d")

dotplot <- DotPlot(revlymphoid, features = genes, cols = "RdYlBu", scale.min = 0, scale.max = 100, dot.scale = 5) +   theme(axis.title.y = element_blank(), 
                                                                                                                            axis.title.x = element_blank(),
                                                                                                                            axis.text.y = element_blank(), 
                                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                                                            legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                                                            legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                                                            legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                                                            legend.title = element_text(size=6), #change legend title font size
                                                                                                                            legend.text = element_text(size=6)) #change legend text font size)
combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -20, unit = "pt")))) + dotplot 
combined_plot
ggsave("SuppFigures/mouselymphoid.svg")


##Prop tables 
lymphoid$group <- NA
lymphoid$group[lymphoid$tissue == "MAT" & lymphoid$condition == "Sham"] <- "Sham MAT"
lymphoid$group[lymphoid$tissue == "MAT" & lymphoid$condition == "Colotomy"] <- "Colotomy MAT"
lymphoid$group[lymphoid$tissue == "Bowel" & lymphoid$condition == "Sham"] <- "Sham Bowel"
lymphoid$group[lymphoid$tissue == "Bowel" & lymphoid$condition == "Colotomy"] <- "Colotomy Bowel"
Idents(lymphoid) <- "group"

#Generate proptable and barplots of fibroblasts in mouse
tab = table(lymphoid$group, lymphoid$immune.annotation)
ptab = prop.table(tab, 1)*100
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
  labs(y = "Percent (%)", 
       fill = "Cluster") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold")) + 
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_discrete(labels=c('Sham', 'Colotomy', 'Sham', 'Colotomy')) + 
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold", margin = margin(0, 0, 5, 0)),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),  # Increasing legend text size
        legend.key.size = unit(1.5, "line")) +  
  theme(panel.spacing = unit(0, "lines"))
ggsave("SuppFigures/proptablemouselymphoid.svg")



#### Myeloid ####

# re-perform BPcells
counts.mat = myeloid[["RNA"]]$counts
#counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as(counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "Myeloid/Myeloid_counts", overwrite = T )
counts.mat <- open_matrix_dir(dir = "Myeloid/Myeloid_counts" )
myeloid = CreateSeuratObject( counts.mat, meta.data = myeloid@meta.data )
saveRDS( myeloid, "Myeloid/Myeloid_counts.rds", compress = F )

# run downstream analysis
myeloid <- NormalizeData(myeloid)
myeloid[["RNA"]] <- split(myeloid[["RNA"]], f = myeloid$orig.ident)
saveRDS( myeloid, "Myeloid/Myeloid_split.rds", compress = F )
myeloid <- FindVariableFeatures(myeloid, verbose = T )
saveRDS( myeloid, "Myeloid/Myeloid_var.rds"  , compress = F )
myeloid <- ScaleData(myeloid, verbose = T)
myeloid <- RunPCA(myeloid, npcs = 30, verbose = T)
saveRDS( myeloid, "Myeloid/Myeloid_pca.rds"  , compress = F )
myeloid <- IntegrateLayers(myeloid, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( myeloid, "Myeloid/Myeloid_harmony.rds", compress = F ) 
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( myeloid, "Myeloid/Myeloid_harmony_umap.rds", compress = F ) 
myeloid <- JoinLayers(myeloid)
saveRDS( myeloid, "Myeloid/Myeloid_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(myeloid, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("Myeloid/Myeloid_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  myeloid <- FindClusters(myeloid, resolution = res )
  DimPlot(myeloid, raster = TRUE, label = T, repel = T)
  ggsave(paste0("Myeloid/Myeloid_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( myeloid, paste0( "Myeloid/Myeloid_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  myeloid[["RNA"]]$data <- as(object = myeloid[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(myeloid, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
  write.csv(markers, paste0("Myeloid/Myeloid_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("Myeloid/Myeloid_harmony_clusters_topmarkers_",res,".csv"))
}
DefaultAssay(myeloid) <- "RNA"

##Round 2 Myeloid 
myeloid <- readRDS("Myeloid/Myeloid_harmony_clusters_0.6.rds")
myeloid <- myeloidset(myeloid, idents = c(6, 11, 14), invert = T)

# re-perform BPcells
counts.mat = myeloid[["RNA"]]$counts
#counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as(counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "Myeloid2/Myeloid2_counts", overwrite = T )
counts.mat <- open_matrix_dir(dir = "Myeloid2/Myeloid2_counts" )
myeloid = CreateSeuratObject( counts.mat, meta.data = myeloid@meta.data )
saveRDS( myeloid, "Myeloid2/Myeloid2_counts.rds", compress = F )

# run downstream analysis
myeloid <- NormalizeData(myeloid)
myeloid[["RNA"]] <- split(myeloid[["RNA"]], f = myeloid$orig.ident)
saveRDS( myeloid, "Myeloid2/Myeloid2_split.rds", compress = F )
myeloid <- FindVariableFeatures(myeloid, verbose = T )
saveRDS( myeloid, "Myeloid2/Myeloid2_var.rds"  , compress = F )
myeloid <- ScaleData(myeloid, verbose = T)
myeloid <- RunPCA(myeloid, npcs = 30, verbose = T)
saveRDS( myeloid, "Myeloid2/Myeloid2_pca.rds"  , compress = F )
myeloid <- IntegrateLayers(myeloid, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( myeloid, "Myeloid2/Myeloid2_harmony.rds", compress = F ) 
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( myeloid, "Myeloid2/Myeloid2_harmony_umap.rds", compress = F ) 
myeloid <- JoinLayers(myeloid)
saveRDS( myeloid, "Myeloid2/Myeloid2_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(myeloid, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("Myeloid2/Myeloid2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  myeloid <- FindClusters(myeloid, resolution = res )
  DimPlot(myeloid, raster = TRUE, label = T, repel = T)
  ggsave(paste0("Myeloid2/Myeloid2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( myeloid, paste0( "Myeloid2/Myeloid2_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  myeloid[["RNA"]]$data <- as(object = myeloid[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(myeloid, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
  write.csv(markers, paste0("Myeloid2/Myeloid2_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("Myeloid2/Myeloid2_harmony_clusters_topmarkers_",res,".csv"))
}
DefaultAssay(myeloid) <- "RNA"
myeloid <- readRDS("Myeloid2/Myeloid2_harmony_clusters_0.5.rds")
myeloid <- subset(myeloid, idents = c(11), invert = T)
myeloid$immune.annotation <- NA
myeloid$immune.annotation[myeloid$seurat_clusters == 0] <- "Macrophages C1qahi"
myeloid$immune.annotation[myeloid$seurat_clusters == 1] <- "Macrophages Fn1+"
myeloid$immune.annotation[myeloid$seurat_clusters == 2] <- "Activated DC Il4i1+"
myeloid$immune.annotation[myeloid$seurat_clusters == 3] <- "Monocytes Thbs1+"
myeloid$immune.annotation[myeloid$seurat_clusters == 4] <- "cDC1 Clec9a+"
myeloid$immune.annotation[myeloid$seurat_clusters == 5] <- "Monocytes Cd209a+"
myeloid$immune.annotation[myeloid$seurat_clusters == 6] <- "Cycling Myeloid Mki67+"
myeloid$immune.annotation[myeloid$seurat_clusters == 7] <- "Neutrophils S100a9+"
myeloid$immune.annotation[myeloid$seurat_clusters == 8] <- "Macrophages Spp1+"
myeloid$immune.annotation[myeloid$seurat_clusters == 9] <- "pDC Siglech+"
myeloid$immune.annotation[myeloid$seurat_clusters == 10] <- "Macrophage Marco+"
Idents(myeloid) <- "immune.annotation"
myeloid$immune.annotation <- factor(myeloid$immune.annotation, levels = c( "Monocytes Thbs1+", "Monocytes Cd209a+",
                                                                           "Macrophages C1qahi", "Macrophages Fn1+", "Macrophages Il4i1+", "Macrophages Spp1+", "Macrophage Marco+", 
                                                                           "Neutrophils S100a9+", 
                                                                           "cDC1 Clec9a+", "pDC Siglech+", "Cycling Myeloid Mki67+"))
Idents(myeloid) <- "immune.annotation"
colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
                    "#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
                    
dimplot <- DimPlot(myeloid, cols = colors)                             

revmyeloid <- myeloid
revmyeloid$immune.annotation <- factor(revmyeloid$immune.annotation, levels = rev(c( "Monocytes Thbs1+", "Monocytes Cd209a+",
                                                                                     "Macrophages C1qahi", "Macrophages Fn1+", "Macrophages Il4i1+", "Macrophages Spp1+", "Macrophage Marco+", 
                                                                                     "Neutrophils S100a9+", 
                                                                                     "cDC1 Clec9a+", "pDC Siglech+", "Cycling Myeloid Mki67+")))
Idents(revmyeloid) <- "immune.annotation"

genes <- c("Thbs1", "Cd209a", "C1qa", "Fn1", "Il4i1", "Spp1", "Marco", 
           "S100a9", "Clec9a", "Siglech", "Mki67")

dotplot <- DotPlot(revmyeloid, features = genes, cols = "RdYlBu", scale.min = 0, scale.max = 100) +   theme(axis.title.y = element_blank(), 
                                                                                                            axis.title.x = element_blank(),
                                                                                                            axis.text.y = element_blank(), 
                                                                                                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                                            legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                                            legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                                            legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                                            legend.title = element_text(size=6), #change legend title font size
                                                                                                            legend.text = element_text(size=6)) #change legend text font size)
combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -20, unit = "pt")))) + dotplot 
combined_plot
ggsave("SuppFigures/mousemyeloid.svg")


##Prop tables 
myeloid$group <- NA
myeloid$group[myeloid$tissue == "MAT" & myeloid$condition == "Sham"] <- "Sham MAT"
myeloid$group[myeloid$tissue == "MAT" & myeloid$condition == "Colotomy"] <- "Colotomy MAT"
myeloid$group[myeloid$tissue == "Bowel" & myeloid$condition == "Sham"] <- "Sham Bowel"
myeloid$group[myeloid$tissue == "Bowel" & myeloid$condition == "Colotomy"] <- "Colotomy Bowel"
Idents(myeloid) <- "group"

#Generate proptable and barplots of fibroblasts in mouse
tab = table(myeloid$group, myeloid$immune.annotation)
ptab = prop.table(tab, 1)*100
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
  labs(y = "Percent (%)", 
       fill = "Cluster") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, face = "bold")) + 
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x", switch = "x") + 
  scale_x_discrete(labels=c('Sham', 'Colotomy', 'Sham', 'Colotomy')) + 
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold", margin = margin(0, 0, 5, 0)),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),  # Increasing legend text size
        legend.key.size = unit(1.5, "line")) +  
  theme(panel.spacing = unit(0, "lines"))
ggsave("SuppFigures/proptablemousemyeloid.svg")


#### All proportions ####
lymphoid$immune <- NA
lymphoid$immune <- "lymphoid"
lymphoid$group <- NA
lymphoid$group[lymphoid$tissue == "MAT" & lymphoid$condition == "Sham"] <- "Sham MAT"
lymphoid$group[lymphoid$tissue == "MAT" & lymphoid$condition == "Colotomy"] <- "Colotomy MAT"
lymphoid$group[lymphoid$tissue == "Bowel" & lymphoid$condition == "Sham"] <- "Sham Bowel"
lymphoid$group[lymphoid$tissue == "Bowel" & lymphoid$condition == "Colotomy"] <- "Colotomy Bowel"
myeloid$immune <- NA
myeloid$immune <- "myeloid"
myeloid$group <- NA
myeloid$group[myeloid$tissue == "MAT" & myeloid$condition == "Sham"] <- "Sham MAT"
myeloid$group[myeloid$tissue == "MAT" & myeloid$condition == "Colotomy"] <- "Colotomy MAT"
myeloid$group[myeloid$tissue == "Bowel" & myeloid$condition == "Sham"] <- "Sham Bowel"
myeloid$group[myeloid$tissue == "Bowel" & myeloid$condition == "Colotomy"] <- "Colotomy Bowel"
immune <- merge(lymphoid, myeloid)
tab = table(immune$group, immune$immune)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Disease State Enrichment", 
       x="Disease State", y = "Percent (%)", fill = "Cluster") + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 



