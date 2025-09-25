setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages",.libPaths()))
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales)
options(Seurat.object.assay.version = "v5")

###GLOBAL ATLAS GENERATION
### QC and pre-processing ---- 
# output file to check if rows = genes; if cols = cells; if there exist non-integer counts

# QC
# WARNING: each sample must be a separate object to properly perform QC! 
obj_list <- list.files(paste0("readyForSeurat/"),"*object.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-11)
write.csv(obj_list, "obj_list_01_preQC.csv")

out <- data.frame(matrix(ncol=6,nrow=length(obj_list)))

for (i in 1:length(obj_list)) {
  ## fails at 413
  objName <- obj_list[i]
  x <- readRDS(paste0("readyForSeurat/",objName,".object.rds"))
  # QC
  x <- subset(x, subset = nFeature_RNA > 200)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-01-pre-filt.jpg"))
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(x$nFeature_RNA, 0.95) & percent.mt < 20)
  VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-02-post-filt.jpg"))
  # update output file
  out[i,] <- c(objName,
               x@assays[["RNA"]]@counts@Dimnames[[1]][1], # confirm rows = genes
               x@assays[["RNA"]]@counts@Dimnames[[2]][1], # confirm cols = cells
               nrow(x),
               ncol(x),
               any(GetAssayData(object = x, slot = "counts")%%1!=0))
  
  dir.create(paste0("processed/",strsplit(objName,"/")[[1]][1]),suppressWarnings(FALSE))
  saveRDS(x, paste0("processed/",objName, "_01_qc.rds"), compress = FALSE)
}
write.csv(out, "qcChecks.csv")
## PLEASE READ QCCHECKS ##
## STOP ##

# pre-processing: format individual post-QC objects (streamline metadata, rename gene names to standard nomenclature)
obj_list <- as.vector(read.csv("obj_list_01_preQC.csv", colClasses=c("NULL",NA))$x)
for (obj in obj_list) {
  x <- readRDS(paste0("processed/",obj,"_01_qc.rds"))
  x <- DietSeurat(x)
  x@meta.data <- x@meta.data %>%
    select(any_of(c("orig.ident", "nCount_RNA","nFeature_RNA", "percent.mt", "study",
                    "hashed","digest","age","sex","platform","treated",
                    "diseasestate","tissue","patientID","specimen_type")))
  x <- RenameCells(x, new.names = paste(x$orig.ident, Cells(x), sep = "_"))
  # deal with duplicate gene names and non-standard gene names
  counts.mat <- GetAssayData( x, slot = 'counts' )
  geneNames <- rownames(counts.mat)
  # download mouse gene names, symbols for GRCm39 using https://www.ensembl.org/biomart/martview/877296f74eb1c2a0b108975726c2386b
  synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
  # download human gene names, symbols for GRCh38.p14
  # synToName <- read.csv("/oak/stanford/groups/longaker/KEBR/IBDAFs/Human/GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
  for (i in 1:length(geneNames)) {
    name <- geneNames[i]
    if (name %in% synToName$Gene.name) {
      next
    } else if (name %in% synToName$Gene.Synonym){
      geneNames[i] <- synToName$Gene.name[which(synToName$Gene.Synonym == name)]
    } 
  }
  rownames(counts.mat) <- geneNames
  # counts mat now has duplicated row names. sum duplicated rows
  geneNamesDup <- unique(geneNames[which(duplicated(geneNames))])
  if (length(geneNamesDup) > 0) {
    counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t" )
    dup.counts.mat <- as(matrix(0,nrow = length(geneNamesDup), ncol = ncol(counts.mat)), Class = "dgCMatrix")
    colnames(dup.counts.mat) <- colnames(counts.mat)
    rownames(dup.counts.mat) <- geneNamesDup
    for (i in 1:length(geneNamesDup)) {
      gene <- geneNamesDup[i]
      rows <- which(rownames(counts.mat) == gene)
      temp <- as(counts.mat[rows,], Class = "dgCMatrix")
      temp <- as(t(colSums(temp)), Class = "dgCMatrix")
      dup.counts.mat[i,] <- temp
    }
    duprows <- (which(duplicated(geneNames) | duplicated(geneNames, fromLast = T)))
    counts.mat <- rbind(
      as(counts.mat[-duprows,], Class = "dgCMatrix"),
      dup.counts.mat
    )
  }
  x <- CreateSeuratObject(counts = counts.mat, meta.data = x@meta.data)
  saveRDS(x, paste0("processed/",obj, "_01b_correctGeneName.rds"), compress = FALSE)
}

# generate list of seurat objects and merge 
# make sure that merge is performed on seurat v3 objects; if using seurat v5 objects,
# there will be multiple layers in the final merged object
obj_list <- as.vector(read.csv("obj_list_01_preQC.csv", colClasses=c("NULL",NA))$x)
ibdaf.list <- list()
for (obj in obj_list) {
  x <- readRDS(paste0("processed/",obj,"_01b_correctGeneName.rds"))
  caf.list <- append(caf.list, x)
}
obj <- merge(caf.list[[1]],caf.list[2:length(caf.list)])
saveRDS(obj, "sketch_v5_all/merged_newNames.rds", compress = F)

# confirm no non-integer values
any(GetAssayData( obj, slot = 'counts' )%%1!=0)

for (cat in colnames(obj@meta.data)) {
  if (cat %in% c("nCount_RNA", "nFeature_RNA")) {
    next
  }
  print(table(obj[[cat]], useNA = "always"))
}

## STOP ##

# remove genes not included in standard nomenclature
counts.mat <- GetAssayData( obj, slot = 'counts' )
counts.mat <- convert_matrix_type( counts.mat, type = "uint32_t" )
# remove genes in less than 0.1% of cells (from MJ)
nCount_Feature <- rowSums( counts.mat > 0)
counts.mat <- counts.mat[nCount_Feature >=0.001*ncol(obj), ]
# download mouse gene names, symbols for GRCm39 using https://www.ensembl.org/biomart/martview/877296f74eb1c2a0b108975726c2386b
synToName <- read.csv("GRCh38.p14_geneNames_geneSynonyms.csv", header = T)
# remove genes not included in standard nomenclature
counts.mat <- counts.mat[rownames(counts.mat) %in% synToName$Gene.name,]
counts.mat <- as(counts.mat, Class = "dgCMatrix")

# switch to seurat v5 for creating final object
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = counts.mat, meta.data = obj@meta.data)
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "sketch_v5_all/unfiltered_counts")
counts.mat <- open_matrix_dir(dir = "sketch_v5_all/unfiltered_counts")
obj[["RNA"]]$counts <- counts.mat
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj.rds", compress = F)

### note: check metadata categories by using table(obj$category, useNA = "always")

### All cells analysis ----
rm(list = ls())
obj <- readRDS("sketch_v5_all/unfiltered_counts_obj.rds")
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study) ###changed from Study
saveRDS(obj, 	"sketch_v5_all/study_split_obj.rds")
obj <- FindVariableFeatures(obj, verbose = FALSE)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar.rds", compress = F)
obj <- SketchData(obj, ncells = 1000, method = "LeverageScore", sketched.assay = "sketch")
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000.rds", compress = F)
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = F)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, npcs = 30, verbose = F)
saveRDS(obj, "sketch_v5_all/unfiltered_counts_obj_findvar_sketch1000_pca.rds", compress = F)

obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, verbose = T)
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony.rds", compress = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap.rds",compress = F)
# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for (i in seq(0.05, 1.0, 0.05)) { 
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find clusters
obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_topmarkers.csv")



#Manual annotation and merging based on top genes
obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65.rds")
obj[["sketch"]] <- split(obj[["sketch"]], f = obj$study)
obj$manual.annotation <- ""
obj$manual.annotation[obj$seurat_clusters == 0] <- "Fibroblasts"
obj$manual.annotation[obj$seurat_clusters == 1] <- "T Cells"
obj$manual.annotation[obj$seurat_clusters == 2] <- "Plasma Cells"
obj$manual.annotation[obj$seurat_clusters == 3] <- "Epithelial Cells"
obj$manual.annotation[obj$seurat_clusters == 4] <- "Vascular ECs"
obj$manual.annotation[obj$seurat_clusters == 5] <- "Epithelial Cells"
obj$manual.annotation[obj$seurat_clusters == 6] <- "T Cells"
obj$manual.annotation[obj$seurat_clusters == 7] <- "Macrophages"
obj$manual.annotation[obj$seurat_clusters == 8] <- "B Cells"
obj$manual.annotation[obj$seurat_clusters == 9] <- "Fibroblasts"
obj$manual.annotation[obj$seurat_clusters == 10] <- "Macrophages"
obj$manual.annotation[obj$seurat_clusters == 11] <- "Epithelial Cells"
obj$manual.annotation[obj$seurat_clusters == 12] <- "SMCs"
obj$manual.annotation[obj$seurat_clusters == 13] <- "Pericytes"
obj$manual.annotation[obj$seurat_clusters == 14] <- "Fibroblasts"
obj$manual.annotation[obj$seurat_clusters == 15] <- "PI Cells"
obj$manual.annotation[obj$seurat_clusters == 16] <- "Glia Cells"
obj$manual.annotation[obj$seurat_clusters == 17] <- "Epithelial Cells"
obj$manual.annotation[obj$seurat_clusters == 18] <- "Mast Cells"
obj$manual.annotation[obj$seurat_clusters == 19] <- "Fibroblasts"
obj$manual.annotation[obj$seurat_clusters == 20] <- "Lymphatic ECs"
obj$manual.annotation[obj$seurat_clusters == 21] <- "Epithelial Cells"
obj$manual.annotation[obj$seurat_clusters == 22] <- "T Cells"
obj$manual.annotation[obj$seurat_clusters == 23] <- "Epithelial Cells"
Idents(obj) <- obj$manual.annotation

# project to full dataset
# obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_reclustered.rds")
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_project_joined.rds", compress = F)


##STOP

### FIBROBLAST ATLAS GENERATION

rm(list = ls())

obj <- readRDS( "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.65_project_joined.rds" )

DefaultAssay(obj) <- "RNA"
obj <- subset(obj, subset = manual.annotation == "Fibroblasts" )
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'Fibroblasts/union_fibro.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "Fibroblasts/union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "Fibroblasts/union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'fibro_normal/union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'Fibroblasts/union_fibro_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'Fibroblasts/union_fibro_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'Fibroblasts/union_fibro_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'Fibroblasts/union_fibro_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'Fibroblasts/union_fibro_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'Fibroblasts/union_fibro_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 1.0, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("Fibroblasts/union_fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibroblasts/union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj <- readRDS("Fibroblasts/union_fibro_harmony_clusters_0.3.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibroblasts/union_fibro_harmony_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibroblasts/union_fibro_harmony_clusters_0.3_topmarkers.csv")

obj <- subset(obj, subset != seurat_clusters %in% c("1", "8"))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'Fibro2/union_fibro2.rds', compress = F )

counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "Fibro2/union_fibro2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "Fibro2/union_fibro2_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'Fibro2/union_fibro2_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'Fibro2/union_fibro2_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'Fibro2/union_fibro2_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'Fibro2/union_fibro2_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'Fibro2/union_fibro2_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'Fibro2/union_fibro2_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'Fibro2/union_fibro2_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("Fibro2/union_fibro2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibro2/union_fibro2_harmony_clusters_', res, '.rds' ), compress = F ) 
}
obj <- readRDS("Fibro2/union_fibro2_harmony_clusters_0.3.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibro2/union_fibro_harmony_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibro2/union_fibro_harmony_clusters_0.3_topmarkers.csv")

obj <- subset(obj, subset != seurat_clusters %in% c("9", "10"))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'Fibro3/union_fibro.rds', compress = F )

 re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "Fibro3/union_fibro3_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "Fibro3/union_fibro3_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'Fibro3/union_fibro3_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'Fibro3/union_fibro3_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'Fibro3/union_fibro3_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'Fibro3/union_fibro3_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'Fibro3/union_fibro3_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'Fibro3/union_fibro3_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'Fibro3/union_fibro3_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("Fibro3/union_fibro3_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibro3/union_fibro3_harmony_clusters_', res, '.rds' ), compress = F ) 
}

obj <- readRDS("Fibro3/union_fibro3_harmony_clusters_0.4.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibro3/union_fibro3_harmony_clusters_0.4_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibro3/union_fibro3_harmony_clusters_0.4_topmarkers.csv")

obj <- subset(obj, subset != seurat_clusters %in% c("5", "10"))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'FibroFinal/union_fibro.rds', compress = F )

 re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "FibroFinal/union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "FibroFinal/union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'FibroFinal/union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'FibroFinal/union_fibro_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'FibroFinal/union_fibro_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'FibroFinal/union_fibro_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'FibroFinal/union_fibro_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'FibroFinal/union_fibro_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'FibroFinal/union_fibro_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("FibroFinal/union_fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'FibroFinal/union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

obj <- readRDS("FibroFinal/union_fibro_harmony_clusters_0.3.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "FibroFinal/union_fibro_harmony_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "FibroFinal/union_fibro_harmony_clusters_0.3_topmarkers.csv")

obj$fibro_annotation[obj$seurat_clusters == 0] <- "ADAMDEC1+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 1] <- "THBS1+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 2] <- "PI16+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 3] <- "F3+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 4] <- "GREM1+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 5] <- "KCNN3+ ssFib"
obj$fibro_annotation[obj$seurat_clusters == 6] <- "LRRC7+ mFib"
obj$fibro_annotation[obj$seurat_clusters == 7] <- "CTHRC1+ mFib"
obj$fibro_annotation[obj$seurat_clusters == 8] <- "MMP3+ iFib"
obj$fibro_annotation[obj$seurat_clusters == 9] <- "CD74+ apFib"
obj$fibro_annotation[obj$seurat_clusters == 10] <- "ADAMDEC1+ ssFib"

Idents(obj) <- obj$fibro_annotation
saveRDS(obj, "Annotated_Fibro_0.3.rds")

#BOWEL FIBROBLAST GENERATION

obj <- readRDS("Annotated_Fibro_0.3.rds")
obj <- subset(obj, tissue %in% c("Bowel"))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS(obj, 'Fibro_Bowel/union_fibro.rds', compress = F)

counts.mat <- obj[['RNA']]$counts
counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t")
counts.mat <- as(counts.mat, Class = 'dgCMatrix')
write_matrix_dir(mat = counts.mat, dir = "Fibro_Bowel/union_fibro_counts", overwrite = T)
counts.mat <- open_matrix_dir(dir = "Fibro_Bowel/union_fibro_counts")
obj <- CreateSeuratObject(counts.mat, meta.data = obj@meta.data)
saveRDS(obj, 'Fibro_Bowel/union_fibro_counts.rds', compress = F)
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS(obj, 'Fibro_Bowel/union_fibro_split.rds', compress = F)

obj <- FindVariableFeatures(obj, verbose = T)
saveRDS(obj, 'Fibro_Bowel/union_fibro_var.rds', compress = F)
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS(obj, 'Fibro_Bowel/union_fibro_pca.rds', compress = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
saveRDS(obj, 'Fibro_Bowel/union_fibro_harmony.rds', compress = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, 'Fibro_Bowel/union_fibro_harmony_umap.rds', compress = F)
obj <- JoinLayers(obj)
saveRDS(obj, 'Fibro_Bowel/union_fibro_harmony_umap_joined.rds', compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("Fibro_Bowel/union_fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibro_Bowel/union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

obj <- readRDS("Fibro_Bowel/union_fibro_harmony_clusters_0.35.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibro_Bowel/union_fibro_harmony_clusters_0.35_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibro_Bowel/union_fibro_harmony_clusters_0.35_topmarkers.csv")

obj$bowel.annotation <- ""
obj$bowel.annotation[obj$seurat_clusters == 0] <- "ADAMDEC1+ lpFib"
obj$bowel.annotation[obj$seurat_clusters == 1] <- "PI16+ ssFRC"
obj$bowel.annotation[obj$seurat_clusters == 2] <- "F3+ Telocytes"
obj$bowel.annotation[obj$seurat_clusters == 3] <- "THBS1+ ssFRC"
obj$bowel.annotation[obj$seurat_clusters == 4] <- "ADAMDEC1+ lpFib"
obj$bowel.annotation[obj$seurat_clusters == 5] <- "KCNN3+ ssFRC"
obj$bowel.annotation[obj$seurat_clusters == 6] <- "GREM1+ ssFRC"
obj$bowel.annotation[obj$seurat_clusters == 7] <- "CTHRC1+ mFib"
obj$bowel.annotation[obj$seurat_clusters == 8] <- "CD74+ apFib"
obj$bowel.annotation[obj$seurat_clusters == 9] <- "LRRC7+ mFib"
obj$bowel.annotation[obj$seurat_clusters == 10] <- "MMP3+ iFib"
obj$bowel.annotation[obj$seurat_clusters == 11] <- "MMP3+ iFib"
obj$bowel.annotation[obj$seurat_clusters == 12] <- "ADAMDEC1+ lpFib"

Idents(obj) <- obj$bowel.annotation
saveRDS(obj, "Fibro_Bowel_0.35.rds", compress = F)

#MAT FIBROBLAST GENERATION

obj <- readRDS("Annotated_Fibro_0.3.rds")
obj <- subset(obj, tissue %in% c("MAT"))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'Fibro_MAT/union_fibro.rds', compress = F )

counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "Fibro_MAT/union_fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "Fibro_MAT/union_fibro_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'Fibro_MAT/union_fibro_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'Fibro_MAT/union_fibro_split.rds', compress = F )
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'Fibro_MAT/union_fibro_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'Fibro_MAT/union_fibro_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'Fibro_MAT/union_fibro_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'Fibro_MAT/union_fibro_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'Fibro_MAT/union_fibro_harmony_umap_joined.rds', compress = F ) 

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("Fibro_MAT/union_fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibro_MAT/union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

obj <- readRDS("Fibro_MAT/union_fibro_harmony_clusters_0.3.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibro_MAT/union_fibro_harmony_clusters_0.3_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibro_MAT/union_fibro_harmony_clusters_0.3_topmarkers.csv")

obj <- subset(obj, subset = !(seurat_clusters %in% c(10, 11, 13, 6)))
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS(obj, 'Fibro_MAT2/union_mat.rds', compress = F)

counts.mat <- obj[['RNA']]$counts
counts.mat <- convert_matrix_type(counts.mat, type = "uint32_t")
counts.mat <- as(counts.mat, Class = 'dgCMatrix')
write_matrix_dir(mat = counts.mat, dir = "Fibro_MAT2/union_mat_counts", overwrite = T)
counts.mat <- open_matrix_dir(dir = "Fibro_MAT2/union_mat_counts")
obj <- CreateSeuratObject(counts.mat, meta.data = obj@meta.data)
saveRDS(obj, 'Fibro_MAT2/union_mat_counts.rds', compress = F)

obj <- NormalizeData(obj)
# Splitting creates separate layers (e.g., for each study)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS(obj, 'Fibro_MAT2/union_mat_split.rds', compress = F)

obj <- FindVariableFeatures(obj, verbose = T)
saveRDS(obj, 'Fibro_MAT2/union_mat_var.rds', compress = F)
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS(obj, 'Fibro_MAT2/union_mat_pca.rds', compress = F)
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
saveRDS(obj, 'Fibro_MAT2/union_mat_harmony.rds', compress = F)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS(obj, 'Fibro_MAT2/union_mat_harmony_umap.rds', compress = F)
obj <- JoinLayers(obj)
saveRDS(obj, 'Fibro_MAT2/union_mat_harmony_umap_joined.rds', compress = F)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.05, 0.5, 0.05) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("Fibro_MAT2/union_fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'Fibro_MAT2/union_fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

obj <- readRDS("Fibro_MAT2/union_fibro_harmony_clusters_0.2.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Fibro_MAT2/union_fibro_harmony_clusters_0.2_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Fibro_MAT2/union_fibro_harmony_clusters_0.2_topmarkers.csv")

obj$mat.annotation <- ""
obj$mat.annotation[obj$seurat_clusters == 0] <- "FMO2+ ssAPC"
obj$mat.annotation[obj$seurat_clusters == 1] <- "DPP4+ ssAPC"
obj$mat.annotation[obj$seurat_clusters == 2] <- "ICAM1+ iAPC"
obj$mat.annotation[obj$seurat_clusters == 3] <- "PPARG+ ssAPC"
obj$mat.annotation[obj$seurat_clusters == 4] <- "CTHRC1+ mFAP"
obj$mat.annotation[obj$seurat_clusters == 5] <- "VIT+ ssAPC"
obj$mat.annotation[obj$seurat_clusters == 6] <- "LRRC7+ mFAP"
obj$mat.annotation[obj$seurat_clusters == 7] <- "NOTCH3+ Pericytes"

Idents(obj) <- obj$mat.annotation
saveRDS(obj, "Fibro_MAT_0.2.rds", compress = F)

##TISSUE ALT
bowel = readRDS("Fibro_Bowel_Annotated_0.35.rds")
mat = readRDS("Fibro_MAT_0.2.rds")
bowel <- FindVariableFeatures(bowel)
mat <- FindVariableFeatures(mat)

anchors <- FindTransferAnchors(reference = bowel, query = mat, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = bowel$bowel.annotation, dims = 1:30)
matquery <- AddMetaData(mat, metadata = predictions)
saveRDS(matquery, "Bowel_MAT_ALT.rds")


### STANFORD DATASET
obj <- readRDS("Annotated_Fibro_0.3.rds")
obj <- subset(global, subset = study %in% c("LongakerBowel", "LongakerMes") )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
new <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( new, 'StanfordSub/newFibro1/fibro1.rds', compress = F )
table(new$study)

# re-perform BPcells
counts.mat = new[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "StanfordSub/newFibro1/fibro1_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "StanfordSub/newFibro1/fibro1_counts" )
new = CreateSeuratObject( counts.mat, meta.data = new@meta.data )
saveRDS( new, 'StanfordSub/newFibro1/fibro1_counts.rds', compress = F )
table(new$study)

rm(list = ls())
fibro <- readRDS('StanfordSub/newFibro1/fibro1_counts.rds')
# run standard anlaysis workflow
fibro <- NormalizeData(fibro)
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- FindVariableFeatures(fibro)
fibro <- ScaleData(fibro)
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
fibro <- JoinLayers(fibro)
fibro <- FindNeighbors(fibro, reduction = "harmony", dims = 1:30)
for (i in seq(0.1, 1.0, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, raster = FALSE, label = TRUE)
  ggsave(paste0("StanfordSub/newFibro1/fibro1_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("StanfordSub/newFibro1/fibro1_harmony_umap_res",i,".rds"), compress = FALSE)
}

fibro <- readRDS("StanfordSub/newFibro1/fibro1_harmony_umap_res_0.5.rds")
fibro <- subset(fibro, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 17))

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
new <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( new, 'StanfordSub/newFibro1B/fibro1B.rds', compress = F )
table(new$study)

# re-perform BPcells
counts.mat = new[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "StanfordSub/newFibro1B/fibro1B_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "StanfordSub/newFibro1B/fibro1B_counts" )
new = CreateSeuratObject( counts.mat, meta.data = new@meta.data )
saveRDS( new, 'StanfordSub/newFibro1B/fibro1B_counts.rds', compress = F )
table(new$study)

rm(list = ls())
fibro <- readRDS('StanfordSub/newFibro1B/fibro1B_counts.rds')
# run standard anlaysis workflow
fibro <- NormalizeData(fibro)
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$orig.ident)
fibro <- FindVariableFeatures(fibro)
fibro <- ScaleData(fibro)
fibro <- RunPCA(fibro, npcs = 30, verbose = T)
fibro <- IntegrateLayers(fibro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
fibro <- RunUMAP(fibro, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
fibro <- JoinLayers(fibro)
fibro <- FindNeighbors(fibro, reduction = "harmony", dims = 1:30)
for (i in seq(0.1, 1.0, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, raster = FALSE, label = TRUE)
  ggsave(paste0("StanfordSub/newFibro1B/fibro1B_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("StanfordSub/newFibro1B/fibro1B_harmony_umap_res",i,".rds"), compress = FALSE)
}

obj <- readRDS("StanfordSub/newFibro1B/fibro1B_harmony_umap_res_0.5.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "StanfordSub/fibro1B_harmony_umap_res0.5_markers.csv")

obj$fibro <- ""
obj$fibro[obj$seurat_clusters == 0] <- "GREM1+ ssFib"
obj$fibro[obj$seurat_clusters == 1] <- "CTHRC1+ mFib"
obj$fibro[obj$seurat_clusters == 2] <- "CD74+ apFib"
obj$fibro[obj$seurat_clusters == 3] <- "FMO2+ ssFib"
obj$fibro[obj$seurat_clusters == 4] <- "ADAMTS16+ ssFib"
obj$fibro[obj$seurat_clusters == 5] <- "GREM1+ ssFib"
obj$fibro[obj$seurat_clusters == 6] <- "LRRC7+ mFib"
obj$fibro[obj$seurat_clusters == 7] <- "PI16+ ssFib"
obj$fibro[obj$seurat_clusters == 8] <- "FMO2+ ssFib"
obj$fibro[obj$seurat_clusters == 9] <- "ADAMDEC1+ ssFib"
obj$fibro[obj$seurat_clusters == 10] <- "MMP3+ iFib"
obj$fibro[obj$seurat_clusters == 11] <- "KCNN3+ ssFib"
obj$fibro[obj$seurat_clusters == 12] <- "NOTCH3+ Pericytes"
obj$fibro[obj$seurat_clusters == 13] <- "VIT+ ssFib"
obj$fibro[obj$seurat_clusters == 14] <- "F3+ ssFib"
obj$fibro[obj$seurat_clusters == 15] <- "MSLN+ Mesothelial Cells"
obj$fibro[obj$seurat_clusters == 16] <- "PPARG+ Preadipocytes"

Idents(obj) <- obj$fibro
saveRDS(obj, "Stanford_Fibro_0.5.rds", compress = F)

##STANFORD ALT
stanford = readRDS("Stanford_Fibro_0.5.rds")
global = readRDS("Annotated_Fibro_0.3.rds")
stanford <- FindVariableFeatures(stanford)
global <- FindVariableFeatures(global)

anchors <- FindTransferAnchors(reference = global, query = stanford, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = global$fibro.annotation, dims = 1:30)
stanfordquery <- AddMetaData(stanford, metadata = predictions)
saveRDS(stanfordquery, "Stanford_Fibro_ALT.rds", compress = F)

