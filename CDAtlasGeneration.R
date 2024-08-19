### META-ANALYSIS CODE
### Written by John Lu
### Date: Dec 28, 2023

### load libraries ----

setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages/updated",.libPaths()))
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales)
options(Seurat.object.assay.version = "v5")

### QC and pre-processing ---- 
# output file to check if rows = genes; if cols = cells; if there exist non-integer counts

# QC
# WARNING: each sample must be a separate object to properly perform QC! 
obj_list <- list.files(paste0("readyForSeurat/"),"*object.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-11)
write.csv(obj_list, "obj_list_01_preQC.csv")

out <- data.frame(matrix(ncol=7,nrow=length(obj_list)))

for (i in 1:length(obj_list)) {
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
caf.list <- list()
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


# remove genes not included in standard nomenclature
counts.mat <- GetAssayData( obj, slot = 'counts' )
counts.mat <- convert_matrix_type( counts.mat, type = "uint32_t" )
# remove genes in less than 0.1% of cells (from MJ)
nCount_Feature <- rowSums( counts.mat > 0)
counts.mat <- counts.mat[nCount_Feature >=0.001*ncol(obj), ]
# download human gene names
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
saveRDS(obj, 	"sketch_v5_all/split_orig.ident/orig.identsplit_obj.rds")
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
FeaturePlot(obj, features = c("PECAM1","EPCAM","RGS5","PDGFRA","PTPRC","CD3D","CD19","IGHG1","ITGAM","CPA3","PLP1","COL1A1","PI16","DPT"), raster = TRUE, reduction = "umap.harmony")
ggsave("sketch_v5_all/unfiltered_obj_sketch1000_harmony_featureplot.jpg", width = 15, height = 15, units = "in", limitsize = FALSE)
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap.rds",compress = F)
# clustering
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
#for (i in seq(0.03, 0.1, 0.01)) { 
for (i in seq(0.1, 0.2, 0.02)) { 
  obj <- FindClusters(obj, resolution = i)
  DimPlot(obj, raster = TRUE, label = TRUE)
  ggsave(paste0("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(obj, paste0("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res",i,".rds"), compress = FALSE)
}

# find cluster markers
obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.5.rds")
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.5_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.5_topmarkers.csv")

# res 0.5:
#Fibroblasts: 0, 11, 14
#T Cells: 1, 7
#Epithelial: 2, 9, 10, 17, 18
#Macrophages: 3
#VECs: 4, 
#LECs: 13
#Plasma cells: 5
#Smooth muscle: 6
#B Cells: 8
#Neutrophils: 12
#Mast cells: 15
#Glia: 16

#Determined that res = 0.2 was best to get all desired clusters

#Manual annotation and merging based on top genes
obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2.rds")
obj[["sketch"]] <- split(obj[["sketch"]], f = obj$study)
# obj <- RenameIdents(obj,  '8' = '2', '9' = '1', '10' = '2', "15" = "2")
# DimPlot(obj, raster = TRUE, label = TRUE)
# ggsave(paste0("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res",0.2,"reclustered.jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
# saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_reclustered.rds")

# project to full dataset
# obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_reclustered.rds")
obj <- ProjectIntegration(object = obj, reduction = "harmony")
options(future.globals.maxSize = 8000 * 1024^2)
obj <- ProjectData(object = obj, sketched.reduction = "harmony.full", full.reduction = "harmony.full", umap.model = "umap.harmony", dims = 1:30, refdata = list(harmony.cluster.full = "seurat_clusters"))
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project.rds",compress = F)
saveRDS(obj$harmony.cluster.full, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_clusterList.rds")

# join layers
DefaultAssay(obj) <- "sketch"
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.rds", compress = F)

# manual annotations for entire dataset
DefaultAssay(obj) <- "RNA"
obj$manual.annotation <- "tbd"
obj$manual.annotation[obj$harmony.cluster.full == 0] <- "T Cells"
obj$manual.annotation[obj$harmony.cluster.full == 1] <- "Fibroblasts"
obj$manual.annotation[obj$harmony.cluster.full == 2] <- "Epithelial Cells"
obj$manual.annotation[obj$harmony.cluster.full == 3] <- "Macrophages"
obj$manual.annotation[obj$harmony.cluster.full == 4] <- "Vascular ECs"
obj$manual.annotation[obj$harmony.cluster.full == 5] <- "Plasma Cells"
obj$manual.annotation[obj$harmony.cluster.full == 6] <- "SMCs"
obj$manual.annotation[obj$harmony.cluster.full == 7] <- "B Cells"
obj$manual.annotation[obj$harmony.cluster.full == 8] <- "Epithelial Cells"
obj$manual.annotation[obj$harmony.cluster.full == 9] <- "Fibroblasts"
obj$manual.annotation[obj$harmony.cluster.full == 10] <- "Epithelial Cells"
obj$manual.annotation[obj$harmony.cluster.full == 11] <- "Myeloid Cells"
obj$manual.annotation[obj$harmony.cluster.full == 12] <- "Lymphatic ECs"
obj$manual.annotation[obj$harmony.cluster.full == 13] <- "Mast Cells"
obj$manual.annotation[obj$harmony.cluster.full == 14] <- "Glia"
obj$manual.annotation[obj$harmony.cluster.full == 15] <- "Epithelial Cells"


#Reorder factor levels 
obj$diseasestate <- factor(x = obj$diseasestate, levels = c("Normal Bowel", "NS Uninflamed", "NS Inflamed", 
                                                                  "Stricture Uninflamed", "Stricture Inflamed", "Stricture", 
                                                                  "Penatrating Uninflamed", "Penatrating Inflamed","Ileostomy Bowel", 
                                                                  "Normal MAT", "Uninvolved", "Involved", "Ileostomy MAT"))
table(obj$diseasestate, obj$tissue)



saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds", compress = F)
saveRDS(obj, "newHuman/lobalcellatlas.annotated.rds", compress = F) 
Idents(obj) <- "manual.annotation"
DimPlot(obj, group.by = "diseasestate", label = TRUE, raster = TRUE, alpha = 0.03)

#Generate proptable and barplots 
tab = table(obj$diseasestate, obj$manual.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Global Cell Type Enrichment", 
       x="Condition", y = "Percent (%)", fill = "Cluster") + 
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))


### fibro subset 1 ----

rm(list = ls())

saveRDS(obj, "sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds", compress = F)
obj <- readRDS("sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds")
DefaultAssay(obj) <- "RNA"
obj <- subset(obj, subset = manual.annotation == "Fibroblasts" )
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'newFibro1/fibro1.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro1/fibro1_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro1/fibro1_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'newFibro1/fibro1_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'newFibro1/fibro1_split.rds', compress = F )


obj<- readRDS('newFibro1/fibro1_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- obj[['RNA']][[Layers(obj)[i]]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(obj)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "newFibro1/Error_list.csv")

error_df <- read.csv("newFibro1/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study

#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'newFibro1/fibro1_split_subset.rds', compress = F)

obj <- obj_subset
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'newFibro1/fibro1_var.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'newFibro1/fibro1_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'newFibro1/fibro1_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'newFibro1/fibro1_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'newFibro1/fibro1_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("newFibro1/fibro1_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.03, 0.1, 0.01) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE)
  ggsave(paste0("newFibro1/fibro1_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'newFibro1/fibro1_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj <- readRDS("newFibro1/fibro1_harmony_clusters_0.07.rds")
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro1/fibro1_harmony_clusters_0.07_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro3/fibro3_harmony_clusters_0.2_topmarkers.csv")

#Clusters to remove: 
#Cluster 1 - low quality cells
#Cluster 6 - SMCs

### fibro subset 2 ----
rm(list = ls())
obj <- readRDS( "newFibro1/fibro1_harmony_clusters_0.07.rds")
DefaultAssay(obj) <- "RNA"
#Remove clusters
obj <- subset(obj, idents = c("1","6"), invert = T )

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'newFibro2/fibro2.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro2/fibro2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro2/fibro2_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'newFibro2/fibro2_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'newFibro2/fibro2_split.rds', compress = F )

obj<- readRDS('newFibro2/fibro2_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- obj[['RNA']][Layers(obj)[i]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(obj)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "newFibro2/Error_list.csv")

error_df <- read.csv("newFibro2/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'newFibro2/fibro2_split_subset.rds', compress = F)

obj <- readRDS('newFibro2/fibro2_split_subset.rds')
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'newFibro2/fibro2_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'newFibro2/fibro2_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'newFibro2/fibro2_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'newFibro2/fibro2_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("newFibro2/fibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.31, 0.35, 0.01) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("newFibro2/fibro2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'newFibro2/fibro2_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro2/fibro2fibro2_harmony_clusters_0.26_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro2/fibro2fibro2_harmony_clusters_0.26_topmarkers.csv")

obj <- subset(obj, idents = "12", invert = T) #remove cluster 12, which is just 2 cells 

#annotation of fibroblasts 
# manual annotations for entire dataset
Idents(obj) <- "seurat_clusters"
obj$fibro.annotation <- "tbd"
obj$fibro.annotation[obj$seurat_clusters == 0] <- "Steady-State PI16+"
obj$fibro.annotation[obj$seurat_clusters == 1] <- "Inflammatory ADAMDEC1+"
obj$fibro.annotation[obj$seurat_clusters == 2] <- "Steady-State FMO2+"
obj$fibro.annotation[obj$seurat_clusters == 3] <- "Steady-State F3+"
obj$fibro.annotation[obj$seurat_clusters == 4] <- "Mechanosensitive LRRC7+"
obj$fibro.annotation[obj$seurat_clusters == 5] <- "Mechanosensitive CTHRC1+"
obj$fibro.annotation[obj$seurat_clusters == 6] <- "Steady-State KCNN3+"
obj$fibro.annotation[obj$seurat_clusters == 7] <- "Steady-State GREM1+"
obj$fibro.annotation[obj$seurat_clusters == 8] <- "Antigen-Presenting CCL19+"
obj$fibro.annotation[obj$seurat_clusters == 9] <- "Inflammatory MMP3+"
obj$fibro.annotation[obj$seurat_clusters == 10] <- "Steady-State KCNN3+"
obj$fibro.annotation[obj$seurat_clusters == 11] <- "Steady-State PI16+"
Idents(obj) <- "fibro.annotation"
#Reorder factor levels 
longsub$fibro.annotation <- factor(x = longsub$fibro.annotation, levels = rev(c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
                                                                    "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
                                                                    "mCDF LRRC7+", "mCDF CTHRC1+")))
Idents(longsub) <- "fibro.annotation"
DimPlot(obj, raster = T)
saveRDS(obj, "finalfibroblastannotations.rds")











#Generate proptable and barplots 
tab = table(obj$diseasestate, obj$fibro.annotation)
ptab = as.data.frame(prop.table(tab, 1)*100)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Enrichment", 
       x="Condition", y = "Percent (%)", fill = "Cluster") + 
  theme(axis.text = element_text(size = 12)) + 
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))

DotPlot(obj, features = c("PI16", "ADAMDEC1", "FMO2", "F3", "LRRC7", "CTHRC1", 
                          "KCNN3", "GREM1", "CCL19", "MMP3")) + 

  
### fibro subset 2 ----
rm(list = ls())
obj <- readRDS( "finalfibroblastannotations.rds")
DefaultAssay(obj) <- "RNA"
#Remove clusters
obj <- subset(obj, diseasestate == "Ileostomy Bowel" | diseasestate == "Ileostomy MAT", invert = T)
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'newFibro3/fibro3.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro3/fibro3_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro3/fibro3_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'newFibro3/fibro3_counts.rds', compress = F )

# run downstream analysis
#obj <- JoinLayers(obj)
#obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
#saveRDS( obj, 'newFibro3/fibro3_split.rds', compress = F )  
  

#######################End of Analysis 

#Additional analyses run to confirm clustering
rm(list = ls())
#obj <- readRDS('newFibro3/fibro3_split.rds')
obj<- readRDS("newFibro3/fibro3_counts.rds")
obj <- NormalizeData(obj)
saveRDS( obj, 'newFibro3/normalized_fibro3_pca.rds'  , compress = F )
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, 'newFibro3/normalized_vfs_fibro3_pca.rds'  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'newFibro3/normalized_vfs_fibro3_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'newFibro3/vf_fibro3_integrated_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'newFibro3/vf_fibro3_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'newFibro3/vf_fibro3_harmony_umap_joined.rds', compress = F ) 

#FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
#ggsave("newFibro2/V3000/V3Kfibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 0.5, 0.1) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("newFibro3/fibro3_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'newFibro3/fibro3_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro2/V5000/V5Kfibro_harmony_clusters_0.20_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro2/V5000/V5Kfibro_harmony_clusters_0.20_topmarkers.csv")

#scale all genes
obj <- readRDS('newFibro2/fibro2_split_subset.rds')
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, features = rownames(obj), verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'newFibro2/Scaled/scaledfibro2_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'newFibro2/Scaled/scaledfibro2_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'newFibro2/Scaled/scaledfibro2_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'newFibro2/Scaled/scaledfibro2_harmony_umap_joined.rds', compress = F ) 

#FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
#ggsave("newFibro2/V3000/V3Kfibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 0.5, 0.1) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = T)
  ggsave(paste0("newFibro2/Scaled/scaledfibro2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( 'newFibro2/Scaled/scaledfibro2_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro2/Scaled/scaledfibro_harmony_clusters_0.20_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro2/Scaled/scaledfibro_harmony_clusters_0.20_topmarkers.csv")

#### Global Cell Atlas UMAP ####
global <- readRDS("globalcellatlas.annotated.rds")
DefaultAssay(global) <- "RNA"
global$manual.annotation <- factor(x = global$manual.annotation, levels = c("Fibroblasts", "SMCs", "Vascular ECs", "Lymphatic ECs", "Glia", 
                                                                            "T Cells", "B Cells", "Plasma Cells", "Macrophages", "Myeloid Cells", 
                                                                            "Mast Cells", "Epithelial Cells"))
colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
                    "#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
                    
# "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", 
# "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", "#5A5156", "#E4E1E3", 
# "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
# "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
# "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6",
# "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79",
# "#66B0FF", "#3B00FB")
Idents(global) <- "manual.annotation"
dimplot <- DimPlot(global, raster = T, cols = colors, alpha = 0.1)
#ggsave("Figures/globalumap.svg", width = 5.51, height = 3.77)

revglobal <- global
revglobal$manual.annotation <- factor(x = revglobal$manual.annotation, levels = rev(c("Fibroblasts", "SMCs", "Vascular ECs", "Lymphatic ECs", "Glia", 
                                                                                      "T Cells", "B Cells", "Plasma Cells", "Macrophages", "Myeloid Cells", 
                                                                                      "Mast Cells", "Epithelial Cells")))
Idents(revglobal) <- "manual.annotation"
dotplot <- DotPlot(revglobal, features = c("PDGFRA", "DES", "PECAM1", "PROX1", "PLP1", "CD3E", 
                                           "BANK1", "JCHAIN", "CD68" , "S100A8", "CPA3", "EPCAM"), cols = "RdYlBu", scale.min = 0, scale.max = 100) + theme(axis.title.y = element_blank(), 
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
ggsave("SuppFigures/globalumapdotplot.svg")

fibro <- readRDS("finalfibroblast.annotations.rds")
fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = rev(c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                            "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                            "Antigen-Presenting CCL19+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")))
Idents(fibro) <- "fibro.annotation"
#Colors for fibroblast studies
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
                     
#### UMAP by study ####
global <- readRDS("globalcellatlas.annotated.rds")
Idents(global) <- "study"
global$study[global$study == "LongakerBowel" ] <- "Stanford_Hospital"
global$study[global$study == "LongakerMes" ] <- "Stanford_Hospital"
global$study[global$study == "Rieder" ] <- "Cleveland_Clinic"
study_idents <- unique(global$study)
cell.list <- WhichCells(global, idents = study_idents, downsample = min(table(global$study)))
global_small <- global[, cell.list]
table(global_small$study)
Idents(global_small) <- "study"
colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
"#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
DimPlot(global_small, cols = colors)
ggsave("SuppFigures/Fibro.studysplit.svg")

#### All Immune Cells ####
global <- readRDS("globalcellatlas.annotated.rds")
global$immune <- NA
global$immune[global$manual.annotation == "T Cells"] <- "Lymphoid"
global$immune[global$manual.annotation == "B Cells"] <- "Lymphoid"
global$immune[global$manual.annotation == "Plasma Cells"] <- "Lymphoid"
global$immune[global$manual.annotation == "Myeloid Cells"] <- "Myeloid"
global$immune[global$manual.annotation == "Macrophages"] <- "Myeloid"
global$immune[global$manual.annotation == "Mast Cells"] <- "Myeloid"
Idents(global) <- "immune"
sub <- subset(global, idents = c("Lymphoid", "Myeloid"))
sub <- subset(sub, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
sub <- subset(sub, specimen_type == "Resection")
sub <- subset(sub, tissue == "Colon", invert = T)
sub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T)
Idents(sub) <- "immune"

Idents(sub) <- "diseasestate"
idents <- c("Normal Bowel", "Stricture Uninflamed", "Stricture", "Normal MAT", "Uninvolved", "Involved")
sub$diseasestate <- factor(sub$diseasestate, levels = idents)
Idents(sub) <- "diseasestate"
tab = table(sub$diseasestate, sub$immune)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Disease State Enrichment", 
       x="Disease State", y = "Percent (%)", fill = "Cluster") + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 

#### Lymphoid ####
global <- readRDS("globalcellatlas.annotated.rds")
Idents(global) <- "manual.annotation"
obj <- subset(global, idents = c("T Cells", "B Cells","Plasma Cells"))

##Regenerate and clean up SO 
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, "Lymphoid/Lymphoid.rds", compress = F )

# re-perform BPcells
counts.mat = obj[["RNA"]]$counts
#counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as(counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "LymphoidNew/LymphoidNew_counts", overwrite = T )
counts.mat <- open_matrix_dir(dir = "LymphoidNew/LymphoidNew_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, "LymphoidNew/LymphoidNew_counts.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, "LymphoidNew/LymphoidNew_split.rds", compress = F )


#obj<- readRDS("LymphoidNew/LymphoidNew_split.rds")
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- obj[["RNA"]][Layers(obj)[i]]
  # Try executing FindVariableFeature
  tryCatch({
    # Find variable features
    print(Layers(obj)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "LymphoidNew/Error_list.csv")

error_df <- read.csv("LymphoidNew/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study

#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, "LymphoidNew/LymphoidNew_split_subset.rds", compress = F)

obj <- obj_subset
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, "LymphoidNew/LymphoidNew_var.rds"  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, "LymphoidNew/LymphoidNew_pca.rds"  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, "LymphoidNew/LymphoidNew_harmony.rds", compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, "LymphoidNew/LymphoidNew_harmony_umap.rds", compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, "LymphoidNew/LymphoidNew_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("LymphoidNew/LymphoidNew_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.1, 0.1) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = T, repel = T)
  ggsave(paste0("LymphoidNew/LymphoidNew_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( "LymphoidNew/LymphoidNew_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
  write.csv(markers, paste0("LymphoidNew/LymphoidNew_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("LymphoidNew/LymphoidNew_harmony_clusters_topmarkers_",res,".csv"))
}

DefaultAssay(obj) <- "RNA"

##Recluster cells
#Remove clusters
obj <- readRDS("LymphoidNew_harmony_clusters_0.2.rds")
obj <- subset(obj, idents = c(0:5, 7:9, 11))

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, "LymphoidNew2/LymphoidNew2_0.2.rds", compress = F )

# re-perform BPcells
counts.mat = obj[["RNA"]]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "LymphoidNew2/LymphoidNew2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "LymphoidNew2/LymphoidNew2_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, "LymphoidNew2/LymphoidNew2_counts.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, "LymphoidNew2/LymphoidNew2_split.rds", compress = F )
obj<- readRDS("LymphoidNew2/LymphoidNew2_split.rds")

# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- obj[["RNA"]][Layers(obj)[i]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(obj)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "LymphoidNew2/Error_list.csv")

error_df <- read.csv("LymphoidNew2/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, "LymphoidNew2/LymphoidNew2_split_subset.rds", compress = F)

obj <- readRDS("LymphoidNew2/LymphoidNew2_split_subset.rds")
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, "LymphoidNew2/LymphoidNew2_pca.rds"  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, "LymphoidNew2/LymphoidNew2_harmony.rds", compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, "LymphoidNew2/LymphoidNew2_harmony_umap.rds", compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, "LymphoidNew2/LymphoidNew2_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("LymphoidNew2/LymphoidNew2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.5, 0.1) ) {
  obj <- FindClusters(obj, resolution = res )
  DimPlot(obj, raster = TRUE, label = T, repel = T)
  ggsave(paste0("LymphoidNew2/LymphoidNew2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( obj, paste0( "LymphoidNew2/LymphoidNew2_harmony_clusters_", res, ".rds" ), compress = F ) 
  # join and find markers
  obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
  markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
  write.csv(markers, paste0("LymphoidNew2/LymphoidNew2_harmony_clusters_markers_",res,".csv"))
  markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
  topmarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 100, order_by = rank)
  write.csv(topmarkers, paste0("LymphoidNew2/LymphoidNew2_harmony_clusters_topmarkers_",res,".csv"))
}


##Annotation of cell types
obj <- readRDS("LymphoidNew2/LymphoidNew2_harmony_clusters_0.31.rds")
obj <- subset(obj, idents = c(16,17), invert = T) #remove low quality
# manual annotations for entire dataset
DefaultAssay(obj) <- "RNA"
obj$immune.annotation <- "tbd"
obj$immune.annotation[obj$seurat_clusters == 0] <- "Plasma Cells IGHA1+"
obj$immune.annotation[obj$seurat_clusters == 1] <- "Follicular B Cells IGHD+"
obj$immune.annotation[obj$seurat_clusters == 2] <- "TH1 T Cells  CCR7+"
obj$immune.annotation[obj$seurat_clusters == 3] <- "CD8 T Cells CD8A+"
obj$immune.annotation[obj$seurat_clusters == 4] <- "TH17 T Cells CCR6+"
obj$immune.annotation[obj$seurat_clusters == 5] <- "Plasma Cells IGHA1+"
obj$immune.annotation[obj$seurat_clusters == 6] <- "Tregs FOXP3+"
obj$immune.annotation[obj$seurat_clusters == 7] <- "NK Cells FCER1G+"
obj$immune.annotation[obj$seurat_clusters == 8] <- "Plasma Cells IGHG1+"
obj$immune.annotation[obj$seurat_clusters == 9] <- "Cycling B Cells MKI67+"
obj$immune.annotation[obj$seurat_clusters == 10] <- "GC B Cells AICDA+"
obj$immune.annotation[obj$seurat_clusters == 11] <- "ILC IL7R+"
obj$immune.annotation[obj$seurat_clusters == 12] <- "Naïve B Cells AFF3+"
obj$immune.annotation[obj$seurat_clusters == 13] <- "Plasma Cells IGHA1+"
obj$immune.annotation[obj$seurat_clusters == 14] <- "Plasma Cells IGHA1+"
obj$immune.annotation[obj$seurat_clusters == 15] <- "Memory B Cells TNFRSF13B+"

obj$immune.annotation <- factor(obj$immune.annotation, levels = c("Naïve B Cells AFF3+", "GC B Cells AICDA+", "Cycling B Cells MKI67+", "Memory B Cells TNFRSF13B+", "Follicular B Cells IGHD+",
                                                                  "Plasma Cells IGHG1+", "Plasma Cells IGHA1+", 
                                                                  "TH1 T Cells  CCR7+", "TH17 T Cells CCR6+", "Tregs FOXP3+", "ILC IL7R+",
                                                                  "CD8 T Cells CD8A+", "NK Cells FCER1G+"))
Idents(obj) <- "immune.annotation"

colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
                    "#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0", "#0072fb", "#de3163", "#3ea5a1")
dimplot <- DimPlot(obj, raster = T, repel = T, cols = colors)

#Correct loss of metadata
global <- readRDS("globalcellatlas.annotated.rds")
Idents(obj) <- "diseasestate"
obj$diseasestate <- global$diseasestate[colnames(global) %in% colnames(obj)]
saveRDS(obj, "lymphoid.annotated.rds", compress = F)
revobj <- obj
revobj$immune.annotation <- factor(revobj$immune.annotation, levels = rev(c("Naïve B Cells AFF3+", "GC B Cells AICDA+", "Cycling B Cells MKI67+", "Memory B Cells TNFRSF13B+", "Follicular B Cells IGHD+",
                                                        "Plasma Cells IGHG1+", "Plasma Cells IGHA1+", 
                                                        "TH1 T Cells  CCR7+", "TH17 T Cells CCR6+", "Tregs FOXP3+", "ILC IL7R+",
                                                        "CD8 T Cells CD8A+", "NK Cells FCER1G+")))
Idents(revobj) <- "immune.annotation"
genes <- c("AFF3", "AICDA", "MKI67", "TNFRSF13B", "IGHD", "IGHG1", "IGHA1", "CCR7", "CCR6", "FOXP3", "IL7R", "CD8A", "FCER1G", "CD3", "CD4", "BANK1", "CD19")
dotplot <- DotPlot(revobj, features = genes, cols = "RdYlBu") +   theme(axis.title.y = element_blank(), 
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
ggsave("SuppFigures/lymphoidallumap.svg")

##Prop tables 
#Generate comparable dataset
sub <- subset(obj, diseasestate == "Ileostomy Bowel" | diseasestate == "Ileostomy MAT", invert = T) #Remove ileostomy samples
sub <- subset(sub, specimen_type == "Resection") #Only compare full thickness due to mucosal bias of endoscopies
sub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T)
sub <- subset(sub, tissue == "Colon", invert = T)

table(sub$diseasestate, sub$tissue)

#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"


sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                                  "Normal MAT", "Uninvolved MAT", "CF"))
Idents(sub) <- "newlabels"
#Equal numbers of cells 
disease_idents <- unique(sub$newlabels)
cell.list <- WhichCells(sub, idents = disease_idents, downsample = min(table(sub$newlabels)))
sub_small <- sub[, cell.list]
table(sub_small$newlabels)
Idents(sub_small) <- "newlabels"

#Generate proptable and barplots of fibroblasts in resection 
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
tab = table(sub_small$diseasestate, sub_small$immune.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
theme_classic() + 
labs(title="Fibroblast Disease State Enrichment", 
x="Disease State", y = "Percent (%)", fill = "Cluster") + 
scale_fill_manual(values = colors) + 
theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) + 
facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") + 
theme(strip.placement = "outside", strip.background = element_blank()) + 
theme(strip.text = element_text(size = 15)) + 
theme(axis.title.x = element_blank(), legend.title = element_blank(), title = element_blank()) + 
scale_x_discrete(labels=c('Normal', 'Uninvolved', 'Stricture', 'Normal', "Uninvolved", "CF")) 
ggsave("SuppFigures/lymphoidpropplot.svg")



#### Myeloid ####
global <- readRDS("globalcellatlas.annotated.rds")
Idents(global) <- "manual.annotation"
obj <- subset(global, idents = c("Myeloid Cells", "Macrophages", "Mast Cells"))

##Regenerate and clean up SO 
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, "MyeloidAll/MyeloidAll.rds", compress = F )

# re-perform BPcells
counts.mat = obj[["RNA"]]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class="dgCMatrix" )
write_matrix_dir(mat = counts.mat, dir = "MyeloidAll/MyeloidAll_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MyeloidAll/MyeloidAll_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, "MyeloidAll/MyeloidAll_counts.rds", compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, "MyeloidAll/MyeloidAllupdated_split.rds", compress = F )
plan("multisession", workers = 8)

#obj<- readRDS("MyeloidAll/MyeloidAll_split.rds")
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
# Access the RNA data
tube <- obj[["RNA"]][Layers(obj)[i]]
# Try executing FindVariableFeatures
tryCatch({
# Find variable features
print(Layers(obj)[i])
tube <- CreateSeuratObject(tube)
tube <- NormalizeData(tube)
tube <- FindVariableFeatures(tube)
}, error = function(e) {
# Record error message
error_messages[[i]] <<- conditionMessage(e)
# Record index of the layer where the error occurred
error_indices[[i]] <<- i
# Create a row for the error dataframe
error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
# Append the error row to the error dataframe
error_df <<- rbind(error_df, error_row)
})
}
print(error_df)
write.csv(error_df, "MyeloidAll/Error_list.csv")

error_df <- read.csv("MyeloidAll/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study

#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, "MyeloidAll/MyeloidAll_split_subset.rds", compress = F)

obj <- obj_subset
obj <- FindVariableFeatures(obj, verbose = T )
saveRDS( obj, "MyeloidAll/MyeloidAll_var.rds"  , compress = F )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, "MyeloidAll/MyeloidAll_pca.rds"  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, "MyeloidAll/MyeloidAll_harmony.rds", compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, "MyeloidAll/MyeloidAll_harmony_umap.rds", compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, "MyeloidAll/MyeloidAll_harmony_umap_joined.rds", compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MyeloidAll/MyeloidAll_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.5, 0.1) ) {
obj <- FindClusters(obj, resolution = res )
DimPlot(obj, raster = TRUE, label = T, repel = T)
ggsave(paste0("MyeloidAll/MyeloidAll_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
saveRDS( obj, paste0( "MyeloidAll/MyeloidAll_harmony_clusters_", res, ".rds" ), compress = F ) 
# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.05)
write.csv(markers, paste0("MyeloidAll/MyeloidAll_harmony_clusters_markers_",res,".csv"))
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
group_by(cluster) %>%
slice_max(n = 100, order_by = rank)
write.csv(topmarkers, paste0("MyeloidAll/MyeloidAll_harmony_clusters_topmarkers_",res,".csv"))
}

DefaultAssay(obj) <- "RNA"

###RECLUSTER MYELOID
#Remove clusters
obj <- readRDS('/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman/MyeloidAll/MyeloidAll_harmony_clusters_0.31.rds')
obj <- subset(obj, idents = c(0, 7, 15), invert = T)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_0.2.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "MyeloidAll2/MyeloidAll2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MyeloidAll2/MyeloidAll2_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_split.rds', compress = F )
obj<- readRDS('MyeloidAll2/MyeloidAll2_split.rds')

# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
# Access the RNA data
tube <- obj[['RNA']][Layers(obj)[i]]
# Try executing FindVariableFeatures
tryCatch({
# Find variable features
print(Layers(obj)[i])
tube <- CreateSeuratObject(tube)
tube <- NormalizeData(tube)
tube <- FindVariableFeatures(tube)
}, error = function(e) {
# Record error message
error_messages[[i]] <<- conditionMessage(e)
# Record index of the layer where the error occurred
error_indices[[i]] <<- i
# Create a row for the error dataframe
error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
# Append the error row to the error dataframe
error_df <<- rbind(error_df, error_row)
})
}
print(error_df)
write.csv(error_df, "MyeloidAll2/Error_list.csv")

error_df <- read.csv("MyeloidAll2/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'MyeloidAll2/MyeloidAll2_split_subset.rds', compress = F)

obj <- readRDS('MyeloidAll2/MyeloidAll2_split_subset.rds')
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'MyeloidAll2/MyeloidAll2_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MyeloidAll2/MyeloidAll2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.5, 0.1) ) {
obj <- FindClusters(obj, resolution = res , verbose = T)
DimPlot(obj, raster = TRUE, label = T, repel = T)
ggsave(paste0("MyeloidAll2/MyeloidAll2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
saveRDS( obj, paste0( 'MyeloidAll2/MyeloidAll2_harmony_clusters_', res, '.rds' ), compress = F ) 
# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, paste0("MyeloidAll2/MyeloidAll2_harmony_clusters_markers_",res,".csv"))
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
group_by(cluster) %>%
slice_max(n = 100, order_by = rank)
write.csv(topmarkers, paste0("MyeloidAll2/MyeloidAll2_harmony_clusters_topmarkers_",res,".csv"))
}

###RECLUSTER MYELOID
#Remove clusters
obj <- readRDS('MyeloidAll2/MyeloidAll2_harmony_clusters_0.11.rds')
obj <- subset(obj, idents = c(8, 10), invert = T)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
obj <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = obj@meta.data)
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_0.2.rds', compress = F )

# re-perform BPcells
counts.mat = obj[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "MyeloidAll3/MyeloidAll3_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MyeloidAll3/MyeloidAll3_counts" )
obj = CreateSeuratObject( counts.mat, meta.data = obj@meta.data )
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_counts.rds', compress = F )

# run downstream analysis
obj <- NormalizeData(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_split.rds', compress = F )
obj<- readRDS('MyeloidAll3/MyeloidAll3_split.rds')

# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(obj)) 
for (i in 1:length) {
# Access the RNA data
tube <- obj[['RNA']][Layers(obj)[i]]
# Try executing FindVariableFeatures
tryCatch({
# Find variable features
print(Layers(obj)[i])
tube <- CreateSeuratObject(tube)
tube <- NormalizeData(tube)
tube <- FindVariableFeatures(tube)
}, error = function(e) {
# Record error message
error_messages[[i]] <<- conditionMessage(e)
# Record index of the layer where the error occurred
error_indices[[i]] <<- i
# Create a row for the error dataframe
error_row <- data.frame(layer_index = i, ident = as.character(obj$orig.ident[i]), layer_name = Layers(obj)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
# Append the error row to the error dataframe
error_df <<- rbind(error_df, error_row)
})
}
print(error_df)
write.csv(error_df, "MyeloidAll3/Error_list.csv")

error_df <- read.csv("MyeloidAll3/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(obj$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(obj$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- obj
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'MyeloidAll3/MyeloidAll3_split_subset.rds', compress = F)

obj <- readRDS('MyeloidAll3/MyeloidAll3_split_subset.rds')
obj <- FindVariableFeatures(obj, verbose = T )
obj <- ScaleData(obj, verbose = T)
obj <- RunPCA(obj, npcs = 30, verbose = T)
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_pca.rds'  , compress = F )
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_harmony.rds', compress = F ) 
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_harmony_umap.rds', compress = F ) 
obj <- JoinLayers(obj)
saveRDS( obj, 'MyeloidAll3/MyeloidAll3_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(obj, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MyeloidAll3/MyeloidAll3_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.5, 0.1) ) {
obj <- FindClusters(obj, resolution = res , verbose = T)
DimPlot(obj, raster = TRUE, label = T, repel = T)
ggsave(paste0("MyeloidAll3/MyeloidAll3_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
saveRDS( obj, paste0( 'MyeloidAll3/MyeloidAll3_harmony_clusters_', res, '.rds' ), compress = F ) 
# join and find markers
obj[["RNA"]]$data <- as(object = obj[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(obj, max.cells.per.ident = 1000, only.pos = T, return.thresh = 0.01)
write.csv(markers, paste0("MyeloidAll3/MyeloidAll3_harmony_clusters_markers_",res,".csv"))
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
group_by(cluster) %>%
slice_max(n = 100, order_by = rank)
write.csv(topmarkers, paste0("MyeloidAll3/MyeloidAll3_harmony_clusters_topmarkers_",res,".csv"))
}

##Annotation of cell types
obj <- readRDS("MyeloidAll13_harmony_clusters_0.15.rds")
old <- obj
obj <- subset(obj, idents = 13, invert = T) #remove low quality
# manual annotations for entire dataset
DefaultAssay(obj) <- "RNA"
obj$immune.annotation <- "tbd"
obj$immune.annotation[obj$seurat_clusters == 0] <- "Mast Cells CPA3+"
obj$immune.annotation[obj$seurat_clusters == 1] <- "Macrophages LYVE1+"
obj$immune.annotation[obj$seurat_clusters == 2] <- "Monocytes CD300E+"
obj$immune.annotation[obj$seurat_clusters == 3] <- "Neutrophils CSF3R+"
obj$immune.annotation[obj$seurat_clusters == 4] <- "Macrophages C1QChi"
obj$immune.annotation[obj$seurat_clusters == 5] <- "cDC2 CD1C+"
obj$immune.annotation[obj$seurat_clusters == 6] <- "mregDC LAMP3+"
obj$immune.annotation[obj$seurat_clusters == 7] <- "Cycling Myeloid MKI67+"
obj$immune.annotation[obj$seurat_clusters == 8] <- "cDC1 CLEC9A+"
obj$immune.annotation[obj$seurat_clusters == 9] <- "Macrophages LYVE1+"
obj$immune.annotation[obj$seurat_clusters == 10] <- "Activated DC IL4I1+"
obj$immune.annotation[obj$seurat_clusters == 11] <- "Macrophages CXCL2hi"
obj$immune.annotation[obj$seurat_clusters == 12] <- "pDC GZMB+"

obj$immune.annotation <- factor(obj$immune.annotation, levels = c("Monocytes CD300E+", 
                                              "Macrophages LYVE1+", "Macrophages C1QChi", "Macrophages IL4I1+", "Macrophages CXCL2hi", 
                                              "cDC1 CLEC9A+", "cDC2 CD1C+", "pDC GZMB+", "mregDC LAMP3+", 
                                              "Mast Cells CPA3+", "Cycling Myeloid",  "Neutrophils CSF3R+"))
Idents(obj) <- "immune.annotation"

colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
"#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
dimplot <- DimPlot(obj, raster = T, repel = T, cols = colors)

#Correct loss of metadata
global <- readRDS("globalcellatlas.annotated.rds")
Idents(obj) <- "diseasestate"
obj$diseasestate <- global$diseasestate[colnames(global) %in% colnames(obj)]
saveRDS(obj, "myeloid.annotated.rds")
revobj <- obj
revobj$immune.annotation <- factor(revobj$immune.annotation, levels = rev(c("Monocytes CD300E+", 
                                                                            "Macrophages LYVE1+", "Macrophages C1QChi", "Macrophages IL4I1+", "Macrophages CXCL2hi", 
                                                                            "cDC1 CLEC9A+", "cDC2 CD1C+", "pDC GZMB+", "mregDC LAMP3+", 
                                                                            "Mast Cells CPA3+", "Cycling Myeloid",  "Neutrophils CSF3R+")))
Idents(revobj) <- "immune.annotation"
genes <- c("CD300E", "LYVE1", "C1QC", "IL4I1", "CXCL2", "CLEC9A", "CD1C", "GZMB", "LAMP3", "CPA3", "MKI67", "CSF3R")
dotplot <- DotPlot(revobj, features = genes, cols = "RdYlBu") +   theme(axis.title.y = element_blank(), 
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
ggsave("SuppFigures/myeloidallumap.svg")

##Prop tables 
#Generate comparable dataset
sub <- subset(obj, diseasestate == "Ileostomy Bowel" | diseasestate == "Ileostomy MAT", invert = T) #Remove ileostomy samples
sub <- subset(sub, specimen_type == "Resection") #Only compare full thickness due to mucosal bias of endoscopies
sub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T)
sub <- subset(sub, tissue == "Colon", invert = T)

table(sub$diseasestate, sub$tissue)

#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"


sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninvolved Bowel", "Stricture", 
                                                      "Normal MAT", "Uninvolved MAT", "CF"))
Idents(sub) <- "newlabels"
#Equal numbers of cells 
disease_idents <- unique(sub$newlabels)
cell.list <- WhichCells(sub, idents = disease_idents, downsample = min(table(sub$newlabels)))
sub_small <- sub[, cell.list]
table(sub_small$newlabels)
Idents(sub_small) <- "newlabels"

#Generate proptable and barplots of fibroblasts in resection 
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
tab = table(sub_small$diseasestate, sub_small$immune.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Disease State Enrichment", 
       x="Disease State", y = "Percent (%)", fill = "Cluster") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) + 
  facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") + 
  theme(strip.placement = "outside", strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 15)) + 
  theme(axis.title.x = element_blank(), legend.title = element_blank(), title = element_blank()) + 
  scale_x_discrete(labels=c('Normal', 'Uninvolved', 'Stricture', 'Normal', "Uninvolved", "CF")) 
ggsave("SuppFigures/myeloidpropplot.svg")

                            

