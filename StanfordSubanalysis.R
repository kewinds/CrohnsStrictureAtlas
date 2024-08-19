setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/LongakerOnly"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers)
library(scales); library(SeuratObject); library(irlba)
options(Seurat.object.assay.version = "v5")

#Import the global object and subset out fibroblasts from Longaker study 
rm(list = ls())
global <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman/sketch_v5_all/humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds")
fibro <- subset(global, subset = (study =="LongakerBowel" | study=="LongakerMes"))
fibro <- subset(fibro, manual.annotation == "Fibroblasts")
table(fibro$study)
new <- JoinLayers(fibro)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
new <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( new, 'newFibro1/fibro1.rds', compress = F )
table(new$study)

# re-perform BPcells
counts.mat = new[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro1/fibro1_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro1/fibro1_counts" )
new = CreateSeuratObject( counts.mat, meta.data = new@meta.data )
saveRDS( new, 'newFibro1/fibro1_counts.rds', compress = F )
table(new$study)

rm(list = ls())
fibro <- readRDS('newFibro1/fibro1_counts.rds')
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
FeaturePlot(fibro, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = FALSE)
ggsave("newFibro1/fibro1_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
for (i in seq(0.1, 1.0, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, raster = FALSE, label = TRUE)
  ggsave(paste0("newFibro1/fibro1_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("newFibro1/fibro1_harmony_umap_res",i,".rds"), compress = FALSE)
}

# join and find markers
fibro[["RNA"]]$data <- as(object = fibro[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro1/fibro1_harmony_umap_feature_0.2_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro1/fibro1_harmony_umap_feature_0.2_topmarkers.csv")

#Cluster 3 has low quality cells, 7, 8, 9, 10 not fibroblasts so we remove and repeat the above pipeline
rm(list = ls())
fibro <- readRDS("newFibro1/fibro1_harmony_umap_res0.2.rds")
fibro <- subset(fibro, idents = c(3,7,8,9,10), invert = T) #remove cluster 3
table(fibro$study)
new <- JoinLayers(fibro)

# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
new <- CreateSeuratObject(counts = fibro[["RNA"]]$counts, meta.data = fibro@meta.data)
saveRDS( new, 'newFibro1B/fibro1B.rds', compress = F )
table(new$study)

# re-perform BPcells
counts.mat = new[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "newFibro1B/fibro1B_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "newFibro1B/fibro1B_counts" )
new = CreateSeuratObject( counts.mat, meta.data = new@meta.data )
saveRDS( new, 'newFibro1B/fibro1B_counts.rds', compress = F )
table(new$study)

rm(list = ls())
fibro <- readRDS('newFibro1B/fibro1B_counts.rds')
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
FeaturePlot(fibro, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = FALSE)
ggsave("newFibro1B/fibro1B_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
for (i in seq(0.1, 1.0, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, raster = FALSE, label = TRUE)
  ggsave(paste0("newFibro1B/fibro1B_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("newFibro1B/fibro1B_harmony_umap_res",i,".rds"), compress = FALSE)
}

# join and find markers
rm(list = ls())
fibro <- readRDS('newFibro1B/fibro1B_harmony_umap_res0.4.rds')
fibro[["RNA"]]$data <- as(object = fibro[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "newFibro1B/fibro1B_harmony_umap_feature_0.4_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "newFibro1B/fibro1B_harmony_umap_feature_0.4_topmarkers.csv")

#Label clusters and generate figures
fibro <- subset(fibro, idents = 11, invert = T) #Remove mesothelial cells
# manual annotations for entire dataset
DefaultAssay(fibro) <- "RNA"
Idents(fibro) <- "seurat_clusters"
fibro$longaker.annotation <- "tbd"
fibro$longaker.annotation[fibro$seurat_clusters == 0] <- "Mechanosensitive LRRC7+"
fibro$longaker.annotation[fibro$seurat_clusters == 1] <- "Mechanosensitive CTHRC1+"
fibro$longaker.annotation[fibro$seurat_clusters == 2] <- "Mechanosensitive LRRC7+"
fibro$longaker.annotation[fibro$seurat_clusters == 3] <- "Antigen-Presenting CD74+"
fibro$longaker.annotation[fibro$seurat_clusters == 4] <- "Steady-State ADAMTS16+"
fibro$longaker.annotation[fibro$seurat_clusters == 5] <- "Steady-State FMO2+"
fibro$longaker.annotation[fibro$seurat_clusters == 6] <- "Steady-State PI16+"
fibro$longaker.annotation[fibro$seurat_clusters == 7] <- "Inflammatory ADAMDEC1+"
fibro$longaker.annotation[fibro$seurat_clusters == 8] <- "Steady-State FMO2+"
fibro$longaker.annotation[fibro$seurat_clusters == 9] <- "Inflammatory MMP3+"
fibro$longaker.annotation[fibro$seurat_clusters == 10] <- "Steady-State VIT+"
Idents(fibro) <- "longaker.annotation"

saveRDS(fibro, "longakerfibro.finalannotations.rds")
#Reorder factor levels 
fibro$longaker.annotation <- factor(x = fibro$longaker.annotation, levels = c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+","Antigen-Presenting CD74+", "Inflammatory ADAMDEC1+", 
                                                                              "Inflammatory MMP3+", "Steady-State ADAMTS16+", "Steady-State VIT+", "Steady-State FMO2+", "Steady-State PI16+"))
Idents(fibro) <- "longaker.annotation"
library(RColorBrewer)
#num_colors <- length(levels(fibro))
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138",
                     "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
dimplot <- DimPlot(fibro, raster = F, cols = colors)

#Create Dotplot
revfibro <- fibro
revfibro$longaker.annotation <- factor(x = revfibro$longaker.annotation, levels = rev(c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+","Antigen-Presenting CD74+", "Inflammatory ADAMDEC1+", 
                                                                                    "Inflammatory MMP3+", "Steady-State ADAMTS16+", "Steady-State VIT+", "Steady-State FMO2+", "Steady-State PI16+")))
Idents(revfibro) = "longaker.annotation"
#levels(revfibro) = c("ssCDF VIT+", "ssCDF ADAMTS16+", "ssCDF FMO2+", "ssCDF PI16+", 
#                             "apCDF CD74+", "iCDF MMP3+", "iCDF ADAMDEC1+", "msCDF LRRC7+", "msCDF POSTN+")
genes <- c("CTHRC1", "LRRC7", "CD74", "ADAMDEC1", "MMP3", "ADAMTS16", "VIT", "FMO2", "PI16")
dotplot<- DotPlot(revfibro, features = genes, cols = "RdYlBu") +   theme(axis.title.y = element_blank(), 
                                                        axis.title.x = element_blank(),
                                                        axis.text.y = element_blank(), 
                                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                        legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                        legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                        legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                        legend.title = element_text(size=6), #change legend title font size
                                                        legend.text = element_text(size=6)) #change legend text font size)
combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -40, unit = "pt")))) + dotplot 
combined_plot
ggsave(paste0("Figures/longakerfibroblastUMAP.jpg"), width = 8.07, height = 3.16, units = "in")

##### Add module scores to get mechanosensitive signatures
#1. YAP gene signature 
library(readxl)
Yapgenes <-as.data.frame(read_excel("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman/AltYAP_list2.xlsx",col_names = "Genes"))
Yapgenes <- Yapgenes[Yapgenes$Genes %in% rownames(revfibro), , drop = FALSE] #Remove genes not present in our dataset
#percent_expressed <- Percent_Expressing(fibro, entire_object = T, features = Yapgenes) #Calculate percent of cells that express
#lowgenes <- percent_expressed[percent_expressed$All_Cells < 1, , drop = FALSE] #Subset genes that are expressed in <10% of fibroblasts
#finalyaplist <- Yapgenes[Yapgenes$Genes %in% rownames(lowgenes), ] 
revfibro <- AddModuleScore(revfibro, features = Yapgenes, name = "Yapsig") #Get score for YAP targets
#2. YAP components signature
components <- list(c("YAP1", "CTNNB1", 
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))
revfibro <- AddModuleScore(revfibro, features = components, name = "Components")
#3. ECM protein signature 
ecmgenes <- list(c("COL1A1", "COL1A2", "COL3A1", "COL4A2", "COL5A1", "COL5A3", 
                   "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", 
                   "LAMB1", "LAMB2", "LAMB3", 
                   "LAMC1", "LAMC2", "LTBP1", 
                   "MFAP2", "MFAP5", "FN1", 
                   "POSTN", "SERPINA1", "SERPINA5", 
                   "TGFBI", "TNC", "VCAN"))
revfibro <- AddModuleScore(revfibro, features = ecmgenes, name = "ECMgenes")
DotPlot(revfibro, features = c("Components1", "Yapsig1", "ECMgenes1"), scale.max = 100, scale.min = 0, cols = "RdBu") + theme(axis.title.x = element_blank(),
                                                                                                                              axis.title.y = element_blank(), 
                                                                                                                              axis.text.x = element_text(size = 15),
                                                                                                                              axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
                                                                                                                              legend.key.size = unit(0.8, 'cm'), #change legend key size
                                                                                                                              legend.key.height = unit(0.4, 'cm'), #change legend key height
                                                                                                                              legend.key.width = unit(0.4, 'cm'), #change legend key width
                                                                                                                              legend.title = element_text(size=9), #change legend title font size
                                                                                                                              legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c('YAP Components', "YAP Signature", "ECM Signature"))
ggsave("Figures/YAPsignatures.jpg" , width = 10, height = 4, units = "in")

#Create new metadata category called general to generalize regions
Idents(fibro) <- "longaker.annotation"
#levels(fibro) <- rev(c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
#                       "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
#                       "mCDF LRRC7+", "mCDF CTHRC1+"))
fibro$general <- "NA"
fibro$general[fibro$tissue == "MAT"] <- "MAT"
fibro$general[fibro$tissue == "SB"] <- "Bowel"
fibro$general[fibro$tissue == "Colon"] <- "Bowel"
sub <- subset(fibro, subset = (diseasestate == "Ileostomy Bowel" | diseasestate == "Ileostomy MAT"), invert = T) #Remove ileostomy samples
sub <- subset(sub, subset = (specimen_type == "Resection")) #Only compare full thickness due to mucosal bias of endoscopies
sub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T)
table(sub$diseasestate, sub$tissue)
DimPlot(sub, cols = colors, split.by = "general") 
ggsave("Figures/BowelvsMATfibro.svg", width = 6.8, height = 3, units = "in")

#Module scores for showing regional YAP differences
sub <- AddModuleScore(sub, features = Yapgenes, name = "Yapsig")
sub <- AddModuleScore(sub, features = components, name = "Components")
sub <- AddModuleScore(sub, features = ecmgenes, name = "ECMgenes")
subsub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T) #Eliminate stricture inflamed for comparable comparisons
DotPlot(subsub, features = c("Components1", "Yapsig1", "ECMgenes1"), group.by = "diseasestate", scale.max = 100, scale.min = 0, cols = "RdBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size))
  scale_x_discrete(labels=c('YAP Components', "YAP Signature", "ECM Signature")) + 
  scale_y_discrete(labels=c("Uninvolved", "Stricture", "Uninvolved", "CF"))
ggsave("Figures/RegionalYAPsignatures.jpg", width = 10, height = 4, units = "in")

#Change labels for clarity
mes <- subset(fibro, tissue == "MAT")
mes <- subset(mes, diseasestate == "Ileostomy MAT", invert = T)
mes$newlabels <- "NA"
mes$newlabels[mes$diseasestate == "Uninvolved"] <- "Uninvolved"
mes$newlabels[mes$diseasestate == "Involved"] <- "CF"
mes$newlabels <- factor(x = mes$newlabels, levels = c("Uninvolved", "CF"))
#num_colors <- length(levels(mes))
#colors <- brewer.pal(num_colors, "RdYlBu")
#Reorder levels
DimPlot(mes, split.by = "newlabels", cols = colors)
ggsave("Figures/MAT_MSexpansion.jpg", width = 8.5, height = 3, units = "in")

#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"


#sub$fibro.annotation <- factor(x = sub$fibro.annotation, levels = rev(c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
 #                                                                       "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
 #                                                                       "mCDF LRRC7+", "mCDF CTHRC1+")))
sub$newlabels <- factor(x = sub$newlabels, levels = c("Uninvolved Bowel", "Stricture", 
                                                      "Uninvolved MAT", "CF"))
#Generate proptable and barplots of fibroblasts in resection 
tab = table(sub$newlabels, sub$longaker.annotation)
ptab = prop.table(tab, 1)*100
tissue <- c("Bowel", "Bowel","MAT", "MAT")
ptab = as.data.frame(ptab)
ptab$tissue <- tissue
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Disease State Enrichment", 
       x="Disease State", y = "Percent (%)", fill = "Cluster") + 
  scale_fill_manual(values = colors) + 
  scale_x_discrete(labels=c('Uninvolved', 'Stricture',"Uninvolved", "CF")) +
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) + 
  facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") + 
  theme(strip.placement = "outside", strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 15)) + 
  theme(axis.title.x = element_blank(), legend.title = element_blank(), title = element_blank()) 
ggsave("Figures/proptablefibroblasts.jpg", width = 6, height = 4)

library(patchwork)
#Anchor label transfer 
rm(list = ls())
ref = readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman/finalfibroblast.annotations.rds")
query = readRDS("longakerfibro.finalannotations.rds")
anchors <- FindTransferAnchors(reference = ref, query = query, dims = 1:30,
                               reference.reduction = "pca")
saveRDS(anchors, "Subanalysisanchors.rds")

predictions <- TransferData(anchorset = anchors, refdata = ref$fibro.annotation, dims = 1:30, prediction.assay = T)
query[["predictions"]] <- predictions
DefaultAssay(query) <- "predictions"
FeaturePlot(query, features = c("Mechanosensitive CTHRC1+"))
query <- AddMetaData(query, metadata = t(GetAssayData(predictions)))
query$predicted.id <- factor(x = query$predicted.id, levels = rev(c("Steady-State PI16+", "Steady-State FMO2+", "Steady-State F3+","Steady-State GREM1+", "Steady-State KCNN3+", 
                                                                                "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", "Antigen-Presenting CCL19+", 
                                                                                "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")))

DimPlot(query, reduction = )

Idents(query) = "predicted.id"
num_colors <- length(levels(query))
colors <- brewer.pal(num_colors, "RdYlBu")
dimplot <- DimPlot(query, cols = colors) + ggtitle("Predicted Metadata Clusters")
featureplot1 <- FeaturePlot(query, features = c("mCDF CTHRC1+")) + ggtitle("mCDF CTHRC1+") + guides(color = guide_legend(title = "Probability")) + 
          theme(legend.title =element_text(size=10)) + 
  scale_colour_gradientn(name = "Score", colors = c("#0d0887", "#cc4778", "#f0f921"))
featureplot2 <-FeaturePlot(query, features = c("mCDF LRRC7+")) + ggtitle("mCDF LRRC7+") + guides(color = guide_legend(title = "Probability")) + 
  theme(legend.title =element_text(size=10)) + 
  scale_colour_gradientn(name = "Score", colors = c("#0d0887", "#cc4778", "#f0f921"))
dimplot + featureplot1 + featureplot2 + plot_layout(guides = "collect")

ggsave("Figures/ALTlongakertometa.jpg", width = 12.5, height = 3)

anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)
saveRDS(anchors, "Subanalysisanchors.rds")

pbmc3k <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = ref,
  refdata = list(fibro.annotation = "fibro.annotation"),
  reference.reduction = "pca", 
  reduction.model = "umap.harmony"
)

Idents(pbmc3k) <- "predicted.fibro.annotation"
levels(pbmc3k) <- c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+","Antigen-Presenting CD74+", "Inflammatory ADAMDEC1+", 
                    "Inflammatory MMP3+", "Steady-State ADAMTS16+", "Steady-State VIT+", "Steady-State FMO2+", "Steady-State PI16+")
# pbmc3k$predicted.fibro.annotation <- factor(pbmc3k$predicted.fibro.annotation, levels = c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+","Antigen-Presenting CD74+", "Inflammatory ADAMDEC1+", 
#                                                                                                     "Inflammatory MMP3+", "Steady-State ADAMTS16+", "Steady-State VIT+", "Steady-State FMO2+", "Steady-State PI16+"))
Idents(pbmc3k) <- "predicted.fibro.annotation"
sub <- subset(pbmc3k, refUMAP_1 < -2)
DimPlot(sub, reduction = "ref.umap", cols = colors)


library(ggplot2)
library(RColorBrewer)
library(patchwork)

fp1 <- FeaturePlot(sub, features = c("Mechanosensitive CTHRC1+"), reduction = "ref.umap") + 
  ggtitle("Mechanosensitive CTHRC1+") + 
  guides(color = guide_legend(title = "Probability")) + 
  theme(title = element_text(size = 10))+ 
  theme(legend.title = element_text(size = 10)) + 
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))
fp2 <- FeaturePlot(sub, features = c("Mechanosensitive LRRC7+"), reduction = "ref.umap") + 
  ggtitle("Mechanosensitive LRRC7+") + 
  guides(color = guide_legend(title = "Probability")) + 
  theme(title = element_text(size = 10))+ 
  theme(legend.title = element_text(size = 10)) + 
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))

fp1 | fp2 & plot_layout(guides = "collect")


#Now we perform the analysis without ileostomy samples 
fibro <- subset(fibro, subset = diseasestate %in% c("Ileostomy MAT", "Ileostomy Bowel"), invert = T) #Remove ileostomy samples 
# re-generate seurat fibroect and remove extraneous data
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

rm(list = ls())
fibro <- readRDS('newFibro2/fibro2_counts.rds')
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
FeaturePlot(fibro, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = FALSE)
ggsave("newFibro2/fibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)
for (i in seq(0.1, 1.0, 0.1)) { 
  fibro <- FindClusters(fibro, resolution = i)
  DimPlot(fibro, raster = FALSE, label = TRUE)
  ggsave(paste0("newFibro2/fibro2_harmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(fibro, paste0("newFibro2/fibro2_harmony_umap_res",i,".rds"), compress = FALSE)
}

