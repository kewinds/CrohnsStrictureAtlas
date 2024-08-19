#KEBR code to integrate visium objects

#Install this if Matrix package fails you: install.packages("SeuratObject", repos = c("https://mojaveazure.r-universe.dev", getOption("repos")))


setwd("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/r_packages",.libPaths()))

library(Seurat); library(ggplot2); library( dplyr ) ;library(Matrix); library(SeuratWrappers); library(scales); library(patchwork); 
library(matrixStats); library(ggplot2)
library(cowplot)
library(grid)
library(readxl)
options(Seurat.object.assay.version = "v5")

rm(list = ls())
#### Integration of Visium objects using Seurat V3 ####
spat9 = Load10X_Spatial('StrictureOutputs/9/outs/', slice = 'slice9'); spat9$orig.ident = 'slice_9'
#spat10 = Load10X_Spatial('StrictureOutputs/10/outs/', slice = 'slice10'); spat10$orig.ident = 'slice_10'
spat19 = Load10X_Spatial('StrictureOutputs/19/outs/', slice = 'slice19'); spat19$orig.ident = 'slice_19'
#spatTO = Load10X_Spatial('StrictureOutputs/TO/outs/', slice = 'sliceTO'); spatTO$orig.ident = 'slice_TO'

# #Manually add data to spat9 due to CellRanger outputs
spat9@images[["slice9"]]@coordinates[["tissue"]] <- as.integer(spat9@images[["slice9"]]@coordinates[["tissue"]])
spat9@images[["slice9"]]@coordinates[["row"]] <- as.integer(spat9@images[["slice9"]]@coordinates[["row"]])
spat9@images[["slice9"]]@coordinates[["col"]] <- as.integer(spat9@images[["slice9"]]@coordinates[["col"]])
spat9@images[["slice9"]]@coordinates[["imagerow"]] <- as.integer(spat9@images[["slice9"]]@coordinates[["imagerow"]])
spat9@images[["slice9"]]@coordinates[["imagecol"]] <- as.integer(spat9@images[["slice9"]]@coordinates[["imagecol"]])

#Manually add data to spat19 due to CellRanger outputs
spat19@images[["slice19"]]@coordinates[["tissue"]] <- as.integer(spat19@images[["slice19"]]@coordinates[["tissue"]])
spat19@images[["slice19"]]@coordinates[["row"]] <- as.integer(spat19@images[["slice19"]]@coordinates[["row"]])
spat19@images[["slice19"]]@coordinates[["col"]] <- as.integer(spat19@images[["slice19"]]@coordinates[["col"]])
spat19@images[["slice19"]]@coordinates[["imagerow"]] <- as.integer(spat19@images[["slice19"]]@coordinates[["imagerow"]])
spat19@images[["slice19"]]@coordinates[["imagecol"]] <- as.integer(spat19@images[["slice19"]]@coordinates[["imagecol"]])

#Add in Global features
spat9_domains <- read.csv(file = "Domains/spat9_domains.csv", header = TRUE, row.names = 1)
#spat10_domains <- read.csv(file = "Domains/spat10_domains.csv", header = TRUE, row.names = 1)
spat19_domains <- read.csv(file = "Domains/spat19_domains.csv", header = TRUE, row.names = 1)
#spatTO_domains <- read.csv(file = "Domains/spatT0_domains.csv", header = TRUE, row.names = 1)

spat9$Domain <- spat9_domains$Domain #This sample came only from the bowel wall. 
#spat10$Domain <- spat10_domains$Domains
spat19$Domain <- spat19_domains$Domains
#spatTO$Domain <- spatTO_domains$Domains

#Import function to rotate plots 
rotate_image <- function(p, rot_angle, scale) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle, width = scale, height = scale)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}

#Visualize H&E sections first 
p1 <- SpatialDimPlot(merged, images = "slice19", alpha = 0) + 
  ggtitle("Outer Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) + NoLegend()

p1 <- rotate_image(p1, 118, 1)
#p2 <- SpatialDimPlot(spat10, group.by = "Domain")
p3 <- SpatialDimPlot(merged, images = "slice9", group.by = "Domain", alpha = 0) 
  ggtitle("Inner Region") + theme(plot.title = element_text(hjust = 0.5, size = 15))  + NoLegend()

p3 <- rotate_image(p3, -153, 0.83)

p3 + p1

#Specify colors manually
Idents(merged) <- "Domain"
colors <- c("#FDAE61", "#D7191C", "#ABD9E9")
names(colors) <- unique(Idents(merged))

Idents(merged) = "stromal"
levels(merged) <- c("Stromal 1", "Stromal 2", "Stromal 3", "Stromal 4")
names(colors) <- unique(Idents(merged))
#Visualize the domains and double check annotations
p1 <- SpatialDimPlot(merged, images = "slice19", group.by = "Domain", image.alpha = 0, cols = colors) + 
ggtitle("Outer Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p1 <- rotate_image(p1, 118, 1)
#p2 <- SpatialDimPlot(spat10, group.by = "Domain")
p3 <- SpatialDimPlot(merged, images = "slice9", group.by = "Domain", image.alpha = 0, cols = colors) + 
  ggtitle("Inner Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p3 <- rotate_image(p3, -153, 0.83)
#p4 <- SpatialDimPlot(spatTO, group.by = "Domain")
#p1 + p2 + p3 + p4 + plot_layout(guides = "collect")
p3 / p1 + plot_layout(guides = "collect")
ggsave("Figures/unprocessed.spatialdimplots.jpg", width = 7, height = 7, units = "in")

#Annotation of stromal areas (generic)
merged$stromal[merged2$seurat_clusters == 0] <- "Stromal 1"
merged$stromal[merged2$seurat_clusters == 1] <- "Stromal 2"
merged$stromal[merged2$seurat_clusters == 2] <- "Stromal 3"
merged$stromal[merged2$seurat_clusters == 3] <- "Stromal 4"
Idents(merged) = "stromal"
levels(merged) <- c("Stromal 1", "Stromal 2", "Stromal 3", "Stromal 4")
library(RColorBrewer)
num_colors <- length(levels(merged))
colors <- brewer.pal(num_colors, "RdYlBu")
colors <- c("#FDAE61" ,"#ABD9E9", "#D7191C", "#2C7BB6")
p1 <- DimPlot(merged, cols = colors)
p2 <- DimPlot(merged, group.by = "Domain", cols =c("#ABD9E9", "#FDAE61", "#D7191C") ) + ggtitle(NULL)
p1 + p2
ggsave("Figures/mergeddimplot.jpg", width = 10, height = 4)
colors <- c("#FDAE61","#D7191C", "#2C7BB6", "#ABD9E9")
names(colors) <- unique(Idents(merged))
p1 <- SpatialDimPlot(merged, images = "slice19", cols = colors, image.alpha = 0)
p1 <- rotate_image(p1, 118, 1)
p2 <- SpatialDimPlot(merged, images = "slice9", cols = colors, image.alpha = 0)
p2 <- rotate_image(p2, -153, 0.83) 
p1 + p2
ggsave("Figures/spatialdimplot.jpg", width = 7, height = 7)
#Save files
saveRDS(spat9, "readyForSeurat/Final/spat9_preprocess.rds")
#saveRDS(spat10, "readyForSeurat/spat10_preprocess.rds")
saveRDS(spat19, "readyForSeurat/Final/spat19_preprocess.rds")
#saveRDS(spatTO, "readyForSeurat/spatTO_preprocess.rds")

#Preprocess the objects and QC 
obj_list <- list.files(paste0("readyForSeurat/Final/"),"*.rds", recursive = TRUE)
obj_list <- substring(obj_list,0,nchar(obj_list)-15)
write.csv(obj_list, "processed/obj_list_01_preQC.csv")
for (i in 1:length(obj_list)) {
  objName <- obj_list[i]
  x <- readRDS(paste0("readyForSeurat/Final/",objName, "_preprocess.rds"))
  # QC
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  VlnPlot(x, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-01-pre-filt.jpg"))
  x <- subset(x, subset = nFeature_Spatial > 200  & percent.mt < 20 & Domain != "Border")
  VlnPlot(x, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
  ggsave(paste0("01_qc/", gsub("/","_",objName), "_qc-02-post-filt.jpg"))
  saveRDS(x, paste0("processed/",objName, "_01_qc.rds"), compress = FALSE)
}
rm(list = ls())
spat9 <- readRDS("processed/spat9_01_qc.rds")
spat19 <- readRDS("processed/spat19_01_qc.rds")
#Visualize the domains and double check annotations
p1 <- SpatialDimPlot(spat9, group.by = "Domain") 
#p2 <- SpatialDimPlot(spat10, group.by = "Domain")
p3 <- SpatialDimPlot(spat19, group.by = "Domain")
#p4 <- SpatialDimPlot(spatTO, group.by = "Domain")
#p1 + p2 + p3 + p4 + plot_layout(guides = "collect")
p1 + p3 + plot_layout(guides = "collect")
ggsave("processed/spatialdimplots.jpg")

#Integration of objects 
# WARNING: each sample must be a separate object to properly perform QC! 
obj_list <- list.files(paste0("processed/"),"*qc.rds", recursive = TRUE)
spatial.list <- list()
for (obj in obj_list) {
  x <- readRDS(paste0("processed/",obj))
  spatial.list <- append(spatial.list, x)
}
merged <- merge(spatial.list[[1]],spatial.list[2:length(spatial.list)]) #Matrix package sensitive. 
merged #check object to make sure layers are not separated
# # switch to seurat v5 for creating final object
# counts.mat <- GetAssayData(merged)
# merged[["Spatial"]]$counts <- counts.mat
# saveRDS(merged, "Integration/spat.merged_postprocess_V3.rds")
# 
# #Split by tube 
# options(Seurat.object.assay.version = "v5")
# merged <- readRDS("Integration/spat.merged_postprocess_V3.rds")
# merged <- UpdateSeuratObject(merged) #update object to be compatible in V5
merged[["Spatial"]] <- split(merged[["Spatial"]], f = merged$orig.ident)
saveRDS(merged, "Integration/spat.merged_postprocess_V5.rds")

####SCT + Harmony Integration ####
rm(list = ls())
merged <- readRDS("Integration/spat.merged_postprocess_V5.rds")
merged <- SCTransform(merged, assay = "Spatial") #Fails if matrixStats is not 1.1.0. Make sure you downgrade. 
merged <- RunPCA(merged)
merged <- IntegrateLayers(merged, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = T)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
#for (i in seq(0.03, 0.1, 0.01)) { 
for (i in seq(0.1, 1.5, 0.1)) { 
  merged <- FindClusters(merged, resolution = i)
  DimPlot(merged, label = TRUE)
  ggsave(paste0("Integration/UMAP_SCT_Harmony/visium_merged_SCTHarmony_umap_res",i,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS(merged, paste0("Integration/UMAP_SCT_Harmony/visium_merged_SCTHarmony_umap_res",i,".rds"), compress = FALSE)
}
# join and find markers
slot(object = merged@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial" #rename umi.assay to "Spatial" 
merged <- PrepSCTFindMarkers(merged)
markers <- FindAllMarkers(merged, only.pos = T, logfc.threshold = 0.25)
write.csv(markers, "Integration/UMAP_SCT_Harmony/SCTharmony_clusters_0.7_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Integration/UMAP_SCT_Harmony/SCTharmony_clusters_0.7_topmarkers.csv")

merged <- FindSpatiallyVariableFeatures(merged, assay = "SCT", features = VariableFeatures(merged)[1:1000],
                                       selection.method = "moransi")



#### Interpret spatial domains ####
#Apply module scores to object 

#1. YAP gene signature 
Yapgenes <- c('ACAT1','ADAM9','CCN1','KLF5','ACTA2','COL4A1','LIF','LACTB','AMOTL2','IL11','SORCS3','ANKRD1','CCN2','C11orf96','LBH','ANLN','HBEGF','EDN1','MAP2K6','AURKB','HSPG2','ATF3','ARHGAP42','BIRC5','LAMC1','F3','ADRA1B','CASP3','LPL','SERPINE1','CRIM1','CASP8','IER3','C18orf63','CAVIN2','THBS1','TRIB1','CDH2','CCNF','TIMP2','FAM178A','CDC42','VCAM1','KLF10','SLC38A2','CDK1','CALM2','SMAD7','ABI3BP','CDK6','COL1A1','FZD8','TBX18','COL1A2','RUNX1','DYNLRB2','COL5A2','FGF2','NRD1','DIAPH3','COL6A2','0','ZRANB1','DICER1','CXCL12','TNS3','E2F1','FBN1','MUSK','INHBA','USP44','EGFR','LAMA2','PMS1','FGF1','LAMB1','PRELID2','FN1','TFPI','JDP2','FOXM1','VEGFA','AL035681.1','GAB1','BDNF','BASP1','IGF1R','CCL2','QSOX1','LATS2','COL4A2','DDHD1','LDHA','CSF1','CITED2','LIFR','SEMA3E','DUSP10','MCM6','ANXA1','OR51B6','MYC','BGN','EFHC1','MYL9','COL11A1','MAGI2','NDUFB3','COL14A1','JARID2','NEGR1','COL5A1','CCDC132','NEK2','COL6A1','ARNT2','NF2','ANGPT1','PELO','PLK1','ADAM12','CENPV','PMAIP1','CCL9','RPL22L1','POLA1','PTN','SYNDIG1','PPP1R3B','SEMA3A','RAD18','PTEN','THBS2','TSPAN3','RAD51','HMGB1','TBC1D2','SLC2A3','IRS1','SNAI2','NUAK2','TAGLN','TRPC3','TEAD4','GAS6','VIM','SNTB1','WWC1','MTMR10','ZEB1','BMP5','FRY','IGFBP3','AC112715.2','SUSD1','RUNX2','ROBO2','AC023590.1','FBXO47','FAM107A','RBMS3','SH3RF1','SYT16','XPNPEP1','CTGF','PIK3C2G','MED13L','RND3','CDH4','TSPAN18','NEDD9','KPNA3','SH2D4A','TRIO','DPYD','LRRTM1','ZNF469','RP11-553A10.1','GLIS3','SLITRK6','IQCJ-SCHIP1','FGF5','NTRK3','ADAMTS8','NAALADL2','PTGER4','ANXA3','ERCC2','FAM172A','GNB1L','PSMG2','NMBR','SLC38A4','ROR1','PARD3B','BET1')
Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(merged)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(merged[["Spatial"]]$counts > 0) >= floor(0.1 * ncol(merged[["Spatial"]]$counts))
counts.sub <- merged[["Spatial"]]$counts[genes.to.keep,]
sub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
merged <- AddModuleScore(merged, features = list(sub), name = "Yapsig") #Get score for YAP targets
#Plot Yap signature
p1 <- SpatialFeaturePlot(merged, features = "Yapsig1", images = "slice19", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = "q95", image.alpha = 0) + 
  theme(legend.position="right") + 
  ggtitle("Outer Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p1 <- rotate_image(p1, 118, 1) 
p2 <- SpatialFeaturePlot(merged, features = "Yapsig1", images = "slice9", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = "q95", image.alpha = 0) + 
  theme(legend.position="right") + 
  ggtitle("Inner Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p2 <- rotate_image(p2, -153, 0.83) 
p2 + p1
ggsave("Figures/SpatialYapsig.svg", width = 7, height = 3, units = "in")

#2. ECM gene signature
ecmgenes <- list(c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL24A1", 
                   "LUM", "DCN", "BGN", "HSPG2", "AGRN", 
                   "FN1", "VTN", "TGFBI", "NID1", "NID2", "LAMC1", "LAMA1", "LAMA2", "SPARC", 
                   "FBLN1", "FBLN2", "FBLN5", "LTBP1", "LTBP5", "EMILIN1", "MFAP4", "EFEMP1", "FBN1", "FBN2", "TNC"))
merged <- AddModuleScore(merged, features = ecmgenes, name = "ECMgenes")
p3 <- SpatialFeaturePlot(merged, features = "ECMgenes1", images = "slice19", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = "q95", image.alpha = 0) + 
   theme(legend.position="right") + 
  ggtitle("Outer Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p3 <- rotate_image(p3, 118, 1)
p4 <- SpatialFeaturePlot(merged, features = "ECMgenes1", images = "slice9", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = "q95", image.alpha = 0) + 
  theme(legend.position="right") + 
  ggtitle("Inner Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p4 <- rotate_image(p4, -153, 0.83)
p4 + p3
ggsave("Figures/SpatialECMsig.svg", width = 7, height = 3, units = "in")

#3 Yap component signature 
components <- list(c("YAP1", "CTNNB1", "MST1",
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))
merged <- AddModuleScore(merged, features = components, name = "Components")
range = range(merged$Components1)
q = "q95" #set quantile to cut off for spatial 
p5 <- SpatialFeaturePlot(merged, features = "Components1", images = "slice19", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = q, image.alpha = 0) + 
  theme(legend.position="right") + 
  ggtitle("Outer Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p5 <- rotate_image(p5, 118, 1)
p6 <- SpatialFeaturePlot(merged, features = "Components1", images = "slice9", pt.size.factor = 2.0, stroke = 0.25, max.cutoff = q, image.alpha = 0) + 
   theme(legend.position="right") + 
  ggtitle("Inner Region") + theme(plot.title = element_text(hjust = 0.5, size = 15)) 
p6 <- rotate_image(p6, -153, 0.83)
p6 + p5
ggsave("Figures/Componentsig.svg", width = 7, height = 3, units = "in")

(p4 | p6 | p2) / (p3 | p5 | p1)   & guides(color = guide_legend(override.aes = list(size = 0.5)))
ggsave("Figures/combospatialplot.svg", height = 5, width = 10)

evels(merged) <- c("Stromal 1", "Stromal 2", "Stromal 3", "Stromal 4")   
DotPlot(merged, features = c("Components1", "Yapsig1", "ECMgenes1"), scale.max = 100, scale.min = 0, cols = "RdBu") + theme(axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
          legend.key.size = unit(0.8, 'cm'), #change legend key size
          legend.key.height = unit(0.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'), #change legend key width
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=9)) +  #change legend text font size)) +
          scale_x_discrete(labels=c('YAP Components', "YAP Signature", "ECM Signature"))
ggsave("Figures/dotplotmerged.jpg", width = 10, height = 4, units = "in")





library(Seurat)
library(ggplot2)
library(patchwork)

# Assuming you have a Seurat object named 'merged'
# and you want to plot features on slice images

# Define a function to create and customize the SpatialFeaturePlot
create_custom_plot <- function(merged, feature, image, rotation_angle) {
  plot <- SpatialFeaturePlot(merged, features = feature, images = image, pt.size.factor = 2.0, 
                             stroke = 0, image.alpha = 0, max.cutoff = "q95", ) + 
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
    ggtitle(feature)
  plot <- rotate_image(plot, rotation_angle, 1)
  return(plot)
}

# Create the plots
p1 <- create_custom_plot(merged, "COL1A1", "slice19", 118)
p2 <- create_custom_plot(merged, "COL1A2", "slice19", 118)
p3 <- create_custom_plot(merged, "COL3A1", "slice19", 118)
p4 <- create_custom_plot(merged, "FN1", "slice19", 118)
p5 <- create_custom_plot(merged, "COL1A1", "slice9", -153)
p6 <- create_custom_plot(merged, "COL1A2", "slice9", -153)
p7 <- create_custom_plot(merged, "COL3A1", "slice9", -153)
p8 <- create_custom_plot(merged, "FN1", "slice9", -153)

# Combine the plots
combined_plot <- (p1 + p2 + p5 + p6) / (p3 + p4 + p7 + p8) + plot_layout(guides = "collect") 
ggsave("Figures/spatialECM.svg")

# Define a function to create and customize the SpatialFeaturePlot
create_custom_plot <- function(merged, feature, image, rotation_angle) {
  plot <- SpatialFeaturePlot(merged, features = feature, images = image, pt.size.factor = 2.0, 
                             stroke = 0, image.alpha = 0, max.cutoff = "q95", ) + 
    theme(legend.position = "right",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
    ggtitle(feature)
  plot <- rotate_image(plot, rotation_angle, 1)
  return(plot)
}

# Create the plots
p1 <- create_custom_plot(merged, "ACTA2", "slice19", 118)
p2 <- create_custom_plot(merged, "TAGLN", "slice19", 118)
p3 <- create_custom_plot(merged, "MYH11", "slice19", 118)
p4 <- create_custom_plot(merged, "DES", "slice19", 118)
p5 <- create_custom_plot(merged, "ACTA2", "slice9", -153)
p6 <- create_custom_plot(merged, "TAGLN", "slice9", -153)
p7 <- create_custom_plot(merged, "MYH11", "slice9", -153)
p8 <- create_custom_plot(merged, "DES", "slice9", -153)

# Combine the plots
combined_plot <- (p1 + p2 + p5 + p6) / (p3 + p4 + p7 + p8) + plot_layout(guides = "collect") 
combined_plot
ggsave("Figures/spatialSMC.svg")

# Create the plots
p1 <- create_custom_plot(merged, "CCN1", "slice19", 118)
p2 <- create_custom_plot(merged, "CCN2", "slice19", 118)
p5 <- create_custom_plot(merged, "CCN1", "slice9", -153)
p6 <- create_custom_plot(merged, "CCN2", "slice9", -153)

# Combine the plots
combined_plot <- (p1 + p5) / (p2+p6) + plot_layout(guides = "collect") 
combined_plot
ggsave("Figures/spatialYAP.svg")