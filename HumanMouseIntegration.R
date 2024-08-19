setwd("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy/r_packages_R4.3.3",.libPaths()))
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers)
library(scales); library(SeuratObject)
library(nichenetr)
options(Seurat.object.assay.version = "v5")

human = readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/LongakerOnly/longakerfibro.finalannotations.rds")
mouse = readRDS("mousefibro.annotated.rds")
mouse_genes <- rownames(mouse)
converted_genes <- convert_mouse_to_human_symbols(mouse_genes)

data = mouse[["RNA"]]$counts
clean_genes = !is.na(converted_genes)  #remove NAs
data_toUse = data[ clean_genes, ]
rownames( data_toUse ) = converted_genes[ clean_genes ]
new_rownames <- rownames(data_toUse)
counts.mat <- data_toUse
mouse[[ 'ortholog_mouse' ]] = CreateAssayObject(data_toUse)
DefaultAssay(mouse) <- "ortholog_mouse"

human <- SCTransform(human, verbose = T)
mouse <- SCTransform(mouse, verbose = T, assay = "ortholog_human")


DefaultAssay(mouse) = "SCT"
DefaultAssay(human) = "SCT"
saveRDS(human, "MHIStanford/SCThuman.rds", compress = F)
saveRDS(mouse, "MHIStanford/SCTmouse.rds", compress = F)

human = readRDS("MHIStanford/SCThuman.rds")
mouse = readRDS("MHIStanford/SCTmouse.rds")

anchors <- FindTransferAnchors(reference = human, query = mouse,  k.anchor = 50, dims = 1:30,  normalization.method='SCT')
saveRDS(anchors, "MHIStanford/anchors_k50_SCT.rds")
predictions.assay <- TransferData(anchorset = anchors, refdata = human$longaker.annotation, prediction.assay = T, 
                                  weight.reduction = mouse@reductions$pca, dims=1:30) #PCA
saveRDS(predictions.assay, "MHIStanford/predictions.assay_k50_SCT.rds")

mouse <- AddMetaData(mouse, metadata = t(GetAssayData(predictions.assay)))
fp <- FeaturePlot(mouse, features = c("Mechanosensitive CTHRC1+")) + ggtitle("CTHRC1+ Mapping") + guides(color = guide_legend(title = "Probability"), order = T) + 
  theme(legend.title =element_text(size=10)) + 
  scale_colour_gradientn(name = "Score", colors = rev(c("#d73027", "#fc8d59", "#fee090", "#ffffbf","#e0f3f8","#91bfdb","#4575b4"))) +  guides(colour = guide_legend(override.aes = list(size=3)))
dp <- DotPlot(mouse, features = c("Mechanosensitive CTHRC1+"), cols = "RdYlBu", col.min = 0, col.max = 100, dot.scale = 2, scale.max = 100) +  theme(axis.title.x = element_blank(),
                                                                                                                                                     axis.title.y = element_blank(), 
                                                                                                                                                     legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                                                                                     legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                                                                                     legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                                                                                     legend.title = element_text(size=6), #change legend title font size                                                                                                                                           legend.text = element_text(size=6)) 
fp | dp
saveRDS(mouse, "mouseintegrated.k50.rds")
VlnPlot(human, features = c("msFibroblast Postn+"), pt.size = 0)

