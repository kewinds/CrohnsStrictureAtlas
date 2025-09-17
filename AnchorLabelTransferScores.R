setwd("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy/r_packages",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers); 
library(scales); library(SeuratObject); library(reshape2); library(nichenetr)
options(Seurat.object.assay.version = "v5")

####Similarity by ALT Mouse vs Human (Figure 5G) ####
human = readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Fibro_Annotated_0.28.rds")
mouse = readRDS("/oak/stanford/groups/longaker/KEBR/Mouse_Colotomy/mousefibro.CCA_0.7.annotated.rds")

human <- FindVariableFeatures(human)
mouse <- JoinLayers(mouse, assay = "RNA")
DefaultAssay(mouse) <- "RNA" #Set the default assay to be "RNA" to change the rownames in entire count matrix
mouse_genes <- rownames(mouse)
converted_genes <- convert_mouse_to_human_symbols(mouse_genes) #convert genes with nichenetr
data = mouse[["RNA"]]$counts 
clean_genes = !is.na(converted_genes)  #remove NAs
data_toUse = data[ clean_genes, ]
rownames( data_toUse ) = converted_genes[ clean_genes ]
new_rownames <- rownames(data_toUse)
counts.mat <- data_toUse
mouse[[ 'ortholog_human' ]] = CreateAssayObject(data_toUse)
DefaultAssay(mouse) <- "ortholog_human"

anchors <- FindTransferAnchors(reference = human, query = mouse, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = human$fibro_annotation, dims = 1:30)
mousequery <- AddMetaData(mouse, metadata = predictions)
table(mousequery$predicted.id)
mousequery <- RunUMAP(mousequery, dims = 1:30, reduction = "cca.integration", return.model = TRUE)
mousequery <- MapQuery(anchorset = anchors, reference = human, query = mousequery,
                       refdata = list(celltype = "fibro_annotation"), reference.reduction = "pca", reduction.model = "umap.harmony")
p1 <- DimPlot(human, reduction = "umap.harmony", group.by = "fibro_annotation", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(mousequery, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

FeaturePlot(mousequery, features = c("prediction.score.MMP3..iFib"))

levels(mousequery) <- rev(c("Cthrc1+ mFib", "Eln+ mFib", "Timp1+ iFib",
                            "Cd74+ apFib", "Adamdec1+ ssFib", "Mgp+ ssFib", 
                            "Grem1+ ssFib", "Abca8a+ ssFib", 
                            "Col15a1+ ssFib", "Pi16+ ssFib"))

DotPlot(mousequery, features = c("prediction.score.CTHRC1..mFib", "prediction.score.GREM1..ssFib", "prediction.score.MMP3..iFib", 
                                 "prediction.score.F3..ssFib", "prediction.score.CD74..apFib", "prediction.score.ADAMDEC1..ssFib", 
                                 "prediction.score.KCNN3..ssFib", "prediction.score.THBS1..ssFib",  "prediction.score.LRRC7..mFib",
                                 "prediction.score.PI16..ssFib"), 
        cols = "RdYlBu", scale = T) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels=c("CTHRC1+ mFib", "GREM1+ ssFib", 
                            "MMP3+ iFib", "F3+ ssFib","CD74+ apFib", "ADAMDEC1+ ssFib", 
                            "KCNN3+ ssFib", "THBS1+ ssFib", "LRRC7+ mFib", "PI16+ ssFib"))

####Similarity by ALT MAT vs Bowel (Figure S3I)####
bowel = readRDS("Fibro_Bowel_Annotated_0.35.rds")
mat = readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/matfibroblast.annotated.rds")
bowel <- FindVariableFeatures(bowel)
mat <- FindVariableFeatures(mat)

anchors <- FindTransferAnchors(reference = bowel, query = mat, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = bowel$bowel.annotation, dims = 1:30)
matquery <- AddMetaData(mat, metadata = predictions)
table(matquery$predicted.id)
matquery <- RunUMAP(matquery, dims = 1:30, reduction = "harmony", return.model = TRUE)
matquery <- MapQuery(anchorset = anchors, reference = bowel, query = matquery,
                     refdata = list(celltype = "bowel.annotation"), reference.reduction = "pca", reduction.model = "umap.harmony")
p1 <- DimPlot(bowel, reduction = "umap.harmony", group.by = "bowel.annotation", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(matquery, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

FeaturePlot(matquery, features = c("prediction.score.MMP3..iFib"))

matquery$mat.annotation <- factor(mat$mat.annotation, levels = rev(c("CTHRC1+ mFAP", "LRRC7+ mFAP", "ICAM1+ iAPC", "PPARG+ ssAPC",
                                                                     "VIT+ ssAPC", "FMO2+ ssAPC", "DPP4+ ssAPC", "NOTCH3+ Pericytes")))

Idents(matquery) <- "mat.annotation"
DotPlot(matquery, features = c("prediction.score.CTHRC1..mFib", "prediction.score.LRRC7..mFib", "prediction.score.MMP3..iFib", "prediction.score.CD74..apFib", 
                               "prediction.score.ADAMDEC1..lpFib", "prediction.score.F3..Telocytes", "prediction.score.KCNN3..ssFRC", 
                               "prediction.score.GREM1..ssFRC", "prediction.score.SOD2..ssFRC" , "prediction.score.PI16..ssFRC"), 
        cols = "RdYlBu", scale = F) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels=c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib", 
                            "CD74+ apFib", "ADAMDEC1+ lpFib","F3+ Telocytes", "KCNN3+ ssFRC", 
                            "GREM1+ ssFRC", "THBS1+ ssFRC", "PI16+ ssFib"))

ggsave("Figures/ALTbowelvsMAT.pdf") #Figure S2

####Similarity by ALT In-House vs Global (Figure S3L)####
global = readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Fibro_Annotated_0.28.rds")
stanford = readRDS("/oak/stanford/groups/longaker/KEBR/stanford_Colotomy/stanfordfibro.CCA_0.7.annotated.rds")
global <- subset(obj, subset = !(study %in% c("LongakerBowel", "LongakerMes"))) #Remove Stanford study
anchors <- FindTransferAnchors(reference = global, query = stanford, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = global$fibro_annotation, dims = 1:30)
stanfordquery <- AddMetaData(stanford, metadata = predictio ( ns)
table(stanfordquery$predicted.id)
stanfordquery <- RunUMAP(stanfordquery, dims = 1:30, reduction = "cca.integration", return.model = TRUE)
stanfordquery <- MapQuery(anchorset = anchors, reference = global, query = stanfordquery,
                       refdata = list(celltype = "fibro_annotation"), reference.reduction = "pca", reduction.model = "umap.harmony")

DotPlot(stanfordquery, features = c("prediction.score.CTHRC1..mFib", "prediction.score.LRRC7..ssFib", "prediction.score.MMP3..iFib", 
                                 "prediction.score.F3..ssFib", "prediction.score.CD74..apFib", "prediction.score.ADAMDEC1..ssFib", 
                                 "prediction.score.KCNN3..ssFib", "prediction.score.THBS1..ssFib",  "prediction.score.GREM1..mFib",
                                 "prediction.score.PI16..ssFib"), 
        cols = "RdYlBu", scale = T) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels=c("CTHRC1+ mFib", "LRRC7+ ssFib", 
                            "MMP3+ iFib","CD74+ apFib", "ADAMDEC1+ ssFib", "F3+ ssFib",
                            "KCNN3+ ssFib", "GREM1+ mFib", "THBS1+ ssFib", "PI16+ ssFib"))

ggsave("Figures/ALTbowelvsMAT.pdf") #Figure S2
