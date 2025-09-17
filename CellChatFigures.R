setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/CellChat"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages_updated",.libPaths()))
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales); library(readxl); library(CellChat); library(future); library(ComplexHeatmap)
options(Seurat.object.assay.version = "v5")
available_cores <- future::availableCores()
workers <- max(1, floor(available_cores * 0.75))
plan("multisession", workers = workers)
options(future.globals.maxSize = 8000 * 1024^2)


#### Generate uninvolved and fibrotic references ####
global <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/Combined_Reference_Adipo_Normalized_Scaled.rds") #load adipocyte-containing reference
Idents(global) <- "subpop"

##Relabel the annotations from specific to general
# Define the general_mapping as a named vector in R 
general_mapping <- list(
  "CTHRC1+ mFib" = "CTHRC1+ mFib",
  "LRRC7+ mFib" = "LRRC7+ mFib",
  "MMP3+ iFib" = "MMP3+ iFib",
  "GREM1+ ssFib" = "GREM1+ ssFib",
  "KCNN3+ ssFib" = "KCNN3+ ssFib",
  "CD74+ apFib" = "CD74+ apFib",
  "F3+ ssFib" = "F3+ ssFib",
  "ADAMDEC1+ ssFib" = "ADAMDEC1+ ssFib",
  "THBS1+ ssFib" = "THBS1+ ssFib",
  "PI16+ ssFib" = "PI16+ ssFib",
  
  "DEShi ssSMCs" = "SMCs",
  
  "EBF1+ Steady-State Peri" = "Pericytes",
  "JUNhi Activated Peri" = "Pericytes",
  
  "MYC+ Synthetic vSMCs" = "vSMCs", 
  "CACNA1C+ Contractile vSMCs" = "vSMCs", 
  
  "STEAP4+ Mural Cells" = "Vascular ECs",
  "ACKR1+ Venous ECs" = "Vascular ECs",
  "CD36+ Capillary ECs" = "Vascular ECs", 
  "EFNB2hi Arterial ECs" = "Vascular ECs", 
  
  "STAB1hi Medullary Sinus LECs" = "Lymphatic ECs",
  "CLDN11+ Valve LECs" = "Lymphatic ECs",
  "CCL21hi Paracortical Sinus LECs" = "Lymphatic ECs",
  "NTShi Ceiling LECs" = "Lymphatic ECs",
  "CXCL3hi Floor LECs" = "Lymphatic ECs",
  
  "NRXN1hi ssGlia" = "Glia",
  "SOCS3hi Submucosal Glia" = "Glia",
  "HLA-DRAhi Activated Glia" = "Glia",
  
  "LYVE1+ Macrophages" = "Macrophages",
  "C1QChi Macrophages" = "Macrophages",
  "CXCL2hi Macrophages" = "Macrophages",
  "CD300E+ Monocytes" = "Macrophages",
  "CPA3+ Mast Cells" = "Mast Cells",
  "CSF3R+ Neutrophils" = "Neutrophils",
  "IL4I1 Activated DC" = "DCs",
  "GZMB+ pDC" = "DCs",
  "LAMP3+ mregDC" = "DCs",
  "CLEC9A+ cDC1" = "DCs",
  "CD1C+ cDC2" = "DCs",
  
  "CD8A+ CD8 T Cells" = "CD8 Cells",
  "CCR6+ TH17 T Cells" = "CD4 Cells",
  "CCR7+ TH1 T Cells" = "CD4 Cells",
  "IL7R+ ILC" = "CD4 Cells",
  "FOXP3+ Tregs" = "CD4 Cells",
  
  "AFF3+ Naïve B Cells" = "B Cells",
  "TNFRSF13B+ Memory B Cells" = "B Cells",
  "IGHD+ Follicular B Cells" = "B Cells",
  "AICDA+ GC B Cells" = "B Cells",
  
  "FCER1G+ NK Cells" = "NK Cells",
  
  "IGHG1+ Plasma Cells" = "Plasma Cells",
  "IGHA1+ Plasma Cells" = "Plasma Cells",
  
  "MKI67+ PI Cells" = "PI Cells",
  
  "LGR5+ ISCs" = "Epithelial Cells",
  "MKI67+ TA Cells" = "Epithelial Cells",
  "DEFA5+ Paneth Cells" = "Epithelial Cells",
  "MUC2+ Goblet Cells" = "Epithelial Cells",
  "CHGB+ Enteroendocrine" = "Epithelial Cells",
  "FABP2+ Enterocytes" = "Epithelial Cells",
  "BEST4+ Enterocytes" = "Epithelial Cells",
  "CEACAM7+ Enterocytes" = "Epithelial Cells",
  "SH2D6+ Tuft Cells" = "Epithelial Cells", 
  
  "PLIN1+ Adipocytes" = "Adipocytes")

fibrotic <- subset(global, diseasestate == "Stricture" | diseasestate == "Involved")
uninvolved <- subset(global, diseasestate == "Stricture Uninflamed" | diseasestate == "Uninvolved" | diseasestate == "Normal Bowel" | diseasestate == "Normal MAT")
uninvolved <- subset(uninvolved, tissue == "Colon", invert = T)

# Convert list to named vector (with original subpop labels as names)
general_mapping_vec <- unlist(general_mapping)

# For fibrosub
fibrotic$general <- unname(general_mapping_vec[as.character(fibrotic$subpop)])

# For uninvolvedsub
uninvolved$general <- unname(general_mapping_vec[as.character(uninvolved$subpop)])

Idents(fibrotic) <- "general"
Idents(uninvolved) <- "general"

#Downsample for faster analysis 
fibrosub <- subset(fibrotic, downsample = 2000)
uninvolvedsub <- subset(uninvolved, downsample = 2000)

# Verify
table(fibrosub$general)  # Should show your general categories
table(uninvolvedsub$general)

saveRDS(fibrosub, "fibrosubset.rds", compress = F)
saveRDS(uninvolvedsub, "uninvolvedsubset.rds", compress = F)

####Processing of fibrotic reference####
fibrosub <- readRDS("fibrosubset.rds")
Idents(fibrosub) <- "general"
#create Cellchat object
data.input <- fibrosub[["RNA"]]$data # normalized data matrix
data.input <- as(object = data.input, Class = "dgCMatrix")
labels <- Idents(fibrosub)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
fibrocellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB)
fibrocellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
fibrocellchat <- subsetData(fibrocellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) # do parallel
fibrocellchat <- identifyOverExpressedGenes(fibrocellchat)
fibrocellchat <- identifyOverExpressedInteractions(fibrocellchat)
fibrocellchat <- computeCommunProb(fibrocellchat, type = "triMean")
fibrocellchat <- filterCommunication(fibrocellchat, min.cells = 10)
fibrocellchat <- computeCommunProbPathway(fibrocellchat)
fibrocellchat <- aggregateNet(fibrocellchat)
groupSize <- as.numeric(table(fibrocellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(fibrocellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(fibrocellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(fibrocellchat,color.heatmap = "Reds", measure = "weight")
netVisual_heatmap(fibrocellchat,color.heatmap = "Reds", measure = "count")

fibrocellchat <- netAnalysis_computeCentrality(fibrocellchat, slot.name = "netP")
saveRDS(fibrocellchat, "fibrocellchat_updated.rds", compress = F)

# Access all the signaling pathways showing significant communications
pathways.show.all <- fibrocellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(fibrocellchat@idents)

niche <- c("CD74+ apFib", "Pericytes", "Macrophages", "Plasma Cells", "THBS1+ ssFib", "Pericytes", "GREM1+ ssFib", "PI16+ ssFib")

# Compute the network centrality scores
fibrocellchat <- netAnalysis_computeCentrality(fibrocellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(fibrocellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Compute the network centrality scores
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
ht1 <- netAnalysis_signalingRole_heatmap(fibrocellchat, pattern = "outgoing", font.size = 5, height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(fibrocellchat, pattern = "incoming", font.size = 5, height = 15)

pdf(file="Figures/netAnalysis_signalfibro.pdf", width=10,height=20)
draw(ht2)
dev.off()

####Processing of uninvolved reference####
uninvolvedsub <- readRDS("uninvolvedsubset.rds")
#create Cellchat object
data.input <- uninvolvedsub[["RNA"]]$data # normalized data matrix
data.input <- as(object = data.input, Class = "dgCMatrix")
labels <- Idents(uninvolvedsub)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
uninvolvedcellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB)
uninvolvedcellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
uninvolvedcellchat <- subsetData(uninvolvedcellchat) # This step is necessary even if using the whole database
#future::plan("multisession", workers = 4) # do parallel
uninvolvedcellchat <- identifyOverExpressedGenes(uninvolvedcellchat)
uninvolvedcellchat <- identifyOverExpressedInteractions(uninvolvedcellchat)
uninvolvedcellchat <- computeCommunProb(uninvolvedcellchat, type = "triMean")
uninvolvedcellchat <- filterCommunication(uninvolvedcellchat, min.cells = 10)
uninvolvedcellchat <- computeCommunProbPathway(uninvolvedcellchat)
uninvolvedcellchat <- aggregateNet(uninvolvedcellchat)
groupSize <- as.numeric(table(uninvolvedcellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(uninvolvedcellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(uninvolvedcellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(uninvolvedcellchat,color.heatmap = "Reds", measure = "count")
netVisual_heatmap(uninvolvedcellchat,color.heatmap = "Reds", measure = "weight")
uninvolvedcellchat <- netAnalysis_computeCentrality(uninvolvedcellchat, slot.name = "netP")

saveRDS(uninvolvedcellchat, "uninvolvedcellchat_updated.rds", compress = F)

#### Comparison of datasets ####
fibrocellchat <- readRDS("fibrocellchat_updated.rds")
uninvolvedcellchat <- readRDS("uninvolvedcellchat_updated.rds")
object.list <- list(Uninvolved = uninvolvedcellchat, Fibrotic = fibrocellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list)) #Merge cellchat objects

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

ht <- netVisual_heatmap(cellchat, measure = "count")
pdf(file="Figures/heatmap.cellchatdiffcount.pdf", width = 10, height = 10) #Figure SF7D
draw(ht)
dev.off()
ht <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file="Figures/heatmap.cellchatdiffweight.pdf", width = 10, height = 10) #Figure SF7E
draw(ht)
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CTHRC1+ mFib")
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
library(RColorBrewer)
pdf(file="Figures/heatmap.cellchatdiffheatmapLRuninvolvedoutgoing.pdf", width = 30, height = 30)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 40, height = 30, color.heatmap = "YlOrRd")
draw(ht1)
dev.off()

pdf(file="Figures/heatmap.cellchatdiffheatmapLRstrictureoutgoing.pdf", width = 30, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 40, height = 30, color.heatmap = "YlOrRd")
draw(ht2)
dev.off()

pdf(file="Figures/heatmap.cellchatdiffheatmapLRuninvolvedincoming.pdf", width = 30, height = 30)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 40, height = 30, color.heatmap = "YlOrRd")
draw(ht1)
dev.off()

pdf(file="Figures/heatmap.cellchatdiffheatmapLRstrictureoutincoming.pdf", width = 30, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 40, height = 30, color.heatmap = "YlOrRd")
draw(ht2)
dev.off()

pdf(file="Figures/heatmap.cellchatdiffheatmapLRuninvolved.pdf", width = 30, height = 30)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], color.heatmap = "YlOrRd", height = 40, width = 30)
draw(ht1)
dev.off()

pdf(file="Figures/heatmap.cellchatdiffheatmapLRstricture.pdf", width = 30, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], color.heatmap = "YlOrRd", height = 40, width = 30)
draw(ht2)
dev.off()

targets = c("CTHRC1+ mFib")
netVisual_diffInteraction(cellchat, weight.scale = T, sources.use = niche, targets.use = targets, remove.isolate = T)
netVisual_diffInteraction(cellchat, weight.scale = T, sources.use = niche, targets.use = targets, measure = "weight", remove.isolate = T)

netVisual_bubble(cellchat, sources.use = niche, 
                 targets.use = c("CTHRC1+ mFib"), comparison = c(1,2), angle.x = 45)
ggsave("Figures/netvisualbubblediffLR.pdf", height = 20)

gg1 <- netVisual_bubble(cellchat, sources.use = niche, 
                        targets.use = c("CTHRC1+ mFib"),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Strictures", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = niche, 
                        targets.use = c("CTHRC1+ mFib"),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Strictures", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

ggsave("Figures/bubbleplotdiffLR")
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = niche, 
                       targets.use = c("CTHRC1+ mFib"), legend.pos.x = 10)
}

for (i in 1:length(object.list)) {
  # Define your sources and targets
  sources <- sources
  targets <- targets
  
  # Get valid cell names
  valid_cells <- dimnames(object.list[[i]]@netP$prob)[[1]]
  
  # Ensure sources and targets exist
  sources <- sources[sources %in% valid_cells]
  targets <- targets[targets %in% valid_cells]
  
  if (length(sources) == 0 | length(targets) == 0) {
    warning("Some sources or targets are missing in dataset. Skipping iteration.")
    next
  }
  
  # Extract pathway data for selected sources and targets
  pathway_probs <- object.list[[i]]@netP$prob[sources, targets, , drop = FALSE]
  
  # Sum pathway activity within selected sources and targets
  pathway_sums <- apply(pathway_probs, 3, sum)  # Sum across cell interactions
  
  # Select top 20 pathways
  top_pathways <- names(sort(pathway_sums, decreasing = TRUE)[1:3])
}
# Visualize only these pathways
pdf(file="Figures/chordLRdiagramgenesuninvolved.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[1]], 
                     sources.use = niche, 
                     targets.use = c("CTHRC1+ mFib"), 
                     #signaling = top_pathways,
                     legend.pos.x = 10, title.name = "Uninvolved", lab.cex = 1.5)
dev.off()

pdf(file="Figures/chordLRdiagramgenesstricture.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = sources, 
                     targets.use = targets, 
                     #signaling = top_pathways,
                     legend.pos.x = 10, title.name = "Stricture", lab.cex = 1.5)
dev.off()

# Visualize only these pathways
pdf(file="Figures/chordLRdiagramgenesuninvolved.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[1]], 
                     sources.use = sources, 
                     targets.use = targets, 
                     signaling = top_pathways,
                     legend.pos.x = 10, title.name = "Uninvolved", lab.cex = 1.5)
dev.off()

pdf(file="Figures/chordLRdiagramgenesstricture.pdf", width = 10, height = 10)

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Fibrotic"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Fibrotic",ligand.logFC = 0.5, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)

# Chord diagram for overall upregulated LR pairs (Figure 3C)
pdf(file="Figures/chordLRdiagramgenesstrictureDEG.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = niche, targets.use = c("CTHRC1+ mFib"), slot.name = 'net', net = net.up, lab.cex = 0.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

#Upregulated LR Pairs in macrophages (Figure 3D)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Fibrotic",ligand.logFC = 0.01, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
pdf(file="Figures/chordLRdiagramgenesstrictureDEG_Macs.pdf", width = 10, height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = c("Macrophages"), targets.use = c("CTHRC1+ mFib"), slot.name = 'net', net = net.up, lab.cex = 0.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

#Upregulated LR pairs in Pericytes (Figure 3D)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Fibrotic",ligand.logFC = 0.01, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
pdf(file="Figures/chordLRdiagramgenesstrictureDEG_Pericytes.pdf", width = 10, height = 10)
sources = c("Pericytes")
netVisual_chord_gene(object.list[[2]], sources.use = sources, targets.use = c("CTHRC1+ mFib"), slot.name = 'net', net = net.up, lab.cex = 0.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2], " ", "(",  sources , ")"))
dev.off()

#Upregulated LR pairs in Fibroblasts (Figure 3D)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Fibrotic",ligand.logFC = 0.25, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
pdf(file="Figures/chordLRdiagramgenesstrictureDEG_Fibroblasts.pdf", width = 10, height = 10)
sources = c("CD74+ apFib" , "THBS1+ ssFib", "GREM1+ ssFib", "PI16+ ssFib" )
netVisual_chord_gene(object.list[[2]], sources.use = sources, targets.use = c("CTHRC1+ mFib"), slot.name = 'net', net = net.up, lab.cex = 0.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2], " ", "(Fibroblasts)"))
dev.off()

#Upregulated LR pairs in Plasma Cells (Figure 3D)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Fibrotic",ligand.logFC = 0.01, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
pdf(file="Figures/chordLRdiagramgenesstrictureDEG_PlasmaCells.pdf", width = 10, height = 10)
sources = c("Plasma Cells")
netVisual_chord_gene(object.list[[2]], sources.use = sources, targets.use = c("CTHRC1+ mFib"), slot.name = 'net', net = net.up, lab.cex = 0.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2], " ", "(",  sources , ")"))
dev.off()