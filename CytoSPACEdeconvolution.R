####Generate Cytospace plots 
####Libraries and functions ####
rm(list = ls())

library(spacexr)
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(BPCells)

.libPaths(c("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/r_packages", .libPaths()))
setwd("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace")

#### Annotate fibrotic spatial objects ####
#load in objects
fibro <- readRDS("finalfibroblastannotations.rds")
global <- readRDS("humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds")
stricture_spatial <- readRDS("Longakerfibroticrefspat.rds")
sc_ref <- readRDS("longakerfineannotated.rds")
#Correct metadata
stricture_spatial$orig.ident[stricture_spatial$orig.ident == "slice_19"] <- "slice19"
stricture_spatial$orig.ident[stricture_spatial$orig.ident == "slice_9"] <- "slice9"
saveRDS(stricture_spatial, "Longakerfibroticrefspat.rds")

#Subset global and fibro objects to only look at our study: 
global_sub <- subset(global, subset = (study == "LongakerBowel" | study == "LongakerMes"))
fibro_sub <- subset(fibro, subset =  (study == "LongakerBowel" | study == "LongakerMes"))

#Add fibroblast subpopulations into Longaker study 
DefaultAssay(global_sub) <- "RNA"
Idents(global_sub) <- "manual.annotation"
DefaultAssay(fibro_sub) <- "RNA"
Idents(fibro_sub) <- "fibro.annotation"
global_sub$subpop.annotation <- "NA"
global_sub$subpop.annotation <- as.character(Idents(global_sub))
global_sub$subpop.annotation <- global_sub$manual.annotation #Add in global populations
global_sub$subpop.annotation[colnames(global_sub) %in% colnames(fibro_sub)] <- "Fibroblasts2"
global_sub$subpop.annotation[Cells(fibro_sub)] <- paste(Idents(fibro_sub))

Idents(global_sub) <- "subpop.annotation"
global_sub$subpop.annotation[global_sub$subpop.annotation == "Fibroblasts2"] <- fibro_sub$fibro.annotation
Idents(global_sub) <- "subpop.annotation"
DimPlot(global_sub)
saveRDS(global_sub, "longakerfineannotated.rds", compress = F)

#### Generate output files for cytospace ####

sc_ref <- readRDS("longakerfineannotated.rds")
Idents(sc_ref) <- "subpop.annotation"
sc_ref <- subset(sc_ref, subset = diseasestate %in% c("Ileostomy Bowel", "Ileostomy MAT"), invert = T) #exclude ileostomy samples 
sc_ref <- subset(sc_ref, idents = "Fibroblasts", invert = T) #Remove ambiguous cells thrown out while subclustering fibroblasts 
table(sc_ref$diseasestate, sc_ref$subpop.annotation)

#Longaker dataset
stricture_spatial[["Spatial"]] <- as(stricture_spatial[["Spatial"]], Class = "Assay")
for (sampleName in unique(stricture_spatial$orig.ident)) {
  stricture_spatial.subset <- subset(stricture_spatial, subset = orig.ident == sampleName)
  imgName <- Images(stricture_spatial)[grepl(sampleName, Images(stricture_spatial))]
  generate_cytospace_from_ST_seurat_object(stricture_spatial.subset, dir_out=paste0('cytospaceInput/',imgName), fout_prefix='', write_sparse=T, slice=imgName)
}

## estimate cell type fractions using RCTD
# generate RCTD ref
counts <- sc_ref[["RNA"]]$counts
cell_types <- sc_ref$manual.annotation
cell_types <- as.factor(cell_types)
nUMI <- sc_ref$nCount_RNA
reference <- Reference(counts, cell_types, nUMI, require_int = F)
table(reference@cell_types) #Mast cells <25 so we combine with Myeloid Cells 
#Create new metadata category 
sc_ref$RCTD.annotation <- "NA"
sc_ref$RCTD.annotation <- sc_ref$manual.annotation
sc_ref$RCTD.annotation[sc_ref$manual.annotation == "Myeloid Cells"] = "Myeloid + Mast Cells"
sc_ref$RCTD.annotation[sc_ref$manual.annotation == "Mast Cells"] = "Myeloid + Mast Cells"
cell_types <- sc_ref$RCTD.annotation
cell_types <- as.factor(cell_types)
nUMI <- sc_ref$nCount_RNA
reference <- Reference(counts, cell_types, nUMI, require_int = F)
table(reference@cell_types) #Mast cells <25 so we combine with Myeloid Cells 
saveRDS(reference, 'cytospaceInput/rctd/SCRefCoarse.rds', compress = F)

#Update RCTD to include fine annotations
sc_ref$fineRCTD <- "NA"
sc_ref$fineRCTD <- sc_ref$subpop.annotation
sc_ref$fineRCTD[sc_ref$subpop.annotation == "Myeloid Cells"] = "Myeloid + Mast Cells"
sc_ref$fineRCTD[sc_ref$subpop.annotation == "Mast Cells"] = "Myeloid + Mast Cells"
table(sc_ref$fineRCTD, sc_ref$diseasestate)
Idents(sc_ref) <- "fineRCTD"
saveRDS(sc_ref, "scREF_RCTDfine_adjusted.rds", compress = F)

#Now export the 
#VariableFeatures(sc_ref) <- readRDS("humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds")
sc_ref[["RNA"]] <- as(sc_ref[["RNA"]], Class = "Assay")
Idents(sc_ref) <- "fineRCTD"
generate_cytospace_from_scRNA_seurat_object(sc_ref, dir_out='cytospaceInput', fout_prefix='', write_sparse=T, rna_assay='RNA')

#Generate cell fractions in stricture spatial objects
rm(list = ls())
reference = readRDS('cytospaceInput/rctd/SCRefCoarse.rds')
query <- readRDS("Longakerfibroticrefspat.rds")
query.list <- unique(query$orig.ident)
for (i in 1:length(query.list)) {
  #for (i in 19:27) {
  sampleName <- query.list[i]
  query.subset <- subset(query, subset = orig.ident == sampleName)
  imgName <- Images(query)[grepl(sampleName, Images(query))]
  coords <- GetTissueCoordinates(query.subset, image = imgName)
  colnames(coords) <- c("x","y")
  counts <- query.subset[["Spatial"]]$counts
  nUMI <- query.subset$nCount_Spatial
  puck <- SpatialRNA(coords, counts, nUMI)
  saveRDS(puck, paste0('cytospaceInput/rctd/',imgName,'.rds'), compress = F)
  myRCTD <- create.RCTD(puck, reference, max_cores = 1, test_mode = F) # here puck is the SpatialRNA object, and reference is the Reference object.
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD, paste0('cytospaceInput/rctd/',imgName,'.RCTDCoarse.rds'), compress = F)
  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights <- normalize_weights(results$weights) 
  col_sums <- colSums(norm_weights)
  cellfrac <- col_sums / sum(col_sums)
  cellfrac <- data.frame(Index = names(cellfrac), Fraction = cellfrac)
  write.table(t(cellfrac), paste0('cytospaceInput/rctd/',imgName,'.RCTDCoarse.cellfrac.txt'), quote = FALSE, sep = "\t", col.names = FALSE)
}

# impute fibroblast subtype proportion for strictures
rm(list = ls())
query <- readRDS("Longakerfibroticrefspat.rds")
refMeta <- readRDS("finalfibroblastannotations.rds") #Determine fibroblast proportions using Longaker dataset subset using labels from meta-analysis as reference
refMeta <- subset(refMeta, subset = (study %in% c("LongakerBowel", "LongakerMes"))) #Extract Longaker study only
refMeta <- subset(refMeta, subset = (diseasestate %in% c("Ileostomy Bowel", "Ileostomy MAT")), invert = T) #Remove Ileostomy samples
refMeta <- subset(refMeta, subset = (diseasestate %in% c("Stricture", "Involved"))) #Keep stricture samples
table(refMeta$diseasestate, refMeta$study)

query.list <- unique(query$orig.ident)
for (i in 1:length(query.list)) {
  sampleName <- query.list[i]
  imgName <- Images(query)[grepl(sampleName, Images(query))]
  cellfrac <- read.table(paste0("cytospaceInput/rctd/",imgName,".RCTDCoarse.cellfrac.txt"), sep = "\t", header = T, check.names = F)
  fib <- cellfrac$Fibroblast
  cellfrac <- data.frame(Index = names(cellfrac), Fraction = t(cellfrac))
  cellfrac <- cellfrac[-1,] #remove dummy first row
  cellfrac <- cellfrac[-which(rownames(cellfrac) == "Fibroblasts"),] #Remove fibroblast row to be replaced with fibroblast cluster specific
  fibDis <- prop.table(table(refMeta$fibro.annotation))
  fibDis <- data.frame(fibDis * fib)
  colnames(fibDis) <- c("Index", "Fraction")
  rownames(fibDis) <- fibDis$Index #Set rownames using index 
  cellfrac <- rbind(cellfrac, fibDis)
  write.table(t(cellfrac), paste0('cytospaceInput/rctd/',imgName,'.RCTDCoarse.cellfrac.fib.txt'), quote = FALSE, sep = "\t", col.names = FALSE)
}

rm(list = ls())

refMeta <- readRDS("finalfibroblastannotations.rds") #Determine fibroblast proportions using Longaker dataset subset using labels from meta-analysis as reference
refMeta <- subset(refMeta, subset = (study %in% c("LongakerBowel", "LongakerMes"))) #Extract Longaker study only
refMeta <- subset(refMeta, subset = (diseasestate %in% c("Ileostomy Bowel", "Ileostomy MAT")), invert = T) #Remove Ileostomy samples
refMeta <- subset(refMeta, subset = (diseasestate %in% c("Stricture Uninflamed", "Uninvolved"))) #Keep uninvolved samples
table(refMeta$diseasestate, refMeta$study)

query.list <- unique(query$orig.ident)
for (i in 1:length(query.list)) {
  sampleName <- query.list[i]
  imgName <- Images(query)[grepl(sampleName, Images(query))]
  cellfrac <- read.table(paste0("cytospaceInput/rctd/",imgName,".RCTDCoarse.cellfrac.txt"), sep = "\t", header = T, check.names = F)
  fib <- cellfrac$Fibroblast
  cellfrac <- data.frame(Index = names(cellfrac), Fraction = t(cellfrac))
  cellfrac <- cellfrac[-1,] #remove dummy first row
  cellfrac <- cellfrac[-which(rownames(cellfrac) == "Fibroblasts"),] #Remove fibroblast row to be replaced with fibroblast cluster specific
  fibDis <- prop.table(table(refMeta$fibro.annotation))
  fibDis <- data.frame(fibDis * fib)
  colnames(fibDis) <- c("Index", "Fraction")
  rownames(fibDis) <- fibDis$Index #Set rownames using index 
  cellfrac <- rbind(cellfrac, fibDis)
  write.table(t(cellfrac), paste0('cytospaceInput/rctd/',imgName,'.RCTDCoarse.cellfrac.fib.txt'), quote = FALSE, sep = "\t", col.names = FALSE)
}

#### Run CytoSPACE ####
#Done in terminal in new tab
#start here 
#module load anaconda; conda env list | grep cuda #activate anaconda environment
cd /oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace/cytospace_v1.0.6/
  conda activate cytospace #activate cytospace environment
#cd /oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace/cytospace_v1.0.6/
conda env create -f environment.yml
cd /oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace/
  
  for SAMPLE in $(ls -d ../cytospaceInput/slice*)
do
SAMPLE=$(echo ${SAMPLE:18})
cytospace \
-sp ../cytospaceInput/scRNA_data.mtx \
-ctp ../cytospaceInput/cell_type_labels.txt \
-stp ../cytospaceInput/$SAMPLE/ST_data.mtx \
-cp ../cytospaceInput/$SAMPLE/Coordinates.txt \
-sm lap_CSPR \
-o ../cytospaceOutput/$SAMPLE \
-ctfep ../cytospaceInput/rctd/$SAMPLE.RCTDCoarse.cellfrac.fib.txt \
-sss -nosss 3000 -nop 16
done

## add cytospace metadata
rm(list=ls())
spat <- readRDS("Longakerfibroticrefspat.rds")
file_list <- list.files(paste0("cytospaceOutput/"),"cell_type_assignments_by_spot.csv", recursive = TRUE)
file_list <- file_list[c(1,4)] #fibrotic slices
for (i in 1:length(file_list)) {
  file <- file_list[i]
  df <- read.csv(paste0("cytospaceOutput/",file), check.names = F)
  rownames(df) <- df[,1]
  df <- df[,-1]
  spat <- AddMetaData(spat, metadata = df)
}

file_list <- list.files(paste0("cytospaceOutput/"),"fractional_abundances_by_spot.csv", recursive = TRUE)
file_list <- file_list[c(1,4)] #fibrotic slices
for (i in 1:length(file_list)) {
  file <- file_list[i]
  df <- read.csv(paste0("cytospaceOutput/",file), check.names = F)
  rownames(df) <- df[,1]
  df <- df[,-1]
  colnames(df) <- paste0("pct_",colnames(df))
  spat <- AddMetaData(spat, metadata = df)
}
spat$keep <- T
spat$keep[is.na(spat$"Total cells")] <- F
spat <- subset(spat, subset = keep)
saveRDS(spat, "LongakerfibroticrefspatCYTOSPACE.rds", compress = F)

#Import function to rotate plots 
rotate_image <- function(p, rot_angle, scale) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle, width = scale, height = scale)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}

#### Plot spatial features numbers per spot and density ####
spat <- readRDS("LongakerfibroticrefspatCYTOSPACE.rds")
q = "q95"
#fill = c("#0d0887", "#cc4778", "#f0f921")
fill = c("grey94","lightyellow", "orange", "red", "darkred")

p1 <- SpatialFeaturePlot(spat, features = "mCDF CTHRC1+", images = "slice19", pt.size.factor = 2.2, stroke = 0, max.cutoff = q, image.alpha = 0) + 
  scale_fill_gradientn(name = "Cell Number", colors = fill) + theme(legend.position="right", legend.key.size = unit(1, 'lines')) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
p1 <- rotate_image(p1, 118, 1)
p2 <- SpatialFeaturePlot(spat, features = "mCDF CTHRC1+", images = "slice9", pt.size.factor = 2.2, stroke = 0, max.cutoff = q, image.alpha = 0) + 
  scale_fill_gradientn(name = "Cell Number", colors = fill) + theme(legend.position="right", legend.key.size = unit(1, 'lines')) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
p2 <- rotate_image(p2, -153, 0.83)
spat$"pct_mCDF CTHRC1+" <- 100*(spat$"pct_mCDF CTHRC1+")
p3 <- SpatialFeaturePlot(spat, features = "pct_mCDF CTHRC1+", images = "slice19", pt.size.factor = 2.2, stroke = 0, max.cutoff = q, image.alpha = 0) + 
  scale_fill_gradientn(name = "Percent", colors = fill) + theme(legend.position="right", legend.key.size = unit(1, 'lines')) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
p3 <- rotate_image(p3, 118, 1)
p4 <- SpatialFeaturePlot(spat, features = "pct_mCDF CTHRC1+", images = "slice9", pt.size.factor = 2.2, stroke = 0, max.cutoff = q, image.alpha = 0) + 
  scale_fill_gradientn(name = "Percent", colors = fill) + theme(legend.position="right", legend.key.size = unit(1, 'lines')) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
p4 <- rotate_image(p4, -153, 0.83)
(p1+p3)/(p2+p4) 
ggsave("Figures/Cytospacepredictionsplot.svg")

####Plot cell locations ####
#Add in assigned locations
rm(list = ls())
file_list <- list.files(paste0("cytospaceOutput/"),"assigned_locations.csv", recursive = TRUE)
file_list <- file_list[c(1,4)] #fibrotic slices
cell_locs <- data.frame()
for (i in 1:length(file_list)) {
  file <- file_list[i]
  df <- read.csv(paste0("cytospaceOutput/",file), check.names = F)
  rownames(df) <- df[,1]
  df <- df[,-1]
  df$Image <- str_extract(file_list[i], "slice\\d+") #Add in image name 
  cell_locs <- rbind(cell_locs, df)
}

cell_locs$UniqueCID <- rownames(cell_locs)

#Add in domains to cell_locs
spat <- readRDS("LongakerfibroticrefspatCYTOSPACE.rds")
SpotID <- colnames(spat)
Domain <- spat$Domain
df <- data.frame(SpotID, Domain)
cell_locs <- full_join(cell_locs, df, by = join_by(SpotID))

saveRDS(cell_locs, "cell_locs_object.rds")

#### Plot all fibroblasts ####
##Slice 19 Outer bowel wall
slice <- cell_locs[cell_locs$Image == "slice19", ]
col = slice$row*0.15885623
row = slice$col*0.15885623
slice$highlight <- slice$CellType
# Sample data
list <- c("SMCs", "B cells", "T cells", "Glia", "Macrophages", "Vascular ECs", "Plasma Cells", "Lymphatic ECs", "Myeloid + Mast Cells", "Epithelial")
for (i in list) {
  slice$highlight[slice$highlight == i] <- "Other"
}
cell.type <- slice$highlight
df <- data.frame(cell.type, row, col)
legend_order <- rev(c("Other", "ssCDF PI16+", "ssCDF FMO2+", "ssCDF GREM1+", 
                  "ssCDF KCNN3+", "ssCDF F3+", "iCDF ADAMDEC1+", "iCDF MMP3+", 
                  "apCDF CCL19+", "mCDF LRRC7+", "mCDF CTHRC1+"))
df$cell.type <- factor(df$cell.type, levels = legend_order)
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394", "#073763", "grey92")
p1 <- ggplot(df) + geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
             data = df[df$cell.type == "Other", ], 
             size = 0.8) +
  geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
             data = df[df$cell.type != "Other", ], 
             size = 0.8) +
  scale_color_manual(values = colors, breaks = legend_order) + 
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  theme(legend.position = "right", legend.key.size = unit(1, 'lines'))
p1 <- rotate_image(p1, 20, 1)
rm(df, highlightdf)

#Slice 9 Inner BW
slice <- cell_locs[cell_locs$Image == "slice9", ]
col = slice$row*0.15885623
row = slice$col*0.15885623
slice$highlight <- slice$CellType
# Sample data
list <- c("SMCs", "B cells", "T cells", "Glia", "Macrophages", "Vascular ECs", "Plasma Cells", "Lymphatic ECs", "Myeloid + Mast Cells", "Epithelial")
for (i in list) {
  slice$highlight[slice$highlight == i] <- "Other"
}
cell.type <- slice$highlight
df <- data.frame(cell.type, row, col)
legend_order <- rev(c("Other", "ssCDF PI16+", "ssCDF FMO2+", "ssCDF GREM1+", 
                      "ssCDF KCNN3+", "ssCDF F3+", "iCDF ADAMDEC1+", "iCDF MMP3+", 
                      "apCDF CCL19+", "mCDF LRRC7+", "mCDF CTHRC1+"))
df$cell.type <- factor(df$cell.type, levels = legend_order)
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394", "#073763", "grey92")
p2 <- ggplot(df) + geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                data = df[df$cell.type == "Other", ], 
                size = 0.8) +
     geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                data = df[df$cell.type != "Other", ], 
                size = 0.8) +
     scale_color_manual(values = colors, breaks = legend_order) + 
     theme_bw() +
     theme(axis.line = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank()) + 
     theme(axis.title = element_blank(),
           axis.text = element_blank(),
           axis.ticks = element_blank()) + 
     theme(legend.position = "right", legend.key.size = unit(1, 'lines'))
p2 <- rotate_image(p2, 120, 0.8)
rm(df, highlightdf)
p1+p2 & plot_layout(guides = "collect")



#### Plot fibrotic niches ####
##Slice 19 Outer bowel wall
slice <- cell_locs[cell_locs$Image == "slice19", ]
col = slice$row*0.15885623
row = slice$col*0.15885623
slice$highlight <- slice$CellType
# Sample data
list <- c("B cells", "T cells", "Glia", "Lymphatic ECs", "Myeloid + Mast Cells", "Epithelial", 
          "ssCDF PI16+", "ssCDF FMO2+", "ssCDF GREM1+", 
          "ssCDF KCNN3+", "ssCDF F3+", "iCDF ADAMDEC1+", "iCDF MMP3+", 
          "apCDF CCL19+", "mCDF LRRC7+")
for (i in list) {
  slice$highlight[slice$highlight == i] <- "Other"
}
cell.type <- slice$highlight
df <- data.frame(cell.type, row, col)
legend_order <- c("mCDF CTHRC1+", "SMCs", "Macrophages", "Vascular ECs", "Plasma Cells", "Other")
df$cell.type <- factor(df$cell.type, levels = legend_order)
colors <- c("#660000", "#ff80be", "#9aceeb", "#f6b26b", "#76a5af", "grey92")
                     p1 <- ggplot(df) + geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                                                   data = df[df$cell.type == "Other", ], 
                                                   size = 0.8) +
                       geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                                  data = df[df$cell.type != "Other", ], 
                                  size = 0.8) +
                       scale_color_manual(values = colors, breaks = legend_order) + 
                       theme_bw() +
                       theme(axis.line = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank()) + 
                       theme(axis.title = element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank()) + 
                       theme(legend.position = "right", legend.key.size = unit(1, 'lines'))
                     p1 <- rotate_image(p1, 20, 1)
                     rm(df, highlightdf)
                     
 #Slice 9 Inner BW
 slice <- cell_locs[cell_locs$Image == "slice9", ]
 col = slice$row*0.15885623
 row = slice$col*0.15885623
 slice$highlight <- slice$CellType
 # Sample data
 list <- c("B cells", "T cells", "Glia", "Lymphatic ECs", "Myeloid + Mast Cells", "Epithelial", 
           "ssCDF PI16+", "ssCDF FMO2+", "ssCDF GREM1+", 
           "ssCDF KCNN3+", "ssCDF F3+", "iCDF ADAMDEC1+", "iCDF MMP3+", 
           "apCDF CCL19+", "mCDF LRRC7+")
 for (i in list) {
   slice$highlight[slice$highlight == i] <- "Other"
 }
 cell.type <- slice$highlight
 df <- data.frame(cell.type, row, col)
 legend_order <- c("mCDF CTHRC1+", "SMCs", "Macrophages", "Vascular ECs", "Plasma Cells", "Other")
 df$cell.type <- factor(df$cell.type, levels = legend_order)
                      p2 <- ggplot(df) + geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                                                    data = df[df$cell.type == "Other", ], 
                                                    size = 0.8) +
                        geom_point(aes(x = jitter(col, factor = 3), y = jitter(row, factor = 3), color = cell.type), 
                                   data = df[df$cell.type != "Other", ], 
                                   size = 0.8) +
                        scale_color_manual(values = colors, breaks = legend_order) + 
                        theme_bw() +
                        theme(axis.line = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank()) + 
                        theme(axis.title = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank()) + 
                        theme(legend.position = "right", legend.key.size = unit(1, 'lines'))
                      p2 <- rotate_image(p2, 120, 0.8)
                      rm(df, highlightdf)
                      p1+p2 & plot_layout(guides = "collect")
ggsave("Figures/Fibroticniches.svg")              
#### Highlight only one fibroblast type ####
#Highlight only fibroblast cell type 
slice <- cell_locs[cell_locs$Image == "slice19", ]
slice$highlight <- slice$CellType
slice$highlight[slice$highlight != "mCDF CTHRC1+"] <- "Other"
row = slice$row*0.15885623
col = slice$col*0.15885623
cell.type = slice$highlight
df = data.frame(cell.type,row,col)
p5 <- ggplot(df, aes(x = jitter(row, factor = 3), y = jitter(col, factor = 3))) + #Add jitter so points don't overlap
  geom_point(data = df[cell.type != "mCDF CTHRC1+", ], size=0.8, aes(x=jitter(row, factor = 3),y=jitter(col, factor = 3), color=cell.type))+ 
  geom_point(data = df[cell.type == "mCDF CTHRC1+", ], size=0.8, aes(x=jitter(row, factor = 3),y=jitter(col, factor = 3), color=cell.type))+ 
  labs(col = "Cell Type") + 
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) + 
  scale_color_manual(values=c("red", "grey90")) + theme(legend.position="right")  + 
  theme(legend.text = element_text(size = "12"), legend.title = element_text(size = "14")) + 
  guides(color = guide_legend(override.aes = list(size = 4)))
p5 <- rotate_image(p5, 23,1.1)
rm(df)

#Highlight only fibroblast cell type 
slice <- cell_locs[cell_locs$Image == "slice9", ]
slice$highlight <- slice$CellType
slice$highlight[slice$highlight != "mCDF CTHRC1+"] <- "Other"
row = slice$row*0.15885623
col = slice$col*0.15885623
cell.type = slice$highlight
df = data.frame(cell.type,row,col)
p6 <- ggplot(df, aes(x = jitter(row, factor = 3), y = jitter(col, factor = 3))) + #Add jitter so points don't overlap
  geom_point(data = df[cell.type != "mCDF CTHRC1+", ], size=0.8, aes(x=jitter(row, factor = 3),y=jitter(col, factor = 3), color=cell.type))+ 
  geom_point(data = df[cell.type == "mCDF CTHRC1+", ], size=0.8, aes(x=jitter(row, factor = 3),y=jitter(col, factor = 3), color=cell.type))+ 
  labs(col = "Cell Type") + 
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) + 
  scale_color_manual(values=c("red", "grey90")) + theme(legend.position="right")  + 
  theme(legend.text = element_text(size = "12"), legend.title = element_text(size = "14")) + 
  guides(color = guide_legend(override.aes = list(size = 4)))
p6 <- rotate_image(p6, 120,0.83)
p6
rm(df)

(p6 + p2 + p4)  / (p5 + p1 + p3) 
ggsave("Figures/CytospaceCTHRC1mappings.jpg", width = 16, height = 7)
ggsave("Figures/CytospaceCTHRC1mappings.svg")
#Plot returned object from CytoSpace
plotnames <- unique(fibro_cell_locs$CellType)
tableau20 <- tableau20(n=20)
my.colors <- tableau20[1:length(plotnames)]
names(my.colors) <- plotnames

#Generate fibroblast proptable

#Generate proptable and barplots of fibroblasts in resection 
tab = table(spat$"mCDF CTHRC1+", spat$Domain)
bowel <- subset(spat, Domain == "BW")
mat <- subset(spat, Domain == "CF")
inter <- subset(spat, Domain == "IF")
#calculate total numbers of CTHRC1+ fibroblasts in each domain 
#calculate total numbers of fibroblasts in each domain using the calculated #predicted cells per spot in the annotated spat object 
btotalcount <- sum(bowel$"ssCDF PI16+", bowel$"ssCDF FMO2+", bowel$"ssCDF F3+", bowel$"ssCDF GREM1+", bowel$"ssCDF KCNN3+", 
                   bowel$"iCDF ADAMDEC1+", bowel$"iCDF MMP3+", bowel$"apCDF CCL19+", 
                   bowel$"mCDF LRRC7+", bowel$"mCDF CTHRC1+")
mtotalcount<- sum(mat$"ssCDF PI16+", mat$"ssCDF FMO2+", mat$"ssCDF F3+", mat$"ssCDF GREM1+", mat$"ssCDF KCNN3+", 
                  mat$"iCDF ADAMDEC1+", mat$"iCDF MMP3+", mat$"apCDF CCL19+", 
                  mat$"mCDF LRRC7+", mat$"mCDF CTHRC1+")
itotalcount <- sum(inter$"ssCDF PI16+", inter$"ssCDF FMO2+", inter$"ssCDF F3+", inter$"ssCDF GREM1+", inter$"ssCDF KCNN3+", 
                   inter$"iCDF ADAMDEC1+", inter$"iCDF MMP3+", inter$"apCDF CCL19+", 
                   inter$"mCDF LRRC7+", inter$"mCDF CTHRC1+")

idents <- c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
            "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
            "mCDF LRRC7+", "mCDF CTHRC1+")
# Create an empty list to store rows
rows <- list()

# Loop over each fibro in idents and add values to a new row each time
for (fibro in idents) {
  new_row <- c(
    sum(bowel[[fibro]]) / btotalcount,
    sum(mat[[fibro]]) / mtotalcount,
    sum(inter[[fibro]]) / itotalcount
  )
  # Append the new row to the list
  rows[[fibro]] <- new_row
}
# Combine the list of rows into a dataframe
df <- do.call(rbind, rows)
# Assign row names
rownames(df) <- idents
# Assign column names
colnames(df) <- c("BW", "CF", "IF")
df <- df*100
print(df)
# Convert matrix to dataframe
df <- as.data.frame(df)

# Load required libraries
library(tidyr)
library(ggplot2)

# Add fibro as a column
df$fibro <- rownames(df)

# Reshape the dataframe into long format
df_long <- gather(data = df, key = "Variable", value = "Value", -fibro)
legend_order <- rev(idents)
df_long$fibro <- factor(df_long$fibro, levels = legend_order)
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
ggplot(df_long,aes(x=Variable,y=Value,fill=fibro)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Disease State Enrichment", 
       x="Disease State", y = "Percent (%)", fill = "fibro") + 
  scale_fill_manual(values = colors) + 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.title.x = element_blank(), legend.title = element_blank())  
ggsave("Figures/Fibroenrichmentspatial.svg")

#Apply CytoSpace metadata to all cells 
global <- readRDS("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace/longakerfineannotated.rds")
barcodes <- cell_locs$OriginalCID
global[["Barcodes"]]<- colnames(global) #Add in barcodes to metadata
global.subset <- subset(global, cells = barcodes)
Idents(global.subset) = "subpop.annotation"
DimPlot(global.subset)

#Add in metadata
global.subset[["SpotID"]] <- cell_locs$SpotID
global.subset[["row"]] <- cell_locs$row
global.subset[["col"]] <- cell_locs$col
global.subset[["Domain"]] <- cell_locs$Domain
global.subset[["Image"]] <- cell_locs$Image
global.subset[["CellType"]] <- cell_locs$CellType
#Plot all fibroblasts 
fibro.subset <- subset(global.subset, manual.annotation == "Fibroblasts") #Subset fibroblasts from object 
slice <- subset(fibro.subset, Image == "slice19")
legend_order = idents <- rev(c("ssCDF PI16+", "ssCDF FMO2+", "ssCDF F3+","ssCDF GREM1+", "ssCDF KCNN3+", 
                               "iCDF ADAMDEC1+", "iCDF MMP3+", "apCDF CCL19+", 
                               "mCDF LRRC7+", "mCDF CTHRC1+"))
