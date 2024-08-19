#KEBR CytoSpace 

## estimate cell type fractions using RCTD

.libPaths(c("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/r_packages", .libPaths()))
setwd("/oak/stanford/groups/longaker/KEBR/Human_Bowel/Visium/CytoSpace")

library(spacexr)
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(BPCells)

#You must run with the most updated SeuratObject and Seurat libraries otherwise the downstream functions fail!
#load in objects
fibro <- readRDS("finalfibroblastannotations.rds")
global <- readRDS("humanMerged_norm_findvar_sketch1000_harmony_umap_res0.2_project_joined.ANNOTATED.rds")
stricture_spatial <- readRDS("Longakerfibroticrefspat.rds")
normal_spatial <- readRDS("Normalbowelrefspat.rds") 
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

##Generate output files for cytospace

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
#Normal reference dataset
normal_spatial[["Spatial"]] <- as(normal_spatial[["Spatial"]], Class = "Assay")
for (sampleName in unique(normal_spatial$orig.ident)) {
  normal_spatial.subset <- subset(normal_spatial, subset = orig.ident == sampleName)
  imgName <- Images(normal_spatial)[grepl(sampleName, Images(normal_spatial))]
  generate_cytospace_from_ST_seurat_object(normal_spatial.subset, dir_out=paste0('cytospaceInput/',imgName), fout_prefix='', write_sparse=T, slice=imgName)
}

## estimate cell type fractions using RCTD
#Generate reference scRNAseq dataset 
library(spacexr)
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

#Generate cell fractions in normal spatial objects
rm(list = ls())
reference = readRDS('cytospaceInput/rctd/SCRefCoarse.rds')
query <- readRDS("Normalbowelrefspat.rds")
#Correct metadata to facilitate downstream analysis
query$orig.ident[query$orig.ident == "P2.noninf"] <- "slice2"
query@images$slice2 <-  query@images$P2.noninf
query$orig.ident[query$orig.ident == "P4.noninf"] <- "slice4"
query@images$slice4 <-  query@images$P4.noninf
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
query <- readRDS("Normalbowelrefspat.rds") 
#Correct metadata to facilitate downstream analysis
query$orig.ident[query$orig.ident == "P2.noninf"] <- "slice2"
query@images$slice2 <-  query@images$P2.noninf
query$orig.ident[query$orig.ident == "P4.noninf"] <- "slice4"
query@images$slice4 <-  query@images$P4.noninf

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

## add cytospace metadata
rm(list=ls())
spat <- readRDS("Normalbowelrefspat.rds")
file_list <- list.files(paste0("cytospaceOutput/"),"cell_type_assignments_by_spot.csv", recursive = TRUE)
file_list <- file_list[c(2,3)] #fibrotic slices
for (i in 1:length(file_list)) {
  file <- file_list[i]
  df <- read.csv(paste0("cytospaceOutput/",file), check.names = F)
  rownames(df) <- df[,1]
  df <- df[,-1]
  spat <- AddMetaData(spat, metadata = df)
}

file_list <- list.files(paste0("cytospaceOutput/"),"fractional_abundances_by_spot.csv", recursive = TRUE)
file_list <- file_list[c(2,3)] #fibrotic slices
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
saveRDS(spat, "NormalbowelrefspatCYTOSPACE.rds", compress = F)
