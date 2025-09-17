setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/SCENIC"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages_updated",.libPaths()))
library(matrixStats); library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(SeuratWrappers)
library(scales); library(SeuratObject); library(readr); library(tidyverse)
options(Seurat.object.assay.version = "v5")

####Global UMAP visualization####
fibro <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/SCENIC/Fibro_SCENIC_UMAP.rds")
#Update labeling 
barcode_to_annotation <- setNames(orig$fibro_annotation, colnames(orig))
fibro$fibro_annotation <- barcode_to_annotation[colnames(fibro)]
levels(fibro) <- "fibro_annotation"
levels = c("CTHRC1+ mFib", "LRRC7+ mFib", "MMP3+ iFib",
           "CD74+ apFib", "F3+ ssFib", "ADAMDEC1+ ssFib" , 
           "KCNN3+ ssFib", "GREM1+ ssFib",
           "THBS1+ ssFib", "PI16+ ssFib")

fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = levels)
Idents(fibro) <- "fibro_annotation"
DimPlot(fibro, cols= palette, raster = T)
fibro$fibro_annotation <- factor(fibro$fibro_annotation, levels = levels)
Idents(fibro) <- "fibro_annotation"
palette = rev(c(
  "#7b3294",  # deeper purple
  "#9467bd",  # original purple
  "#1f77b4",  # original blue
  "#17becf",  # original cyan
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff9900",  # more vivid orange
  "#ff7f0e",  # original orange
  "#db4d6d",  # new warm rose/red (fills the gap)
  #  "#e377c2",  # original pink
  "#8b0000"   # original dark red
))

DefaultAssay(fibro) <- "AUC"
DimPlot(fibro, cols= palette, raster = T) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave("Figures/globalumap.pdf")

#### MAT vs Bowel ####
#Generate controlled dataset
sub <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
sub <- subset(sub, specimen_type == "Resection")
sub <- subset(sub, tissue == "Colon", invert = T)
Idents(sub) <- "fibro_annotation"

#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"
Idents(sub) <- "newlabels"
#Bowel vs MAT
Idents(sub) <- "tissue"
mins <- min(table(sub$tissue))
downsample <- subset(sub, downsample = mins)
Idents(downsample) <- "fibro_annotation"
DimPlot(downsample, split.by = "tissue", cols = palette, pt.size = 0.05, raster =F) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave("Figures/matvsbowelumap.pdf")

#Disease State
sub$newlabels <- factor(x = sub$newlabels, levels = c("Normal Bowel", "Uninflamed Bowel","Inflamed Bowel" , "Stricture", 
                                                      "Normal MAT", "Uninvolved MAT", "CF"))
Idents(sub) <- "newlabels"
Idents(sub) <- "fibro_annotation"
mins <- min(table(sub$diseasestate))
downsample <- subset(sub, downsample = mins)
Idents(downsample) <- "fibro_annotation"
DimPlot(downsample, split.by = "diseasestate", cols = palette)
ggsave("Figures/diseasestateumap.pdf")

##highlight CTHRC1+ fibroblasts in the UMAP
# Identify cells labelled as "CTHRC1+ mFib"
cthrc1_cells <- WhichCells(fibro, expression = fibro_annotation == "CTHRC1+ mFib")
# Create a new metadata column for cell highlight
fibro$highlight <- NA
fibro$highlight <- ifelse(colnames(fibro) %in% cthrc1_cells, "highlight", "Other")

# Create a UMAP plot with all cells grey
DimPlot(fibro, group.by = "highlight") + 
  scale_color_manual(values = c("highlight" = "#8b0000", "other" = "grey92")) +
  theme(legend.position = "none") + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave("Figures/cthrc1highlight.pdf")

#### Signature Analysis ####
DefaultAssay(sub) <- "RNA"
Yapgenes <- c('ACAT1','ADAM9','CCN1','KLF5','ACTA2','COL4A1','LIF','LACTB',
              'AMOTL2','IL11','SORCS3','ANKRD1','CCN2','C11orf96','LBH','ANLN',
              'HBEGF','EDN1','MAP2K6','AURKB','HSPG2','ATF3','ARHGAP42','BIRC5',
              'LAMC1','F3','ADRA1B','CASP3','LPL','SERPINE1','CRIM1','CASP8','IER3',
              'C18orf63','CAVIN2','THBS1','TRIB1','CDH2','CCNF','TIMP2','FAM178A','CDC42',
              'VCAM1','KLF10','SLC38A2','CDK1','CALM2','SMAD7','ABI3BP','CDK6','COL1A1','FZD8',
              'TBX18','COL1A2','RUNX1','DYNLRB2','COL5A2','FGF2','NRD1','DIAPH3','COL6A2',
              'ZRANB1','DICER1','CXCL12','TNS3','E2F1','FBN1','MUSK','INHBA','USP44','EGFR',
              'LAMA2','PMS1','FGF1','LAMB1','PRELID2','FN1','TFPI','JDP2','FOXM1','VEGFA','AL035681.1',
              'GAB1','BDNF','BASP1','IGF1R','CCL2','QSOX1','LATS2','COL4A2','DDHD1','LDHA','CSF1','CITED2',
              'LIFR','SEMA3E','DUSP10','MCM6','ANXA1','OR51B6','MYC','BGN','EFHC1','MYL9','COL11A1','MAGI2',
              'NDUFB3','COL14A1','JARID2','NEGR1','COL5A1','CCDC132','NEK2','COL6A1','ARNT2','NF2','ANGPT1',
              'PELO','PLK1','ADAM12','CENPV','PMAIP1','CCL9','RPL22L1','POLA1','PTN','SYNDIG1','PPP1R3B','SEMA3A',
              'RAD18','PTEN','THBS2','TSPAN3','RAD51','HMGB1','TBC1D2','SLC2A3','IRS1','SNAI2','NUAK2','TAGLN','TRPC3',
              'TEAD4','GAS6','VIM','SNTB1','WWC1','MTMR10','ZEB1','BMP5','FRY','IGFBP3','AC112715.2','SUSD1','RUNX2',
              'ROBO2','AC023590.1','FBXO47','FAM107A','RBMS3','SH3RF1','SYT16','XPNPEP1','PIK3C2G',
              'MED13L','RND3','CDH4','TSPAN18','NEDD9','KPNA3','SH2D4A','TRIO','DPYD','LRRTM1','ZNF469',
              'RP11-553A10.1','GLIS3','SLITRK6','IQCJ-SCHIP1','FGF5','NTRK3','ADAMTS8','NAALADL2',
              'PTGER4','ANXA3','ERCC2','FAM172A','GNB1L','PSMG2','NMBR','SLC38A4','ROR1','PARD3B','BET1', 
              "FOXF2", "CCDC80", "MYOF","FJX1","PTPN14","ARHGEF17", "DOCK5","NT5E", "ASAP1", "GADD45A", "TGFB2", "AXL",
              "YAP1", "CTNNB1", "MST1",
              "PTK2", "ITGA5", "ITGAV", "ITGB1",  
              "SAV1", "STK3", "LATS1", "LATS2", 
              "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
              "VIM", "VCL", "VASP")

Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(sub)] #retain YAP genes only in rownames of sub 
genes.to.keep <- Matrix::rowSums(sub[["RNA"]]$counts > 0) >= floor(0.1 * ncol(sub[["RNA"]]$counts))
counts.sub <- sub[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts

ecmgenes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL24A1", 
              "LUM", "DCN", "BGN", "HSPG2", "AGRN", 
              "FN1", "VTN", "TGFBI", "NID1", "NID2", "LAMC1", "LAMA1", "LAMA2", "SPARC", 
              "FBLN1", "FBLN2", "FBLN5", "LTBP1", "LTBP5", "EMILIN1", "MFAP4", "EFEMP1", 
              "FBN1", "FBN2", "TNC")

sub <- AddModuleScore(sub, features = list(Yapsub), name = "Mechano")
sub <- AddModuleScore(sub, features = list(ecmgenes), name = "ECMgen")

DefaultAssay(sub) <- "AUC"
library(RColorBrewer)
FeaturePlot(sub, features = c("Mechano1")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("Figures/mechanosig.pdf")

FeaturePlot(sub, features = c("Mechano1"), split.by = "diseasestate") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("Figures/mechanosigbyDState.pdf")

FeaturePlot(sub, features = c("ECMgen1")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("Figures/ECMsig.pdf")

FeaturePlot(sub, features = c("ECMgen1"), split.by = "diseasestate") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("Figures/ECMsigbyDState.pdf")

DotPlot(fibro, features = c("ECMgen1", "Mechano1"), group.by = "diseasestate", scale.max = 100, scale.min = 0, cols = "RdBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size))
  scale_x_discrete(labels=c("ECM", "Mechanical")) + 
  scale_y_discrete(labels=c('Normal', "Uninvolved", "Stricture", "Normal", "Uninvolved", "CF"))

#### Heatmap of regulons ####
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
options(Seurat.object.assay.version = "v5")

fibro <- readRDS("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/SCENIC/Fibro_SCENIC_UMAP.rds")
# Step 1: Extract regulon activity matrix
regulon_matrix <- fibro[["AUC"]]$counts

# Step 2: Clean gene names FIRST
original_names <- rownames(regulon_matrix)
clean_names <- gsub("\\.{3}$", "", original_names)
rownames(regulon_matrix) <- clean_names  # Update in original matrix

# Step 3: Set identities
Idents(fibro) <- "fibro_annotation"

# Step 4: Compute average activity with CLEAN names
avg_activity <- sapply(levels(fibro), function(cluster) {
  cells <- WhichCells(fibro, idents = cluster)
  rowMeans(regulon_matrix[, cells, drop = FALSE])
})

# Step 5: Z-score normalization
z_scores <- t(scale(t(avg_activity)))
sorted_zscores <- t(apply(z_scores, 1, function(x) sort(x, decreasing = TRUE)))

#Find top z scores
top_50 <- z_scores[order(z_scores[,"CTHRC1+ mFib"], decreasing = TRUE), "CTHRC1+  mFib"]
head(top_50, 50)

# Step 6: Define genes using CLEANED names
genes_to_label <- c("WT1", "HOXC10", "CREB3L1", "TWIST1", "HOXC4", "EGR3" ,
                    "PRRX1", "GATA6", "SIRT6", "GLI3", "BACH1", "BACH2", 
                    "KLF16", "CEBPG")
  
  # c("WT1", "HOXC10", "CREB3L1", "TWIST1", "HOXC4", 
  #                   "PRRX1", "GATA6",  "SIRT6", "ZNF274", "GLI3", 
  #                   "KLF16", "NFATC1", "CEBPG")

# Step 7: Create SINGLE heatmap with annotations
library(ComplexHeatmap)
library(circlize)

ha = rowAnnotation(foo = anno_mark(at = which(rownames(z_scores) %in% genes_to_label),
                                   labels = rownames(z_scores)[rownames(z_scores)%in%genes_to_label]))
ht <- Heatmap(
  z_scores,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE, 
  right_annotation = ha
)
#### Binary heatmap of regulon contribution to signatures ####
# Load required libraries
library(ComplexHeatmap)
library(readr)
library(dplyr)

# Define regulon names and target genes (in specified order)
regulon_names <- c("WT1", "HOXC10", "CREB3L1", "TWIST1", "HOXC4", 
                   "PRRX1", "GATA6", "SIRT6", "GLI3", "BACH1", "BACH2", 
                   "KLF16", "CEBPG")
target_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL24A1", "LAMA1", "SPARC", "BGN", "CNN3","PTK2", "ITGB1", 
                  "ITGAV", "SDC1", "YAP1", "AMOTL1", "SMAD7", "IL1RL1",  "WWTR1", "POSTN", "PDGFD" ,  "CTHRC1", "FAP", "CDH11")

# Path to regulon files
regulon_dir <- "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/SCENIC/Regulons"

# Initialize list to store regulon genes
regulon_gene_lists <- list()

# Read each regulon file
for (reg in regulon_names) {
  file_name <- paste0(reg, "_+_regulon_genes.tsv")
  file_path <- file.path(regulon_dir, file_name)
  
  if (file.exists(file_path)) {
    regulon_genes <- read_tsv(file_path, col_names = FALSE, show_col_types = FALSE) %>%
      pull(X1)  # Assuming genes are in the first column
    regulon_gene_lists[[reg]] <- regulon_genes
  } else {
    stop(paste("Regulon file not found:", file_path))
  }
}

# Create binary matrix
binary_matrix <- matrix(0, nrow = length(target_genes), ncol = length(regulon_names),
                        dimnames = list(target_genes, regulon_names))

for (tg in target_genes) {
  for (reg in regulon_names) {
    if (tg %in% regulon_gene_lists[[reg]]) {
      binary_matrix[tg, reg] <- 1
    }
  }
}

# Create heatmap
Heatmap(binary_matrix,
        name = "Presence",
        col = c("white", "black"),
        border = TRUE,
        
        # Axis labels positioning
        column_title = "Regulons",
        row_title = "Target Genes",
        column_title_side = "bottom",  # Moved to bottom
        row_title_side = "left",
        
        # Font size adjustments
        row_names_gp = gpar(fontsize = 10),  # Row gene names
        column_names_gp = gpar(fontsize = 10),  # Column regulon names
        
        # Formatting
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "left",
        column_names_rot = 45,  # Angled column names
        rect_gp = gpar(col = "gray80", lwd = 0.5),
        
        heatmap_legend_param = list(
          title = "Present",
          at = c(0, 1),
          labels = c("No", "Yes"),
          border = "black"
        ))


#### MAT subanalysis ####
setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/SCENIC/MAT")
fibro <- readRDS("new_MAT_Fibro_SCENIC_UMAP.rds")
Idents(fibro) <- "fibro.annotation"
palette = rev(c(
  "#7b3294",  # deeper purple
  "#9467bd",  # original purple
  "#1f77b4",  # original blue
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff7f0e",  # original orange
  "#8b0000"   # original dark red
))
ggsave("Figures/globalumap.pdf")

#### MAT only analysis ####
#Motivated by the difference in TF level expression btwn bowel and MAT
mat <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/SCENIC/new_MAT_Fibro_SCENIC_UMAP.rds")
palette = rev(c(
  "#7b3294",  # deeper purple
  "#1f77b4",  # original blue
  "#2ca02c",  # original green
  "#bcbd22",  # original olive
  "#f5b800",  # original yellow
  "#ff7f0e",  # original orange
  "#e377c2",  # original pink
  "#8b0000"   # original dark red
))

DimPlot(mat, cols= palette, raster = F, pt.size = 0.01) + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave("Figures/matumaponly.pdf")

#highlight CTHRC1
# Identify cells labelled as "CTHRC1+ mFib"
cthrc1_cells <- WhichCells(mat, expression = mat.annotation == "NOTCH3+ Pericytes")
# Create a new metadata column for cell highlight
mat$highlight <- NA
mat$highlight <- ifelse(colnames(mat) %in% cthrc1_cells, "highlight", "Other")

DimPlot(mat, group.by = "highlight") + 
  scale_color_manual(values = c("highlight" = "#8b0000", "other" = "grey92")) +
  theme(legend.position = "none") + 
  coord_fixed() +  # Make the UMAP square
  theme_void() +   # Removes all gridlines, axis lines, labels, and ticks
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
ggsave("Figures/matonlycthrc1.pdf")

#Find regulons that are significantly upregulated in fibroblasts
mat <- JoinLayers(mat)
mat[["AUC"]]$data <- as(object = mat[["AUC"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(mat, max.cells.per.ident = 1000, only.pos = T, group.by = "mat.annotation", 
                          assay = "AUC", test = "MAST")
write.csv(markers, "mat_clusterregulons.csv")


#Heatmap
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
options(Seurat.object.assay.version = "v5")

# Step 1: Extract regulon activity matrix
regulon_matrix <- mat[["AUC"]]$counts

# Step 2: Clean gene names FIRST
rownames(regulon_matrix) <- paste0(rownames(regulon_matrix), '(+)')

# Step 3: Set identities
Idents(mat) <- "mat.annotation"

# Step 4: Compute average activity with CLEAN names
avg_activity <- sapply(levels(mat), function(cluster) {
  cells <- WhichCells(mat, idents = cluster)
  rowMeans(regulon_matrix[, cells, drop = FALSE])
})

# Step 5: Z-score normalization
z_scores <- t(scale(t(avg_activity)))
sorted_zscores <- t(apply(z_scores, 1, function(x) sort(x, decreasing = TRUE)))

#Find top z scores
top_50 <- z_scores[order(z_scores[,"CTHRC1+ mFAP"], decreasing = TRUE), "CTHRC1+ mFAP"]
head(top_50, 50)

# Step 6: Define genes using CLEANED names
genes_to_label <- c("HES1(+)",  "CREB3L1(+)", "STAT2(+)", "NFKB2(+)", 
                    "REL(+)", "JUNB(+)", "SOX4(+)", "EGR3(+)", "KLF6(+)", "TWIST1(+)")

# c("WT1", "HOXC10", "CREB3L1", "TWIST1", "HOXC4", 
#                   "PRRX1", "GATA6",  "SIRT6", "ZNF274", "GLI3", 
#                   "KLF16", "NFATC1", "CEBPG")

# Step 7: Create SINGLE heatmap with annotations

ha = rowAnnotation(foo = anno_mark(at = which(rownames(z_scores) %in% genes_to_label),
                                   labels = rownames(z_scores)[rownames(z_scores)%in%genes_to_label]))
ht <- Heatmap(
  z_scores,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), brewer.pal(11, "RdYlBu")[c(11, 6, 1)]),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE, 
  right_annotation = ha
)

pdf(file="Figures/regulonheatmap.pdf", width = 3, height = 5)
draw(ht)
dev.off()

# Load libraries
library(igraph)
library(ggraph)
library(tidyverse)

# Read the tsv file

df<- read.csv("Regulons/CREB3L1_+_regulon_genes.tsv" ) # specify the path to your excel file

# Assuming your Excel file has two columns: 'TF' and 'Target'
# Check the first few rows of the dataframe
head(df)
df$TF <- "CREB3L1"
df$Target <- df$gene
# Create a graph object in igraph
# Assuming `df` has columns 'TF' for transcription factors and 'Target' for target genes
edges <- df %>% rename(from = TF, to = Target)
graph <- graph_from_data_frame(edges, directed = TRUE)

# Visualize using ggraph
ggraph(graph, layout = 'fr') +   # 'fr' refers to the Fruchterman-Reingold layout algorithm
  geom_edge_link(aes(edge_alpha = 0.8)) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

#Comparison across conditions
auc_mat <- as.matrix(mat[["AUC"]]$counts)

gene_metadata <- as.data.frame(t(auc_mat))

# Make sure row names of gene_metadata are cell barcodes
rownames(gene_metadata) <- colnames(auc_mat)
colnames(gene_metadata) <- paste0(colnames(gene_metadata), '(+)')


# Add this new metadata to the Seurat object
# Merge with the existing metadata
new <- AddMetaData(mat, metadata = gene_metadata)

new <- subset(new, diseasestate == "Ileostomy MAT", invert = T)

# Verify that the metadata has been added correctly
head(new[[]])
#Create a barplot showing the percentage changes 
new$newlabels <- "NA"
new$newlabels[new$diseasestate == "Normal MAT"] <- "Normal"
new$newlabels[new$diseasestate == "Uninvolved"] <- "Uninvolved"
new$newlabels[new$diseasestate == "Involved"] <- "Creeping Fat"
Idents(new) <- "newlabels"

new$newlabels <- factor(new$newlabels, levels = c("Normal", "Uninvolved", "Creeping Fat"))
Idents(new) <- "newlabels"
Idents(new) <- "mat.annotation"

cthrc1 <- subset(new, idents = "CTHRC1+ mFAP")
Idents(cthrc1) <- "newlabels"

DotPlot(cthrc1, features = genes_to_label, cols = "RdYlBu", scale.min = 0, scale.max = 100) + coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("Figures/regulondotplot.pdf")
  
#### Construct TF network analysis ####

# Define regulon names and target genes (in specified order)
regulon_names <- c("REL", "KLF6", "HES1", "EGR3", "CREB3L1", 
                   "MYB", "TWIST1", "NFKB2", "STAT2", "SOX4")

# Path to regulon files
regulon_dir <- "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/SCENIC/Regulons"

# Initialize list to store regulon genes
regulon_gene_lists <- list()

# Read each regulon file
for (reg in regulon_names) {
  file_name <- paste0(reg, "(+)_regulon_genes.tsv")
  file_path <- file.path(regulon_dir, file_name)
  
  if (file.exists(file_path)) {
    regulon_genes <- readr::read_tsv(file_path, show_col_types = FALSE)
    regulon_gene_lists[[reg]] <- regulon_genes
  } else {
    stop(paste("Regulon file not found:", file_path))
  }
}

fibro <- readRDS("/home/kbauerro/Oak/IBDAFs/MetaUpdate/globalfibroannotated_0.3.rds")
fibro[["RNA"]]$data <- as(object = fibro[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(fibro, max.cells.per.ident = 1000, only.pos = T)

# Filter markers for the "CTHRC1+ mFib" population
fib_markers <- markers[markers$cluster == "CTHRC1+ mFib", ]

# Extract top 200 genes based on adjusted p-values for "CTHRC1+ mFib"
top200_genes <- fib_markers[order(fib_markers$p_val_adj), "gene"][1:200]
top200_genes <- as.vector(top200_genes)

# Print relevant variables for debugging
print("Top 200 genes:")
print(top200_genes)

# Find intersection between top 200 genes and genes in regulon_gene_lists
all_regulated_genes <- unlist(regulon_gene_lists)
intersection_genes <- intersect(as.vector(all_regulated_genes), top200_genes)
intersection_genes <- as.vector(intersection_genes)

print("Intersection genes:")
print(intersection_genes)

# Assuming each element in regulon_gene_lists is a data frame with a 'gene' column
filtered_regulon_gene_lists <- lapply(regulon_gene_lists, function(regulon) {
  # Extract gene names from the regulon's data frame
  regulon_genes <- regulon$gene
  # Find intersection with the top genes
  intersect(regulon_genes, intersection_genes)
})

# Filter out regulons with no overlapping genes
filtered_regulon_gene_lists <- filtered_regulon_gene_lists[sapply(filtered_regulon_gene_lists, length) > 0]

# Convert the list to an edge data frame
edges <- do.call(rbind, lapply(names(filtered_regulon_gene_lists), function(tf) {
  data.frame(from = tf, to = filtered_regulon_gene_lists[[tf]])
}))

# Create graph from data frame
g <- graph_from_data_frame(d = edges, directed = TRUE)

# Define color palette
color_palette <- rainbow(length(names(filtered_regulon_gene_lists)))

# Map colors to TFs
tf_colors <- setNames(color_palette, names(filtered_regulon_gene_lists))

# Filter nodes with degree greater than threshold
degree_threshold <- 0
high_degree_nodes <- V(g)[degree(g) > degree_threshold]

# Create subgraph with nodes having degree greater than threshold
sub_g <- induced_subgraph(g, high_degree_nodes)

# Assign colors to nodes based on TFs
V(sub_g)$color <- sapply(V(sub_g)$name, function(node) {
  if (node %in% names(filtered_regulon_gene_lists)) {
    tf_colors[node] # Color for TF
  } else {
    "white" # Default color for target genes
  }
})

# Keep the outline color to distinguish nodes
V(sub_g)$frame.color <- "black"
V(sub_g)$label.color <- "black" # Color for text labels
V(sub_g)$label.cex <- 0.7       # Smaller font size for labels

# Assign edge colors based on source (TF)
E(sub_g)$color <- sapply(E(sub_g), function(edge) {
  from <- ends(sub_g, edge)[1] # Get the source node of the edge
  if (from %in% names(filtered_regulon_gene_lists)) {
    tf_colors[from] # Color for edge if it originates from a TF
  } else {
    "grey" # Default color for edges not originating from TFs
  }
})

# Ensure that vertices and labels are displayed correctly, with filled colors and outlines
V(sub_g)$shape <- "circle"
V(sub_g)$size <- 15  # Larger size for clarity

# Open PDF device
# Open PDF device
pdf("Filtered_Colored_TF_Network.pdf", width = 12, height = 12)  # Slightly larger canvas

# Use Fruchterman-Reingold layout with stronger repulsion
coords <- layout_with_fr(
  sub_g,
  niter = 10000,
  repulserad = vcount(sub_g)^3 * 5,  # Increased repulsion radius
  start.temp = 0.1,                  # Lower initial temperature
  coolexp = 1.5                      # Slower cooling to allow more movement
)

# Plot the filtered graph with smaller elements
plot(sub_g, 
     vertex.label = V(sub_g)$name,
     vertex.label.color = V(sub_g)$label.color, 
     vertex.label.cex = 0.5,          # Smaller text (reduced from 0.7)
     vertex.shape = V(sub_g)$shape,
     vertex.color = V(sub_g)$color,
     vertex.size = 8,                 # Smaller nodes (reduced from 15)
     vertex.frame.color = V(sub_g)$frame.color,
     edge.width = 0.3,                # Thinner edges
     edge.arrow.size = 0.2,           # Smaller arrows
     layout = coords,
     main = "Filtered Colored Transcription Factor Network")

# Close PDF device
dev.off()

# Close PDF device
dev.off()

#### Feature plots of regulons ####
regulons <- c("KLF6", "TWIST1", "CREB3L1", "REL")
# Create individual FeaturePlots for each regulon feature
FeaturePlot(fibro, features = regulons, order = FALSE, raster = FALSE, ncol = 2, max.cutoff = "q95") &
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu")), guide = "none") & 
  theme_void() & 
  theme(
    plot.title = element_text(face = "plain", size = rel(1), hjust = 0.5, vjust = 1), # Centered title
    plot.background = element_rect(fill = NA, colour = "black"),
    aspect.ratio = 1) # Ensure square aspect ratio

ggsave("Figures/regulonfeatureplots.pdf")