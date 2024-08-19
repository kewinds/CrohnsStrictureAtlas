setwd("/oak/stanford/groups/longaker/KEBR/IBDAFs/newHuman"); set.seed(54000)
.libPaths(c("/oak/stanford/groups/longaker/KEBR/IBDAFs/r_packages_updated",.libPaths()))
library(Seurat); library(ggplot2); library( dplyr )
library(BPCells); library(Matrix); library(SeuratWrappers)
library(scales); library(readxl); 
options(Seurat.object.assay.version = "v5")

#Downstream analysis of fibroblasts 
#Reorder factor levels 
##### Annotate fibroblasts ####
fibro <- readRDS("finalfibroblastannotations.rds")
Idents(fibro) <- "seurat_clusters"
fibro$fibro.annotation <- "tbd"
fibro$fibro.annotation[fibro$seurat_clusters == 0] <- "Steady-State PI16+"
fibro$fibro.annotation[fibro$seurat_clusters == 1] <- "Inflammatory ADAMDEC1+"
fibro$fibro.annotation[fibro$seurat_clusters == 2] <- "Steady-State FMO2+"
fibro$fibro.annotation[fibro$seurat_clusters == 3] <- "Steady-State F3+"
fibro$fibro.annotation[fibro$seurat_clusters == 4] <- "Mechanosensitive LRRC7+"
fibro$fibro.annotation[fibro$seurat_clusters == 5] <- "Mechanosensitive CTHRC1+"
fibro$fibro.annotation[fibro$seurat_clusters == 6] <- "Steady-State KCNN3+"
fibro$fibro.annotation[fibro$seurat_clusters == 7] <- "Steady-State GREM1+"
fibro$fibro.annotation[fibro$seurat_clusters == 8] <- "Antigen-Presenting CD74+"
fibro$fibro.annotation[fibro$seurat_clusters == 9] <- "Inflammatory MMP3+"
fibro$fibro.annotation[fibro$seurat_clusters == 10] <- "Steady-State KCNN3+"
fibro$fibro.annotation[fibro$seurat_clusters == 11] <- "Steady-State PI16+"
Idents(fibro) <- "fibro.annotation"

fibro$fibro.annotation <- factor(fibro$fibro.annotation, levels = rev(c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                        "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                        "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")))
Idents(fibro) <- "fibro.annotation"
#DotPlot of relevant genes
library(RColorBrewer)
#num_colors <- length(levels(fibro))
#colors <- brewer.pal(num_colors, "RdYlBu")
# colors <- c("#700d2b", "#c0242c", "#E34234", "#f47422", "#f2af1e", 
#                      "#8da13b", "#35bdcc", "#0a95ce","#09415f", "#471466")
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
dimplot <- DimPlot(fibro, raster = T, cols = colors)
revfibro <- fibro
revfibro$fibro.annotation <- factor(revfibro$fibro.annotation, levels = c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                        "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                        "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+"))
Idents(revfibro) <- "fibro.annotation"
genes <- rev(c("PI16", "FMO2", "GREM1","KCNN3", "F3", "ADAMDEC1", "MMP3", "CD74", "LRRC7", "CTHRC1"))
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
ggsave("Figures/fibroblastUMAP.svg")

##### Add module scores to get mechanosensitive signatures ####
#1. YAP gene signature 
#Yapgenes <-as.data.frame(read_excel("AltYAP_list2.xlsx",col_names = "Genes"))
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
              "FOXF2", "CCDC80", "MYOF","FJX1","PTPN14","ARHGEF17", "DOCK5","NT5E", "ASAP1", "GADD45A", "TGFB2", "AXL")
Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- fibro[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
revfibro <- AddModuleScore(revfibro, features = list(Yapsub), name = "Yapsig") #Get score for YAP targets
#2. YAP components signature
components <- list(c("YAP1", "CTNNB1", "MST1",
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))
revfibro <- AddModuleScore(revfibro, features = components, name = "Components")
#3. ECM protein signature 
ecmgenes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL24A1", 
              "LUM", "DCN", "BGN", "HSPG2", "AGRN", 
              "FN1", "VTN", "TGFBI", "NID1", "NID2", "LAMC1", "LAMA1", "LAMA2", "SPARC", 
              "FBLN1", "FBLN2", "FBLN5", "LTBP1", "LTBP5", "EMILIN1", "MFAP4", "EFEMP1", "FBN1", "FBN2", "TNC")
revfibro <- AddModuleScore(revfibro, features = list(ecmgenes), name = "ECMgenes")
DotPlot(revfibro, features = c("ECMgenes1", "Components1", "Yapsig1"), scale.max = 100, scale.min = 0, cols = "RdYlBu") + theme(axis.title.x = element_blank(),
                                                      axis.title.y = element_blank(), 
                                                      axis.text.x = element_text(size = 15),
                                                      axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
                                                      legend.key.size = unit(0.8, 'cm'), #change legend key size
                                                      legend.key.height = unit(0.4, 'cm'), #change legend key height
                                                      legend.key.width = unit(0.4, 'cm'), #change legend key width
                                                      legend.title = element_text(size=9), #change legend title font size
                                                      legend.text = element_text(size=9)) +  #change legend text font size)) +
                                                      scale_x_discrete(labels=c('ECM Genes', "YAP Target Genes", "YAP Pathway Genes"))
ggsave("Figures/YAPsignatures.svg" , width = 10, height = 4, units = "in")
sub <- subset(revfibro, idents = c("Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")) 
sub <- subset(revfibro, specimen_type == "Resection")
DotPlot(sub, features = "Yapsig1", group.by = "diseasestate", scale.max = 100, scale.min = 0) #MS fibroblasts only

#### Regional analysis of bowel and MAT fibroblasts ####
#Create new metadata category called general to generalize regions
Idents(fibro) <- "fibro.annotation"
fibro$general[fibro$tissue == "MAT"] <- "MAT"
fibro$general[fibro$tissue == "SB"] <- "Bowel"
fibro$general[fibro$tissue == "Colon"] <- "Bowel"
sub <- subset(fibro, subset = (diseasestate == "Ileostomy Bowel" | diseasestate == "Ileostomy MAT"), invert = T) #Remove ileostomy samples
sub <- subset(sub, specimen_type == "Resection") #Only compare full thickness due to mucosal bias of endoscopies
sub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T)
table(sub$diseasestate, sub$tissue)
DimPlot(sub, raster = T, cols = colors, split.by = "general") 
ggsave("Figures/BowelvsMATfibro.svg")

bowel <- subset(fibro, general == "Bowel")
DimPlot(bowel, raster = T, cols = colors, split.by = "specimen_type") 
ggsave("SuppFigures/ResectionvsFT.svg")
####Activation of mechanosensitive fibroblasts ####
#Do this analysis with just the controlled dataset
sub <- readRDS("fibrocontrolled_dataset.rds")
Yapgenes <- unique(Yapgenes)
Yapgenes <- Yapgenes[Yapgenes %in% rownames(fibro)] #retain YAP genes only in rownames of fibro 
genes.to.keep <- Matrix::rowSums(fibro[["RNA"]]$counts > 0) >= floor(0.1 * ncol(fibro[["RNA"]]$counts))
counts.sub <- sub[["RNA"]]$counts[genes.to.keep,]
Yapsub <- Yapgenes[Yapgenes %in% rownames(counts.sub)] #Remove genes not in 10% of fibroblasts
sub <- AddModuleScore(sub, features = list(Yapsub), name = "Yapsig") #Get score for YAP targets
sub <- AddModuleScore(sub, features = components, name = "Components")
sub <- AddModuleScore(sub, features = list(ecmgenes), name = "ECMgenes")
Idents(sub) <- "fibro.annotation"
CP <- subset(sub, idents = c( "Steady-State PI16+", "Steady-State FMO2+", "Mechanosensitive CTHRC1+" , "Mechanosensitive LRRC7+"))
#Create a barplot showing the percentage changes 
CP$newlabels <- "NA"
CP$newlabels[CP$diseasestate == "Normal MAT"] <- "Normal"
CP$newlabels[CP$diseasestate == "Uninvolved"] <- "Uninvolved"
CP$newlabels[CP$diseasestate == "Involved"] <- "Stricture"
CP$newlabels[CP$diseasestate == "Normal Bowel"] <- "Normal"
CP$newlabels[CP$diseasestate == "Stricture Uninflamed"] <- "Uninvolved"
CP$newlabels[CP$diseasestate == "Stricture"] <- "Stricture"
Idents(CP) <- "newlabels"
dotplot <- DotPlot(CP, features = c("Components1", "Yapsig1", "ECMgenes1"), 
        cols = "RdYlBu", split.by = "fibro.annotation", scale.min = 0, scale.max = 100) 
#Rearrange the data 
dotplot$data$id <- factor(dotplot$data$id, levels = c("Normal_Steady-State PI16+", "Uninvolved_Steady-State PI16+", "Stricture_Steady-State PI16+", 
                                                      "Normal_Steady-State FMO2+", "Uninvolved_Steady-State FMO2+", "Stricture_Steady-State FMO2+", 
                                                      "Normal_Mechanosensitive LRRC7+", "Uninvolved_Mechanosensitive LRRC7+", "Stricture_Mechanosensitive LRRC7+", 
                                                      "Normal_Mechanosensitive CTHRC1+", "Uninvolved_Mechanosensitive CTHRC1+", "Stricture_Mechanosensitive CTHRC1+"))
                                                      
dotplot + scale_y_discrete(labels=c('Normal', "Uninvolved", "Stricture",
                                    'Normal', "Uninvolved", "Stricture", 
                                    'Normal', "Uninvolved", "Stricture", 
                                    'Normal', "Uninvolved", "Stricture")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("SuppFigures/CTHRC1activationsigsall.svg")

#Myofibroblast markers 
fibro$fibro.annotation <- factor(x = fibro$fibro.annotation, levels = c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                        "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                        "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+"))
Idents(fibro) <- "fibro.annotation"
DotPlot(fibro, features = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "CDH11", "POSTN", "FAP", "SPARC"), 
        cols = "RdYlBu", scale.min = 0, scale.max = 100) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
ggsave("SuppFigures/pathogenicfibroblasts.svg")

#Module scores for showing regional YAP differences
sub <- AddModuleScore(sub, features = list(Yapgenes), name = "Yapsig")
sub <- AddModuleScore(sub, features = components, name = "Components")
sub <- AddModuleScore(sub, features = ecmgenes, name = "ECMgenes")
subsub <- subset(sub, diseasestate == "Stricture Inflamed", invert = T) #Eliminate stricture inflamed for comparable comparisons
DotPlot(subsub, features = c("ECMgenes1", "Yapsig1", "Components1"), group.by = "diseasestate", scale.max = 100, scale.min = 0, cols = "RdBu") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)) +  #change legend text font size))
        scale_x_discrete(labels=c('ECM Genes', "YAP Target Genes", "YAP Pathway Genes")) + 
        scale_y_discrete(labels=c('Normal', "Uninvolved", "Stricture", "Normal", "Uninvolved", "CF"))
ggsave("Figures/RegionalYAPsignatures.svg", width = 10, height = 4, units = "in")
#DimPlot MAT cells only
#Change labels for clarity
mes <- subset(fibro, tissue == "MAT")
mes <- subset(mes, diseasestate == "Ileostomy MAT", invert = T)
mes$newlabels <- "NA"
mes$newlabels[mes$diseasestate == "Normal MAT"] <- "Normal"
mes$newlabels[mes$diseasestate == "Uninvolved"] <- "Uninvolved"
mes$newlabels[mes$diseasestate == "Involved"] <- "Stricture"
mes$newlabels <- factor(x = mes$newlabels, levels = c("Normal", "Uninvolved", "Stricture"))
#Reorder levels
levels(mes) <- rev(c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                     "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                     "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")) 
Idents(mes) <- "newlabels"
#Equal numbers of cells 
disease_idents <- unique(mes$newlabels)
cell.list <- WhichCells(mes, idents = disease_idents, downsample = min(table(mes$newlabels)))
mes_small <- mes[, cell.list]
table(mes_small$newlabels, mes_small$fibro.annotation)

submes <- subset(mes_small, idents = c("Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                       "Antigen-Presenting CD74+"), invert = T) #remove bowel fibroblasts
submes <- subset(submes, subset = (umapharmony_1 < 0 & umapharmony_2 > -5)) #remove extra umap space
colors <- c("#660000", "#cc0000", "#cfe2f3", "#0b5394","#073763")

DimPlot(submes, split.by = "newlabels", cols = colors, raster = F)
#Eliminate contaminating bowel fibroblasts on corners of plot 
ggsave("Figures/MAT_MSexpansion.svg", width = 8.5, height = 3, units = "in")



#Create a barplot showing the percentage changes 
sub$newlabels <- "NA"
sub$newlabels[sub$diseasestate == "Normal MAT"] <- "Normal MAT"
sub$newlabels[sub$diseasestate == "Uninvolved"] <- "Uninvolved MAT"
sub$newlabels[sub$diseasestate == "Involved"] <- "CF"
sub$newlabels[sub$diseasestate == "Normal Bowel"] <- "Normal Bowel"
sub$newlabels[sub$diseasestate == "Stricture Uninflamed"] <- "Uninvolved Bowel"
sub$newlabels[sub$diseasestate == "Stricture"] <- "Stricture"


sub$fibro.annotation <- factor(x = sub$fibro.annotation, levels = rev(c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                        "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                        "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+")))
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
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
                     "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
tab = table(sub$newlabels, sub$fibro.annotation)
ptab = prop.table(tab, 1)*100
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
ptab = as.data.frame(ptab)
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
ggsave("Figures/proptablefibroblasts.svg", width = 8, height = 4, units = "in" )

#### Representative signature gene plots ####

genes <- c("YAP1", "PTK2", "TEAD1", "LATS1", "LATS2", 
           "CCN1", "CCN2", "SERPINE1", "ADAM12", "RUNX1", 
           "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2")
DotPlot(revfibro, features = genes, cols = "RdYlBu", scale.min = 0, scale.max = 100) +   theme(axis.title.x = element_blank(),
                                                                                               axis.title.y = element_blank(),
                                                                                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                               legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                               legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                               legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                               legend.title = element_text(size=6), #change legend title font size
                                                                                               legend.text = element_text(size=6)) 

ggsave("SuppFigures/YAPsig_repgenes.dotplot.svg")

#### Migration Genes ####
migration <- c("CTHRC1", "TNC", "CEMIP", "CDH2", "CD44", "ARID5A", "DDR2", "IQGAP1", 
               "SLC8A1", "TNS1", "FGF2", "PLAU", "PRKCE", "PODXL", "ID1", 
               "CNN1", "SNAI2", "ACTC1", "MEOX2", "BMP2", "LUM", "FGF18", 
               "PIK3CA", "EDN1", "GRB2", "THBS1", "ITGB1", "ITGA3", 
               "ZEB1", "ZEB2", "TWIST1", "TWIST2")
revfibro <- AddModuleScore(revfibro, features = list(migration), name = "Migration")
DotPlot(revfibro, features = "Migration1", cols = "RdYlBu", scale.min = 0, scale.max = 100) + theme(axis.title.x = element_blank(),
                                                                                                    axis.title.y = element_blank(), 
                                                                                                    axis.text.x = element_text(size = 15),
                                                                                                    axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15), 
                                                                                                    legend.key.size = unit(0.8, 'cm'), #change legend key size
                                                                                                    legend.key.height = unit(0.4, 'cm'), #change legend key height
                                                                                                    legend.key.width = unit(0.4, 'cm'), #change legend key width
                                                                                                    legend.title = element_text(size=9), #change legend title font size
                                                                                                    legend.text = element_text(size=9)) +  #change legend text font size)) +
  scale_x_discrete(labels=c("Migration Signature"))
ggsave("SuppFigures/Fibromigrationsig.svg")

#####Display our study on meta-analysis UMAP ####
fibro <- readRDS("finalfibroblastannotations.rds")
##write this as a for loop 
fibro$longaker <- "NA"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Mechanosensitive CTHRC1+"] <- "Mechanosensitive CTHRC1+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Mechanosensitive CTHRC1+"] <- "Mechanosensitive CTHRC1+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Mechanosensitive LRRC7+"] <- "Mechanosensitive LRRC7+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Mechanosensitive LRRC7+"] <- "Mechanosensitive LRRC7+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Antigen-Presenting CD74+"] <- "Antigen-Presenting CD74+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Antigen-Presenting CD74+"] <- "Antigen-Presenting CD74+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Inflammatory MMP3+"] <- "Inflammatory MMP3+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Inflammatory MMP3+"] <- "Inflammatory MMP3+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Inflammatory ADAMDEC1+"] <- "Inflammatory ADAMDEC1+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Inflammatory ADAMDEC1+"] <- "Inflammatory ADAMDEC1+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State KCNN3+"] <- "Steady-State KCNN3+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State KCNN3+"] <- "Steady-State KCNN3+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State GREM1+"] <- "Steady-State GREM1+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State GREM1+"] <- "Steady-State GREM1+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State F3+"] <- "Steady-State F3+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State F3+"] <- "Steady-State F3+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State FMO2+"] <- "Steady-State FMO2+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State FMO2+"] <- "Steady-State FMO2+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State FMO2+"] <- "Steady-State FMO2+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State FMO2+"] <- "Steady-State FMO2+"
fibro$longaker[fibro$study == "LongakerBowel" & fibro$fibro.annotation == "Steady-State PI16+"] <- "Steady-State PI16+"
fibro$longaker[fibro$study == "LongakerMes" & fibro$fibro.annotation == "Steady-State PI16+"] <- "Steady-State PI16+"
fibro$longaker[fibro$longaker == "NA"] <- "Meta-analysis"
#Reorder factor levels 
fibro$longaker <- factor(x = fibro$longaker, levels = rev(c("Steady-State PI16+", "Steady-State FMO2+", "Steady-State F3+","Steady-State GREM1+", "Steady-State KCNN3+", 
                                                            "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", "Antigen-Presenting CD74+", 
                                                            "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+", "Meta-analysis")))
Idents(fibro) <- "longaker"
table(fibro$diseasestate, fibro$tissue)
colors <- c("grey94", colors) #Add in light grey color for the meta-analysis data
DimPlot(fibro, cols = colors, order = T)
ggsave("Figures/Longakermaptometa.svg")

#### Fibroblast characteristics ####
##Markers

revfibro <- fibro
revfibro$fibro.annotation <- factor(revfibro$fibro.annotation, levels = c("Steady-State PI16+", "Steady-State FMO2+","Steady-State GREM1+", 
                                                                          "Steady-State KCNN3+", "Steady-State F3+", "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", 
                                                                          "Antigen-Presenting CD74+", "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+"))
Idents(revfibro) <- "fibro.annotation"

genes <- c("CTHRC1", "POSTN", "FAP", "SERPINE1", "SFRP4", 
           "LRRC7", "PIEZO2", "NLGN1", "SDK1", "ROBO2",
           "CD74", "HLA-DRA", "CD74", "CCL21", "TNFSF13B", 
           "MMP3", "MMP1", "CHI3L1", "IL32", "INHBA",
           "ADAMDEC1", "HAPLN1", "MYLK", "ADAM28", "CCL13",
           "F3", "BMP5", "EDNRB", "TRPA1", "FRZB", 
           "KCNN3", "THBS4", "LY6H", "EPHA3", "P2RY1", 
           "GREM1", "IGSF10", "MGP", "DEPP1", "RSPO3", 
           "FMO2", "MYC", "APOD", "MT2A", "IGF1", 
           "PI16", "MFAP5", "SEMA3C", "CD55", "CD248")

DotPlot(revfibro, features = genes, cols = "RdYlBu", scale.min = 0, scale.max = 100) +   theme(axis.title.x = element_blank(),
                                                                                               axis.title.y = element_blank(),
                                                                                               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                               legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                               legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                               legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                               legend.title = element_text(size=6), #change legend title font size
                                                                                               legend.text = element_text(size=6)) 
ggsave("SuppFigures/Fibrorep.genes.svg")




#### UMAP by study (fibroblasts) ####
fibro <- readRDS("fibrocellatlas.annotated.rds")
Idents(fibro) <- "study"
fibro$study[fibro$study == "LongakerBowel" ] <- "Stanford_Hospital"
fibro$study[fibro$study == "LongakerMes" ] <- "Stanford_Hospital"
fibro$study[fibro$study == "Rieder" ] <- "Cleveland_Clinic"
study_idents <- unique(fibro$study)
cell.list <- WhichCells(fibro, idents = study_idents, downsample = min(table(fibro$study)))
fibro_small <- fibro[, cell.list]
table(fibro_small$study)
Idents(fibro_small) <- "study"
colors = c("#e06666", "#ff80be", "#f6b26b", "#fbcd07","#ffeb84", "#93c47d", "#6fcc9f", 
"#76a5af", "#9aceeb", "#6fa8dc", "#8e7cc3", "#c27ba0")
DimPlot(fibro_small, cols = colors)
ggsave("SuppFigures/Fibrointegrationbystudy.svg")

####Integration QC Plots####
fibro[["percent.mt"]] <- PercentageFeatureSet(fibro, pattern = "^MT-")
# Create individual violin plots
plot1 <- VlnPlot(fibro, features = "nCount_RNA", pt.size = 0, group.by = "fibro.annotation", cols = colors) & coord_flip() & theme(legend.position = "none")
plot2 <- VlnPlot(fibro, features = "nFeature_RNA", pt.size = 0, group.by = "fibro.annotation", cols = colors) & coord_flip() & theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
plot3 <- VlnPlot(fibro, features = "percent.mt", pt.size = 0, group.by = "fibro.annotation", cols = colors) & coord_flip() & theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
# Combine the plots
combined_plot <- plot1 | plot2 | plot3
combined_plot
ggsave("SuppFigures/IntegrationQC.svg")

#### Disease State UMAPS ####
#Split by diseasestate
bowel <- subset(fibro, tissue == "SB" | tissue == "Colon")
bowel$diseasestate <- factor(bowel$diseasestate, levels = c("Normal Bowel", "Ileostomy Bowel", "NS Uninflamed", "NS Inflamed", 
                  "Stricture Uninflamed", "Stricture Inflamed", "Stricture", 
                  "Penatrating Uninflamed", "Penatrating Inflamed"))
DimPlot(bowel, split.by = "diseasestate", cols = colors) + theme(strip.text.x = element_text(size = 10, face = "bold") )
ggsave("SuppFigures/bowelfibro_diseasestate_split.svg", width = 21, height = 2.93)

mat <- subset(fibro, tissue == "MAT")
mat$diseasestate <- factor(mat$diseasestate, levels = c("Normal MAT", "Ileostomy MAT", "Uninvolved", "Involved"))
DimPlot(mat, split.by = "diseasestate", cols = colors) + theme(strip.text.x = element_text(size = 10, face = "bold") )
ggsave("SuppFigures/MATfibro_diseasestate_split.svg", width = 12, height = 2.93)

mat <- subset(fibro, tissue == "MAT")
DimPlot(mat, split.by = "diseasestate", cols = "RdBu") + theme(strip.text.x = element_text(size = 10, face = "bold") )
ggsave("SuppFigures/MATfibro_diseasestate_split.jpg", width = 9.5, height = 2.93)

#Generate proptable of CTHRC1 population by patient 
fibro_subset <- subset(fibro, subset = (diseasestate == c("Normal Bowel", "Stricture Uninflamed", "Stricture", 
                "Normal MAT", "Uninvolved", "Involved") & specimen_type == "Resection") )
tab1 = table(fibro_subset$fibro.annotation, fibro_subset$patientID)
ptab = prop.table(tab1, 1)*100
pdf1 <- as.data.frame(ptab)
colnames(pdf1) <- c("Fibroblast", "patientID", "Percent")
pdf1$diseasestate <- NA
metadata <- fibro_subset[[]]
for (ident in pdf1$patientID){
pdf1$diseasestate[pdf1[["patientID"]] == ident] <- toString(metadata$diseasestate[metadata[["patientID"]] == ident][1])
}
CTHRC1 <- pdf1[pdf1$Fibroblast == 'mCDF CTHRC1+', ]
diseasestate <- CTHRC1$diseasestate
Percent <- CTHRC1$percent
color = CTHRC1$PatientID
ggplot(CTHRC1, aes(x = diseasestate, y = Percent, fill=color)) + #Add jitter so points don't overlap
geom_point(alpha = 0.8, position = position_jitter(width = .2))

ggplot(data = CTHRC1, aes(x=diseasestate, y=Percent)) +
geom_boxplot(mapping = aes(fill=diseasestate))

#### Generate a controlled dataset ####
fibro_small <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
fibro_small <- subset(fibro_small, specimen_type == "Resection")
fibro_small <- subset(fibro_small, tissue == "Colon", invert = T)
fibro_small <- subset(fibro_small, diseasestate == "Stricture Inflamed", invert = T)
Idents(fibro_small) <- "fibro.annotation"
saveRDS(fibro_small, "fibrocontrolled_dataset.rds", compress = F)
#### Patient comparison ####
#Separate bowel and MAT
fibro_small <- readRDS("fibrocontrolled_dataset.rds")
mat <- subset(fibro_small, tissue == "MAT")
bowel <- subset(fibro_small, tissue == "SB")
table(mat$diseasestate, mat$age)
table(bowel$diseasestate, bowel$age)
#Find minimum number of cells for any condition
minimum = min(table(bowel$age)[table(bowel$age) > 0], 
table(mat$age)[table(mat$age) > 0])
Idents(bowel) <- "age"
Idents(mat) <- "age"
idents <- c("Adult", "Pediatric")
cell.list <- WhichCells(bowel, idents = idents, downsample = 2*minimum)
sub <- bowel[, cell.list]
table(sub$diseasestate, sub$age)

Idents(bowel) <- "fibro.annotation"
table(bowel$diseasestate, bowel$age)








Idents(bowel) <- "diseasestate"
Idents(mat) <- "diseasestate"
idents <- c("Normal Bowel", "Stricture Uninflamed", "Stricture")
cell.list <- WhichCells(bowel, idents = idents, downsample = min(table(bowel$diseasestate)[table(bowel$diseasestate) > 0]))
bowel <- bowel[, cell.list]
table(bowel$age)
Idents(bowel) <- "fibro.annotation"
table(bowel$diseasestate, bowel$age)

idents <- c("Normal MAT", "Uninvolved", "Invovled")
cell.list <- WhichCells(mat, idents = idents, downsample = min(table(mat$diseasestate)[table(mat$diseasestate) > 0]))
mat <- mat[, cell.list]
table(mat$age)
Idents(mat) <- "fibro.annotation"
table(mat$diseasestate, mat$age)

#Dimplot of age
DimPlot(fibro_small, split.by = "age", cols = "RdBu")
ggsave("SuppFigures/AdultPeddimplot.jpg", width = 7, height = 3.33 )

ped <- subset(fibro_small, age == "Pediatric")
adult <- subset(fibro_small, age == "Adult")

#Generate proptable and barplots of fibroblasts in resection 
tab = table(fibro_small$age, fibro_small$fibro.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
theme_classic() + 
labs(title="Fibroblast Enrichment", 
x="Age", y = "Percent (%)", fill = "Cluster") + 
scale_fill_brewer(palette = "RdBu") 
theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 
ggsave("SuppFigures/agebarplot.jpg", width = 3.73, height.= 3.69)

#Biopsy vs Resection in bowel fibroblasts to distinguish regional separation
fibro_small <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
fibro_small <- subset(fibro_small, general == "Bowel")
Idents(fibro_small) <- "specimen_type"
cell.list <- WhichCells(fibro_small, idents = c("Resection", "Biopsy"), downsample = min(table(fibro_small$specimen_type)))
fibro_small <- fibro_small[, cell.list]
table(fibro_small$specimen_type)
Idents(fibro_small) <- "fibro.annotation"
DimPlot(fibro_small, split.by = "specimen_type", cols = "RdBu")
ggsave("SuppFigures/biopsyvsresection.jpg", width = 8, height = 3.65)

#Drug treatments
fibro_small <- subset(fibro, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
fibro_small <- subset(fibro_small, specimen_type == "Resection")
#Idents(fibro_small) <- "treated"
#cell.list <- WhichCells(fibro_small, downsample = min(table(fibro_small$treated)))
#fibro_small <- fibro_small[, cell.list]
#table(fibro_small$treated)
Idents(fibro_small) <- "fibro.annotation"
dimplot <- DimPlot(fibro_small, split.by = "treated", cols = "RdBu")
#ggsave("SuppFigures/treatmentumaps.jpg", width = 13.5, height = 3.21)

#Generate proptable and barplots of fibroblasts in drug treatments 
treat <- subset(fibro, specimen_type == "Resection")
Idents(treat) <- "diseasestate"
treat <- subset(treat, idents = c("Normal Bowel", "Stricture", "Normal MAT", "Involved"))
# Change labels of normal to "None"
treat$treated[treat$diseasestate == "Normal Bowel"] <- "Non-CD"
treat$treated[treat$diseasestate == "Normal MAT"] <- "Non-CD"
treat <- subset(treat, treated == "NA", invert = T) #Remove cells where treatment unknown
treat$treated <- factor(treat$treated, levels = c("Non-CD", "IM", "Biologic", "IM + Biologic"))
tab = table(treat$treated, treat$fibro.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
theme_classic() + 
labs(title="Fibroblast Enrichment", 
x="Treatment", y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) + 
scale_fill_brewer(palette = "RdBu") +
theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) 
ggsave("SuppFigures/treatmentpropplot.jpg")

#### Adult vs Pediatric ####
# Define colors
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", 
"#cfe2f3", "#6fa8dc", "#0b5394","#073763")
names(colors) <- levels(fibro_small$fibro.annotation)

# Custom x-axis labels
x_labels <- c('Normal', 'Uninvolved', 'Stricture', 'Normal', 'Uninvolved', 'Stricture')

# Adult plot
adult <- subset(fibro_small, age == "Adult")
tab <- table(adult$diseasestate, adult$fibro.annotation)
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT", "MAT")
ptab <- prop.table(tab, 1) * 100
ptab <- as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue

adult_plot <- ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) +
geom_col() +
theme_classic() +
labs(title = "Adult",
y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) +
theme(axis.text.x = element_text(size = 12, colour = "black"),
axis.text.y = element_text(size = 12, colour = "black"),
axis.title.x = element_blank(),
strip.placement = "outside", 
strip.background = element_blank(), axis.title.y = element_text(size = 15),
strip.text = element_text(size = 15, hjust = 0.5)) +
facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") +
scale_fill_manual(values = colors) +
scale_x_discrete(labels = x_labels)

# Pediatric plot
ped <- subset(fibro_small, age == "Pediatric")
tab <- table(ped$diseasestate, ped$fibro.annotation)
tissue <- c("Bowel", "Bowel", "MAT", "MAT")
ptab <- prop.table(tab, 1) * 100
ptab <- as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue

ped_plot <- ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) +
geom_col() +
theme_classic() +
labs(title = "Pediatric",
y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) +
theme(axis.text.x = element_text(size = 12, colour = "black"),
axis.text.y = element_text(size = 12, colour = "black"),
axis.title.x = element_blank(),
strip.placement = "outside", 
strip.background = element_blank(), axis.title.y = element_text(size = 15),
strip.text = element_text(size = 15, hjust = 0.5)) +
facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") +
scale_fill_manual(values = colors) +
scale_x_discrete(labels = c("Uninvolved", "Stricture"))

# Combine plots
library(patchwork)
combined_plot <- adult_plot + ped_plot + plot_layout(guides = "collect", widths = c(1.5, 1))

combined_plot
ggsave("SuppFigures/agepropplot.svg") 

#### Sex differences ####
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
names(colors) <- levels(fibro_small$fibro.annotation)

# Custom x-axis labels
x_labels <- c('Normal', 'Uninvolved', 'Stricture', 'Normal', 'Uninvolved', 'Stricture')

# Male plot
male <- subset(fibro_small, sex == "Male")
tab <- table(male$diseasestate, male$fibro.annotation)
tissue <- c("Bowel", "Bowel", "Bowel", "MAT", "MAT")
ptab <- prop.table(tab, 1) * 100
ptab <- as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue

male_plot <- ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) +
geom_col() +
theme_classic() +
labs(title = "Male",
y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) +
theme(axis.text.x = element_text(size = 12, colour = "black"),
axis.text.y = element_text(size = 12, colour = "black"),
axis.title.x = element_blank(),
strip.placement = "outside", 
strip.background = element_blank(), axis.title.y = element_text(size = 15),
strip.text = element_text(size = 15, hjust = 0.5)) +
facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") +
scale_fill_manual(values = colors) +
scale_x_discrete(labels = x_labels)

# Female plot
female <- subset(fibro_small, sex == "Female")
tab <- table(female$diseasestate, female$fibro.annotation)
tissue <- c("Bowel", "Bowel", "MAT", "MAT", "MAT")
ptab <- prop.table(tab, 1) * 100
ptab <- as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ptab$tissue <- tissue

female_plot <- ggplot(ptab, aes(x = Var1, y = Freq, fill = Var2)) +
geom_col() +
theme_classic() +
labs(title = "Female",
y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) +
theme(axis.text.x = element_text(size = 12, colour = "black"),
axis.text.y = element_text(size = 12, colour = "black"),
axis.title.x = element_blank(),
strip.placement = "outside", 
strip.background = element_blank(), axis.title.y = element_text(size = 15),
strip.text = element_text(size = 15, hjust = 0.5)) +
facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") +
scale_fill_manual(values = colors) +
scale_x_discrete(labels = x_labels)

# Combine plots
library(patchwork)
combined_plot <- male_plot + female_plot + plot_layout(guides = "collect", widths = c(1, 1))

combined_plot
ggsave("SuppFigures/sexpropplot.svg") 




#### Drug treatments ####
colors <- c("#660000", "#cc0000", "#8e7cc3", "#b45f06", "#e69138", "#6aa84f", "#cfe2f3", "#6fa8dc", "#0b5394","#073763")
names(colors) <- levels(fibro_small$fibro.annotation)
#Generate proptable and barplots of fibroblasts in drug treatments 
treat <- fibro_small
Idents(treat) <- "diseasestate"
treat <- subset(treat, idents = c("Normal Bowel", "Stricture", "Normal MAT", "Involved"))
# Change labels of normal to "None"
treat$treated[treat$diseasestate == "Normal Bowel"] <- "Non-CD"
treat$treated[treat$diseasestate == "Normal MAT"] <- "Non-CD"
treat <- subset(treat, treated == "NA", invert = T) #Remove cells where treatment unknown
treat$treated <- factor(treat$treated, levels = c("Non-CD", "IM", "Biologic", "IM + Biologic"))
tab = table(treat$treated, treat$fibro.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
theme_classic() + 
labs(title="Fibroblast Enrichment", 
x="Treatment", y = "Percent (%)", fill = "Cluster", element_text(colour = "black", face = "bold")) + 
scale_fill_manual(values = colors) + 
theme(axis.text.x = element_text(size = 12, colour = "black"), 
axis.text.y = element_text(size = 12, colour = "black"), 
axis.title.x = element_text(size = 15, colour = "black"),
axis.title.y = element_text(size = 15, colour = "black")) 
ggsave("SuppFigures/treatmentpropplot.svg")

#### Subanalyses ####

#Fibroblast subanalyses for bowel
rm(list = ls())
fibro <- readRDS("finalfibroblastannotations.rds")
bowel <- subset(fibro, tissue == "SB" | tissue == "Colon")
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
bowel <- CreateSeuratObject(counts = bowel[["RNA"]]$counts, meta.data = bowel@meta.data)
saveRDS( bowel, 'Bowel/fibro.rds', compress = F )

# re-perform BPcells
counts.mat = bowel[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "Bowel/fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "Bowel/fibro_counts" )
bowel = CreateSeuratObject( counts.mat, meta.data = bowel@meta.data )
saveRDS( bowel, 'Bowel/fibro_counts.rds', compress = F )

# run downstream analysis
bowel <- NormalizeData(bowel)
bowel[["RNA"]] <- split(bowel[["RNA"]], f = bowel$orig.ident)
saveRDS( bowel, 'Bowel/fibro_split.rds', compress = F )

bowel<- readRDS('Bowel/fibro_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(bowel)) 
for (i in 1:length) {
# Access the RNA data
tube <- bowel[['RNA']][[Layers(bowel)[i]]]
# Try executing FindVariableFeatures
tryCatch({
# Find variable features
print(Layers(bowel)[i])
tube <- CreateSeuratObject(tube)
tube <- NormalizeData(tube)
tube <- FindVariableFeatures(tube)
}, error = function(e) {
# Record error message
error_messages[[i]] <<- conditionMessage(e)
# Record index of the layer where the error occurred
error_indices[[i]] <<- i
# Create a row for the error dataframe
error_row <- data.frame(layer_index = i, ident = as.character(bowel$orig.ident[i]), layer_name = Layers(bowel)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
# Append the error row to the error dataframe
error_df <<- rbind(error_df, error_row)
})
}
print(error_df)
write.csv(error_df, "Bowel/Error_list.csv")

error_df <- read.csv("Bowel/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(bowel$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(bowel$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- bowel
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'Bowel/fibro_split_subset.rds', compress = F)

bowel <- readRDS('Bowel/fibro_split_subset.rds')
bowel <- FindVariableFeatures(bowel, verbose = T )
bowel <- ScaleData(bowel, verbose = T)
bowel <- RunPCA(bowel, npcs = 30, verbose = T)
saveRDS( bowel, 'Bowel/fibro_pca.rds'  , compress = F )
bowel <- IntegrateLayers(bowel, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( bowel, 'Bowel/fibro_harmony.rds', compress = F ) 
bowel <- RunUMAP(bowel, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( bowel, 'Bowel/fibro_harmony_umap.rds', compress = F ) 
bowel <- JoinLayers(bowel)
saveRDS( bowel, 'Bowel/fibro_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(bowel, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("Bowel/fibro_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

bowel <- FindNeighbors(bowel, reduction = "harmony", dims = 1:30)
for( res in seq(0.01, 0.1, 0.05) ) {
bowel <- FindClusters(bowel, resolution = res )
DimPlot(bowel, raster = TRUE, label = T)
ggsave(paste0("Bowel/fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
saveRDS( bowel, paste0( 'Bowel/fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
bowel <- readRDS("Bowel/fibro_harmony_clusters_0.31.rds")
bowel[["RNA"]]$data <- as(object = bowel[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(bowel, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "Bowel/fibrofibro_harmony_clusters_0.31_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
group_by(cluster) %>%
slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "Bowel/fibrofibro_harmony_clusters_0.31_topmarkers.csv")

# manual annotations for entire dataset
Idents(bowel) <- "seurat_clusters"
bowel$bowel.annotation <- "tbd"
bowel$bowel.annotation[bowel$seurat_clusters == 0] <- "Inflammatory ADAMDEC1+" #LP
bowel$bowel.annotation[bowel$seurat_clusters == 1] <- "Steady-State PI16+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 2] <- "Telocytes F3+" #Villus
bowel$bowel.annotation[bowel$seurat_clusters == 3] <- "Steady-State GREM1+" #Crypt
bowel$bowel.annotation[bowel$seurat_clusters == 4] <- "Mechanosensitive LRRC7+" #Submucosa
bowel$bowel.annotation[bowel$seurat_clusters == 5] <- "Steady-State KCNN3+"#Villus
bowel$bowel.annotation[bowel$seurat_clusters == 6] <- "Inflammatory ADAMDEC1+" #LP
bowel$bowel.annotation[bowel$seurat_clusters == 7] <- "Mechanosensitive CTHRC1+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 8] <- "Antigen-Presenting CD74+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 9] <- "Inflammatory ADAMDEC1+" #LP
bowel$bowel.annotation[bowel$seurat_clusters == 10] <- "Inflammatory MMP3+"#LP
bowel$bowel.annotation[bowel$seurat_clusters == 11] <- "Steady-State PI16+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 12] <- "Steady-State PI16+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 13] <- "Inflammatory ADAMDEC1+" #Adventitial
bowel$bowel.annotation[bowel$seurat_clusters == 14] <- "Inflammatory ADAMDEC1+" #LP

bowel$bowel.annotation <- factor(bowel$bowel.annotation, levels = c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+", 
                                               "Inflammatory MMP3+", "Inflammatory ADAMDEC1+", "Antigen-Presenting CD74+", 
                                               "Telocytes F3+", "Steady-State KCNN3+", "Steady-State GREM1+", "Steady-State PI16+"))

Idents(bowel) <- "bowel.annotation"
colors <- c("#660000", "#cc0000", "#b45f06", "#e69138", "#8e7cc3", "#6aa84f", 
"#cfe2f3", "#6fa8dc", "#0b5394","#073763")
dimplot <- DimPlot(bowel, cols = colors) 
revbowel <- bowel 
levels(revbowel) <- rev(c("Mechanosensitive CTHRC1+", "Mechanosensitive LRRC7+", 
                          "Inflammatory MMP3+", "Inflammatory ADAMDEC1+", "Antigen-Presenting CD74+", 
                          "Telocytes F3+", "Steady-State KCNN3+", "Steady-State GREM1+", "Steady-State PI16+"))
dotplot <- DotPlot(revbowel, features = c("CTHRC1", "POSTN", "LRRC7", "NLGN1", "MMP3", "CHI3L1",  "ADAMDEC1", "HAPLN1",
                                          "CD74", "HLA-DRB1", "F3", "BMP5", "KCNN3", "C7", "GREM1",  "C3", "PI16", "MFAP5" ), 
                   cols = "RdYlBu", col.min = 0, col.max = 100, dot.scale = 3) +  
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) #change legend text font size)

combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -30, unit = "pt")))) + dotplot 
combined_plot

ggsave("Bowel/Figures/bowelfibroumap.svg")
saveRDS(bowel, "Bowel/bowelfibroblasts.annotated.rds", compress = F)

Idents(bowel) <- "diseasestate"
bsub <- subset(bowel, idents = c("Normal Bowel", "Stricture Uninflamed", "Stricture"))
Idents(bsub) <- "bowel.annotation"
DimPlot(bsub, split.by = "diseasestate", cols = colors)
ggsave("Bowel/Figures/bowelfibroblast_umap_sdiseasestate.svg")

#Organize fibroblasts by location 
# manual annotations for entire dataset
Idents(bowel) <- "bowel.annotation"
bowel$location <- "tbd"
bowel$location[bowel$bowel.annotation == "iCDF ADAMDEC1+"] <- "Mucosa" #LP
bowel$location[bowel$bowel.annotation == "ssCDF PI16+"] <- "SubSer" #Adventitial
bowel$location[bowel$bowel.annotation == "Telocytes F3+"] <- "Mucosa" #Villus
bowel$location[bowel$bowel.annotation == "ssCDF GREM1+"] <- "SubSer" #Crypt
bowel$location[bowel$bowel.annotation == "msCDF LRRC7+"] <- "SubSer" #Submucosa
bowel$location[bowel$bowel.annotation == "ssCDF KCNN3+"] <- "SubSer"  #Villus
bowel$location[bowel$bowel.annotation == "msCDF CTHRC1+"] <-  "SubSer" #LP
bowel$location[bowel$bowel.annotation == "apCDF CD74+"] <- "LP" #Adventitial
bowel$location[bowel$bowel.annotation == "iCDF MMP3+"] <- "LP" #Adventitial
Idents(bowel) <- "location"
DimPlot(bowel)
library(readxl)
Yapgenes <-as.data.frame(read_excel("AltYAP_list2.xlsx",col_names = "Genes"))
Yapgenes <- Yapgenes[Yapgenes$Genes %in% rownames(bowel), ,drop = FALSE] #Remove genes not present in our dataset
bowel <- AddModuleScore(bowel, features = Yapgenes, name = "YAP")
DotPlot(bowel, features = "YAP1")

ecmgenes <- list(c("COL1A1", "COL1A2", "COL3A1", "COL4A2", "COL5A1", "COL5A3", 
                   "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", 
                   "LAMB1", "LAMB2", "LAMB3", 
                   "LAMC1", "LAMC2", "LTBP1", 
                   "MFAP2", "MFAP5", "FN1", 
                   "POSTN", "SERPINA1", "SERPINA5", 
                   "TGFBI", "TNC", "VCAN"))

bowel <- AddModuleScore(bowel, features = ecmgenes, name = "ECM")

components <- list(c("YAP1", "CTNNB1", "MST1",
                     "PTK2", "ITGA5", "ITGAV", "ITGB1",  
                     "SAV1", "STK3", "LATS1", "LATS2", 
                     "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", 
                     "VIM", "VCL", "VASP"))

bowel <- AddModuleScore(bowel, features = ecmgenes, name = "components")

Idents(bowel) <- "bowel.annotation"
DotPlot(bowel, features = c("components1", "ECM1", "YAP1"), scale.min = 0, scale.max = 100, split.by = "specimen_type")

#Compare conditions
Idents(bowel) <- "specimen_type"
cell.list <- WhichCells(bowel, idents = c("Biopsy", "Resection"), downsample = min(table(bowel$specimen_type)))
bowel_sub <- bowel[, cell.list]
table(bowel_sub$specimen_type)
Idents(bowel_sub) <- "bowel.annotation"
DimPlot(bowel_sub)


Idents(bowel) <- "diseasestate"
bowel_sub <- subset(bowel, idents = c("Normal Bowel", "Stricture Uninflamed", "Stricture"))
tab = table(bowel_sub$diseasestate, bowel_sub$bowel.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(y = "Percent (%)", fill = "Cluster")+ 
  scale_fill_manual(values = colors) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) + 
  scale_x_discrete(labels=c('Normal', 'Uninvolved', 'Stricture'))
ggsave("Bowel/Figures/Bowelpropplot.svg")



#Now for MAT 
#Fibroblast subanalyses for mat
rm(list = ls())
fibro <- readRDS("finalfibroblastannotations.rds")
mat <- subset(fibro, tissue == "MAT")
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
mat <- CreateSeuratObject(counts = mat[["RNA"]]$counts, meta.data = mat@meta.data)
saveRDS( mat, 'MAT/fibro.rds', compress = F )

# re-perform BPcells
counts.mat = mat[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "MAT/fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MAT/fibro_counts" )
mat = CreateSeuratObject( counts.mat, meta.data = mat@meta.data )
saveRDS( mat, 'MAT/fibro_counts.rds', compress = F )

# run downstream analysis
mat <- NormalizeData(mat)
mat[["RNA"]] <- split(mat[["RNA"]], f = mat$orig.ident)
saveRDS( mat, 'MAT/fibro_split.rds', compress = F )

mat<- readRDS('MAT/fibro_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(mat)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- mat[['RNA']][[Layers(mat)[i]]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(mat)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(mat$orig.ident[i]), layer_name = Layers(mat)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "MAT/Error_list.csv")

error_df <- read.csv("MAT/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(mat$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(mat$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- mat
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'MAT/fibro_split_subset.rds', compress = F)

mat <- readRDS('MAT/fibro_split_subset.rds')
mat <- FindVariableFeatures(mat, verbose = T )
mat <- ScaleData(mat, verbose = T)
mat <- RunPCA(mat, npcs = 30, verbose = T)
saveRDS( mat, 'MAT/fibro_pca.rds'  , compress = F )
mat <- IntegrateLayers(mat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( mat, 'MAT/fibro_harmony.rds', compress = F ) 
mat <- RunUMAP(mat, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( mat, 'MAT/fibro_harmony_umap.rds', compress = F ) 
mat <- JoinLayers(mat)
saveRDS( mat, 'MAT/fibro_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(mat, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MAT/fibro_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

mat <- FindNeighbors(mat, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  mat <- FindClusters(mat, resolution = res )
  DimPlot(mat, raster = TRUE, label = T)
  ggsave(paste0("MAT/fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( mat, paste0( 'MAT/fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
mat[["RNA"]]$data <- as(object = mat[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(mat, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "MAT/fibrofibro_harmony_clusters_0.5_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "MAT/fibrofibro_harmony_clusters_0.5_topmarkers.csv")

#Now for MAT 
#Fibroblast subanalyses for mat
rm(list = ls())
fibro <- readRDS("finalfibroblastannotations.rds")
mat <- subset(fibro, tissue == "MAT")
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
mat <- CreateSeuratObject(counts = mat[["RNA"]]$counts, meta.data = mat@meta.data)
saveRDS( mat, 'MAT/fibro.rds', compress = F )

# re-perform BPcells
counts.mat = mat[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "MAT/fibro_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MAT/fibro_counts" )
mat = CreateSeuratObject( counts.mat, meta.data = mat@meta.data )
saveRDS( mat, 'MAT/fibro_counts.rds', compress = F )

# run downstream analysis
mat <- NormalizeData(mat)
mat[["RNA"]] <- split(mat[["RNA"]], f = mat$orig.ident)
saveRDS( mat, 'MAT/fibro_split.rds', compress = F )

mat<- readRDS('MAT/fibro_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(mat)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- mat[['RNA']][[Layers(mat)[i]]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(mat)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(mat$orig.ident[i]), layer_name = Layers(mat)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "MAT/Error_list.csv")

error_df <- read.csv("MAT/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(mat$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(mat$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- mat
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'MAT/fibro_split_subset.rds', compress = F)

mat <- readRDS('MAT/fibro_split_subset.rds')
mat <- FindVariableFeatures(mat, verbose = T )
mat <- ScaleData(mat, verbose = T)
mat <- RunPCA(mat, npcs = 30, verbose = T)
saveRDS( mat, 'MAT/fibro_pca.rds'  , compress = F )
mat <- IntegrateLayers(mat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( mat, 'MAT/fibro_harmony.rds', compress = F ) 
mat <- RunUMAP(mat, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( mat, 'MAT/fibro_harmony_umap.rds', compress = F ) 
mat <- JoinLayers(mat)
saveRDS( mat, 'MAT/fibro_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(mat, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MAT/fibro_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

mat <- FindNeighbors(mat, reduction = "harmony", dims = 1:30)
for( res in seq(0.1, 1, 0.1) ) {
  mat <- FindClusters(mat, resolution = res )
  DimPlot(mat, raster = TRUE, label = T)
  ggsave(paste0("MAT/fibro_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( mat, paste0( 'MAT/fibro_harmony_clusters_', res, '.rds' ), compress = F ) 
}

# join and find markers
mat[["RNA"]]$data <- as(object = mat[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(mat, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "MAT/fibrofibro_harmony_clusters_0.5_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "MAT/fibrofibro_harmony_clusters_0.5_topmarkers.csv")

#Remove non-fibroblast populations based on markers 
#Fibroblast subanalyses for mat
rm(list = ls())
mat <- readRDS("MAT/newFibro1/fibro_harmony_clusters_0.4.rds")
mat <- subset(mat, idents = c(7, 9, 10,11), invert = T)
DimPlot(mat)
# re-generate seurat object and remove extraneous data
options(Seurat.object.assay.version = "v5")
mat <- CreateSeuratObject(counts = mat[["RNA"]]$counts, meta.data = mat@meta.data)
saveRDS( mat, 'MAT/newFibro2/fibro2.rds', compress = F )

# re-perform BPcells
counts.mat = mat[['RNA']]$counts
counts.mat = convert_matrix_type( counts.mat, type="uint32_t" )
counts.mat = as( counts.mat, Class='dgCMatrix' )
write_matrix_dir(mat = counts.mat, dir = "MAT/newFibro2/fibro2_counts", overwrite = T )
counts.mat <- open_matrix_dir(     dir = "MAT/newFibro2/fibro2_counts" )
mat = CreateSeuratObject( counts.mat, meta.data = mat@meta.data )
saveRDS( mat, 'MAT/newFibro2/fibro2_counts.rds', compress = F )

# run downstream analysis
mat <- NormalizeData(mat)
mat[["RNA"]] <- split(mat[["RNA"]], f = mat$orig.ident)
saveRDS( mat, 'MAT/newFibro2/fibro2_split.rds', compress = F )

mat<- readRDS('MAT/newFibro2/fibro2_split.rds')
# Create an empty list to store error messages
error_messages <- list()
# Create an empty list to store the layer indices where errors occurred
error_indices <- list()
# Create an empty dataframe to store information about errors
error_df <- data.frame(layer_index = integer(), ident = character(), layer_name = character(), error_message = character(), stringsAsFactors = FALSE)

# Loop through layers
length = length(Layers(mat)) 
for (i in 1:length) {
  # Access the RNA data
  tube <- mat[['RNA']][[Layers(mat)[i]]]
  # Try executing FindVariableFeatures
  tryCatch({
    # Find variable features
    print(Layers(mat)[i])
    tube <- CreateSeuratObject(tube)
    tube <- NormalizeData(tube)
    tube <- FindVariableFeatures(tube)
  }, error = function(e) {
    # Record error message
    error_messages[[i]] <<- conditionMessage(e)
    # Record index of the layer where the error occurred
    error_indices[[i]] <<- i
    # Create a row for the error dataframe
    error_row <- data.frame(layer_index = i, ident = as.character(mat$orig.ident[i]), layer_name = Layers(mat)[i], error_message = conditionMessage(e), stringsAsFactors = FALSE)
    # Append the error row to the error dataframe
    error_df <<- rbind(error_df, error_row)
  })
}
print(error_df)
write.csv(error_df, "MAT/newFibro2/Error_list.csv")

error_df <- read.csv("MAT/newFibro2/Error_list.csv")
idents <- error_df$layer_name
remove <- sub("^[^.]+\\.", "", idents)
tab <- as.data.frame(table(mat$orig.ident))
tab$Names <- as.character(tab$Var1)
tab <- tab[tab$Names %in% remove, ] #subset the bad studies 
sum(tab$Freq) #Number of cells that will be omitted
print(sum(tab$Freq)/sum(table(mat$orig.ident))) #Percentage of cells removed from study
###STOPP
#Remove the objects: 
obj_subset <- mat
obj_subset$keep <- T
obj_subset$keep[obj_subset$orig.ident %in% remove] <- F
obj_subset <- subset(obj_subset, subset = keep == T)

saveRDS(obj_subset, 'MAT/newFibro2/fibro2_split_subset.rds', compress = F)

mat <- readRDS('MAT/newFibro2/fibro2_split_subset.rds')
mat <- FindVariableFeatures(mat, verbose = T )
mat <- ScaleData(mat, verbose = T)
mat <- RunPCA(mat, npcs = 30, verbose = T)
saveRDS( mat, 'MAT/newFibro2/fibro2_pca.rds'  , compress = F )
mat <- IntegrateLayers(mat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony" )
saveRDS( mat, 'MAT/newFibro2/fibro2_harmony.rds', compress = F ) 
mat <- RunUMAP(mat, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = T, min.dist = 0.001)
saveRDS( mat, 'MAT/newFibro2/fibro2_harmony_umap.rds', compress = F ) 
mat <- JoinLayers(mat)
saveRDS( mat, 'MAT/newFibro2/fibro2_harmony_umap_joined.rds', compress = F ) 

FeaturePlot(mat, features = c("PECAM1", "MKI67", "PTPRC", "PLP1","IGF2", "DACH1", "PDGFRA", "DPT", "EPCAM", "KRT19", "KCNJ8", "RGS5", "MYH11", "LRAT"), raster = TRUE)
ggsave("MAT/newFibro2/fibro2_harmony_umap_feature.jpg", width = 8, height = 10, units = "in", limitsize = FALSE)

mat <- FindNeighbors(mat, reduction = "harmony", dims = 1:30)
for( res in seq(0.2, 0.25, 0.01) ) {
  mat <- FindClusters(mat, resolution = res )
  DimPlot(mat, raster = TRUE, label = T)
  ggsave(paste0("MAT/newFibro2/fibro2_harmony_clusters_",res,".jpg"), width = 5, height = 5, units = "in", limitsize = FALSE)
  saveRDS( mat, paste0( 'MAT/newFibro2/fibro2_harmony_clusters_', res, '.rds' ), compress = F ) 
}

mat <- readRDS("MAT/newFibro2/fibro2_harmony_clusters_0.25.rds")
# join and find markers
mat[["RNA"]]$data <- as(object = mat[["RNA"]]$data, Class = "dgCMatrix")
markers <- FindAllMarkers(mat, max.cells.per.ident = 1000, only.pos = T)
write.csv(markers, "MAT/newFibro2/MATfibro2_harmony_clusters_0.25_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "MAT/newFibro2/MATfibro2_harmony_clusters_0.25_topmarkers.csv")

# manual annotations for entire dataset
Idents(mat) <- "seurat_clusters"
mat$mat.annotation <- "tbd"
mat$mat.annotation[mat$seurat_clusters == 0] <- "Steady-State APC FMO2+"
mat$mat.annotation[mat$seurat_clusters == 1] <- "Steady-State APC DPP4+"
mat$mat.annotation[mat$seurat_clusters == 2] <- "Steady-State APC ICAM1+"
mat$mat.annotation[mat$seurat_clusters == 3] <- "Inflammatory APC EBF2+"
mat$mat.annotation[mat$seurat_clusters == 4] <- "Mechanosenitive FAP CTHRC1+"
mat$mat.annotation[mat$seurat_clusters == 5] <- "Steady-State APC F3+"
mat$mat.annotation[mat$seurat_clusters == 6] <- "Mechanosensitive FAP LRRC7+"
mat$mat.annotation[mat$seurat_clusters == 7] <- "Inflammatory FAP ADAMDEC1+"

levels <- c("Mechanosenitive FAP CTHRC1+", "Mechanosensitive FAP LRRC7+", "Steady-State APC EBF2+",
            "Steady-State APC F3+", "Steady-State APC ICAM1+", "Steady-State APC FMO2+", "Steady-State APC DPP4+", 
            "Inflammatory FAP ADAMDEC1+")
mat$mat.annotation <- factor(x = mat$mat.annotation, levels = levels)
Idents(mat) <- "mat.annotation"
colors <- c("#660000", "#cc0000", "#6aa84f", "#cfe2f3", "#6fa8dc", "#0b5394","#073763", "#e69138")
dimplot <- DimPlot(mat, cols = colors)
revmat <- mat
levels(revmat) <- rev(c("Mechanosenitive FAP CTHRC1+", "Mechanosensitive FAP LRRC7+", "Steady-State APC EBF2+",
                        "Steady-State APC F3+", "Steady-State APC ICAM1+", "Steady-State APC FMO2+", "Steady-State APC DPP4+", 
                        "Inflammatory FAP ADAMDEC1+"))
dotplot <- DotPlot(revmat, features = c("CTHRC1", "POSTN", "LRRC7", "NLGN1", "EBF2","ABCA9-AS1", "F3", "VIT", "ICAM1", "CXCL2", 
                                        "FMO2", "DEPP1", "DPP4", "PI16", "ADAMDEC1", "HAPLN1"), cols = "RdYlBu", col.min = 0, col.max = 100, dot.scale = 2, scale.max = 100) +  theme(axis.title.x = element_blank(),
                                                                                                                                                                                      axis.title.y = element_blank(), axis.text.y = element_blank(),
                                                                                                                                                                                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                                                                                                                      legend.key.size = unit(0.4, 'cm'), #change legend key size
                                                                                                                                                                                      legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                                                                                                                                      legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                                                                                                                                      legend.title = element_text(size=6), #change legend title font size
                                                                                                                                                                                      legend.text = element_text(size=6)) 

combined_plot <- (dimplot + theme(axis.title.x = element_text(margin = margin(t = -40, unit = "pt")))) + dotplot 
combined_plot
ggsave("MAT/Figures/matfibro.umap.svg")

#Generate proptable and barplots of fibroblasts in resection 
mat_sub <- subset(mat, diseasestate == "Ileostomy MAT", invert = T)
tab = table(mat_sub$diseasestate, mat_sub$mat.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
ptab <- ptab[!is.na(ptab$Freq), ]
ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(y = "Percent (%)", fill = "Cluster")+ 
  scale_fill_manual(values = colors) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 
ggsave("MAT/Figures/MATpropplot.svg")


Idents(mat_sub) <- "diseasestate"
idents <- c("Normal MAT", "Uninvolved", "Involved")
cell.list <- WhichCells(mat_sub, idents = ident, downsample = min(table(mat_sub$diseasestate)))
mat_sub <- mat_sub[, cell.list]
table(mat_sub$diseasestate)
DimPlot(mat_sub, cols = colors, split.by = "diseasestate")
ggsave("MAT/Figures/matfibro.umapdiseasestate.svg")

saveRDS(mat, "MAT/matfibroblasts.annotated.rds", compress = F)


#Find markers between bowel and MAT fibroblasts 
markers <- FindMarkers(fibro, ident.1 = "MAT", ident.2 = "Bowel", assay = "RNA")
write.csv(markers, "BowelvsFibro_markers.csv")
markers$rank <- markers$pct.1 / (markers$pct.2 + 0.001) * markers$avg_log2FC
topmarkers <- markers %>%
  slice_max(n = 100, order_by = rank)
write.csv(topmarkers, "BowelvsFibro_topmarkers.csv")
Idents(fibro) <- "fibro.annotation"

p <- DotPlot(fibro, features = c("EBF2", "CLU", "FMO2", "SFRP4", "PDGFRL", "ABCA10")) 
data <- as.data.frame(p$data)
ggplot(p$data, aes(x =id, y = pct.exp, fill = "ave.exp.scaled")) + 
  geom_dotplot() + 
  facet_wrap(~unique(fibro$general))


fibro$alt.annotation <- "tbd"
fibro$alt.annotation[fibro$seurat_clusters == 0] <- "Steady-State PI16+"
fibro$alt.annotation[fibro$seurat_clusters == 1] <- "Inflammatory ADAMDEC1+"
fibro$alt.annotation[fibro$seurat_clusters == 2] <- "Steady-State FMO2+"
fibro$alt.annotation[fibro$seurat_clusters == 3] <- "Steady-State F3+"
fibro$alt.annotation[fibro$seurat_clusters == 4] <- "Mechanosensitive LRRC7+"
fibro$alt.annotation[fibro$seurat_clusters == 5] <- "Mechanosensitive CTHRC1+"
fibro$alt.annotation[fibro$seurat_clusters == 6] <- "Steady-State KCNN3+"
fibro$alt.annotation[fibro$seurat_clusters == 7] <- "Steady-State GREM1+"
fibro$alt.annotation[fibro$seurat_clusters == 8] <- "Antigen-Presenting CD74+"
fibro$alt.annotation[fibro$seurat_clusters == 9] <- "Inflammatory MMP3+"
fibro$alt.annotation[fibro$seurat_clusters == 10] <- "Steady-State KCNN3+"
fibro$alt.annotation[fibro$seurat_clusters == 11] <- "Steady-State PI16+"
Idents(fibro) <- "alt.annotation"

fibro$alt.annotation <- factor(x = fibro$alt.annotation, levels = c("Steady-State PI16+", "Steady-State FMO2+", "Steady-State F3+","Steady-State GREM1+", "Steady-State KCNN3+", 
                                                                    "Inflammatory ADAMDEC1+", "Inflammatory MMP3+", "Antigen-Presenting CD74+", 
                                                                    "Mechanosensitive LRRC7+", "Mechanosensitive CTHRC1+"))
Idents(fibro) <- "alt.annotation"
palette <- c("#b2182b", "#ef8a62", "#fddbc7", "grey92", "#d1e5f0", "#67a9cf", "#2166ac")

DotPlot(fibro, features = c("EBF2", "CLU", "FMO2", "SFRP4", "PDGFRL", "ABCA10"),
        scale.min = 0, scale.max = 100, cols = "RdYlBu", split.by = "general") + theme(axis.title.x = element_blank(),
                                                                                       axis.title.y = element_blank(), 
                                                                                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                                                                       axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 10), 
                                                                                       legend.key.size = unit(0.8, 'cm'), #change legend key size
                                                                                       legend.key.height = unit(0.4, 'cm'), #change legend key height
                                                                                       legend.key.width = unit(0.4, 'cm'), #change legend key width
                                                                                       legend.title = element_text(size=9), #change legend title font size
                                                                                       legend.text = element_text(size=9)) 


#Immune Population Changes
global_small <- subset(global, diseasestate == "Ileostomy MAT" | diseasestate == "Ileostomy Bowel", invert = T)
global_small <- subset(global_small, specimen_type == "Resection")
Idents(global_small) <- "manual.annotation"
myeloid <- subset(global_small, idents = c("Myeloid Cells", "Macrophages", "Mast Cells"))
lymphoid <- subset(global_small, idents = c("B Cells", "Plasma Cells", "T Cells"))
immune <- subset(global_small, idents = c("Myeloid Cells", "Macrophages", "Mast Cells", "B Cells", "Plasma Cells", "T Cells"))
#Generate proptable and barplots of fibroblasts in resection 
tab = table(myeloid$diseasestate, myeloid$manual.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
p1 <- ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Enrichment", 
       x="Age", y = "Percent (%)", fill = "Cluster") + 
  # scale_fill_brewer(palette = "RdBu") 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 

#Generate proptable and barplots of fibroblasts in resection 
tab = table(lymphoid$diseasestate, lymphoid$manual.annotation)
ptab = prop.table(tab, 1)*100
ptab = as.data.frame(ptab)
p2 <- ggplot(ptab,aes(x=Var1,y=Freq,fill=Var2)) + geom_col() +
  theme_classic() + 
  labs(title="Fibroblast Enrichment", 
       x="Age", y = "Percent (%)", fill = "Cluster") + 
  # scale_fill_brewer(palette = "RdBu") 
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) 






Idents(global_small) <- "diseasestate"
idents <- c("Normal", "Uninvolved", "Involved")
cell.list <- WhichCells(global_small, idents = idents, downsample = min(table(global_small$diseasestate)))
global_small <- global_small[, cell.list]
table(global_small$age)
