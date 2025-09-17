#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:16:15 2025

@author: kebr
"""
/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/python3.11

import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix
import os

#set wd 
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Inputs/")

os.chdir("/oak/stanford/groups/longaker/BPP/BA_Xenium/")

# 1️ Load the Counts Matrix (MTX format)
counts = mmread("matrix.mtx").tocsc()  # Load as sparse matrix

# 2️ Load Barcodes (Cell Names)
barcodes = pd.read_csv("barcodes.tsv", header=None, sep="\t").squeeze("columns").tolist()

# 3️ Load Gene Names
genes = pd.read_csv("features.tsv", header=None, sep="\t").squeeze("columns").tolist()

# 4️ Load Metadata
metadata = pd.read_csv("metadata.csv", index_col=0)
metadata = metadata.loc[metadata.index.intersection(barcodes)]  # Filter metadata to match cells

# 5️ Transpose Counts Matrix (if needed)
if counts.shape[1] == len(barcodes):  # Check if cells are in columns instead of rows
    counts = counts.T  # Transpose to (cells, genes)

# 6️ Create AnnData Object
adata = sc.AnnData(X=counts, obs=metadata, var=pd.DataFrame(index=genes))

# 8️ Save AnnData Object (Optional)
adata.write("global_reference.h5ad")  # Save in H5AD format for future use

for col in adata.obs.columns:
    adata.obs[col] = adata.obs[col].astype(str)

# ✅ Check Final Structure
print(adata)

def downsample_anndata(ad, column, n_cells):
    ad.obs['subpop_downsample'] = ad.obs[column]
    downsampled_data = ad.obs.groupby('subpop_downsample').sample(n=n_cells, replace=False, random_state=1)
    ad_downsampled = ad[downsampled_data.index]
    return ad_downsampled

# Downsample the AnnData object
downsampled_scdata = downsample_anndata(adata, 'cellTypeFine', 1000)

# Verify the shape of the downsampled data
print(downsampled_scdata)



import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import zipfile
import .wpy as sq
import spatialdata as sd
from spatialdata_io import xenium
import anndata as ad
import nico
from nico import Annotations as sann
import warnings
import time
warnings.filterwarnings('ignore')
#export PYTHONWARNINGS='ignore:Multiprocessing-backed parallel loops:UserWarning'
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
os.chdir("/oak/stanford/groups/longaker/BPP/New_NiCo")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
file_path = "/oak/stanford/groups/longaker/BPP/New_NiCo/Input/global_reference.h5ad"
# Load AnnData object
scdata = sc.read_h5ad(file_path)
scdata

#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/BPP/New_NiCo/Input/Region_1_0038414/outs/"
zarr_path = "/oak/stanford/groups/longaker/BPP/New_NiCo/Input/Region_1_0038414/outs/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)
sc.pp.filter_cells(ad_seq_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)
sc.pp.filter_genes(ad_seq_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)

#Process data
Original_counts=ad_seq_ori.copy()
Original_counts.raw=Original_counts.copy()

# Standard scanpy analysis

sc.pp.normalize_total(Original_counts)
sc.pp.log1p(Original_counts)

sc.tl.pca(Original_counts)
sc.pp.neighbors(Original_counts)
sc.tl.umap(Original_counts)
sc.pl.umap(Original_counts, color = "subpop")

scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
Original_counts.write_h5ad('/oak/stanford/groups/longaker/BPP/New_NiCo/Input/Original_counts.h5ad')

#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data
sc.experimental.pp.normalize_pearson_residuals(ad_seq_common,inplace=True) #ad_seq_common.X[ad_seq_common.X<0]=0

ad_seq_common.write_h5ad('/oak/stanford/groups/longaker/BPP/New_NiCo/Input/sct_singleCell.h5ad')

ad_spatial_common.raw=ad_spatial_common.copy()
sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='Region_1_0038555_res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/BPP/New_NiCo/Output/"

ad_spatial_common.write_h5ad(spdatapath+'Region_1_0038555_sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/BPP/New_NiCo/Output/Region_2_0038556/")
scdatapath = "/oak/stanford/groups/longaker/BPP/New_NiCo/Input/"
spdatapath = "/oak/stanford/groups/longaker/BPP/New_NiCo/Output/Region_2_0038556/"
ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='./nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_2_0038556_nico_celltype_annotation.h5ad'
inputRadius=0

annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib', 'STEAP4+ Mural Cells', 'CD300E+ Monocytes']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))

#Covariation analysis
os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
# #Remove the unmapped cells
# spdata_filtered = spdata[spdata.obs["nico_ct"] != "NM"].copy()
# spdata_filtered.write_h5ad("./nico_out/nico_celltype_annotation.h5ad")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
ref_cluster_tag = "subpop"
cov_out=scov.gene_covariation_analysis(iNMFmode=True,
        Radius=inputRadius,
        no_of_factors=3,
        refpath=ref_datapath,
        quepath=query_datapath,
        spatial_integration_modality='double',
        output_niche_prediction_dir=output_nico_dir,
        ref_cluster_tag=ref_cluster_tag) #LRdbFilename='NiCoLRdb.txt'

scov.plot_cosine_and_spearman_correlation_to_factors(cov_out,
choose_celltypes=['CTHRC1+ mFib'],
NOG_Fa=30,saveas=saveas,transparent_mode=transparent_mode,
figsize=(15,10))

dataFrame=scov.extract_and_plot_top_genes_from_chosen_factor_in_celltype(
cov_out,
choose_celltype='CTHRC1+ mFib',
choose_factor_id=1,
top_NOG=20,
correlation_with_spearman=True,
positively_correlated=False,
saveas=saveas,transparent_mode=transparent_mode )

dataFrame

scov.make_excel_sheet_for_gene_correlation(cov_out)

choose_celltypes=[]
scov.plot_significant_regression_covariations_as_circleplot(cov_out,
choose_celltypes=choose_celltypes,
mention_pvalue=True,
saveas=saveas,transparent_mode=transparent_mode,
figsize=(6,1.25))

scov.plot_significant_regression_covariations_as_heatmap(cov_out,
choose_celltypes=[],
saveas=saveas,transparent_mode=transparent_mode, figsize=(6,1.25))

scov.save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(cov_out,
pvalueCutoff=0.05,correlation_with_spearman=True,
LR_plot_NMF_Fa_thres=0.1,LR_plot_Exp_thres=0.1,number_of_top_genes_to_print=50)

scov.find_LR_interactions_in_interacting_cell_types(cov_out,
choose_interacting_celltype_pair=['CTHRC1+ mFib'],
choose_factors_id=[1,1],
pvalueCutoff=0.05,
LR_plot_NMF_Fa_thres=0.3,
LR_plot_Exp_thres=0.2,
saveas=saveas,transparent_mode=transparent_mode,figsize=(12, 10))

scov.pathway_analysis(cov_out,
choose_celltypes=['CTHRC1+ mFib'],
NOG_pathway=50,
choose_factors_id=[3],
savefigure=True,
positively_correlated=True,
saveas='pdf',
rps_rpl_mt_genes_included=False,
display_plot_as='dotplot',
correlation_with_spearman=True,
circlesize=12,
database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
object_for_color='Adjusted P-value',
object_for_xaxis='Combined Score',
fontsize=12,
showit=True,
input_colormap='plasma')

scov.pathway_analysis(cov_out,
choose_celltypes=['CTHRC1+ mFib'],
NOG_pathway=50,
choose_factors_id=[],
positively_correlated=True,
savefigure=True,
rps_rpl_mt_genes_included=False,
#database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
display_plot_as='barplot',
object_for_color='Adjusted P-value',
object_for_xaxis='Odds Ratio',
showit=True,
input_colormap='plasma')


scov.plot_top_genes_for_a_given_celltype_from_all_factors(
cov_out,choose_celltypes=['CTHRC1+ mFib'],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors(cov_out,
choose_interacting_celltype_pair=['CTHRC1+ mFib'],
visualize_factors_id=[3,3],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.visualize_factors_in_scRNAseq_umap(cov_out,
choose_interacting_celltype_pair=['CTHRC1+ mFib'],
visualize_factors_id=[1,3],
saveas=saveas,transparent_mode=transparent_mode,figsize=(8,3.5))

from numpy import load 
data = load('/oak/stanford/groups/longaker/KEBR/Xenium/NiCo/nico_out/niche_prediction_linear/classifier_matrices_0.npz')
lst = data.files
for item in lst:
    print(item)
    print(data[item])
    
    
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

coefficients = data['coef']
# Create a heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(coefficients, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5)
plt.title("Heatmap of Coefficients")
plt.show()




#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038414/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")

####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038414/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038414__Region_2__20241127__220634"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038414__Region_2__20241127__220634/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data

sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038414/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038414/")

scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038414/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='./nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_2_0038414_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib', 'STEAP4+ Mural Cells']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))

#Covariation analysis
os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
# #Remove the unmapped cells
# spdata_filtered = spdata[spdata.obs["nico_ct"] != "NM"].copy()
# spdata_filtered.write_h5ad("./nico_out/nico_celltype_annotation.h5ad")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
ref_cluster_tag = "subpop"
cov_out=scov.gene_covariation_analysis(iNMFmode=True,
        Radius=inputRadius,
        no_of_factors=3,
        refpath=ref_datapath,
        quepath=query_datapath,
        spatial_integration_modality='double',
        output_niche_prediction_dir=output_nico_dir,
        ref_cluster_tag=ref_cluster_tag) #LRdbFilename='NiCoLRdb.txt'

scov.plot_cosine_and_spearman_correlation_to_factors(cov_out,
choose_celltypes=['mCDF CTHRC1+'],
NOG_Fa=30,saveas=saveas,transparent_mode=transparent_mode,
figsize=(15,10))

dataFrame=scov.extract_and_plot_top_genes_from_chosen_factor_in_celltype(
cov_out,
choose_celltype='mCDF CTHRC1+',
choose_factor_id=1,
top_NOG=20,
correlation_with_spearman=True,
positively_correlated=False,
saveas=saveas,transparent_mode=transparent_mode )

dataFrame

scov.make_excel_sheet_for_gene_correlation(cov_out)

choose_celltypes=[]
scov.plot_significant_regression_covariations_as_circleplot(cov_out,
choose_celltypes=choose_celltypes,
mention_pvalue=True,
saveas=saveas,transparent_mode=transparent_mode,
figsize=(6,1.25))

scov.plot_significant_regression_covariations_as_heatmap(cov_out,
choose_celltypes=[],
saveas=saveas,transparent_mode=transparent_mode, figsize=(6,1.25))

scov.save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(cov_out,
pvalueCutoff=0.05,correlation_with_spearman=True,
LR_plot_NMF_Fa_thres=0.1,LR_plot_Exp_thres=0.1,number_of_top_genes_to_print=50)

scov.find_LR_interactions_in_interacting_cell_types(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Macrophages'],
choose_factors_id=[1,1],
pvalueCutoff=0.05,
LR_plot_NMF_Fa_thres=0.3,
LR_plot_Exp_thres=0.2,
saveas=saveas,transparent_mode=transparent_mode,figsize=(12, 10))

scov.pathway_analysis(cov_out,
choose_celltypes=['mCDF CTHRC1+', 'Macrophages'],
NOG_pathway=50,
choose_factors_id=[3],
savefigure=True,
positively_correlated=True,
saveas='pdf',
rps_rpl_mt_genes_included=False,
display_plot_as='dotplot',
correlation_with_spearman=True,
circlesize=12,
database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
object_for_color='Adjusted P-value',
object_for_xaxis='Combined Score',
fontsize=12,
showit=True,
input_colormap='plasma')

scov.pathway_analysis(cov_out,
choose_celltypes=['mCDF CTHRC1+'],
NOG_pathway=50,
choose_factors_id=[],
positively_correlated=True,
savefigure=True,
rps_rpl_mt_genes_included=False,
#database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
display_plot_as='barplot',
object_for_color='Adjusted P-value',
object_for_xaxis='Odds Ratio',
showit=True,
input_colormap='plasma')


scov.plot_top_genes_for_a_given_celltype_from_all_factors(
cov_out,choose_celltypes=['mCDF CTHRC1+','Macrophages', 'Vascular ECs'],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Vascular ECs'],
visualize_factors_id=[3,3],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.visualize_factors_in_scRNAseq_umap(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Vascular ECs'],
visualize_factors_id=[1,3],
saveas=saveas,transparent_mode=transparent_mode,figsize=(8,3.5))

from numpy import load 
data = load('/oak/stanford/groups/longaker/KEBR/Xenium/NiCo/nico_out/niche_prediction_linear/classifier_matrices_0.npz')
lst = data.files
for item in lst:
    print(item)
    print(data[item])
    
    
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

coefficients = data['coef']
# Create a heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(coefficients, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5)
plt.title("Heatmap of Coefficients")
plt.show()



#Region_1_0038438

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038438__Region_1__20241127__220634"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038438__Region_1__20241127__220634/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/")

scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038438/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='./nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_1_0038438_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib', 'STEAP4+ Mural Cells']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))

#Covariation analysis
os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
# #Remove the unmapped cells
# spdata_filtered = spdata[spdata.obs["nico_ct"] != "NM"].copy()
# spdata_filtered.write_h5ad("./nico_out/nico_celltype_annotation.h5ad")

# spdata = sc.read_h5ad("./nico_out/nico_celltype_annotation.h5ad")
# sc.pl.umap(spdata, color=["nico_ct"], 
#            title=["NiCo Clusters"],wspace=0.4,show=True)
ref_cluster_tag = "subpop"
cov_out=scov.gene_covariation_analysis(iNMFmode=True,
        Radius=inputRadius,
        no_of_factors=3,
        refpath=ref_datapath,
        quepath=query_datapath,
        spatial_integration_modality='double',
        output_niche_prediction_dir=output_nico_dir,
        ref_cluster_tag=ref_cluster_tag) #LRdbFilename='NiCoLRdb.txt'

scov.plot_cosine_and_spearman_correlation_to_factors(cov_out,
choose_celltypes=['mCDF CTHRC1+'],
NOG_Fa=30,saveas=saveas,transparent_mode=transparent_mode,
figsize=(15,10))

dataFrame=scov.extract_and_plot_top_genes_from_chosen_factor_in_celltype(
cov_out,
choose_celltype='mCDF CTHRC1+',
choose_factor_id=1,
top_NOG=20,
correlation_with_spearman=True,
positively_correlated=False,
saveas=saveas,transparent_mode=transparent_mode )

dataFrame

scov.make_excel_sheet_for_gene_correlation(cov_out)

choose_celltypes=[]
scov.plot_significant_regression_covariations_as_circleplot(cov_out,
choose_celltypes=choose_celltypes,
mention_pvalue=True,
saveas=saveas,transparent_mode=transparent_mode,
figsize=(6,1.25))

scov.plot_significant_regression_covariations_as_heatmap(cov_out,
choose_celltypes=[],
saveas=saveas,transparent_mode=transparent_mode, figsize=(6,1.25))

scov.save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(cov_out,
pvalueCutoff=0.05,correlation_with_spearman=True,
LR_plot_NMF_Fa_thres=0.1,LR_plot_Exp_thres=0.1,number_of_top_genes_to_print=50)

scov.find_LR_interactions_in_interacting_cell_types(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Macrophages'],
choose_factors_id=[1,1],
pvalueCutoff=0.05,
LR_plot_NMF_Fa_thres=0.3,
LR_plot_Exp_thres=0.2,
saveas=saveas,transparent_mode=transparent_mode,figsize=(12, 10))

scov.pathway_analysis(cov_out,
choose_celltypes=['mCDF CTHRC1+', 'Macrophages'],
NOG_pathway=50,
choose_factors_id=[3],
savefigure=True,
positively_correlated=True,
saveas='pdf',
rps_rpl_mt_genes_included=False,
display_plot_as='dotplot',
correlation_with_spearman=True,
circlesize=12,
database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
object_for_color='Adjusted P-value',
object_for_xaxis='Combined Score',
fontsize=12,
showit=True,
input_colormap='plasma')

scov.pathway_analysis(cov_out,
choose_celltypes=['mCDF CTHRC1+'],
NOG_pathway=50,
choose_factors_id=[],
positively_correlated=True,
savefigure=True,
rps_rpl_mt_genes_included=False,
#database=['GO_Biological_Process_2021'], #database=['BioPlanet_2019'],
display_plot_as='barplot',
object_for_color='Adjusted P-value',
object_for_xaxis='Odds Ratio',
showit=True,
input_colormap='plasma')


scov.plot_top_genes_for_a_given_celltype_from_all_factors(
cov_out,choose_celltypes=['mCDF CTHRC1+','Macrophages', 'Vascular ECs'],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Vascular ECs'],
visualize_factors_id=[3,3],
top_NOG=20,saveas=saveas,transparent_mode=transparent_mode)

scov.visualize_factors_in_scRNAseq_umap(cov_out,
choose_interacting_celltype_pair=['mCDF CTHRC1+','Vascular ECs'],
visualize_factors_id=[1,3],
saveas=saveas,transparent_mode=transparent_mode,figsize=(8,3.5))

from numpy import load 
data = load('/oak/stanford/groups/longaker/KEBR/Xenium/NiCo/nico_out/niche_prediction_linear/classifier_matrices_0.npz')
lst = data.files
for item in lst:
    print(item)
    print(data[item])
    
    
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

coefficients = data['coef']
# Create a heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(coefficients, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5)
plt.title("Heatmap of Coefficients")
plt.show()




#Region_2_0038438



os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038438/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038438/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038438__Region_2__20241127__220634"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038438__Region_2__20241127__220634/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038438/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038438/")

scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038438/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='./nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_2_0038438_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib', 'STEAP4+ Mural Cells']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))




#Region_1_0038555
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038555__Region_1__20241202__215325"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038555__Region_1__20241202__215325/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/")

scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038555/nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_1_0038555_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib', 'LRRC7+ mFib', 'PLIN1+ Adipocytes']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))



#Region_2_0038555
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038555__Region_2__20241202__215325"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038555__Region_2__20241202__215325/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')



####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/")
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/"
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038555/nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_2_0038555_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))






#Region_2_0038556
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038556__Region_1__20241202__215324"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038556__Region_1__20241202__215324/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')





####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/")
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_1_0038556/nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_1_0038556_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))






#Region_2_0038556
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/")
#os.chdir("/labs/longaker/USR/KE_Bauer-Rowe/MetaUpdate")


####Part 1: Prepare data 
##import datasets
# Load AnnData object
scdata = sc.read_h5ad("sct_singleCell.h5ad")
scdata
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/"
#generate object using python3.11
xenium_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038556__Region_2__20241202__215324"
zarr_path = "/oak/stanford/groups/longaker/KEBR/Xenium/Output/output-XETG00102__0038556__Region_2__20241202__215324/Xenium.zarr"

sdata = xenium(xenium_path)
spdata = sdata.tables["table"]
#Convert the integers to strings
spdata.obs_names = ["cell" + str(x) for x in map(str, spdata.obs_names)]
spdata

#Copy variables
ad_seq_ori = scdata
ad_spatial_ori = spdata

# Filter the cells and genes
sc.pp.filter_cells(ad_spatial_ori, min_counts=5)

sc.pp.filter_genes(ad_spatial_ori, min_cells=1)

print(ad_spatial_ori)
print(ad_seq_ori)


#Find shared genes 

sp_genename=ad_spatial_ori.var_names.to_numpy()
sc_genename=ad_seq_ori.var_names.to_numpy()

index_sp,index_sc=sann.find_index(sp_genename,sc_genename)
#ad_seq_common=ad_seq_ori[:,index_sc].copy()
ad_seq_common=ad_seq_ori.copy()
ad_spatial_common=ad_spatial_ori[:,index_sp].copy()

print(len(sp_genename[index_sp]))
print(len(sc_genename[index_sc]))

# Alternative 1
# The sctransform normalization function from scanpy

ad_seq_common.raw=ad_seq_common.copy()
ad_spatial_common.raw=ad_spatial_common.copy()
# perform scTranform normalization common gene space for spatial data and scRNAseq data


sc.experimental.pp.normalize_pearson_residuals(ad_spatial_common,inplace=True) #ad_spatial_common.X[ad_spatial_common.X<0]=0
#print(ad_spatial_common.X.toarray()

# standard scanpy analysis
sc.pp.pca(ad_spatial_common)
sc.pp.neighbors(ad_spatial_common,n_pcs=30)
sc.tl.umap(ad_spatial_common)

# visualize umap
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(ad_spatial_common, title=["Spatial umap on common gene space"],wspace=0.4,show=True)
sc.tl.leiden(ad_spatial_common, resolution=0.7,key_added="leiden0.7")

# Visualize your initial spatial clustering in the umap
# A good resolution parameter should yield clusters corresponding to major cell types.

sc.pl.umap(ad_spatial_common, color=["leiden0.7"], title=["Spatial umap"],wspace=0.4,
           show=True, save='res0.7_spatial_umap.png')
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/"

ad_spatial_common.write_h5ad(spdatapath+'sct_spatial.h5ad')





####Part 2
# if you installed the nico package

from nico import Annotations as sann
from nico import Interactions as sint
from nico import Covariations as scov

import matplotlib as plt


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.linewidth'] = 0.1 #set the value globally

# please use Helvetica font according to your OS to ensure compatibility with Adobe Illustrator.
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Use the default font for all the figures
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans','Lucida Grande', 'Verdana']

import warnings
warnings.filterwarnings("ignore")

#parameters for saving plots
saveas='png'
transparent_mode=False

#os.chdir("/oak/stanford/groups/longaker/KEBR/Xenium/Python/bin/nicoUser")
#ref_datapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/"
os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/")
scdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/"
spdatapath = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/"

ref_datapath = scdatapath
#query_datapath = "/oak/stanford/groups/longaker/KEBR/Xenium/"
query_datapath = spdatapath
output_nico_dir='/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/nico_out/'
output_annotation_dir=None #uses default location
#output_annotation_dir=output_nico_dir+'annotations/'
annotation_save_fname= 'Region_2_0038556_nico_celltype_annotation.h5ad'
inputRadius=0

ref_cluster_tag='subpop' #scRNAseq cell type slot
annotation_slot='leiden0.7' #spatial cell type slot

anchors_and_neighbors_info= sann.find_anchor_cells_between_ref_and_query(
refpath=ref_datapath,
quepath=query_datapath,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir)


output_info=sann.nico_based_annotation(anchors_and_neighbors_info,
guiding_spatial_cluster_resolution_tag='leiden0.7',
across_spatial_clusters_dispersion_cutoff=0.15,
ref_cluster_tag="subpop",
resolved_tie_issue_with_weighted_nearest_neighbor='No')

[]

sann.save_annotations_in_spatial_object(output_info,
anndata_object_name=annotation_save_fname)

sann.visualize_umap_and_cell_coordinates_with_all_celltypes(
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
#spatial_cluster_tag='nico_ct',
spatial_cluster_tag='nico_ct',
spatial_coordinate_tag='spatial',
umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

choose_celltypes=[['CTHRC1+ mFib']]

sann.visualize_umap_and_cell_coordinates_with_selected_celltypes(
choose_celltypes=choose_celltypes,
output_nico_dir=output_nico_dir,
output_annotation_dir=output_annotation_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',spatial_coordinate_tag='spatial',umap_tag='X_umap',
saveas=saveas,transparent_mode=transparent_mode)

#do_not_use_following_CT_in_niche=['Basophils','Cycling/GC B cell','pDC']

niche_pred_output=sint.spatial_neighborhood_analysis(
Radius=inputRadius,
output_nico_dir=output_nico_dir,
anndata_object_name=annotation_save_fname,
spatial_cluster_tag='nico_ct',
#removed_CTs_before_finding_CT_CT_interactions=do_not_use_following_CT_in_niche
)

celltype_niche_interaction_cutoff=0.1


import pygraphviz
a=pygraphviz.AGraph()
a._get_prog('neato')

# import os
# if  not '/home/[username]/miniforge3/envs/SC/bin/' in os.environ["PATH"]:
#     os.environ["PATH"] += os.pathsep + '/home/[username]/miniforge3/envs/SC/bin/'
dpi = 300
sint.plot_niche_interactions_without_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,                #Resolution in dots per inch for saving the figure.
input_colormap='jet',   #Colormap for node colors, from matplotlib colormaps.
with_labels=True,       #Display cell type labels on the nodes, if True.
node_size=500,          #Size of the nodes.
linewidths=0.5,         #Width of the node border lines.
node_font_size=6,       #Font size for node labels.
alpha=0.5,              #Opacity level for nodes and edges. 1 is fully opaque, and 0 is fully transparent.
font_weight='bold'      #Font weight for node labels; 'bold' for emphasis, 'normal' otherwise.
)

sint.plot_niche_interactions_with_edge_weight(niche_pred_output,
niche_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,
transparent_mode=transparent_mode,
showit=True,
figsize=(10,7),
dpi=dpi,
input_colormap='jet',
with_labels=True,
node_size=500,
linewidths=1,
node_font_size=8,
alpha=0.5,
font_weight='normal',
edge_label_pos=0.35,   #Relative position of the weight label along the edge.
edge_font_size=3       #Font size for edge labels.
)

#Select niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,
choose_celltypes=['CTHRC1+ mFib'],
celltype_niche_interaction_cutoff=celltype_niche_interaction_cutoff,
saveas=saveas,transparent_mode=transparent_mode,figsize=(4.0,2.0))

#All niche interaction scores
sint.find_interacting_cell_types(niche_pred_output,choose_celltypes=[])

#Check the classifier metrics
sint.plot_confusion_matrix(niche_pred_output,
saveas=saveas,transparent_mode=transparent_mode)

sint.plot_evaluation_scores(niche_pred_output,
saveas=saveas, transparent_mode=transparent_mode,
figsize=(4,3))




import os
import anndata as ad
import pandas as pd

# Define the base directory
base_dir = "/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo"

# List of region directories
region_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

# Loop through each region directory
for region in region_dirs:
    nico_out_path = os.path.join(base_dir, region, "nico_out")
    
    # Check if nico_out exists
    if os.path.isdir(nico_out_path):
        # Find the h5ad file inside nico_out
        h5ad_files = [f for f in os.listdir(nico_out_path) if f.endswith(".h5ad")]
        
        if h5ad_files:
            h5ad_path = os.path.join(nico_out_path, h5ad_files[0])  # Assuming one file per directory
            
            # Read h5ad file
            nico = ad.read_h5ad(h5ad_path)
            
            # Extract relevant columns
            df = nico.obs[['cell_id', 'nico_ct']].copy()
            
            # Save to CSV
            csv_path = os.path.join(nico_out_path, f"{region}_nico_cell_assignments.csv")
            df.to_csv(csv_path, index=False)
            
            print(f"CSV file saved at: {csv_path}")


# Read the CSV file into a DataFrame. Adjust the file path as necessary.
csv_file_path = 'Region_2_0038556_2.csv'
df = pd.read_csv(csv_file_path)

csv_cell_ids = df['cell_id']
# Find matching cell_ids in the AnnData object
matching_indices = adata.obs[adata.obs['cell_id'].isin(csv_cell_ids)].index

# Ensure the new category is added to the 'nico_ct' Categorical column
if 'nico_ct' in adata.obs:
    if not 'DEShi ssSMCs' in adata.obs['nico_ct'].cat.categories:
        adata.obs['nico_ct'] = adata.obs['nico_ct'].cat.add_categories(['DEShi ssSMCs'])

# Update the 'nico_ct' metadata column for the matching indices
adata.obs.loc[matching_indices, 'nico_ct'] = 'DEShi ssSMCs'

os.chdir("/oak/stanford/groups/longaker/KEBR/IBDAFs/MetaUpdate/NiCo/Region_2_0038556/")
file_path = 'nico_out/Region_2_0038556_nico_celltype_annotation.h5ad'
adata = ad.read_h5ad(file_path)

replacements = {
    "JUNBhi Stress-Responsive SMCs": "MYC+ Synthetic vSMCs",
    "LUM+ Fibroblast-Like SMCs": "DEShi ssSMCs",
    "PRKG1+ ssSMCs": "CACNA1C+ Contractile vSMCs"
}

# Ensure they exist in the categories if the column is categorical
if pd.api.types.is_categorical_dtype(adata.obs['nico_ct']):
    new_categories = adata.obs['nico_ct'].cat.categories.tolist()
    for old_value, new_value in replacements.items():
        if old_value in new_categories and new_value not in new_categories:
            new_categories.append(new_value)
    adata.obs['nico_ct'] = adata.obs['nico_ct'].cat.set_categories(new_categories)

# Apply the replacements
adata.obs['nico_ct'] = adata.obs['nico_ct'].replace(replacements)

# Verify the unique values after replacement
unique_nico_ct = adata.obs['nico_ct'].unique()
print(unique_nico_ct)
adata.write_h5ad(file_path)


csv_path = os.path.join(nico_out_path, f"Region_2_0034556_nico_cell_assignments.csv")
df.to_csv(csv_path, index=False)


 