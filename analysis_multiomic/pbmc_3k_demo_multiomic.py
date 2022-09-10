import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from natsort import natsorted

sc.settings.autoshow = False

multiomic_figure_dir = 'analysis_multiomic/figures/'
if not os.path.isdir(multiomic_figure_dir):
    os.mkdir(multiomic_figure_dir)

# Load QC'd h5ad AnnData object
pbmc_qc = sc.read_h5ad('analysis_gexp/pbmc_gexp_qc.h5ad')

# Load ATAC data
atac_col = pd.read_csv('analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_colData.csv', index_col=0);
atac_emd = pd.read_csv('analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_elementMetadata.csv', index_col=0);
atac_gsm = pd.read_csv('analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_gene_score_matrix.csv', index_col=0);
atac_umap_coords = pd.read_csv('analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_lsi_umap_coords.csv', index_col=0);
atac_imp = pd.read_csv('analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_imputed_matrix.csv', index_col=0);

# The column indicating ATAC-based cluster membership has an unfortunate name - let's rename it
atac_col.columns = ['BlacklistRatio', 'DoubletEnrichment', 'DoubletScore', 'nDiFrags', 'nFrags',
                    'nMonoFrags', 'nMultiFrags', 'NucleosomeRatio', 'PassQC', 'PromoterRatio',
                    'ReadsInBlacklist', 'ReadsInPromoter', 'ReadsInTSS', 'Sample', 'TSSEnrichment',
                    'ATAC_Cluster'];

# Plot ATAC-based QC data on the ATAC-based UMAP projection
atac_summary_data = pd.concat([atac_umap_coords,atac_col],axis=1);
for q in atac_col.columns[:-1]:
    f,ax = plt.subplots();
    sns.scatterplot(data=atac_summary_data,
                    y='IterativeLSI#UMAP_Dimension_2', 
                    x='IterativeLSI#UMAP_Dimension_1',
                    hue=q, palette='viridis');
    plt.legend(bbox_to_anchor=(1,1));
    plt.title(q);
    plt.savefig(multiomic_figure_dir+'/atac_umap_qc_'+q+'.png');
    
# Same UMAP, colored by ATAC cluster
sorted_clust_names = natsorted(atac_summary_data.loc[:,'ATAC_Cluster'].drop_duplicates().values)
q = 'ATAC_Cluster';
f,ax = plt.subplots();
sns.scatterplot(data=atac_summary_data,
                y='IterativeLSI#UMAP_Dimension_2',
                x='IterativeLSI#UMAP_Dimension_1',
                hue=q, palette='Set3', hue_order=sorted_clust_names);
plt.legend(bbox_to_anchor=(1,1));
plt.title(q);
plt.savefig(multiomic_figure_dir+'/atac_umap_qc_'+q+'.png');

# Marker genes for immune cells from https://www.archrproject.com/articles/Articles/tutorial.html
marker_genes = [
    "CD34",  # Early Progenitor
    "GATA1", # Erythroid
    "PAX5", "MS4A1", "MME", # B-Cell Trajectory
    "CD14", "MPO", # Monocytes
    "CD3D", "CD8A" # TCells
];

# New directory for gene score UMAPs
gscore_umap_dir = 'analysis_multiomic/gscore_figures/'
if not os.path.isdir(gscore_umap_dir):
    os.mkdir(gscore_umap_dir)

# Raw gene scores without imputation
atac_gscore_data = pd.concat([atac_umap_coords,atac_gsm.T],axis=1);
for q in marker_genes:
    f,ax = plt.subplots();
    sns.scatterplot(data=atac_gscore_data,
                    y='IterativeLSI#UMAP_Dimension_2',
                    x='IterativeLSI#UMAP_Dimension_1', hue=q, palette='viridis');
    plt.legend(bbox_to_anchor=(1,1));
    plt.title(q);
    plt.savefig(gscore_umap_dir+'/atac_umap_gscore_'+q+'.png');
    
# Gene scores with imputation
atac_imp_data = pd.concat([atac_umap_coords,atac_imp.T],axis=1);
for q in marker_genes:
    f,ax = plt.subplots();
    sns.scatterplot(data=atac_imp_data,
                    y='IterativeLSI#UMAP_Dimension_2',
                    x='IterativeLSI#UMAP_Dimension_1', hue=q, palette='viridis');
    plt.legend(bbox_to_anchor=(1,1));
    plt.title(q);
    plt.savefig(gscore_umap_dir+'/atac_umap_imputedgscore_'+q+'.png');

#### Merge ATAC and gene expression data for comparison

# We make some UMAPs colored by various summary statistics
other_umap_figure_dir = 'analysis_multiomic/other_umap_figures/'
if not os.path.isdir(other_umap_figure_dir):
    os.mkdir(other_umap_figure_dir)


# Rename cells in the gene expression object to have experimental rep name prepended
# This will make cell IDs match ArchR output
gexp_index = pd.Index(['pbmc_3k#'+q for q in pbmc_qc.obs.index]);
pbmc_qc.obs.index = gexp_index

keep_cells = atac_col.index.intersection(gexp_index);

pbmc_filt = pbmc_qc[keep_cells];

gexp_umap_coords_filt = pd.DataFrame(data=pbmc_filt.obsm['X_umap'],
                                     index=pbmc_filt.obs.index,
                                     columns=['GEXP_UMAP1','GEXP_UMAP2']);
atac_umap_coords_filt = atac_umap_coords.loc[keep_cells];

atac_col_filt = atac_col.loc[keep_cells];
atac_emd_filt = atac_emd;
atac_gsm_filt = atac_gsm.loc[:,keep_cells];
atac_imp_filt = atac_imp.loc[:,keep_cells];

keep_genes = pd.Index(atac_emd_filt.loc[:,'name']).intersection(pbmc_filt.var.index);

pbmc_filt = pbmc_filt[:,keep_genes];

atac_emd_filt = atac_emd_filt.where(atac_emd_filt.loc[:,'name'].isin(keep_genes)).dropna();
atac_gsm_filt = atac_gsm_filt.loc[keep_genes];
atac_imp_filt = atac_imp_filt.loc[keep_genes];

sorted_clust_names = natsorted(atac_summary_data.loc[:,'ATAC_Cluster'].drop_duplicates().values);

######### Helper functions for plotting comparisons

def modal_rel_scatterplot(marker_case, x_scores, y_scores, hue_scores, hue_order=sorted_clust_names, x_name='Imputed ATAC Score', y_name='Normed Expression Level', hue_name='ATAC_Cluster', palette='Set3', point_size=1.5):
    ''' This automates plotting across modes in single-cell gene expression data,
        for example ATAC versus gene expression. x_scores & hue_scores should be
        Pandas DataFrames; y_scores should be an AnnData object. '''
    
    gene_score_case = pd.concat(
        [hue_scores.loc[:,hue_name],
         x_scores.loc[marker_case],
         y_scores[:,marker_case].to_df()],axis=1);
    gene_score_case.columns = [hue_name, x_name, y_name];
    
    f,ax = plt.subplots(figsize=(5,3));
    sns.scatterplot(data=gene_score_case,
                    x=x_name, y=y_name,
                    hue=hue_name, palette=palette, hue_order=sorted_clust_names);
    plt.legend(bbox_to_anchor=(1,1), title=hue_name);
    plt.title(marker_case);
    plt.savefig(gscore_umap_dir+'/atac_umap_imputedgscore_'+marker_case+'.png');

def color_umap(umap_coords, hue_data, title=None, palette='viridis', point_size=1.5):
    ''' Take a dataframe with UMAP coordinates and concat it with a dataframe
        having the same row names containing data for coloring the UMAP. '''
    
    plot_dat = pd.concat([umap_coords,hue_data],axis=1);
    
    f,ax = plt.subplots();
    sns.scatterplot(data=plot_dat,
                    x=plot_dat.columns[0], y=plot_dat.columns[1],
                    hue=plot_dat.columns[2], palette=palette, size=point_size);
    plt.title(title);
    if not title==None:
        plt.savefig(other_umap_figure_dir+'/atac_umap_'+title+'.png');
    else:
        plt.savefig(other_umap_figure_dir+'/atac_umap_.png');
    
def scan_markers_scatter(marker_list, x_scores=atac_imp_filt, y_scores=pbmc_filt, hue_scores=atac_col_filt):
    marker_not_present = [];
    for marker_case in marker_list:
        if np.any(x_scores.index.values == marker_case):
            modal_rel_scatterplot(marker_case, x_scores=x_scores, y_scores=y_scores, hue_scores=hue_scores);
        else:
            marker_not_present.append(marker_case);
    if len(marker_not_present) > 0:
        print('Markers not in dataset:\t', marker_not_present)

        
######## Plot comparisons using helper functions

# Myeloid lineage markers
markers_myeloid = ['ITGAM']; # aka CD11b
markers_hla = ['HLA-DRA','HLA-DRB1', 'HLA-DRB2', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQB1'];

# Note the high number of zeroes in the gene expression data - this is dropouts
# scATAC-seq has dropouts as well, but this is showing the IMPUTED gene scores
# That means the ATAC-based gene scores are statistical estimates based on properties
# of the surrounding genome compared to similar cells, not raw scores.
marker_list = markers_myeloid;
scan_markers_scatter(marker_list);

marker_list = markers_hla;
scan_markers_scatter(marker_list);

# MDSC markers
marker_list = [
    'FUT4', #CD15
    'CD14'];

scan_markers_scatter(marker_list);

marker_case = 'CD14';
color_umap(atac_umap_coords_filt, atac_imp_filt.loc[marker_case], title='ATAC Score');

marker_case = 'CD14';
color_umap(atac_umap_coords_filt, pbmc_filt[:,marker_case].to_df(), title='Normalized log1p Gene Expression');

# T Cell Markers
marker_list = [
    'CD3D',
    'CD3E',
    'CD3G',
    'CD8A',
    'CD8B',
    'CD4',
    'NCAM1', # 'CD56',
    'CD19',
    'TNFRSF17'
];

scan_markers_scatter(marker_list);

# First, let's revisit and annotate the putative major
# lineages of the ATAC clusters based on individual markers
# C1, C2, C3, and C4 are generally high in HLA and lower in T cell markers - possibly myeloid.
# C5 is similar to the low-quality cells filtered out. They may be approaching apoptosis.
# C6, C7, C8, C9, C10, and C11 have T cell markers but lower scores for HLA. Possibly lymphoid.
atac_summary_data = pd.concat([atac_umap_coords_filt,atac_col_filt],axis=1);
sorted_clust_names = natsorted(atac_summary_data.loc[:,'ATAC_Cluster'].drop_duplicates().values)
q = 'ATAC_Cluster';
f,ax = plt.subplots();
sns.scatterplot(data=atac_summary_data,
                y='IterativeLSI#UMAP_Dimension_2',
                x='IterativeLSI#UMAP_Dimension_1',
                hue=q, palette='Set3', hue_order=sorted_clust_names);
plt.legend(bbox_to_anchor=(1,1));
plt.title(q);

# This just draws the annotations
plt.text(-1, -3, 'Apoptotic?');
plt.arrow(-1.25,-2.75,-0.7,1.6, linewidth = 0.25);
plt.arrow(-1.25,-2.75,-1.2,-0.75, linewidth = 0.25);
plt.text(0, -5, 'Myeloid');
plt.arrow(-0.25,-4.75,-2,-1, linewidth = 0.25);
plt.arrow(2.75,-4.75,2,1, linewidth = 0.25);
plt.text(3, 3, 'Lymphoid');
plt.arrow(2.75,3.25,-2,0, linewidth = 0.25);
plt.savefig(other_umap_figure_dir+'/atac_umap_annotated.png');

# Within the lymphoid lineage, there are clear subpopulations of 
# CD4+/CD8- T helper cells (members of C7)
# CD4-/CD8+ cytotoxic T cells (members of C10)
# and other cells types that will require different marker combinations to identify.
marker_x = 'CD4';
marker_y = 'CD8A';
hue_scores = atac_col_filt;
hue_name = 'ATAC_Cluster';
gene_score_case = pd.concat([hue_scores.loc[:,hue_name], atac_imp_filt.loc[marker_x], atac_imp_filt.loc[marker_y]],axis=1);
gene_score_case.columns = [hue_name, marker_x, marker_y];

f,ax = plt.subplots(figsize=(5,3));
sns.scatterplot(data=gene_score_case, x=marker_x, y=marker_y, hue=hue_name, palette='Set3', hue_order=sorted_clust_names);
plt.legend(bbox_to_anchor=(1,1), title=hue_name);
plt.savefig(multiomic_figure_dir+'/scatter_CD4_CD8.png');
