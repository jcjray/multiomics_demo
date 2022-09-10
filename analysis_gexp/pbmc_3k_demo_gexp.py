import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
# from natsort import natsorted

sc.settings.autoshow = False
gexp_figure_dir = 'analysis_gexp/figures/'
sc._settings.ScanpyConfig.figdir = pathlib.Path(gexp_figure_dir)

if not os.path.isdir(gexp_figure_dir):
    os.mkdir(gexp_figure_dir)

# Load 10X dataset into Scanpy as an AnnData object
pbmc = sc.read_10x_h5('pbmc_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5');
pbmc.var_names_make_unique()

# Calculate mitochondrial fraction for QC filtering
mt_counts = pbmc[:,pbmc.var.filter(regex='^MT-',axis=0).index].to_df().sum(axis=1);
total_counts = pbmc.to_df().sum(axis=1);
pct_mt = mt_counts.div(total_counts).mul(100);
pct_mt.name = 'pct_mt'

# Run preprocessing
sc.pp.calculate_qc_metrics(pbmc, inplace=True)
sc.pp.filter_cells(pbmc, min_genes=1)
sc.pp.filter_genes(pbmc, min_counts=1)
# Note that while many of the downstream analyses in Scanpy depend on using log1p data,
# statistical models are often dependent on linear data.
sc.pp.log1p(pbmc)
sc.pp.highly_variable_genes(pbmc)
sc.pp.pca(pbmc)
sc.pp.normalize_total(pbmc,exclude_highly_expressed=True)
sc.pp.neighbors(pbmc)
pbmc.obs = pd.concat([pbmc.obs, pct_mt], axis=1)

# Unsurprisingly, the most highly expressed genes are MALAT1 and mitochondrial
sc.pl.highest_expr_genes(pbmc, save='.png')
sc.pl.highly_variable_genes(pbmc,  save='.png')

# We run the UMAP algorithm for visualization but be cautious interpreting the distances between points
# See, for example: https://www.biorxiv.org/content/10.1101/689851v6
sc.tl.umap(pbmc)

def check_marker_genes_in_adata(adata, gene_list):
    return [x for x in gene_list if x in adata.var_names]

# The first 43 genes are S phase the others are markers of G2M.
# See here: https://github.com/scverse/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
cc_genes = pd.read_csv('./ref_data/regev_lab_cell_cycle_genes.txt', header=None).T.values[0];
s_genes = cc_genes[:43];
g2m_genes = cc_genes[43:];

cc_genes = check_marker_genes_in_adata(pbmc, cc_genes);
# s_genes = check_marker_genes_in_adata(pbmc, s_genes);
# g2m_genes = check_marker_genes_in_adata(pbmc, g2m_genes);

sc.tl.score_genes_cell_cycle(pbmc, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pl.umap(pbmc, color = 'phase', save='_cell_cycle.png')

pbmc_regressed = sc.pp.regress_out(pbmc, ['S_score','G2M_score'], copy=True);

# Quality control metrics suggest that two clusters contain lower quality cells than others
sc.pl.umap(pbmc_regressed, color = 'log1p_total_counts', save='_log1p_cts.png')
sc.pl.umap(pbmc_regressed, color = 'pct_counts_in_top_500_genes', save='_pct_cts_500.png')
sc.pl.umap(pbmc_regressed, color = 'pct_mt', save='_pct_mt.png')

plot_dat = pbmc_regressed.obs.loc[:,['log1p_n_genes_by_counts','log1p_total_counts']];
plot_dat_melt = plot_dat.melt();

f,ax=plt.subplots();
sns.violinplot(data=plot_dat_melt, y='value', x='variable', palette='Set2');
sns.stripplot(data=plot_dat_melt, y='value', x='variable', palette='Dark2', alpha=0.5, size=1.5);
plt.savefig(gexp_figure_dir+'/violin_count_summary.png');

pct_summary_vars = ['pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'pct_mt'];
plot_dat = pbmc_regressed.obs.loc[:,pct_summary_vars];
plot_dat_melt = plot_dat.melt();

f,ax=plt.subplots();
sns.violinplot(data=plot_dat_melt, y='value', x='variable', palette='Set2');
sns.stripplot(data=plot_dat_melt, y='value', x='variable', palette='Dark2', alpha=0.5, size=1.5);
plt.xticks(rotation=45, ha='right');
plt.savefig(gexp_figure_dir+'/violin_pct_summary.png');

# Applying thresholds in two of the dimensions appears to rid the dataset of the most extreme outliers
pct_summary_vars = ['pct_counts_in_top_100_genes', 'pct_mt', 'log1p_total_counts'];
plot_dat = pbmc_regressed.obs.loc[:,pct_summary_vars];

f,ax=plt.subplots();
sns.scatterplot(data=plot_dat, y='pct_counts_in_top_100_genes', x='pct_mt', hue='log1p_total_counts', palette='viridis', size=1.5);
top_100_threshold = 57;
pct_mt_threshold = 32;
plt.hlines(top_100_threshold, 0, 80, linestyles='dashed', linewidths=0.5, colors='k');
plt.vlines(pct_mt_threshold, 0, 100, linestyles='dashed', linewidths=0.5, colors='k');
plt.savefig(gexp_figure_dir+'/scatter_qc.png');

pbmc_obs_filt = pbmc_regressed.obs.where(pbmc_regressed.obs.loc[:,'pct_mt'].lt(pct_mt_threshold));
pbmc_obs_filt = pbmc_obs_filt.where(pbmc_obs_filt.loc[:,'pct_counts_in_top_100_genes'].lt(top_100_threshold)).dropna();

pbmc_qc = pbmc_regressed[pbmc_obs_filt.index].copy();

pbmc_qc.write('analysis_gexp/pbmc_gexp_qc.h5ad', compression="gzip")
