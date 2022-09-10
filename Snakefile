rule all:
    input:
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_colData.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_elementMetadata.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_gene_score_matrix.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_imputed_matrix.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_lsi_umap_coords.csv",
        "analysis_gexp/pbmc_gexp_qc.h5ad",
        "analysis_multiomic/figures/scatter_CD4_CD8.png"

rule atac_archr:
    input:
        "analysis_atac/pbmc_3k_demo.R"
    output:
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_colData.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_elementMetadata.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_gene_score_matrix.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_imputed_matrix.csv",
        "analysis_atac/pbmc_3k_archr_csv/pbmc_3k_atac_lsi_umap_coords.csv",
    shell:
        "R -s -f {input}"

rule gexp_scanpy:
    input:
        "analysis_gexp/pbmc_3k_demo_gexp.py"
    output:
        "analysis_gexp/pbmc_gexp_qc.h5ad"
    shell:
        "python {input}"

rule multiomic_scanpy:
    input:
        "analysis_multiomic/pbmc_3k_demo_multiomic.py"
    output:
        "analysis_multiomic/figures/scatter_CD4_CD8.png"
    shell:
        "python {input}"
