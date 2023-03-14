### Author: Komal S. Rathi
 
### Purpose

This module performs the following analyses:

1) Perform clustering on input data matrix using `ConsensusClusterPlus` and parameters chosen by user.
2) Identify differentially expressed genes per cluster of interest and perform pre-ranked pathway enrichment using `fgsea` on those genes.
3) Identify differentially regulated pathways per cluster of interest using `GSVA`.

The module can be run using a single bash command:

```
bash run.sh
```

### Input

Inputs used for all scripts in this module:

```
input
├── kegg_geneset_mrna.rds # kegg gene set for enrichment analyses
├── non_expr_pan_cancer_splice_subset.rds # pbta splicing data matrix. Rows are unique features and Columns are Kids_First_Biospecimen_Identifier
└── raw_counts_pbta_subset.rds # pbta mRNA data matrix subsetted to samples in the splicing matrix. Rows are unique features and Columns are Kids_First_Biospecimen_Identifier
```

### 01-get-clustering-output.R

This script performs clustering on input data matrix using `ConsensusClusterPlus` and parameters chosen by user.

#### Inputs

Input for this script is a data matrix with unique features as rows and sample identifiers as columns.

```
input
├── non_expr_pan_cancer_splice_subset.rds # pbta splicing data matrix.
└── raw_counts_pbta_subset.rds # pbta mRNA expression data matrix
```

#### Outputs

There are three outputs to this script per input data matrix:

1. `*_ccp.rds`: Full CCP output. 
2. `*_consensus.pdf`: Consensus PDF output.
3. `*_matrix.rds`: Filtered and normalized data matrix for downstream analyses.

```
# pbta splicing output
output/ccp_output
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_consensus.pdf
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds

# pbta mRNA output
output/ccp_output
├── raw_counts_pbta_subset_km_euclidean_0_ccp.rds
├── raw_counts_pbta_subset_km_euclidean_0_consensus.pdf
└── raw_counts_pbta_subset_km_euclidean_0_matrix.rds
```

**NOTE**: we will not be using the CCP clustering output of PBTA mRNA expression data matrix but only use the filtered/normalized data matrix for downstream plotting.

### 02-diff-genes-per-clusters.R

This script identifies differentially expressed genes per cluster of interest (`FDR < 0.05`) and perform pre-ranked pathway enrichment using `fgsea` on those genes (`FDR < 0.05`). 

#### Inputs

```
# KEGG geneset
input
└── kegg_geneset_mrna.rds # kegg gene set for enrichment analyses

# pbta splicing inputs
output/ccp_output
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds # full CCP output for PBTA splicing dataset
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds # filtered and normalized PBTA splicing data matrix for downstream analyses

# pbta mRNA inputs
output/ccp_output
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds # full CCP output for PBTA splicing dataset
└── raw_counts_pbta_subset_km_euclidean_0_matrix # filtered and normalized PBTA mRNA data matrix for downstream analyses 
```

#### Outputs

There are 2 types of outputs to this script per input data matrix:

1. `*_k_limma_output.tsv`: Differentially expressed genes in a specific cluster `k` obtained after running `limma`. 
2. `*_k_limma_output_pathway_ranked.tsv`: Differentially regulated pathways corresponding to the differential genes in a specific cluster `k` obtained after running pre-ranked GSEA (`fGSEA`).

In our case, we are interested in k = 3, so there are 2 output files per cluster i.e., a total of 6 output files per input data matrix.

```
# pbta splicing outputs
output/diff_genes
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_1_limma_output.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_1_limma_output_pathway_ranked.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_2_limma_output.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_2_limma_output_pathway_ranked.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_3_limma_output.tsv
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_3_limma_output_pathway_ranked.tsv

# pbta mRNA outputs
output/diff_genes
├── raw_counts_pbta_subset_km_euclidean_0_cluster_1_limma_output.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_1_limma_output_pathway_ranked.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_2_limma_output.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_2_limma_output_pathway_ranked.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_3_limma_output.tsv
└── raw_counts_pbta_subset_km_euclidean_0_cluster_3_limma_output_pathway_ranked.tsv
```

### 03-diff-pathways-per-clusters.R

This script identifies differentially regulated pathways per cluster of interest using `GSVA` (`FDR < 0.05`).

#### Inputs

```
# KEGG geneset
input
└── kegg_geneset_mrna.rds # kegg gene set for enrichment analyses

# pbta splicing inputs
output/ccp_output
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds # full CCP output for PBTA splicing dataset
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds # filtered and normalized PBTA splicing data matrix for downstream analyses

# pbta mRNA inputs
output/ccp_output
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds # full CCP output for PBTA splicing dataset
└── raw_counts_pbta_subset_km_euclidean_0_matrix # filtered and normalized PBTA mRNA data matrix for downstream analyses 
```

#### Outputs

There are 3 types of outputs to this script per input data matrix:

1. `*_gsva_output.tsv`: Full GSVA output scores per sample using `GSVA`. 
2. `*_k_pathway.tsv`: Differentially regulated pathways in a specific cluster `k` using `GSVA`.
3. `*_top20_pathways.pdf`: Heatmap of top 20 differentially regulated pathways per cluster `k`.

```
# pbta splicing outputs
output/diff_pathways
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_1_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_2_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_3_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_gsva_output.tsv
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_top20_pathways.pdf

# pbta mRNA outputs
output/diff_pathways
├── raw_counts_pbta_subset_km_euclidean_0_cluster_1_pathway.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_2_pathway.tsv
├── raw_counts_pbta_subset_km_euclidean_0_cluster_3_pathway.tsv
├── raw_counts_pbta_subset_km_euclidean_0_gsva_output.tsv
└── raw_counts_pbta_subset_km_euclidean_0_top20_pathways.pdf
```