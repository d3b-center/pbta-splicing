### Author: Komal S. Rathi, Ammar S. Naqvi, Jo Lynne Rokita

### Purpose

This module performs the following analyses:
1) Pre-process splicing data to create matrix and input
2) Optimal cluster selection by evaluating multiple combinations of clustering algorithm, distances and k values using [ClusTarIDseq](https://github.com/d3b-center/ClusTarIDseq)
3) Perform clustering on input data matrix using `ConsensusClusterPlus` and parameters chosen by user.
4) Identify differentially expressed genes per cluster of interest and perform pre-ranked pathway enrichment using `fgsea` on those genes.
5) Identify differentially regulated pathways per cluster of interest using `GSVA`.
6) Visualize cluster members categorized by SBI high vs low
7) Visualize histology membership across clusters

The module can be run using a single bash command:

```
bash run.sh
```

### Input

Inputs used for all scripts are created using the pre-processing scripts:

```
perl code/00-create_matrix_of_PSI_SE_gene.pl ../../data/histologies.tsv ../../data/splice-events-rmats.tsv.gz ../../data/independent-specimens.rnaseqpanel.primary.tsv input/pan_cancer_splicing_SE.gene.tsv
Rscript --vanilla code/01-convert_to_rds.R
Rscript --vanilla code/02-create-inputs.R
```

These scripts will generate the following input files in the `input` directory:

```
input
├── pan_cancer_splicing_SE.gene.rds #  pbta splicing data matrix. Rows are unique features and Columns are Kids_First_Biospecimen_Identifier
├── kegg_geneset_mrna.rds # kegg gene set for enrichment analyses
└── raw_counts_pbta_subset.rds # pbta mRNA data matrix subsetted to samples in the splicing matrix. Rows are unique features and Columns are Kids_First_Biospecimen_Identifier
```

### 03-optimal-clustering.R

This script runs the `lspline_clustering` function from [ClusTarIDseq](https://github.com/d3b-center/ClusTarIDseq) on the input matrix. 

No transformation was applied to the input matrix and PSI values were used as is. To reduce the dimensionality of the input matrix, there are two methods that were used for feature selection: 1. `Variance-based` filtering and 2. `Hartigans' dip test`. The former identifies top n % (user-defined) variable features and the latter identifies dips in the distribution of input features and selects features that have a bi or multi-modal distribution across the input samples. The idea behind this strategy is that these “dips” in the distribution may correspond to differences within underlying clinical variables of interest.

The `lspline_clustering` function was applied on the non-transformed, feature selected input matrix. The output of the function is a tsv file with the following metrics for each input `clustering algorithm`, `distance metric`, `k-value` and `feature selection` combination:

* `Stretch`: Calculation of the most linear/stable k from the consensus CDF curve plot using lspline package
* `P-value`: Calculate p-value (KS test) between adjacent k curves
* `Delta AUC`: Calculate increase in consensus (Delta k/Delta AUC) between adjacent k clustering

Cluster compactness and separation coefficients using `fpc::cluster.stats`:

* `Average.between`: average between cluster dissimilarity or separation
* `Average.within`: average within cluster dissimilarity or separation
* `Within.cluster.ss`: Silhouette score; -1 to 1. Close to 1 is good. Means that data point is close to the points in the assigned cluster and farther from the points in other clusters.
* `Dunn`: dunn index: a ratio of the smallest inter-cluster distance and the largest intra-cluster distance. A higher Dunn Index will indicate compact, well-separated clusters, while a lower index will indicate less compact or less well-separated clusters.
* `Entropy`: entropy (entropy should be less i.e. amount of misclassification, purity should be more i.e. amount of true classification)
* `Avg_sil`: Average silhouette width (average silhouette score across all clusters).

Next, it uses the [COINr](https://github.com/bluefoxr/COINr) R package to assign weights to the above metrics in order to generate a `composite score` that represents  `cluster quality` and a corresponding `rank` to all the combinations. Largest composite score gets the top rank and so on.

Following are the output files created:

```
# feature_selection can be either dip_test or if variance is selected, it will be of the format "var_{selected_proportion}"
output/optimal_clustering
├── {hc, km, pam}_{pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski}_{feature_selection}.pdf # consensus PDF file
├── {hc, km, pam}_{pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski}_{feature_selection}.rds # CCP output
├── {hc, km, pam}_{pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski}_{feature_selection}.tsv # clustering metrics
├── lspline_output.tsv # clustering and ranking metrics for all evaluated input combinations
└── ccp_optimal_clusters.tsv # classification of top ranking cluster
```

### 04-get-clustering-output.R

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

### 05-diff-genes-per-clusters.R

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

### 06-diff-pathways-per-clusters.R

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
3. `*_top*n*_pathways.pdf`: Heatmap of top n differentially regulated pathways per cluster `k`.

```
# pbta splicing outputs
output/diff_pathways
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_1_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_2_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_cluster_3_pathway.tsv
├── non_expr_pan_cancer_splice_subset_km_euclidean_0_gsva_output.tsv
└── non_expr_pan_cancer_splice_subset_km_euclidean_0_top*n*_pathways.pdf


### 07-plot-clustering-heatmap.R

This script generates heatmap based on the CCP matrix from above

#### Inputs

```
input
├── ../cohort_summary/results/histologies-plot-group.tsv ## histology file with plot groups
└── output/ccp_output/non_expr_pan_cancer_splice_subset_pam_canberra_0_ccp.rds ## results from CCP run
```

#### Outputs

There are 3 types of outputs to this script per input data matrix:

1. `non_expr_pan_cancer_splice_subset_ccp_heatmap.tiff`: Heatmap of CCP matrix colored by plot groups`


### 08-plot-sbi_with_cluster-mem.R

This script generates plot of cluster members categorized by splicing burden

#### Inputs

```
input
├── lspline_output.tsv ## lspline table from CCP script
└── input/splicing_index.total.txt ## results SBI module
```

#### Outputs

There are 3 types of outputs to this script per input data matrix:

1. `cluster_by_sbi.pdf`: Bar plots of clusters stratified by high vs low SBI


### 09-plot-histology-distr-across-clusters.R

This script generates plot histology membership across the different clusters

#### Inputs

```
input
└── ccp_optimal_clusters.tsv ## output from optimal clustering script
```

#### Outputs

There are 3 types of outputs to this script per input data matrix:

1. `cluster_membership.pdf`: Bar plots of cluster membership categorized by histologies/plot groups
2. ` cluster_members_by_cancer_group_subtype.tsv`: Table of results
