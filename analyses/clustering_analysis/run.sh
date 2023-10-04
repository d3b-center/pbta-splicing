# 0) run script for optimal clustering
Rscript code/00-optimal-clustering.R \
--input_mat "input/non_expr_pan_cancer_splice_subset.rds" \
--cluster_algorithm "hc, km, pam" \
--cluster_distance "pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski" \
--filter_expr FALSE \
--protein_coding_only FALSE \
--feature_selection "dip.test" \
--transformation_type "none" \
--max_k 10

# 1) PBTA splicing data
# get ccp clustering output for a specific combination of distance + algorithm + % variable genes
# for the splicing dataset, we chose pam, pearson at 0% (based on Ammar's manual inspection)
Rscript code/01-get-clustering-output.R \
--input_mat "input/non_expr_pan_cancer_splice_subset.rds" \
--data_type "non_expr" \
--var_genes "0" \
--cluster_algorithm "pam" \
--cluster_distance "spearman" \
--prefix "non_expr_pan_cancer_splice_subset"

# get differential genes per cluster and perform pre-ranked gsea using those genes
# this was done for k = 3 with pam  + pearson + 0% genes
Rscript code/02-diff-genes-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_ccp.rds" \
--n_cluster "8" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_pam_spearman_0" \
--output_dir "output/diff_genes"

# get differential pathways per cluster using GSVA
# this was done for k = 3 with pam  + pearson + 0% genes
Rscript code/03-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_matrix.rds" \
--input_clin "../../data/histologies.tsv" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_ccp.rds" \
--n_cluster "8" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_pam_spearman_0" \
--output_dir "output/diff_pathways"

# get heatmap of the CCP matrix
Rscript code/04-plot-clustering-heatmap.R \
--ccp_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_ccp.rds" \
--input_clin "../../data/histologies.tsv" \
--n_cluster "8" \
--prefix "non_expr_pan_cancer_splice_subset"

# 2) PBTA mRNA expression data (subsetted to samples in the splicing data)
# get clustering output with same clustering parameters used for pbta splicing dataset
# we will not be using the CCP clustering output but only the filtered/normalized data matrix for downstream plotting
# we will be using the clustering output obtained from the splicing dataset
Rscript code/01-get-clustering-output.R \
--input_mat "input/raw_counts_pbta_subset.rds" \
--data_type "raw_counts" \
--var_genes "0" \
--cluster_algorithm "pam" \
--cluster_distance "spearman" \
--prefix "raw_counts_pbta_subset"

# get differential genes per cluster and perform pre-ranked gsea using those genes
# this was done for k = 3 with pam  + pearson + 0% genes
Rscript code/02-diff-genes-per-clusters.R \
--input_mat "output/ccp_output/raw_counts_pbta_subset_pam_spearman_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_ccp.rds" \
--n_cluster "8" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "raw_counts_pbta_subset_pam_spearman_0" \
--output_dir "output/diff_genes"

# get differential pathways per cluster using GSVA
# this was done for k = 3 with pam  + pearson + 0% genes
Rscript code/03-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/raw_counts_pbta_subset_pam_spearman_0_matrix.rds" \
--input_clin "../../data/histologies.tsv" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_pam_spearman_0_ccp.rds" \
--n_cluster "8" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "raw_counts_pbta_subset_pam_spearman_0" \
--output_dir "output/diff_pathways"

##plot cluster members categorized by SBI high vs low
Rscript code/05-plot-sbi_with_cluster-mem.R

## plot histologies across clusters
Rscript code/06-plot-histology-distr-across-clusters.R
