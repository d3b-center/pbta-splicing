# 1) PBTA splicing data
# get ccp clustering output for a specific combination of distance + algorithm + % variable genes 
# for the splicing dataset, we chose km, euclidean at 0% (based on Ammar's manual inspection)
Rscript code/01-get-clustering-output.R \
--input_mat "input/non_expr_pan_cancer_splice_subset.rds" \
--data_type "non_expr" \
--var_genes "0" \
--cluster_algorithm "km" \
--cluster_distance "euclidean" \
--prefix "non_expr_pan_cancer_splice_subset"

# get differential genes per cluster and perform pre-ranked gsea using those genes
# this was done for k = 3 with km + euclidean + 0% genes
Rscript code/02-diff-genes-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds" \
--n_cluster "3" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_km_euclidean_0" \
--output_dir "output/diff_genes"

# get differential pathways per cluster using GSVA
# this was done for k = 3 with km + euclidean + 0% genes
Rscript code/03-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds" \
--n_cluster "3" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "non_expr_pan_cancer_splice_subset_km_euclidean_0" \
--output_dir "output/diff_pathways"

# get heatmap of the CCP matrix
Rscript code/04-plot-clustering-heatmap.R \
--ccp_output "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds" \
--n_cluster "3" \
--prefix "non_expr_pan_cancer_splice_subset"

# 2) PBTA mRNA expression data (subsetted to samples in the splicing data)
# get clustering output with same clustering parameters used for pbta splicing dataset
# we will not be using the CCP clustering output but only the filtered/normalized data matrix for downstream plotting 
# we will be using the clustering output obtained from the splicing dataset 
Rscript code/01-get-clustering-output.R \
--input_mat "input/raw_counts_pbta_subset.rds" \
--data_type "raw_counts" \
--var_genes "0" \
--cluster_algorithm "km" \
--cluster_distance "euclidean" \
--prefix "raw_counts_pbta_subset"

# get differential genes per cluster and perform pre-ranked gsea using those genes
# this was done for k = 3 with km + euclidean + 0% genes
Rscript code/02-diff-genes-per-clusters.R \
--input_mat "output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds" \
--n_cluster "3" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "raw_counts_pbta_subset_km_euclidean_0" \
--output_dir "output/diff_genes"

# get differential pathways per cluster using GSVA
# this was done for k = 3 with km + euclidean + 0% genes
Rscript code/03-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/raw_counts_pbta_subset_km_euclidean_0_matrix.rds" \
--cluster_output "output/ccp_output/non_expr_pan_cancer_splice_subset_km_euclidean_0_ccp.rds" \
--n_cluster "3" \
--gene_set "input/kegg_geneset_mrna.rds" \
--prefix "raw_counts_pbta_subset_km_euclidean_0" \
--output_dir "output/diff_pathways"

##plot cluster members categorized by SBI high vs low
Rscript 05-plot-sbi_with_cluster-mem.R
