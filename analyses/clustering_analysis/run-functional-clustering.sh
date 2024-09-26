# 1) scripts to create matrix and input
Rscript --vanilla code/01a-convert_func_to_rds.R

# 2) run script for optimal clustering
Rscript code/03-optimal-clustering.R \
--input_mat "input/functional-sites/pan_cancer_splicing_SE.gene.rds" \
--cluster_algorithm "hc, km, pam" \
--cluster_distance "pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski" \
--filter_expr FALSE \
--protein_coding_only FALSE \
--feature_selection "dip.test" \
--transformation_type "none" \
--max_k 15 \
--output_dir "output/functional-sites/optimal_clustering"

# get ccp clustering output for a specific combination of distance + algorithm + % variable genes
# for the splicing dataset, pam + canberra k = 14, 0% based on optimal clustering
Rscript code/04-get-clustering-output.R \
--input_mat "input/functional-sites/pan_cancer_splicing_SE.gene.rds" \
--data_type "non_expr" \
--var_genes "0" \
--cluster_algorithm "pam" \
--cluster_distance "canberra" \
--prefix "functional_non_expr_pan_cancer_splice_subset"

# get differential pathways per cluster using GSVA
# this was done for k = 14 with pam  + canberra + 0% genes
Rscript code/06-diff-pathways-per-clusters.R \
--input_mat "output/ccp_output/functional_non_expr_pan_cancer_splice_subset_pam_canberra_0_matrix.rds" \
--input_clin "../../data/histologies-plot-group.tsv" \
--cluster_output "output/ccp_output/functional_non_expr_pan_cancer_splice_subset_pam_canberra_0_ccp.rds" \
--n_cluster "14" \
--gene_set "input/hallmark_splice_geneset_mrna.rds" \
--prefix "functional_non_expr_pan_cancer_splice_subset_pam_canberra_0" \
--output_dir "output/functional-sites/diff_pathways"

# get heatmap of the CCP matrix
Rscript code/07-plot-clustering-heatmap.R \
--ccp_output "output/ccp_output/functional_non_expr_pan_cancer_splice_subset_pam_canberra_0_ccp.rds" \
--input_clin "../../data/histologies-plot-group.tsv" \
--n_cluster "14" \
--prefix "functional_non_expr_pan_cancer_splice_subset"

## plot histologies across clusters
Rscript code/09a-plot-histology-distr-across-clusters.R
