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
--max_k 17 \
--output_dir "output/functional-sites/optimal_clustering"
