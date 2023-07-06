# run script for optimal clustering
Rscript code/00-optimal-clustering.R \
--input_mat "input/non_expr_pan_cancer_splice_subset.rds" \
--cluster_algorithm "hc, km, pam" \
--cluster_distance "pearson, spearman, euclidean, manhattan, binary, maximum, canberra, minkowski" \
--filter_expr FALSE \
--protein_coding_only FALSE \
--feature_selection "dip.test" \
--transformation_type "none" \
--max_k 10
