# pbta splicing ranked
Rscript code/diff-expr-with-clusters.R \
--mat '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset.rds' \
--type 'non_expr' \
--var_genes 0 \
--cluster_output '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds' \
--cluster 3 \
--pathways '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds' \
--prefix 'non_expr_pan_cancer_splice_subset_km_euclidean_0_3' \
--output_dir '~/Projects/Test/pbta_splicing/output'

# pbta splicing unbiased
Rscript code/diff-pathways-with-clusters.R \
--mat '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset.rds' \
--type 'non_expr' \
--var_genes 0 \
--cluster_output '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds' \
--cluster 3 \
--pathways '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds' \
--prefix 'non_expr_pan_cancer_splice_subset_km_euclidean_0_3' \
--output_dir '~/Projects/Test/pbta_splicing/output'

# pbta gene expression ranked
Rscript code/diff-expr-with-clusters.R \
--mat '~/Projects/Test/pbta_splicing/data/raw_counts_pbta_subset.rds' \
--type 'raw_counts' \
--var_genes 0 \
--cluster_output '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds' \
--cluster 3 \
--pathways '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds' \
--prefix 'raw_counts_pbta_subset_km_euclidean_0_3' \
--output_dir '~/Projects/Test/pbta_splicing/output'

# pbta splicing unbiased
Rscript code/diff-pathways-with-clusters.R \
--mat '~/Projects/Test/pbta_splicing/data/raw_counts_pbta_subset.rds' \
--type 'raw_counts' \
--var_genes 0 \
--cluster_output '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds' \
--cluster 3 \
--pathways '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds' \
--prefix 'raw_counts_pbta_subset_km_euclidean_0_3' \
--output_dir '~/Projects/Test/pbta_splicing/output'
