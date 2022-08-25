setwd('~/Projects/d3b-miRNA-analysis/')

# unbiased analysis using pbta splicing data
mat = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset.rds'
type = "non_expr"
var_genes = 0
cluster_output = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds'
cluster = 3
pathways = '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds'
prefix = "non_expr_pan_cancer_splice_subset_km_euclidean_0_3"
output_dir = "~/Projects/Test/pbta_splicing/output"

# pre-ranked analysis using pbta splicing data
mat = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset.rds'
type = "non_expr"
var_genes = 0
cluster_output = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds'
cluster = 3
pathways = '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds'
prefix = "non_expr_pan_cancer_splice_subset_km_euclidean_0_3"
output_dir = "~/Projects/Test/pbta_splicing/output"

# subset gene-expression to pbta splicing 
splice_dat = readRDS('~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset.rds')
dat = readRDS('../OpenPedCan-analysis/data/gene-counts-rsem-expected_count-collapsed.rds')
dat <- dat %>%
  dplyr::select(colnames(splice_dat))
common_genes <- intersect(rownames(splice_dat), rownames(dat))
dat_subset <- dat[rownames(dat) %in% common_genes,]
identical(colnames(splice_dat), colnames(dat_subset))
saveRDS(dat_subset, file = '../Test/pbta_splicing/data/raw_counts_pbta_subset.rds')

# pre-ranked analysis using pbta gene expression data
mat = '~/Projects/Test/pbta_splicing/data/raw_counts_pbta_subset.rds'
type = "raw_counts"
var_genes = 0
cluster_output = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds'
cluster = 3
pathways = '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds'
prefix = "raw_counts_pbta_subset_km_euclidean_0_3"
output_dir = "~/Projects/Test/pbta_splicing/output"

# unbiased analysis using pbta gene expression data
mat = '~/Projects/Test/pbta_splicing/data/raw_counts_pbta_subset.rds'
type = "raw_counts"
var_genes = 0
cluster_output = '~/Projects/Test/pbta_splicing/data/non_expr_pan_cancer_splice_subset_km_euclidean_0.rds'
cluster = 3
pathways = '~/Projects/Test/pbta_splicing/data/kegg_geneset_mrna.rds'
prefix = "raw_counts_pbta_subset_km_euclidean_0_3"
output_dir = "~/Projects/Test/pbta_splicing/output"
