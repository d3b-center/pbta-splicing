## subset rMATs for CLK1
Rscript --vanilla 00-subset-CLK1-from-rMATs.R

## plot HGGs with high vs low exon 4
Rscript --vanilla 01-plot_highExon4_vs_lowExon4_and_SBI.R

## plot correlations of splicing vs expr
Rscript --vanilla 02-plot_splicing_vs_expr.R
Rscript --vanilla 03-plot_SR-phosp_vs_CLK1-RNA.R

## plot CLK1 stacked barplots
Rscript --vanilla 04-CLK1_PSI_plots.R

# run CLK1-SRSF RNA expression correlation script
Rscript --vanilla 05-CLK-SRSF-expr-correlations.R

# run CLK1-SRSF protein/phosphoprotein expression correlation script
Rscript --vanilla 06-CLK1-psi-expr-SRSF-expr-prot-phospho-heatmap.R