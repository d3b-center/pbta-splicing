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

# plot CLK1-201 expression in normals
Rscript --vanilla 07-plot-clk1ex4-hgg-normals.R   

## run correlation analyses for CLK1 and NF1 transcripts
#Rscript -e "rmarkdown::render('08-CLK1-impact-NF1-splicing.Rmd', clean = TRUE)" 

# run CLK1-NF1 protein/phosphoprotein expression correlation script
#Rscript --vanilla 09-clk1-nf1-protein-correlations.R

# plot DMG clk1/nf1 rna, protein, splicing z-scores
#Rscript --vanilla 10-clk1-nf1-single-sample-heatmap.R

rm ./Rplots.pdf
