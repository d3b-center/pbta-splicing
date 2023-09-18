## Format DESeq2 results for GSEA script

echo -e "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tSignificant" > results/ctrl_vs_treated.de.formatted.tsv
cat results/ctrl_vs_treated.de.tsv|  awk -F "\t" '{print $8"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | awk -F "_" '{print $2}'  >> results/ctrl_vs_treated.de.formatted.tsv
