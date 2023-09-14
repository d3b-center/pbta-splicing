## Format DESeq2 results for GSEA script

echo "gene      baseMean        log2FoldChange  lfcSE   stat    pvalue  padj    Significant" > results/ctrl_vs_treated.de.formatted.ts
cat results/ctrl_vs_treated.de.tsv|  awk -F "\t" '{print $8"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | awk -F "_" '{print $2}'  >> results/ctrl_vs_treated.de.formatted.tsv
