cat input/diffSplicing_cand.filterDiffExpr.txt | awk '{print $2}' | perl extract_psi_from_list.pl
ls results/results_diff/*txt | xargs -n 1 echo "perl quant_expr_vs_psi_PCT.pl ~/Desktop/AS-DMG/data/stranded_trans_rsem_counts_tab.tsv ~/Desktop/AS-DMG/data/stranded_trans_rsem_counts_tab.tsv " | bash
perl corr_calc.pl | grep PSI | grep -v ^PSI > results/expr_vs_psi_corr_res.txt
