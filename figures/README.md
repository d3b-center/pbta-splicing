# Figures

## CLK1 Exon 4 splicing
### Input
gencode.v39.primary_assembly.annotation.CLK1-select.gtf

### Samples:
BS_HRJ9145M (high Ex4 inclusion)
BS_XM1AHBDJ (low Ex4 inclusion)

### Generate sashimi plot

```./ggsashimi.py -b examples/input_bams.tsv -g ../gencode.v39.primary_assembly.annotation.CLK1-select.gtf -c chr2:200853009-200864691 -M 10 -C 3 -O 3 --alpha 1 --base-size=10 -R 300 --height=1.5 --width=8 --ann-height 1 -P examples/palette.txt -o sashimi-CLK1-v2-plot --fix-y-scal --shrink```
