# CLK1 Specific

Module authors: Ammar Naqvi (@naqvia)

The purpose of this module is to identify differential splicing and expression as a result from increased CLK1 Exon 4 inclusion (high exon 4 vs low exon 4 tumors)

## Usage
### Run scripts on rMATS output to generate differential splicing and gene expression volcano plots:
<br>**Run shell script to make plots**
```
./run_module.sh
```

Input files:
```
input/dca735c2-6e0e-4239-8a68-10c6d2aa9015.CLK1_EI_vs_CLK1_ES.non_denovo.SE.MATS.JC.txt
```


## Folder content
* `run_module.sh` takes the files from above and generates differential splicing and expression plots
* `01-plot_CLK1_EI_vs_ES_PSI_volcano.R` processing output from rMATS and generates volcano plot deposited in `plots/`
* `02-plot_highExon4_vs_lowExon4_diffExpr_volcano.R` generates volcano plot from count table between low vs high exon 4 inclusion tumors and outputs plot in `plots/`
