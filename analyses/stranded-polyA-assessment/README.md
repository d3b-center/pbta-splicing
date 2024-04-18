# Assess impact of PSI between polyA and stranded library types

This module correlates polyA and stranded PSIs from patient samples known to have both library types.

## Usage
`bash run_module`

## Folder content 
`01-plot-str-vs-polyA.Rmd` runs correlation analysis on samples with both stranded and polyA RNA-seq and generates scatter plots.

```
.
├── 01-plot-str-vs-polyA.Rmd
├── 01-plot-str-vs-polyA.html
├── README.md
├── plots
│   ├── PT_RYMG3M91_polyA_v_stranded_psi.pdf
│   └── PT_W5GP3F6B_polyA_v_stranded_psi.pdf
└── run_module.sh
```