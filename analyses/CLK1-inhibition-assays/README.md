# CLK1 inhibition assays
This module plots CLK1 inhibition assays (cell proliferation and viability) using the Cirtuvivint inhibitor. This will allow us to assess CLK1 impact on proliferation and viability in KNS42 cells.
## Usage
To run analysis module from the command line:
```
bash run_module.sh
```

## Folder Content
```01-plot_cell-proliferation-assay.R```  plots the results from the cell proliferation assay comparing ctrl vs inhibition <br>
```02-plot_cell-proliferation-assay-res.R``` plots the results from the cell viability assay (CTG) comparing ctrl vs inhibition <br>

## Directory Structure
```
.
├── 01-plot_cell-proliferation-assay.R
├── 02-plot_ctg-assay.R
├── input
│   ├── 2024.05.30_KNS-42_3D_cirtuvivint.txt
│   ├── 2024.05.30_KNS-42_3D_cirtuvivint_nogroups.txt
│   ├── 2024.05.30_KNS-42_6D_cirtuvivint.txt
│   └── 2024.05.30_KNS-42_6D_cirtuvivint_nogroups.txt
├── plots
│   ├── cell-prolif.pdf
│   └── cell_viability-barplot.pdf
└── run_module.sh
```
