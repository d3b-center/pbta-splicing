# pbta-splicing
Aberrant splicing in pediatric brain tumors

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the manuscript noted above.

1. Clone the repository
```
git clone https://github.com/d3b-center/pbta-splicing.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/naqvia/pbta-splicing:latest
```

3. Start the docker container, from the `pbta-splicing` folder, run:
```
docker run --platform=linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/naqvia/pbta-splicing:latest
```

4. To execute shell within the docker image, from the `pbta-splicing` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the download_data.sh shell script to obtain latest data files:
```
bash download_data.sh
```

## Directory structure
```
.
├── analyses
│   ├── CLK1-splicing-impact
│   ├── CLK1_splicing_correlations
│   │   ├── input
│   │   ├── plots
│   │   └── util
│   ├── CLK1_splicing_impact
│   │   ├── input
│   │   ├── plots
│   │   └── results
│   ├── KNS42_cell-line
│   │   ├── input
│   │   └── plots
│   ├── clustering_analysis
│   │   ├── code
│   │   ├── input
│   │   ├── output
│   │   │   ├── ccp_output
│   │   │   ├── diff_genes
│   │   │   ├── diff_pathways
│   │   │   └── optimal_clustering
│   │   ├── plots
│   │   └── util
│   ├── cohort_summary
│   │   ├── input
│   │   ├── plots
│   │   └── results
│   ├── histology-specific-splicing
│   │   ├── plots
│   │   └── results
│   ├── long-read-CLK1-validation
│   ├── proteomics_correlation
│   │   └── plots
│   ├── splicing-factor_dysregulation
│   │   ├── input
│   │   ├── plots
│   │   └── results
│   ├── splicing_events_functional_sites
│   │   ├── input
│   │   ├── plots
│   │   └── results
│   │       └── archive
│   ├── splicing_index
│   │   ├── plots
│   │   └── results
│   └── survival
│       ├── input
│       ├── plots
│       │   ├── ATRT
│       │   ├── CPG
│       │   ├── EPN
│       │   ├── GNG
│       │   ├── HGG
│       │   ├── LGG
│       │   └── MB
│       ├── results
│       │   ├── ATRT
│       │   ├── CPG
│       │   ├── EPN
│       │   ├── GNG
│       │   ├── HGG
│       │   ├── LGG
│       │   └── MB
│       └── util
├── data
│   ├── v2
│   ├── v3
│   ├── v4
│   ├── v5
│   └── v6
├── doc
├── figures
├── palettes
├── scripts
├── tools
├── util
└── workflows
```

## Code Authors
Ammar Naqvi ([@naqvia](https://github.com/naqvia)) Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)) and
