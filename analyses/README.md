## Analysis Modules
This directory contains various analysis modules in the pbta-splicing-hgat project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance
The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
This is in service of documenting interdependent analyses.

Note that _nearly all_ modules use the harmonized clinical data file (`pbta-histologies.tsv`) even when it is not explicitly included in the table below.

| Module | Input Files | Brief Description | Output Files Consumed by Other Analyses |
|--------|-------|-------------------|--------------|
| [`pan_cancer`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/pan_cancer) | `pbta-histologies.RNA-Seq.initial.tsv` <br> `rmats_output.gz` <br> | Consensus clustering of samples into cluster for downstream analyses (eg. survival). <br> Generates splicing index table and plot to assess aberrant splicing (compared to n=9 healthy samples). <br> Onoplot generation to help visualize mutations, fusions and splicing. | N/A
| [`cnv-chrom-plot`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-chrom) | `pbta-cnv-consensus-gistic.zip` <br> `analyses/copy_number_consensus_call/results/pbta-cnv-consensus.seg` | Plots genome wide visualizations relating to copy number results | N/A
| [`cnv-comparison`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-comparison) | Earlier version of SEG files | *Deprecated*; compared earlier version of the CNV methods. | N/A
