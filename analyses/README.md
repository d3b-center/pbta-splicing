### Get rMATS data files
<br>**Run shell script to get merged rMATS result tables**
```
./download_data.sh
```
## Analysis Modules
This directory contains various analysis modules in the pbta-splicing-hgat project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance
The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
This is in service of documenting interdependent analyses.

Note that _nearly all_ modules use the harmonized clinical data file (`pbta-histologies.tsv`) even when it is not explicitly included in the table below.

| Module | Link| Brief Description |
|--------|-------|-------------------|
| [`psi_clustering`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/psi_clustering) | Consensus clustering of tumor splicing quantifications
| [`splicing burden index`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/splicing_index) | Compute splicing burden indices for each tumor sample across pediatric brain histologies
| [`splicing_events_functional_sites`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/splicing_events_functional_sites) | Identify signficiant aberrant splicing events that result in loss/gain of functional sites
