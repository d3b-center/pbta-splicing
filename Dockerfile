FROM rocker/tidyverse:4.4.0
LABEL maintainer = "Ammar S. Naqvi (naqvia@chop.edu)"
WORKDIR /rocker-build/


### Install apt-getable packages to start
#########################################
RUN apt-get -y update && apt-get install -y --no-install-recommends

# Install dev libraries and curl
RUN apt install -y  \
  bedtools \
	build-essential \
	bzip2 \
	cpanminus \
	cmake \
	curl \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgdal-dev \
  ibglpk40 \
  libglpk-dev \
	libgmp-dev \
	liblzma-dev \
	libmpfr-dev \
	libncurses5-dev \
	libproj-dev \
	libreadline-dev \
	libssl-dev \
	libv8-dev \
	libxt-dev \
	libudunits2-dev \
	zlib1g-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

# Install packages
RUN R -e 'BiocManager::install(c( \
  "Biobase", \
  "broom", \
  "circlize", \
  "COINr", \
  "clusterProfiler", \
  "ComplexHeatmap", \
  "ConsensusClusterPlus", \
  "corrplot", \
  "cowplot", \
  "DGCA", \
  "DESeq2", \
  "DOSE", \
  "diptest", \
  "edgeR", \
  "EnhancedVolcano", \
  "factoextra", \
  "fgsea", \
  "fpc", \
  "ggpubr", \
  "ggstatsplot", \
  "ggthemes", \
  "ggVennDiagram", \
  "GSVA", \
  "gridExtra", \
  "Hmisc", \
  "hrbrthemes", \
  "limma", \
  "lspline", \
  "maftools",\
  "msigdbr", \
  "optparse", \
  "org.Hs.eg.db", \
  "PMCMRplus", \
  "pheatmap", \
  "reshape2", \
  "rstatix", \
  "rtracklayer", \
  "R.utils", \
  "sva", \
  "survival", \
  "survminer", \
  "UpSetR" \
  "VennDiagram" \
))'


## install GitHub packages
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('d3b-center/annoFuseData', ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66', dependencies = TRUE)"

# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
