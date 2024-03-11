FROM --platform=linux/amd64 rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/

RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY scripts/install_bioc.r .
COPY scripts/install_github.r .

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
	libgmp-dev \
	liblzma-dev \
	libmpfr-dev \
	libncurses5-dev \
	libproj-dev \
	libreadline-dev \
	libssl-dev \
	libv8-dev \
	libxt-dev \
	libudunits2-dev\
	zlib1g-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6


## install annoFuse
RUN ./install_github.r 'd3b-center/annoFuseData' --ref '321bc4f6db6e9a21358f0d09297142f6029ac7aa'

# install R packages
RUN ./install_bioc.r \
	Biobase \
	BiocManager \
	broom \
  circlize \
	COINr \
  clusterProfiler \
  ComplexHeatmap \
	ConsensusClusterPlus \
	corrplot \
  cowplot \
  DGCA \
	DESeq2 \
  DOSE \
	diptest \
	edgeR \
	EnhancedVolcano \
	factoextra \
	fgsea \
	fpc \
	ggpubr \
	ggstatsplot \
  ggthemes \
  ggVennDiagram \
  gridExtra \
	GSVA \
	hrbrthemes \
	limma \
	lspline \
  msigdbr \
	optparse \
  org.Hs.eg.db \
  PMCMRplus \
	patchwork \
	pheatmap \
  reshape2 \
  rstatix \
  rtracklayer \
  sva \
  survival \
  survminer \
  UpSetR

# install R packages from GitHub
RUN ./install_github.r \
	PoisonAlien/maftools

# Patchwork for plot compositions
RUN ./install_github.r  'thomasp85/patchwork' --ref 'c67c6603ba59dd46899f17197f9858bc5672e9f4'
RUN ./install_github.r 'clauswilke/colorblindr' --ref '90d64f8fc50bee7060be577f180ae019a9bbbb84'

# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
