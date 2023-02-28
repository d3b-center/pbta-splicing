FROM rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/

RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY scripts/install_bioc.r .

COPY scripts/install_github.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install dev libraries and curl
RUN apt update && apt install -y \
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
	zlib1g-dev \

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6

# install R packages
RUN ./install_bioc.r \
	Biobase \
	BiocManager \
	broom \
	ConsensusClusterPlus \
	corrplot \
  cowplot \
  DGCA \
	DESeq2 \
	diptest \
	EnhancedVolcano \
	fgsea \
	ggpubr \
	ggstatsplot \
  ggthemes \
  gridExtra \
	GSVA \
	hrbrthemes \
	limma \
	optparse \
	pheatmap \
  reshape2 \
  sva \
  survminer \
  UpSetR

# install R packages from GitHub
RUN ./install_github.r \
	PoisonAlien/maftools

# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
