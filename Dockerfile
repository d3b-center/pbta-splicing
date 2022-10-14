FROM rocker/tidyverse:4.2
MAINTAINER rokita@chop.edu
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
	build-essential \
	bzip2 \
	cpanminus \
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
	zlib1g-dev 

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6

# install R packages from CRAN
RUN install2.r \
	BiocManager \
  corrplot \
  cowplot \
	DCGA \
	diptest \
	ggpubr \
  ggstatsplot \
	ggthemes \
  grid \
  gridExtra \
  hrbrthemes \
	optparse \
	pheatmap \
  reshape2 \
  sva \
  UpSetR \
	viridis 

# install R packages from GitHub
RUN ./install_github.r \
	PoisonAlien/maftools

# install R packages from BioC
RUN ./install_bioc.r \
	Biobase \
	ConsensusClusterPlus \
	DESeq2 \
	EnhancedVolcano \
	fgsea \
	GSVA \
	limma 

# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
