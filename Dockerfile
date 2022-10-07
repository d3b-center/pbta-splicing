FROM rocker/tidyverse:4.2
MAINTAINER="Jo Lynne Rokita (rokita@chop.edu)"
WORKDIR /rocker-build/

RUN RSPM="https://packagemanager.rstudio.com/cran/2022-09-06" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site
  
COPY scripts/install_bioc.r .

COPY scripts/install_github.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	curl \
	cpanminus \
	bzip2 \
	zlib1g \
	libreadline-dev \
  	build-essential \
	libxt-dev \
	libproj-dev \
	libv8-dev \
	cpanminus \
	libgdal-dev \
	libgmp-dev \
	libmpfr-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk \
   libxt6

# install R packages from CRAN
RUN install2.r \
	BiocManager \
	pheatmap \
	optparse \
	hrbrthemes \
	viridis \
	plyr \
	ggstatsplot \
	diptest \
	ggthemes \
	UpSetR \
	cowplot \
	grid
	
# install R packages from GitHub
RUN ./install_github.r \
	PoisonAlien/maftools

# install R packages from BioC
RUN ./install_bioc.r \
	ConsensusClusterPlus \
	sva \
	EnhancedVolcano \
	DESeq2 \
	fgsea \
	GSVA \
	limma \
	pheatmap \
	Biobase \
	DGCA

# install perl packages
RUN cpanm install Statistics::Lite

ADD Dockerfile .
