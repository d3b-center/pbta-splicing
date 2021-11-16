# pbta-splicing
Aberrant splicing in brain tumors

## Docker set-up

### docker pull and run
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:0.1
docker run --name pbta-splicing -d -e PASSWORD=pass -p 8787:8787 -v "$PWD":/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:0.1

### docker execute
docker exec -ti pbta-splicing bash
