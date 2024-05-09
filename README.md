# pbta-splicing
Aberrant splicing in brain tumors

## Docker set-up

### docker pull and run
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:latest
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:latest
```
### docker execute
```
docker exec -ti pbta-splicing bash
```
### Get rMATS data files

<br>**Run shell script to get merged rMATS result tables**
```
bash download_data.sh
```
