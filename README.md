# pbta-splicing
Aberrant splicing in brain tumors

## Docker set-up

### docker pull and run
```
docker pull pgc-images.sbgenomics.com/naqvia/pbta-splicing:pbta_splicing
docker run --platform linux/amd64 --name splicing -d -v $PWD:/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/naqvia/pbta-splicing:pbta_splicing
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
