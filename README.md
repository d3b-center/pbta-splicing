# pbta-splicing
Aberrant splicing in brain tumors

## Docker set-up

### docker pull and run
```
docker pull jrokita1/pbta-splicing:latest
docker run --name test -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing jrokita1/pbta-splicing:latest
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
