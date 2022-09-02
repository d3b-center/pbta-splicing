# pbta-splicing
Aberrant splicing in brain tumors

## Docker set-up

### docker pull and run
```
docker pull jrokita1/pbta-splicing:version1.1
docker run --name test -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing pbta-splicing:version1.1
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
