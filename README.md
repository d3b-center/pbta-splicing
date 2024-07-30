# pbta-splicing
The splicing modulator CLK1 is a transcriptional dependency in pediatric high-grade gliomas


## Docker set-up
### docker pull
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:v1.0.0
```
### docker run
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:v1.0.0
```
### docker execute
```
docker exec -ti pbta-splicing bash
```
### Get project data files
**Run shell script to get merged rMATS result tables**
```
bash download_data.sh
```
### Generate paper figures
```
bash scripts/run_code.sh
```
### Contact
For questions, please submit an issue or send an email to Ammar S. Naqvi (@naqvia): naqvia@chop.edu
