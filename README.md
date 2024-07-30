# The splicing modulator CLK1 is a transcriptional dependency in pediatric high-grade gliomas
Ammar S. Naqvi, Ryan J. Corbett, Priyanka Seghal, Karina L. Conkrite, Komal S. Rathi, Brian M. Ennis, Katharina E Hayer, Bo Zhang, Miguel A. Brown, Daniel P. Miller, Adam A. Kraya, Joseph M. Dybas, Zhuangzhuang Geng, Christopher Blackden, Shebheel Arif, Antonia Chroni, Aditya Lahiri, Madison L. Hollawell, Phillip B. Storm, Jessica B. Foster, Matuesz Koptyra, Peter J. Madsen, Sharon J. Diskin, Andrei Thomas Tikhonenko, Adam C. Resnick, Jo Lynne Rokita 

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
```
bash download_data.sh
```
### Generate paper figures
```
bash scripts/run_code.sh
```
### Contact
For questions, please submit an issue or send an email to Ammar S. Naqvi (@naqvia): naqvia@chop.edu
