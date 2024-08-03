# The splicing modulator CLK1 is a candidate oncogenic dependency in pediatric high-grade gliomas

Ammar S. Naqvi, Ryan J. Corbett, Priyanka Seghal, Karina L. Conkrite, Komal S. Rathi, Brian M. Ennis, Katharina E Hayer, Bo Zhang, Miguel A. Brown, Daniel P. Miller, Adam A. Kraya, Joseph M. Dybas, Zhuangzhuang Geng, Christopher Blackden, Shehbeel Arif, Antonia Chroni, Aditya Lahiri, Madison L. Hollawell, Phillip B. Storm, Jessica B. Foster, Mateusz Koptyra, Peter J. Madsen, Sharon J. Diskin, Andrei Thomas-Tikhonenko, Adam C. Resnick, Jo Lynne Rokita

## Docker set-up

### docker pull and run
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:v1.0.0
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-splicing pgc-images.sbgenomics.com/d3b-bixu/pbta-splicing:v1.0.0
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

### Contact
For questions, please submit an issue or send an email to Jo Lynne Rokita (@jharenza): rokita@chop.edu and Ammar S. Naqvi (@naqvia): naqvia@chop.edu.
