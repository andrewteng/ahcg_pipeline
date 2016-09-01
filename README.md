### AppliedHumanCompGenomics

Code for BIOL8803F-Applied Human Computational Genomics

###Progress Updates 

####08/30/16
See HW1 README

####09/01/16
Created fasta file for BRCA1 NM_007294
```{sh}
sudo apt-get install samtools
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
grep -w "NM_007294" hg19_refGenome.txt > NM_007294.txt 
````
Create a bed file from the reference (not the best way)
```{sh}
awk '{print $10}' NM_007294.txt > col2.txt
awk '{print $11}' NM_007294.txt > col3.txt
```
Extracting start and stop positions
```{sh}
tr , '\n' < col2.txt > col2-1.txt
tr , '\n' < col3.txt > col3-1.txt

printf 'chr17\n%.0s' {1..23} > col1.txt

paste col1.txt col2-1.txt col3-1.txt > NM007294_exome.txt
```
