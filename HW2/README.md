###HW2
####Andrew Teng
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
printf 'NM_007294\n%.0s' {1..23} > col4.txt
printf '0\n%.0s' {1..23} > col5.txt
printf ' -\n%.0s' {1..23} > col6.txt
awk '{$1=$1}1' col6.txt > col6-1.txt

paste col1.txt col2-1.txt col3-1.txt col4.txt col5.txt col6-1.txt> NM007294_exome.bed
```
Create FASTA file
```{sh}
bedtools getfasta -fi resources/genome/hg19.fa -bed NM007294_exome.bed -fo NM007294.fa
```
