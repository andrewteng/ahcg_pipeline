##AHCG Pipeline
Andrew Teng-BIOL 8803F    
Variant calling pipeline for genomic data analysis  
  

###Software Requirements

1. Python3 - version 3.4.1
2. Trimmomatic - version 0.36
3. Bowtie2 - version 2.2.9
4. Picard tools - version 2.6.0
5. GATK - version 3.4
6. Samtools- version 1.x
7. Java- version 1.8

####Reference sequence: UCSC hg19.fa

####Test data: NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

###Add to path 
```{sh}
export PATH="$PATH:/home/basespace/projects/ahcg_pipeline/lib/bowtie2-2.2.9"
export PATH="$PATH:/home/basespace/projects/ahcg_pipeline/lib/Trimmomatic-0.36"
```

###Build bowtie index 
```{sh}
bowtie2-build -f resources/genome/hg19.fa bt2_base
```

###Build fasta index
```{sh}
samtools faidx resources/genome/ref.fa
```

###Build genome dict file 
```{sh}
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict
```

###Run script
```{sh}
python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p  
lib/picard.jar -g lib/GenomeAnalysisTK.jar -i fastq/*fastq -w bowtie_index/hg19 -d resources/dbsnp/dbsnp_138.hg19.vcf  
-a lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -r resources/genome/hg19.fa -o ./
```
###----
###Progress Updates
####08/30/16
Added directories to .gitignore and changed path in config file
```{sh}
.gitignore
vi .git/config
```
Pushed output vcf file, AHCG pipeline, and README file
```{sh}
git add <FILENAME>
git commit <FILENAME> -m <MESSAGE>
git push origin master
```
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
