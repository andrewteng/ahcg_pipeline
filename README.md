##AHCG Pipeline
Andrew Teng || BIOL 8803F || Variant calling pipeline for genomic data analysis    

---
####Software Requirements
```{sh}
1. Python3 - version 3.4.1
2. Trimmomatic - version 0.36
3. Bowtie2 - version 2.2.9
4. Picard tools - version 2.6.0
5. GATK - version 3.4
6. Samtools- version 0.1.19
7. Java- version 1.8
```
####Reference sequence: UCSC (hg19)
```{sh}
wget ftp://ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz  
tar -xvzf www.prism.gatech.edu/~sravishankar9/resources.tar.gz
```

####Test data: NIST NA12878 sample

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

####Add to path 
```{sh}
export PATH="$PATH:/home/basespace/projects/ahcg_pipeline/lib/bowtie2-2.2.9"
export PATH="$PATH:/home/basespace/projects/ahcg_pipeline/lib/Trimmomatic-0.36"
```

####Build bowtie index 
```{sh}
bowtie2-build -f resources/genome/hg19.fa bt2_base
```

####Build fasta index
```{sh}
samtools faidx resources/genome/ref.fa
```

####Build genome dict file 
```{sh}
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict
```

####Run script
```{sh}
python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p  
lib/picard.jar -g lib/GenomeAnalysisTK.jar -i fastq/*fastq -w bowtie_index/hg19 -d resources/dbsnp/dbsnp_138.hg19.vcf  
-a lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -r resources/genome/hg19.fa -o ./
```
---
###Progress Updates
---
####08/30/2016: Setting up github repository.
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
---
####09/01/2016: Extraction of desired gene and BED file formation
Downloaded gene coordinates file for hg19  
Used `grep` to extract specific gene
```{sh}
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
grep -w "BRCA1" hg19_refGenome.txt > NM_007294.txt 
grep -w "NM_007294" NM_007294.txt > NM_007294.txt
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
Create FASTA file to extract sequences from original files
```{sh}
bedtools getfasta -fi resources/genome/hg19.fa -bed NM007294_exome.bed -fo NM007294.fa
```
---
####09/06/2016: Variant Calling
Download BAM files for NUST NA12878 genome
```{sh}
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bamedb7bba8479cf224bf3015fdfda44f39ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.baieaaad4ad3400ab03cb54fa1f898134de
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam90d7a35bd59971c44f528427a0b2da45ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bai3937b9d067979cfa74f1f8dd717e52b5
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bam1246c31ecfe53e9f55bb4890d16ebb9aftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bai1c9437d4ada3a5c8278c46cc2654b354
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.bam08f63aba86cad1cde5ace41b602cb347ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.baiacb06b877735a4bed4b310d7f08eecfa
```
Merge BAM files together
```{sh}
samtools merge <OUTPUT BAM FILE> *.bam
```
Subset NA12878 sample using samtools
```{sh}
samtools view -L <BED FILE> -b -o <OUTPUT BAM FILE> <INPUT BAM FILE>
```
Convert BAM to FASTQ
```{sh}
bedtools bamtofastq -i <BAM FILE> -fq <FASTQ R1> -fq2 <FASTQ R2>
```
---
####09/08/2016: re-ran previous commands (yielded an empty VCF file) 
---
####09/13/2016: re-ran pipeline using new reference file
Download new reference file 
```{sh}
wget http://vannberg.biology.gatech.edu/data/chr17.fa
```
Re-ran entire pipeline  
Download gold standard VCF files from Illumina and Genome in a Bottle (GIAB)
```{sh}
Illumina: ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz   
Genome in a Bottle: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
```
Compare variants in the VCF file
```{sh}
vcftools --vcf <file> --diff <file>
```
---
####09/15/2016: Finding variants
```{sh}
bgzip <VCF file>
tabix -p vcf <VCF file>
```
TODO:  
- Sliding window read-depth coverage.  
	```{sh}
	java -jar GenomeAnalysisTK.jar -T DepthOfCoverage -R reference.fasta -I input_bams.list
	```
- For list of bam files, a text file needs to be creates with each line being a path to a unique bam file.  
- R script for venn diagram visualization.  

---
####09/20/2016: Gene List for Breast and Ovarian Cancers
Common genes in both the Color Genetics paper and Otogenetics Breast Cancer Panel Gene List 
```{sh}
Color Genetics Paper: https://s3.amazonaws.com/color-static-prod/pdfs/validationWhitePaper.pdf
Otogenetics Gene List: http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf

```

|Gene|NCBI Reference|OMIM|# of Exons|Source|
|----|--------------|----|----|---|
|AR|NM_000044.3|313700|8|OTO|
|ATM|NM_000051.3|607585|68|Both|
|BARD1|NM_000465.3|601593|11|Both|
|BRCA1|NM_007298.3|113705|22|Both|
|BRCA2|NM_000059.3|600185|27|Both|
|BRIP1|NM_032043.2|605882|20|Both|
|CASP8|NM_001080124.1|601763|9|OTO|
|CDH1|NM_004360.3|192090|16|Both|
|CHEK2|NM_001005735.1|604373|16|Both|
|DIRAS3|NM_004675.2|605193|2|OTO|
|ERBB2|NM_001005862.1|164870|30|OTO|
|EPCAM|NM_002354.2|185535|9|Color|
|MLH1|NM_000249.3|120436|21|Color|
|MSH2|NM_000251.2|609309|21|Color|
|MSH6|NM_000179.2|600678|12|Color|
|NBN|NM_002485.4|602667|16|Both|
|PALB2|NM_024675.3|601355|13|Both|
|PTEN|NM_000314.4|601728|9|Both|
|PMS2|NM_000535.6|600259|16|Color|
|RAD50|NM_005732.3|604040|25|OTO|
|RAD51|NM_001164269.1|179617|10|Both|
|STK11|NM_000455.4|602216|10|Both|
|TGFFB1|NM_000660.4|190180|7|OTO|
|TP53|NM_000546.5|191770|11|Both|
---
####09/22/2016:

---
