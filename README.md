##AHCG Pipeline
Andrew Teng || BIOL 8803F || Variant calling pipeline for genomic data analysis    

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
####Reference genome: UCSC (hg19)
```{sh}
wget ftp://ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz  
*used curated file due to size limits
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
Added directories to `.gitignore` and changed path in `config` file.
```{sh}
.gitignore
vi .git/config
```
Pushed output vcf file, AHCG pipeline, and README file.
```{sh}
git add <FILENAME>
git commit <FILENAME> -m <MESSAGE>
git push origin master
```
---
####09/01/2016: Extraction of desired gene and BED file formation.
Downloaded gene coordinates file for hg19.  
Used `grep` to extract specific gene.
```{sh}
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
grep -w "BRCA1" hg19_refGenome.txt > NM_007294.txt 
grep -w "NM_007294" NM_007294.txt > NM_007294.txt
````
Created a bed file from the reference (not the best way).
```{sh}
awk '{print $10}' NM_007294.txt > col2.txt
awk '{print $11}' NM_007294.txt > col3.txt
```
Extracted start and stop positions.
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
Created FASTA file to extract sequences from original files.
```{sh}
bedtools getfasta -fi resources/genome/hg19.fa -bed NM007294_exome.bed -fo NM007294.fa
```
---
####09/06/2016: Variant Calling.
Downloaded BAM files for NIST NA12878 genome.
```{sh}
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bamedb7bba8479cf224bf3015fdfda44f39ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.baieaaad4ad3400ab03cb54fa1f898134de
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam90d7a35bd59971c44f528427a0b2da45ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bai3937b9d067979cfa74f1f8dd717e52b5
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bam1246c31ecfe53e9f55bb4890d16ebb9aftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bai1c9437d4ada3a5c8278c46cc2654b354
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.bam08f63aba86cad1cde5ace41b602cb347ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.baiacb06b877735a4bed4b310d7f08eecfa
```
Merged BAM files together.
```{sh}
samtools merge <OUTPUT BAM FILE> *.bam
```
Subset NA12878 sample using samtools.
```{sh}
samtools view -L <BED FILE> -b -o <OUTPUT BAM FILE> <INPUT BAM FILE>
```
Converted BAM to FASTQ.
```{sh}
bedtools bamtofastq -i <BAM FILE> -fq <FASTQ R1> -fq2 <FASTQ R2>
```
---
####09/08/2016: Re-ran previous commands (yielded an empty VCF file). 
---
####09/13/2016: Re-ran pipeline using new reference file.
Downloaded new reference file. 
```{sh}
wget http://vannberg.biology.gatech.edu/data/chr17.fa
```
Re-ran entire pipeline.  
Downloaded gold standard VCF files from Illumina and Genome in a Bottle (GIAB).
```{sh}
Illumina: ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz   
Genome in a Bottle: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
```
Compared variants in the VCF file.
```{sh}
vcftools --vcf <FILE1> --diff <FILE2>
```
---
####09/15/2016: Finding variants.
```{sh}
bgzip <VCF file>
tabix -p vcf <VCF file>
```
Sliding window read-depth coverage.  
```{sh}
java -jar GenomeAnalysisTK.jar -T DepthOfCoverage -R reference.fasta -I input_bams.list
```
---
####09/20/2016: Generation of common genes in breast and ovarian cancers.
Common genes in both the Color Genetics paper and Otogenetics Breast Cancer Panel Gene List. 
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
|RAD51C|NM_058216.1|602774|NA|OTO|
|RAD51D|NM_001142571.1|602954|NA|Color
|STK11|NM_000455.4|602216|10|Both|
|TGFB1|NM_000660.4|190180|7|OTO|
|TP53|NM_000546.5|191770|11|Both|


---
####09/22/2016: Generation of BED file for new gene list.
Used `grep` to extract.
```{sh}
grep -Ff gene_list_breast_ovarian.txt hg19_refGene.txt > breast_ovarian.txt
grep -Ff omim_breast_ovarian.txt breast_ovarian.txt > breast_ovarian.txt
``` 
Created BED file and added 20 nucleotides (capture step) to each side for each exon.
```{sh}
perl generatebed.pl > exomes_breast_ovarian.bed
```
---
####09/27/2016: Creation of VCF file and extraction of variants.
Used commands from [09/06](https://github.com/andrewteng/ahcg_pipeline#09062016-variant-calling) to create two FASTQ read files.

Re-ran [pipeline](https://github.com/andrewteng/ahcg_pipeline#build-bowtie-index) to generate a final VCF file.

Extracted variants from VCF file using coordinates from BED file.
```{sh}
bedtools intersect -a variants.vcf -b exomes_breast_ovarian.bed > extractedVariants.vcf
```

---
####09/29/2016: Comparison of extracted variants with GIAB variants.
Obtained Genome in a Bottle (GIAB) variants.
```{sh}
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
gunzip *.gz
```
First column was incorrect, added 'chr' column.
```{sh}
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf > GIAB.vcf
```
Compared VCF files.
```{sh}
bedtools intersect -a <FILE1> -b <FILE2>
```
---
####10/04/2016: Used GATK VariantRecalibrator and ApplyRecalibration on variants file.
Downloaded new VCF file.
```{sh}
wget http://vannberglab.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf
```
Extracted variants and compared. 

```{sh}
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf.gz  
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz  
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz  
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz  
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz  
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz  
  
java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R resources/genome/hg19.fa -input variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 resources/hapmap_3.3.hg19.sites.vcf.gz -resource:omni,known=false,training=true,truth=false,prior=12.0 resources/1000G_omni2.5.hg19.sites.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10.0 resources/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 resources/dbsnp/dbsnp_138.hg19.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP -recalFile output.recal -tranchesFile output.tranches -rscriptFile output.plots.R
```

---
####10/06/2016: Reran previous commands, previously errored out. 
Created a `tabix` indexed VCF file.
```{sh}
bgzip <VCF FILE>  
tabix -p <GZ VCF FILE>  
```
Reran GATK `VariantRecalibrator` command.
Ran `ApplyRecalibration` command.

```{sh}
java -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R resources/genome/hg19.fa -input variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile output.recal -tranchesFile output.tranches -o recalibrated_variants.vcf
```

---

####10/11/2016: Fall Break :leaves: :fallen_leaf:

---

####10/13/2016: BRCA Exchange database comparison.

Utilized the data for BRCA1 and BRCA2 from the BRCA Exchange database.
```{sh}
wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv
wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv
```
The database uses hg36/37/38, while the reference previously used was hg19; therefore extraction based on hg36. Also, because database files include small indels, this function should be turned on.

Downloaded `compare.py` script that compares vcf and csv files.
```{sh}
pip install PyVCF
python compare.py <VCF FILE> <CSV FILE>
```

---
####10/18/2016: Coverage calculation and report.
Bedtools genome coverage (genomecov)

```{sh}
bedtools genomecov [OPTIONS] [-i|-ibam] -g (iff. -i) 
bedtooms intersect
```
---
####10/20/2016: Attended ASHG Conference.
Excused Absence.

---
####10/25-27/2016: Coverage (cont.)
Matched variants from VCF with clinical risk.
```{sh}
python compare_clin_with_vcf.py vcf/NA12878_variants.vcf BRCA1_brca_exchange_variants.csv BRCA2_brca_exchange_variants.csv > brca_clinical_xref.txt

grep -vi "Benign" brca_clinical_xref.txt > brca_clinical_nonbenign_xref.txt

python convert_clintxt_to_vcf.py brca_clinical_nonbenign_xref.txt brca_clinical_nonbenign_xref.vcf
```

Coverage calculator. 
```{sh}
grep 'NM_007298' bcoc_padded.bed > brca1.bed

samtools view -L brca1.bed data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > new.bam

bedtools genomecov -ibam new.bam -bga na12878.bga.bed

bedtools intersect -split -a brca1.bed -b na12878.bga.bed -bed > brca1.final.bed

awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$10,$6)}' brca1.final.bed > brca1.coverage_final.bed

bedtools intersect -a brca1.final.bed -b brca_clinical_nonbenign_xref.bed -wo > brca_clinical_nonbenign_final.bed

cat brca_clinical_nonbenign_final.bed | cut -f4,5,7,8,10
```
---
####11/01/2016
Repeated the same procedure for BRCA2.

---
####11/03/2016
Extended the pipeline for all the genes, as seen in the gene list.

---
####11/08/2016: New gene list for dilated cardiomyopathy (DCM).
|Gene|NCBI Reference|
|----|----|
|LMNA|NM_170707|
|MYBPC3|NM_000256|
|MYH7|NM_000257|
|MYH6|NM_002471|
|SCN5A|NM_198056|
|TNNT2|NM_001001430|

---
####11/10/2016: Lecture.
---
####11/15/2016: Generation of BED file for DCM.
```
grep -Ff dcm_gene_list.txt hg19_refGene.txt > dcm.txt
grep -Ff dcm_omim.txt dcm.txt > dcm1.txt
```
Create BED file and add 20 nucleotides to each side. 
```
perl generatebed.pl > exomes_dcm.bed
```

Download BAM files using `wget`.
```
Control1:
http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bai

Control2:
http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bai

Patient1:
http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai

Patient2:
http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bai

Patient3:
http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bai

Patient4:
http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bam
http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bai
```
Following previous steps from 09/06, merge the BAM files, and convert to 2 FASTQ files.

Using the 2 FASTQ files, rerun the entire pipeline and get a new variants VCF file. 

---
####11/17/2016

Using the generated VCF files from the controls and patients, GATK `VariantRecalibrator` was used. 

Reran pipeline on Patient 1.

---
####11/22/2016: Lecture.

---
####11/24/2016: Thanksgiving Break :poultry_leg: :meat_on_bone:

---
####11/29-12/01/2016: Created Master Script.
To use:
```
./wrapper.sh -b [BAM FILE] -g [GENE LIST] -n [OMIM FILE]
```

---
####12/06/2016: Finishing Touches.
