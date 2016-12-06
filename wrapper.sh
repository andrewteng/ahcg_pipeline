#!/bin/bash

while getopts "b:g:n:" options
do
	case $options in
		b) bamFile=$OPTARG;;
		g) geneList=$OPTARG;;
		n) nm_omim=$OPTARG;;
	esac
done

if [ ! -e "$bamFile" ]
	then printf "ERROR: Missing BAM file (-b flag)\n\n";
	exit;
fi
if [ ! -e "$geneList" ]
	then printf "ERROR: Missing gene list file (-g flag)\n\n";
	exit;
fi
if [ ! -e "$nm_omim" ]
	then printf "ERROR: Missing omim file with NM #'s (-n flag)\n\n";
	exit;
fi


printf "=============================================\n"
echo "AHCG Pipeline"
printf "=============================================\n\n"

echo -n "Starting"
sleep .5
echo -n "."
sleep .5
echo -n "."
sleep .5
echo -n "."
sleep .5 
echo -n "."
sleep .5

printf "\n(1/8) Converting BAM to FASTQ\n\n"
bedtools bamtofastq -i "$bamFile" -fq fastq/read1.fq -fq2 fastq/read2.fq

printf "(2/8) Finding variants\n\n" 
python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i fastq/*fq -w bowtie_index/hg19 -d resources/dbsnp/dbsnp_138.hg19.vcf -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -r resources/genome/hg19.fa -o .

printf "(3/8) Recalibrating variant file\n"
java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R resources/genome/hg19.fa -input variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 resources/hapmap_3.3.hg19.sites.vcf.gz -resource:omni,known=false,training=true,truth=false,prior=12.0 resources/1000G_omni2.5.hg19.sites.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10.0 resources/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 resources/dbsnp/dbsnp_138.hg19.vcf.gz -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -recalFile output.recal -tranchesFile output.tranches -rscriptFile output.plots.R

java -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R resources/genome/hg19.fa --input variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile output.recal -tranchesFile output.tranches -o recalibrated_variants.vcf

printf "\n(4/8) Creating BED file\n\n"
grep -Ff "$geneList" hg19_refGene.txt > temp.txt
grep -Ff "$nm_omim" temp.txt > temp1.txt
perl generatebed.pl > exomes.bed

printf "(5/8) Extracting variants\n\n"
bedtools intersect -a variants.vcf -b exomes.bed -header > extracted_variants.vcf

printf "(6/8) Matching variants with clinical risk\n\n"
bedtools intersect -a extracted_variants.vcf -b new_clinvar.vcf -header > extracted_clinvar.vcf

printf "(7/8) Calculating coverage\n\n"
samtools view -L dcm_gene_list.bed $bamFile -b > new.bam
bedtools genomecov -ibam new.bam -bga > coverage_output.bed
bedtools intersect -loj -split -a genes.bed -b coverage_output.bed >  cov.bed
awk '{printf("%s\t%s\t\%s\t%s\t%s\n", $1,$2,$3,$4,$10,$6)}' cov.bed > final_cov.bed

printf "(7/8) Generating report\n\n"
python3 parse_clnsig.py -i extracted_clinvar.vcf.gz 2>&1 | tee report.txt 
cut -c 24- report.txt

printf "DONE.'\n\n"

