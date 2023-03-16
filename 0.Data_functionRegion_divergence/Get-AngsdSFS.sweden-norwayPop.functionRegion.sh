#!/bin/bash
#SBATCH -A naiss2023-5-36
#SBATCH -o %JSFS.out
#SBATCH -e %JSFS.err
#SBATCH -J %JSFS.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL

batch=$1 ###functional region bed

module load bioinfo-tools
module load vcftools/0.1.15
module load perl
angsd=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/angsd
realSFS=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/realSFS

All_data=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/all.35inds.vcf.gz
pop=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/25inds.sweden-norway.txt
fasta=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai

output=/home/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data

# may need to change how many cores to use regarding size of file
###get sweden-norway pop data from all inds
vcftools --gzvcf $All_data --keep $pop --recode --out $output/sweden_norway.noPab002Pab034
mv $output/sweden_norway.noPab002Pab034.recode.vcf $output/sweden_norway.noPab002Pab034.vcf
bgzip $output/sweden_norway.noPab002Pab034.vcf && tabix -p vcf $output/sweden_norway.noPab002Pab034.vcf.gz

#### generate function region:4fold, 0fold, introns, intergenic, conserved, promoters
vcftools --gzvcf $output/sweden_norway.noPab002Pab034.vcf.gz --bed $batch --recode --out $output/sweden_norway.noPab002Pab034.$batch
mv $output/sweden_norway.noPab002Pab034.$batch.recode.vcf $output/sweden_norway.noPab002Pab034.$batch.vcf
bgzip $output/sweden_norway.noPab002Pab034.$batch.vcf && tabix -p vcf $output/sweden_norway.noPab002Pab034.$batch.vcf.gz

#########little test if we want to downsample snps by randomly selecting snp
#perl $output/thin_vcf.pl $output/sweden_norway.noPab002Pab034.$batch.vcf.gz 0.5
#gunzip -d $output/sweden_norway.noPab002Pab034.$batch.thin0.5.vcf.gz
#bgzip $output/sweden_norway.noPab002Pab034.$batch.thin0.5.vcf && tabix -p vcf $output/sweden_norway.noPab002Pab034.$batch.thin0.5.vcf.gz

### get sfs for each function_region
$angsd -vcf-PL $output/sweden_norway.noPab002Pab034.$batch.vcf.gz -out $output/sweden_norway.noPab002Pab034.$batch -doSaf 1 -doMajorMinor 1 -anc $fasta -nInd 25 -minQ 30
$realSFS $output/sweden_norway.noPab002Pab034.$batch.saf.idx -cores 1 -fold 1 > $output/sweden_norway.noPab002Pab034.$batch.saf.idx.sfs



