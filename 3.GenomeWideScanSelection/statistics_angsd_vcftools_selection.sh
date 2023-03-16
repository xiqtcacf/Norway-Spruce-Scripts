#!/bin/bash
#SBATCH -A naiss2023-5-36
#SBATCH -o %J_Selection.out
#SBATCH -e %J_Selection.err
#SBATCH -J %J_Selection.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load vcftools/0.1.15
module load bcftools/1.9
module load python/2.7.15
module load gsl
angsd=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/angsd
realSFS=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/realSFS
thetaStat=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/thetaStat

batch=$1 ###scaffolds under positive and balancing selection

input=/proj/uppstore2017066/processedData/analysis/3_natural_selection/vcf_selection
output=/proj/uppstore2017066/processedData/analysis/3_natural_selection/vcf_selection
fasta=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai

### generate vcf under selection for sweden-norway population
vcftools --gzvcf $input/sweden_norway.vcf.gz --bed $bed --recode --out $output/top0.1.$batch.sweden_norway.bed
mv $output/top0.1.$batch.sweden_norway.bed.recode.vcf $output/top0.1.$batch.sweden_norway.bed.vcf
bgzip $output/top0.1.$batch.sweden_norway.bed.vcf && tabix -p vcf $output/top0.1.$batch.sweden_norway.bed.vcf.gz

### calculate statistics  by ANGSD 25inds sweden-norway population
$angsd -vcf-PL $output/top0.1.$batch.sweden_norway.bed.vcf.gz -out $output/top0.1.$batch.sweden_norway.bed -doSaf 1 -doMajorMinor 1 -anc $fasta -nInd 25 -minQ 30
$realSFS $output/top0.1.$batch.sweden_norway.bed.saf.idx -fold 1 > $output/top0.1.$batch.sweden_norway.bed.sfs
$angsd -vcf-PL $input/top0.1.$batch.sweden_norway.bed.vcf.gz -out $output/top0.1.$batch.sweden_norway.bed -doThetas 1 -doSaf 1 -pest $output/top0.1.$batch.sweden_norway.bed.sfs -anc $fasta -nInd 25 -minQ 30
$thetaStat do_stat $output/top0.1.$batch.sweden_norway.bed.thetas.idx -outnames $output/top0.1.$batch.sweden_norway.bed.theta.gz
