#!/bin/bash
#SBATCH -A naiss2023-5-36
#SBATCH -o %J_NoSelection.out
#SBATCH -e %J_NoSelection.err
#SBATCH -J %J_NoSelection.job
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

batch=$1 ###neutral randomly scaffolds

input=/proj/uppstore2017066/processedData/analysis/3_natural_selection/vcf_selection
output=/proj/uppstore2017066/processedData/analysis/3_natural_selection/vcf_selection

fasta=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai
scaffold_no_selection=/proj/uppstore2017066/processedData/analysis/3_natural_selection/vcf_no_selection/background_5000Scaff.$batch.bed

### generate vcf no selection for sweden-norway population
vcftools --gzvcf $input/sweden_norway.vcf.gz --bed $scaffold_no_selection --recode --out $output/background.5000Scaff.$batch.sweden_norway
mv $output/background.5000Scaff.$batch.sweden_norway.recode.vcf $output/background.5000Scaff.$batch.sweden_norway.vcf
bgzip $output/background.5000Scaff.$batch.sweden_norway.vcf && tabix -p vcf $output/background.5000Scaff.$batch.sweden_norway.vcf.gz

### calculate statistics for neutral regions
$angsd -vcf-PL $output/background.5000Scaff.$batch.sweden_norway.vcf.gz -out $output/background.5000Scaff.$batch.sweden_norway -doSaf 1 -doMajorMinor 1 -anc $fasta -nInd 25 -minQ 30
$realSFS $output/background.5000Scaff.$batch.sweden_norway.saf.idx -fold 1 > $output/background.5000Scaff.$batch.sweden_norway.sfs
$angsd -vcf-PL $output/background.5000Scaff.$batch.sweden_norway.vcf.gz -out $output/background.5000Scaff.$batch.sweden_norway -doThetas 1 -doSaf 1 -pest $output/background.5000Scaff.$batch.sweden_norway.sfs -anc $fasta -nInd 25  -minQ 30
$thetaStat do_stat $output/background.5000Scaff.$batch.sweden_norway.thetas.idx -outnames $output/background.5000Scaff.$batch.sweden_norway.theta.gz
