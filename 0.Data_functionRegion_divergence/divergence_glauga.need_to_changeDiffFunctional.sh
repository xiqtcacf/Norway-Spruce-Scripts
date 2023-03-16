#!/bin/bash
#SBATCH -A snic2022-5-68
#SBATCH -o %JfixDiff.out
#SBATCH -e %JfixDiff.err
#SBATCH -J %JfixDiff.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL


module load bioinfo-tools
module load htslib
module load vcftools

batch=$1 ###functional region

input=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/sweden_norway.noPab002Pab034.$batch.vcf.gz
output=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/divergence_glauga

###step0.generate GT:1/1 from P.glauga and further generate chr position for each functional region
#0/0 : homozygous for the reference allele
#0/1 : heterozygous (one ref allele, one alt allele)
#1/1:  homozygous for the alternate allele

###step1.generate chr position from P.abies for each functional region
zless $input | grep -v '#' | cut -f1,2 > $output/chr.position.25inds.sweden_norway.noPab002Pab034.$batch
cat header1.txt $output/chr.position.25inds.sweden_norway.noPab002Pab034.$batch > t
mv t $output/chr.position.25inds.sweden_norway.noPab002Pab034.$batch

###step2.get the fixed differences 1/1 for glauga - (1/1 for glauga overlaping with abies)
Rscript divergence_glauga.R $output/function_abies $output/function_glauga
mv fixed_difference_function.glauga.txt fixed_difference_${batch}.glauga.txt
cut -f1 fixed_difference_${batch}.glauga.txt | sort | uniq -c > fixed_difference_${batch}.glauga.scaffolds
sed 's/      //g' fixed_difference_${batch}.glauga.scaffolds -i
cat header2.txt fixed_difference_${batch}.glauga.scaffolds > tt
mv tt fixed_difference_${batch}.glauga.scaffolds
