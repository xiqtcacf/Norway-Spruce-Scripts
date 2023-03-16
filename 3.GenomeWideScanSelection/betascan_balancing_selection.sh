#!/bin/bash

#SBATCH -A snic2019-3-555
#SBATCH -o %J_betascan.out
#SBATCH -e %J_betascan.err
#SBATCH -J %J_betascan.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type ALL

set -e
set -x

module load bioinfo-tools
module load vcftools
module load htslib
module load BEDTools/2.27.1
module load python
module load samtools


inputFile=/proj/uppstore2017066/processedData/analysis/3_natural_selection/betascan/sweden_norway.noPab002Pab034.vcf.gz
output=/proj/uppstore2017066/processedData/analysis/3_natural_selection/betascan/subset_ref
glactools=/proj/uppstore2017066/tools/glactools/glactools
BetaScan=/proj/uppstore2017066/scripts/Analysis_whole_genome/3_natural-selection/betascan/BetaScan.py

batch=$1

file="/proj/uppstore2017066/processedData/analysis/3_natural_selection/betascan/subset_ref/scaffold_name/$batch"
name=$(cat "$file")

for i in $name
do

tabix -h $inputFile $i > $output/temp/$i.vcf
var=$(grep -v "##" $output/temp/$i.vcf | wc -l)

	if [ $var -gt 2 ];
	then
		$glactools vcfm2acf --onlyGT --fai $output/fasta/$batch.ref.fa.fai $output/temp/$i.vcf > $output/temp/$i.acf.gz
		$glactools acf2betascan --fold $output/temp/$i.acf.gz | gzip > $output/temp/$i.beta.txt.gz
		python2 $BetaScan -i $output/temp/$i.beta.txt.gz -fold -m 0.2 -o $output/betascore/$batch/$i.swe_nor25.betascores.txt
		rm $output/temp/$i.vcf $output/temp/$i.acf.gz $output/temp/$i.beta.txt.gz;
	else
                rm $output/temp/$i.vcf
	fi

done
