#!/bin/bash

#SBATCH -A snic2019-3-555
#SBATCH -o %J_RAiSD_subset.out
#SBATCH -e %J_RAiSD_subset.err
#SBATCH -J %J_RAiSD_subset.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6-00:00:00
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

batch=$1
inputFile=/proj/uppstore2017066/processedData/analysis/3_natural_selection/RAiSD/sweden_norway.noPab002Pab034.fill_tags.vcf.gz
output=/proj/uppstore2017066/processedData/analysis/3_natural_selection/RAiSD/subset
RAiSD=/proj/uppstore2017066/tools/RAiSD/raisd-master/RAiSD

export DIR=/proj/uppstore2017066/processedData/analysis/3_natural_selection/RAiSD/subset/RAiSD_result/$batch
cd $DIR

file="/proj/uppstore2017066/processedData/analysis/3_natural_selection/RAiSD/subset/scaffold_name/$batch"
name=$(cat "$file")

for i in $name
do

tabix -h $inputFile $i > $output/temp/$i.vcf
var=$(grep -v "#" $output/temp/$i.vcf | wc -l)

	if [ $var -gt 49 ];
	then
		$RAiSD -n swe_nor25.$i -I $output/temp/$i.vcf -t -R -y 2
		rm $output/temp/$i.vcf
		rm RAiSD_Info.swe_nor25.$i;
	else
                rm $output/temp/$i.vcf
	fi

done
