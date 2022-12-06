#!/bin/bash
#SBATCH -A snic2022-5-68
#SBATCH -o %J_AdaptiveRate.out
#SBATCH -e %J_AdaptiveRate.err
#SBATCH -J %J_AdaptiveRate.job
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 05:00:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load ANGSD/0.921
module load htslib
module load vcftools
module load gsl

batch=$1 #### recombination or gene_density

input_0fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/sweden_norway.noPab002Pab034.0fold.vcf.gz
input_4fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/sweden_norway.noPab002Pab034.4fold.vcf.gz
output=/home/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/2.adaptive_rate/bins10.$batch
fasta=/domus/h1/xiwa/xi/reference_genome/1k_reference/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/reference_genome/1k_reference/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai
dfe_fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16

export DIR=/home/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/2.adaptive_rate/bins10.$batch
cd $DIR

### first step bootstraps bed files and calculate sfs by ANGSD

for i in {1..10}
do
 vcftools --gzvcf $input_0fold --bed bins_bed/bin${i}.10.$batch.bed --recode --out bins_vcf_sfs/$i.$batch
 vcftools --gzvcf $input_4fold --bed bins_bed/bin${i}.10.$batch.bed --recode --out bins_vcf_sfs/$i.4fold
 angsd -vcf-gl bins_vcf_sfs/$i.$batch.recode.vcf -out bins_vcf_sfs/$i.$batch -doSaf 1 -anc $fasta -fai $fai -nInd 25 -fold 1 -P 2
 angsd -vcf-gl bins_vcf_sfs/$i.4fold.recode.vcf -out bins_vcf_sfs/$i.4fold -doSaf 1 -anc $fasta -fai $fai -nInd 25 -fold 1 -P 2
 realSFS bins_vcf_sfs/$i.$batch.saf.idx -P 2 > bins_vcf_sfs/$i.$batch.sfs
 realSFS bins_vcf_sfs/$i.4fold.saf.idx -P 2 > bins_vcf_sfs/$i.4fold.sfs 
# rm bins_vcf_sfs/$i.$batch.recode.vcf
# rm bootstrap_vcf_sfs/$i.4fold.recode.vcf
 rm bins_vcf_sfs/$i.$batch.arg
 rm bins_vcf_sfs/$i.$batch.saf.gz
 rm bins_vcf_sfs/$i.$batch.saf.idx
 rm bins_vcf_sfs/$i.$batch.saf.pos.gz
 rm bins_vcf_sfs/$i.4fold.arg
 rm bins_vcf_sfs/$i.4fold.saf.gz
 rm bins_vcf_sfs/$i.4fold.saf.idx
 rm bins_vcf_sfs/$i.4fold.saf.pos.gz

### second step, calculate dfe by dfe-alpha
input_dfe=/home/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/2.adaptive_rate/bins10.$batch
#output_dfe=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/$batch.4fold

est_dfe=/crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_dfe
est_alpha_omega=/crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_alpha_omega
prop_muts_in_s_ranges=/crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/prop_muts_in_s_ranges

######Step0:calculate dfe input sfs from angsd sfs
cp all_info_bed_origional/bin${i}.10.$batch.divergence.0fold.bed $batch.txt
cp all_info_bed_origional/bin${i}.10.$batch.divergence.4fold.bed 4fold.txt
cp bins_vcf_sfs/$i.$batch.sfs $batch.sfs
cp bins_vcf_sfs/$i.4fold.sfs 4fold.sfs
Rscript ../adaptive_rate.getDivergence.R $input_dfe $batch.txt 4fold.txt $batch.sfs 4fold.sfs
awk 'NF' divergence1 divergence2 >  input/$i.divergence
paste sfs3 zero.txt > sfs33
paste sfs4 zero.txt > sfs44
awk 'NF' sfs1 sfs2 sfs33 sfs44 > input/$batch.$i.sfs
rm $batch.txt 4fold.txt $batch.sfs 4fold.sfs divergence1 divergence2 sfs1 sfs2 sfs3 sfs33 sfs4 sfs44

######Step1:site_class_0
echo  "data_path_1    /crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo  "data_path_2    data-three-epoch/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "sfs_input_file  input/$batch.$i.sfs" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "est_dfe_results_dir   result/$i.results_dir_neut/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "site_class  0" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "fold   1">> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "epochs  2" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "search_n2  1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2_variable   1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2 50" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt

$est_dfe -c result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt

######Step2:site_class_1
echo  "data_path_1    /crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo  "data_path_2    data-three-epoch/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "sfs_input_file  input/$batch.$i.sfs" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "est_dfe_results_dir  result/$i.results_dir_sel/ " >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "est_dfe_demography_results_file   result/$i.results_dir_neut/est_dfe.out" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "site_class  1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "fold   1">> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "epochs  2" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "mean_s_variable  1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "mean_s   -0.1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "beta_variable  1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "beta   0.5" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt

$est_dfe -c result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
$prop_muts_in_s_ranges -c result/$i.results_dir_sel/est_dfe.out -o result/$i.NeS.txt

#####Step3: est_alpha_omega
echo "data_path_1    /crex/proj/uppstore2017145/V3/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "divergence_file  input/$i.divergence" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "est_alpha_omega_results_file   result/$i.$batch.4fold.25inds_est_alpha_omega.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "est_dfe_results_file  result/$i.results_dir_sel/est_dfe.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "neut_egf_file  result/$i.results_dir_neut/neut_egf.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "sel_egf_file  result/$i.results_dir_sel/sel_egf.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "do_jukes_cantor  1" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "remove_poly  0" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt

$est_alpha_omega -c result/$i.$batch.4fold.25inds_est_alpha_omega.txt

done
