#!/bin/bash
#SBATCH -A naiss2023-5-36
#SBATCH -o %J_dfe_bootstrap.out
#SBATCH -e %J_dfe_bootstrap.err
#SBATCH -J %J_dfe_bootstrap.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load vcftools/0.1.15
module load gsl
angsd=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/angsd
realSFS=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/realSFS

batch=$1 #### functional region

input=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/sweden_norway.noPab002Pab034.$batch.vcf.gz
input_4fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/sweden_norway.noPab002Pab034.4fold.vcf.gz
output=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/$batch.4fold
fasta=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai
dfe_fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16

export DIR=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/$batch.4fold
cd $DIR

#need to change cores relate to different function region
### first step bootstraps bed files and calculate sfs by ANGSD

for i in {1..100}
do
 Rscript ../1.bootstraps.dfe.$batch.vcf.R
 mv $batch.bed bootstrap_bed/$i.$batch.txt
 mv 4fold.bed bootstrap_bed/$i.4fold.txt
 cut -f1,2,3 bootstrap_bed/$i.$batch.txt > bootstrap_bed/$i.$batch.bed
 cut -f1,2,3 bootstrap_bed/$i.4fold.txt > bootstrap_bed/$i.4fold.bed
 vcftools --gzvcf $input --bed bootstrap_bed/$i.$batch.bed --recode --out bootstrap_vcf_sfs/$i.$batch
 vcftools --gzvcf $input_4fold --bed bootstrap_bed/$i.4fold.bed --recode --out bootstrap_vcf_sfs/$i.4fold
 $angsd -vcf-PL bootstrap_vcf_sfs/$i.$batch.recode.vcf -out bootstrap_vcf_sfs/$i.$batch -doSaf 1 -doMajorMinor 1 -anc $fasta -nInd 25 -minQ 30
 $angsd -vcf-PL bootstrap_vcf_sfs/$i.4fold.recode.vcf -out bootstrap_vcf_sfs/$i.4fold -doSaf 1 -doMajorMinor 1 -anc $fasta -nInd 25 -minQ 30
 $realSFS bootstrap_vcf_sfs/$i.$batch.saf.idx -cores 1 -fold 1 > bootstrap_vcf_sfs/$i.$batch.sfs
 $realSFS bootstrap_vcf_sfs/$i.4fold.saf.idx -cores 1 -fold 1 > bootstrap_vcf_sfs/$i.4fold.sfs 
 rm bootstrap_vcf_sfs/$i.$batch.recode.vcf
 rm bootstrap_vcf_sfs/$i.4fold.recode.vcf
 rm bootstrap_vcf_sfs/$i.$batch.arg
 rm bootstrap_vcf_sfs/$i.$batch.saf.gz
 rm bootstrap_vcf_sfs/$i.$batch.saf.idx
 rm bootstrap_vcf_sfs/$i.$batch.saf.pos.gz
 rm bootstrap_vcf_sfs/$i.4fold.arg
 rm bootstrap_vcf_sfs/$i.4fold.saf.gz
 rm bootstrap_vcf_sfs/$i.4fold.saf.idx
 rm bootstrap_vcf_sfs/$i.4fold.saf.pos.gz

### second step, calculate dfe by dfe-alpha
input_dfe=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/$batch.4fold
output_dfe=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/$batch.4fold
est_dfe=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_dfe
est_alpha_omega=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_alpha_omega
prop_muts_in_s_ranges=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/prop_muts_in_s_ranges

######Step0:calculate dfe input sfs from angsd sfs
cp bootstrap_bed/$i.$batch.txt input/$batch.txt
cp bootstrap_bed/$i.4fold.txt input/4fold.txt
cp bootstrap_vcf_sfs/$i.$batch.sfs input/$batch.sfs
cp bootstrap_vcf_sfs/$i.4fold.sfs input/4fold.sfs
Rscript ../2.bootstraps.dfe.divergence_data.R $input_dfe input/$batch.txt input/4fold.txt input/$batch.sfs input/4fold.sfs
awk 'NF' divergence1 divergence2 > input/bootstrap.$i.divergence
#paste sfs3 zero.txt > sfs33
#paste sfs4 zero.txt > sfs44
awk 'NF' sfs1 sfs2 sfs3 sfs4 > input/bootstrap.$i.sfs
rm input/$batch.txt input/4fold.txt input/$batch.sfs input/4fold.sfs divergence1 divergence2 sfs1 sfs2 sfs3 sfs4

######Step1:site_class_0
echo  "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo  "data_path_2    data-three-epoch/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "sfs_input_file  input/bootstrap.$i.sfs" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "est_dfe_results_dir   result/$i.results_dir_neut/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "site_class  0" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "fold   1">> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "epochs  2" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "search_n2  1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2_variable   1" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2 50" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt

$est_dfe -c result/$i.$batch.4fold.25inds_est_dfe-site_class-0.txt

######Step2:site_class_1
echo  "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo  "data_path_2    data-three-epoch/" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "sfs_input_file  input/bootstrap.$i.sfs" >> result/$i.$batch.4fold.25inds_est_dfe-site_class-1.txt
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
echo "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "divergence_file  input/bootstrap.$i.divergence" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "est_alpha_omega_results_file   result/$i.$batch.4fold.25inds_est_alpha_omega.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "est_dfe_results_file  result/$i.results_dir_sel/est_dfe.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "neut_egf_file  result/$i.results_dir_neut/neut_egf.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "sel_egf_file  result/$i.results_dir_sel/sel_egf.out" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "do_jukes_cantor  1" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt
echo "remove_poly  0" >> result/$i.$batch.4fold.25inds_est_alpha_omega.txt

$est_alpha_omega -c result/$i.$batch.4fold.25inds_est_alpha_omega.txt

done
