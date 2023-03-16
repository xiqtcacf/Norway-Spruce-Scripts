#!/bin/bash
#SBATCH -A naiss2023-5-36
#SBATCH -o %J_dfe_coding.out
#SBATCH -e %J_dfe_coding.err
#SBATCH -J %J_dfe_coding.job
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-user wxi12345.mail@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load vcftools
module load gsl
angsd=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/angsd
realSFS=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/realSFS

batch=$1 ### functional region

input_path=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/
output=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE
fasta=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa
fai=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/data/1k_picea_abies.master-rna-scaff.nov2012_sorted.fa.fai
dfe_fold=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16
input_dfe=/crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE
output_dfe=/crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE

export DIR=/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE
cd $DIR

est_dfe=/crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_dfe
est_alpha_omega=/crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/est_alpha_omega
prop_muts_in_s_ranges=/crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/prop_muts_in_s_ranges

######Step0.1:calculate dfe input sfs from angsd sfs
Rscript $output/dfe.angsdSFS_to_dfeSFS.R $input_path/sweden_norway.noPab002Pab034.$batch.saf.idx.sfs $output/allsites.$batch $input_path/sweden_norway.noPab002Pab034.4fold.saf.idx.sfs $output/allsites.4fold
awk 'NF' $output/sfs1 $output/sfs2 $output/sfs3 $output/sfs4 > $output/$batch.4fold.sfs
rm $output/sfs1 $output/sfs2 $output/sfs3 $output/sfs4

######Step1:site_class_0
echo  "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo  "data_path_2    data-three-epoch/" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "sfs_input_file  $output/$batch.4fold.sfs" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "est_dfe_results_dir   $output/$batch.4fold.results_dir_neut/" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "site_class  0" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "fold   1">> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "epochs  2" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "search_n2  1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2_variable   1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt
echo "t2 50" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt

$est_dfe -c $input_dfe/$batch.4fold.25inds_est_dfe-site_class-0.txt

######Step2:site_class_1
echo  "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo  "data_path_2    data-three-epoch/" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "sfs_input_file  $output/$batch.4fold.sfs" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "est_dfe_results_dir  $output/$batch.4fold.results_dir_sel/ " >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "est_dfe_demography_results_file   $output/$batch.4fold.results_dir_neut/est_dfe.out" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "site_class  1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "fold   1">> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "epochs  2" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "mean_s_variable  1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "mean_s   -0.1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "beta_variable  1" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
echo "beta   0.5" >> $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt

$est_dfe -c $input_dfe/$batch.4fold.25inds_est_dfe-site_class-1.txt
$prop_muts_in_s_ranges -c $input_dfe/$batch.4fold.results_dir_sel/est_dfe.out -o $input_dfe/$batch.NeS

#####Step3: est_alpha_omega
echo "data_path_1    /crex/proj/uppstore2017145/V2/users/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/software/dfe-alpha-release-2.16/data/" > $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "divergence_file  $output/divergence.allVSfixed.glauga.$batch.4fold.txt" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "est_alpha_omega_results_file   $output/$batch.4fold.25inds_est_alpha_omega.out" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "est_dfe_results_file  $output/$batch.4fold.results_dir_sel/est_dfe.out" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "neut_egf_file  $output/$batch.4fold.results_dir_neut/neut_egf.out" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "sel_egf_file  $output/$batch.4fold.results_dir_sel/sel_egf.out" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "do_jukes_cantor  1" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
echo "remove_poly  0" >> $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt

$est_alpha_omega -c $input_dfe/$batch.4fold.25inds_est_alpha_omega.txt
