# Scripts for 
## Wang X, Ingvarsson PK. Quantifying adaptive evolution and the effects of natural selection across the Norway spruce genome. bioRxiv. 2020 Jun 27:2020-06.

### 0.get functional regions of data and divergence
Get-AngsdSFS.sweden-norwayPop.functionRegion.sh - Script to get 'Sweden-Norway' population and each functional region, e.g. 0fold; 4fold; introns; conserved; promoters and intergenic regions.

divergence_glauga.R
divergence_glauga.need_to_changeDiffFunctional.sh -Scripts to get fix difference sites for each functional regions.

### 1.dfe-alpha
dfe.angsdSFS_to_dfeSFS.R
dfe.sh - Scripts to calculate distribution of deleterious fitness effects and the fraction of adaptive substitutions for both genic and non-genic regions.

### 1.dfe-alpha bootstrap
0.get_SitesPerFunctionalScaffold.R
2.bootstraps.dfe.divergence_data.R
bootstraps.dfe.sh -Scripts to calculate 100 bootstraps for distribution of deleterious fitness effects and the fraction of adaptive substitutions for both genic and non-genic regions.

Plot.dfe_95CI.R -Script to plot NeS and adaptive substitutions for each functional regions.

### 2.adaptive rate
AdaptiveRate.GetStatistic.GeneDensity.R -Script to get bed file and calculate statistics for 10 bins split by gene density.
AdaptiveRate.GetStatistic.recombination.R -Script to get bed file and calculate statistics for 10 bins split by recombination rate (length > 10k).

adaptive_rate.getDivergence.R
adaptive_rate.recombination_geneDensity.sh - Scripts to calculate fraction of adaptive substitutions and adaptive rate for 10 bins.

Plot.Adaptive_rate.AndModelChoce.R -Script to plot correlation between adaptive rate and gene density/recombination rate for 10 bins. 

### 3.Genome-wide scan for Selection
RAiSD_positive_selection.sh -Script to perform genome-wide scan for regions under positive selection.

betascan_balancing_selection.sh -Script to perform genome-wide scan for regions under balancing selection.

statistics_angsd_vcftools_Noselection.sh -Script to calculate pairwise nucleotide diversity and Tajima’s D for regions with randomly neutral scaffolds.
statistics_angsd_vcftools_selection.sh -Script to calculate pairwise nucleotide diversity and Tajima’s D for regions under both positive selection and balancing selection.
