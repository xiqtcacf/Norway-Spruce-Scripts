library(ggplot2)
library(reshape2)
library("gridExtra")
library("ggpubr")
library(car)
library(reshape)
library(dplyr)
setwd("/home/xi/Desktop/ME_Revisions_Spruce/MERevision-Xi-Spruce-manuscript-2022-09-28/adaptive_rate_ME")

##### bin genes for recombination rates

lengthAbove10k <- read.table('lengthAbove10k.1k_reference_0fixed.bed',header = T)
recom_swe_nor <- read.table('only.sweden-norway.recombination_rate.txt',header = T)
recom <- merge(lengthAbove10k,recom_swe_nor, by='Chr' ) ### get recombination rate focus on scaffold length above 10k
gene <- read.table('gene.scaffolds',header = T)

recom_gene <- merge(gene,recom, by='Chr' )

quantile(recom_gene$recombination_rate,probs=c(0.10))
10% 
0.0002902048  
quantile(recom_gene$recombination_rate,probs=c(0.20))
20% 
0.0006386398  
quantile(recom_gene$recombination_rate,probs=c(0.30))
30% 
0.001032138  
quantile(recom_gene$recombination_rate,probs=c(0.40))
40% 
0.001474118 
quantile(recom_gene$recombination_rate,probs=c(0.50))
50% 
0.001969355 
quantile(recom_gene$recombination_rate,probs=c(0.60))
60% 
0.002513759  
quantile(recom_gene$recombination_rate,probs=c(0.70))
70% 
0.00313483  
quantile(recom_gene$recombination_rate,probs=c(0.80))
80% 
0.003981203 
quantile(recom_gene$recombination_rate,probs=c(0.90))
90% 
0.005318895 

bin1 <- subset(recom_gene, recom_gene$recombination_rate < 0.0002902048)
bin2 <- subset(recom_gene, recom_gene$recombination_rate > 0.0002902048 & recom_gene$recombination_rate < 0.0006386398)
bin3 <- subset(recom_gene, recom_gene$recombination_rate > 0.0006386398 & recom_gene$recombination_rate < 0.001032138)
bin4 <- subset(recom_gene, recom_gene$recombination_rate > 0.001032138 & recom_gene$recombination_rate < 0.001474118)
bin5 <- subset(recom_gene, recom_gene$recombination_rate > 0.001474118 & recom_gene$recombination_rate < 0.001969355)
bin6 <- subset(recom_gene, recom_gene$recombination_rate > 0.001969355 & recom_gene$recombination_rate < 0.002513759)
bin7 <- subset(recom_gene, recom_gene$recombination_rate > 0.002513759 & recom_gene$recombination_rate < 0.00313483)
bin8 <- subset(recom_gene, recom_gene$recombination_rate > 0.00313483 & recom_gene$recombination_rate < 0.003981203)
bin9 <- subset(recom_gene, recom_gene$recombination_rate > 0.003981203 & recom_gene$recombination_rate < 0.005318895)
bin10 <- subset(recom_gene, recom_gene$recombination_rate > 0.005318895)

#### summary statistics for each bins
sum(bin1$Num_genes)
sum(bin2$Num_genes)
sum(bin3$Num_genes)
sum(bin4$Num_genes)
sum(bin5$Num_genes)
sum(bin6$Num_genes)
sum(bin7$Num_genes)
sum(bin8$Num_genes)
sum(bin9$Num_genes)
sum(bin10$Num_genes)

min(bin1$recombination_rate)
max(bin1$recombination_rate)
mean(bin1$recombination_rate)
min(bin2$recombination_rate)
max(bin2$recombination_rate)
mean(bin2$recombination_rate)
min(bin3$recombination_rate)
max(bin3$recombination_rate)
mean(bin3$recombination_rate)
min(bin4$recombination_rate)
max(bin4$recombination_rate)
mean(bin4$recombination_rate)
min(bin5$recombination_rate)
max(bin5$recombination_rate)
mean(bin5$recombination_rate)
min(bin6$recombination_rate)
max(bin6$recombination_rate)
mean(bin6$recombination_rate)
min(bin7$recombination_rate)
max(bin7$recombination_rate)
mean(bin7$recombination_rate)
min(bin8$recombination_rate)
max(bin8$recombination_rate)
mean(bin8$recombination_rate)
min(bin9$recombination_rate)
max(bin9$recombination_rate)
mean(bin9$recombination_rate)
min(bin10$recombination_rate)
max(bin10$recombination_rate)
mean(bin10$recombination_rate)

##### get bed files for each bins
#all_scaffolds <- read.table('1k_reference_0fixed.bed',header = T)
#bin1_bed <- merge(bin1,all_scaffolds, by='Chr')
#bin2_bed <- merge(bin2,all_scaffolds, by='Chr')
#bin3_bed <- merge(bin3,all_scaffolds, by='Chr')
#bin4_bed <- merge(bin4,all_scaffolds, by='Chr')
#bin5_bed <- merge(bin5,all_scaffolds, by='Chr')
#bin6_bed <- merge(bin6,all_scaffolds, by='Chr')
#bin7_bed <- merge(bin7,all_scaffolds, by='Chr')
#bin8_bed <- merge(bin8,all_scaffolds, by='Chr')
#bin9_bed <- merge(bin9,all_scaffolds, by='Chr')
#bin10_bed <- merge(bin10,all_scaffolds, by='Chr')
write.table(bin1, "bins10_recombination_bed/bin1.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin2, "bins10_recombination_bed/bin2.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin3, "bins10_recombination_bed/bin3.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin4, "bins10_recombination_bed/bin4.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin5, "bins10_recombination_bed/bin5.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin6, "bins10_recombination_bed/bin6.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin7, "bins10_recombination_bed/bin7.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin8, "bins10_recombination_bed/bin8.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin9, "bins10_recombination_bed/bin9.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin10, "bins10_recombination_bed/bin10.10.recombination.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)

### how many 0fold / 4fold sites (all sites instead of polymorphic sites!!!!!!) and fix difference
fold0 <- read.table('FinalUsing.0fold.all_sites.fixed_diff.bed',header = T)
fold4 <- read.table('FinalUsing.4fold.all_sites.fixed_diff.bed',header = T)
sum(fold4$all_sites)

sites_fold0_bin1 <- merge(fold0, bin1, by="Chr")
sites_fold4_bin1 <- merge(fold4, bin1, by="Chr")
sum(sites_fold0_bin1$all_sites)
sum(sites_fold4_bin1$all_sites)

sites_fold0_bin2 <- merge(fold0, bin2, by="Chr")
sites_fold4_bin2 <- merge(fold4, bin2, by="Chr")
sum(sites_fold0_bin2$all_sites)
sum(sites_fold4_bin2$all_sites)

sites_fold0_bin3 <- merge(fold0, bin3, by="Chr")
sites_fold4_bin3 <- merge(fold4, bin3, by="Chr")
sum(sites_fold0_bin3$all_sites)
sum(sites_fold4_bin3$all_sites)

sites_fold0_bin4 <- merge(fold0, bin4, by="Chr")
sites_fold4_bin4 <- merge(fold4, bin4, by="Chr")
sum(sites_fold0_bin4$all_sites)
sum(sites_fold4_bin4$all_sites)

sites_fold0_bin5 <- merge(fold0, bin5, by="Chr")
sites_fold4_bin5 <- merge(fold4, bin5, by="Chr")
sum(sites_fold0_bin5$all_sites)
sum(sites_fold4_bin5$all_sites)

sites_fold0_bin6 <- merge(fold0, bin6, by="Chr")
sites_fold4_bin6 <- merge(fold4, bin6, by="Chr")
sum(sites_fold0_bin6$all_sites)
sum(sites_fold4_bin6$all_sites)

sites_fold0_bin7 <- merge(fold0, bin7, by="Chr")
sites_fold4_bin7 <- merge(fold4, bin7, by="Chr")
sum(sites_fold0_bin7$all_sites)
sum(sites_fold4_bin7$all_sites)

sites_fold0_bin8 <- merge(fold0, bin8, by="Chr")
sites_fold4_bin8 <- merge(fold4, bin8, by="Chr")
sum(sites_fold0_bin8$all_sites)
sum(sites_fold4_bin8$all_sites)

sites_fold0_bin9 <- merge(fold0, bin9, by="Chr")
sites_fold4_bin9 <- merge(fold4, bin9, by="Chr")
sum(sites_fold0_bin9$all_sites)
sum(sites_fold4_bin9$all_sites)

sites_fold0_bin10 <- merge(fold0, bin10, by="Chr")
sites_fold4_bin10 <- merge(fold4, bin10, by="Chr")
sum(sites_fold0_bin10$all_sites)
sum(sites_fold4_bin10$all_sites)


##### get bed files for each bins
write.table(sites_fold0_bin1, "bins10_recombination_bed/bin1.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin1, "bins10_recombination_bed/bin1.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin2, "bins10_recombination_bed/bin2.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin2, "bins10_recombination_bed/bin2.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin3, "bins10_recombination_bed/bin3.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin3, "bins10_recombination_bed/bin3.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin4, "bins10_recombination_bed/bin4.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin4, "bins10_recombination_bed/bin4.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin5, "bins10_recombination_bed/bin5.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin5, "bins10_recombination_bed/bin5.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin6, "bins10_recombination_bed/bin6.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin6, "bins10_recombination_bed/bin6.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin7, "bins10_recombination_bed/bin7.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin7, "bins10_recombination_bed/bin7.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin8, "bins10_recombination_bed/bin8.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin8, "bins10_recombination_bed/bin8.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin9, "bins10_recombination_bed/bin9.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin9, "bins10_recombination_bed/bin9.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin10, "bins10_recombination_bed/bin10.10.recombination.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin10, "bins10_recombination_bed/bin10.10.recombination.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)










