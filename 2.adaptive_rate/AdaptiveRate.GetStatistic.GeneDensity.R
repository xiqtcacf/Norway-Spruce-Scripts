library(ggplot2)
library(reshape2)
library("gridExtra")
library("ggpubr")
library(car)
library(reshape)
library(dplyr)

setwd("/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/adaptive_rate_ME")
##### bin genes for gene density
all_scaffolds <- read.table('1k_reference_0fixed.bed',header = T)
#lengthAbove10k <- read.table('lengthAbove10k.1k_reference_0fixed.bed',header = T)
gene_density <- read.table('factor/Pab01b_gene_density_genebase_per_scaffoldlength.txt',header = T)
geneD_gene1 <- merge(all_scaffolds,gene_density, by='Chr') 
gene <- read.table('gene.scaffolds',header = T)
geneD_gene <- merge(gene, geneD_gene1, by='Chr' )
#hist(geneD_gene$gene_Den, breaks=20)

quantile(geneD_gene$gene_Den,probs=c(0.10))
10% 
0.01862171
quantile(geneD_gene$gene_Den,probs=c(0.20))
20% 
0.03564025 
quantile(geneD_gene$gene_Den,probs=c(0.30))
30% 
0.05758504 
quantile(geneD_gene$gene_Den,probs=c(0.40))
40% 
0.08827245
quantile(geneD_gene$gene_Den,probs=c(0.50))
50% 
0.1318233 
quantile(geneD_gene$gene_Den,probs=c(0.60))
60% 
0.1937714 
quantile(geneD_gene$gene_Den,probs=c(0.70))
70% 
0.2844086  
quantile(geneD_gene$gene_Den,probs=c(0.80))
80% 
0.4203016 
quantile(geneD_gene$gene_Den,probs=c(0.90))
90% 
0.6393386 


bin1 <- subset(geneD_gene, geneD_gene$gene_Den < 0.01862171)
bin2 <- subset(geneD_gene, geneD_gene$gene_Den > 0.01862171 & geneD_gene$gene_Den < 0.03564025)
bin3 <- subset(geneD_gene, geneD_gene$gene_Den > 0.03564025 & geneD_gene$gene_Den < 0.05758504)
bin4 <- subset(geneD_gene, geneD_gene$gene_Den > 0.05758504 & geneD_gene$gene_Den < 0.08827245)
bin5 <- subset(geneD_gene, geneD_gene$gene_Den > 0.08827245 & geneD_gene$gene_Den < 0.1318233)
bin6 <- subset(geneD_gene, geneD_gene$gene_Den > 0.1318233 & geneD_gene$gene_Den < 0.1937714)
bin7 <- subset(geneD_gene, geneD_gene$gene_Den > 0.1937714 & geneD_gene$gene_Den < 0.2844086)
bin8 <- subset(geneD_gene, geneD_gene$gene_Den > 0.2844086 & geneD_gene$gene_Den < 0.4203016)
bin9 <- subset(geneD_gene, geneD_gene$gene_Den > 0.4203016 & geneD_gene$gene_Den < 0.6393386)
bin10 <- subset(geneD_gene, geneD_gene$gene_Den > 0.6393386)


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


min(bin1$gene_Den)
max(bin1$gene_Den)
mean(bin1$gene_Den)
min(bin2$gene_Den)
max(bin2$gene_Den)
mean(bin2$gene_Den)
min(bin3$gene_Den)
max(bin3$gene_Den)
mean(bin3$gene_Den)
min(bin4$gene_Den)
max(bin4$gene_Den)
mean(bin4$gene_Den)
min(bin5$gene_Den)
max(bin5$gene_Den)
mean(bin5$gene_Den)
min(bin6$gene_Den)
max(bin6$gene_Den)
mean(bin6$gene_Den)
min(bin7$gene_Den)
max(bin7$gene_Den)
mean(bin7$gene_Den)
min(bin8$gene_Den)
max(bin8$gene_Den)
mean(bin8$gene_Den)
min(bin9$gene_Den)
max(bin9$gene_Den)
mean(bin9$gene_Den)
min(bin10$gene_Den)
max(bin10$gene_Den)
mean(bin10$gene_Den)


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

write.table(bin1, "bins10_gene_bed/bin1.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin2, "bins10_gene_bed/bin2.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin3, "bins10_gene_bed/bin3.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin4, "bins10_gene_bed/bin4.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin5, "bins10_gene_bed/bin5.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin6, "bins10_gene_bed/bin6.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin7, "bins10_gene_bed/bin7.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin8, "bins10_gene_bed/bin8.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin9, "bins10_gene_bed/bin9.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bin10, "bins10_gene_bed/bin10.10.gene_density.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)


### how many 0fold / 4fold sites (all sites instead of polymorphic sites!!!!!!) and all fix difference
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
write.table(sites_fold0_bin1, "bins10_gene_bed//bin1.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin1, "bins10_gene_bed//bin1.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin2, "bins10_gene_bed//bin2.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin2, "bins10_gene_bed//bin2.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin3, "bins10_gene_bed//bin3.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin3, "bins10_gene_bed//bin3.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin4, "bins10_gene_bed//bin4.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin4, "bins10_gene_bed//bin4.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin5, "bins10_gene_bed//bin5.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin5, "bins10_gene_bed//bin5.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin6, "bins10_gene_bed//bin6.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin6, "bins10_gene_bed//bin6.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin7, "bins10_gene_bed//bin7.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin7, "bins10_gene_bed//bin7.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin8, "bins10_gene_bed//bin8.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin8, "bins10_gene_bed//bin8.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin9, "bins10_gene_bed//bin9.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin9, "bins10_gene_bed//bin9.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold0_bin10, "bins10_gene_bed//bin10.10.gene_density.divergence.0fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(sites_fold4_bin10, "bins10_gene_bed//bin10.10.gene_density.divergence.4fold.bed", sep="\t", quote=FALSE, col.names = T, row.names = FALSE)



