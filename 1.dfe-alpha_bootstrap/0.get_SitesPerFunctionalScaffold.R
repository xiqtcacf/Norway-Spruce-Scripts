library(dplyr)
#### get all sites for each scaffold and each functional region
#####for 0fold and 4 fold, they are point distribution, only calculate duplicates for each scaffolds:
#step1.in linux 
cut -f1 0fold.bed | sort | uniq -c > 0fold.perScaffold.allSite
cut -f1 4fold.bed | sort | uniq -c > 4fold.perScaffold.allSite
###step2.adding name chr and all_sites
#####for intron and intergenic, they are continusously distribution, calculate length for each duplicates for each scaffolds and group same scaffolds:
intron <- read.table('intron.bed',header = F)
intron$diff <- intron$V3 - intron$V2
intron_allsites <-
  intron %>%
    group_by(V1) %>%
    summarise(all_sites=(sum(diff)))

#### put all information together: all sites, chr position, and fixed difference
### 0fold
bed_1k <- read.table('1k_reference_0fixed.bed',header = T)
t <- read.table('0fold.perScaffold.allSite',header = T)
all_sites_0fold <- merge(bed_1k, t, by= 'Chr')
fixed_0fold <- read.table('fixed_difference_0fold.glauga.scaffolds',header = T)
fixed <- merge(all_sites_0fold,fixed_0fold,by="Chr")
no_fixed <- anti_join(all_sites_0fold, fixed_0fold, by = "Chr")
no_fixed$Num_fixed_diff_0fold <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

### 4fold
t <- read.table('4fold.perScaffold.allSite',header = T)
all_sites_4fold <- merge(bed_1k, t, by= 'Chr')
fixed_4fold <- read.table('fixed_difference_4fold.glauga.scaffolds',header = T)
fixed <- merge(all_sites_4fold,fixed_4fold,by="Chr")
no_fixed <- anti_join(all_sites_4fold, fixed_4fold, by = "Chr")
no_fixed$Num_fixed_diff_4fold <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

### intron
t <- read.table('intron.perScaffold.allSite',header = T)
all_sites_intron <- merge(bed_1k, t, by= 'Chr')
fixed_intron <- read.table('fixed_difference_intron.glauga.scaffolds',header = T)
fixed <- merge(all_sites_intron,fixed_intron,by="Chr")
no_fixed <- anti_join(all_sites_intron, fixed_intron, by = "Chr")
no_fixed$Num_fixed_diff_intron <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

### conserved
bed_1k <- read.table('1k_reference_0fixed.bed',header = T)
t <- read.table('conserved.perScaffold.allSite',header = T)
all_sites_conserved <- merge(bed_1k, t, by= 'Chr')
fixed_conserved <- read.table('fixed_difference_conserved.glauga.scaffolds',header = T)
fixed <- merge(all_sites_conserved,fixed_conserved,by="Chr")
no_fixed <- anti_join(all_sites_conserved, fixed_conserved, by = "Chr")
no_fixed$Num_fixed_diff_conserved <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

### promoters
bed_1k <- read.table('1k_reference_0fixed.bed',header = T)
t <- read.table('promoters.perScaffold.allSite',header = T)
all_sites_promoters <- merge(bed_1k, t, by= 'Chr')
fixed_promoters <- read.table('fixed_difference_promoters.glauga.scaffolds',header = T)
fixed <- merge(all_sites_promoters, fixed_promoters, by="Chr")
no_fixed <- anti_join(all_sites_promoters, fixed_promoters, by = "Chr")
no_fixed$Num_fixed_diff_promoters <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

####intergenic same with others
####linux cat fixed.txt no_fixed.txt > FinalUsing.function.all_sites.fixed_diff.bed
### intergenic
fixed_intergenic<- read.table('fixed_difference_intergenic.glauga.scaffolds',header = T)
sum(fixed_intergenic$Num_fixed_diff_intergenic)
all_sites_intergenic <- read.table('All_sites_per_scaffold.intergenic',header = T)
all_sites_intergenic <- 
  all_sites_intergenic %>%
    group_by(Chr) %>%
    summarise(all_sites_intergenic=(sum(all_sites)))
sum(as.numeric(all_sites_intergenic$all_sites_intergenic))
fixed <- merge(all_sites_intergenic,fixed_intergenic,by="Chr")
sum(fixed$Num_fixed_diff_intergenic)
no_fixed <- anti_join(all_sites_intergenic, fixed_intergenic, by = "Chr")
no_fixed$Num_fixed_diff_intergenic <- 0
write.table(fixed, "fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(no_fixed, "no_fixed.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
