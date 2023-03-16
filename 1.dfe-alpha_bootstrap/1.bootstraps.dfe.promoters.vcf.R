library(dplyr)
data <- read.table('/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/FinalUsing.promoters.all_sites.fixed_diff.bed',header = T)
fold4 <- read.table('/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/FinalUsing.4fold.all_sites.fixed_diff.bed',header = T)

bootstrap_data <- sample_n(data, 5434, replace = T)
bootstrap_4fold <- sample_n(fold4, 5472, replace = T)

write.table(bootstrap_data, "/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/promoters.4fold/promoters.bed",sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
write.table(bootstrap_4fold, "/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.3.DFE.bootstraps/promoters.4fold/4fold.bed",sep="\t", quote=FALSE, col.names = T, row.names = FALSE)
