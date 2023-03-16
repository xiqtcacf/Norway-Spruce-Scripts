args <- commandArgs(T)

#### get sfs input file for dfe-alpha from angsd output SFS
data_sfs <- read.table(args[1],header = F)
data_allsites <- read.table(args[2],header = F)
fold4_sfs <- read.table(args[3],header = F)
fold4_allsites <- read.table(args[4],header = F)

data_sfs$V1 <- NULL
data_sfs$V1 <- data_allsites$V1 - rowSums(data_sfs)
data_sfs_new <- data_sfs[,c(ncol(data_sfs),1:(ncol(data_sfs)-1))]

fold4_sfs$V1 <- NULL
fold4_sfs$V1 <- fold4_allsites$V1 - rowSums(fold4_sfs)
fold4_sfs_new <- fold4_sfs[,c(ncol(fold4_sfs),1:(ncol(fold4_sfs)-1))]

cat("1",file="/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE/sfs1")
cat("50",file="/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE/sfs2")
write.table(data_sfs_new, "/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE/sfs3",sep="\t", quote=FALSE, col.names = F, row.names = FALSE)
write.table(fold4_sfs_new, "/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/1.1.DFE/sfs4",sep="\t", quote=FALSE, col.names = F, row.names = FALSE)
