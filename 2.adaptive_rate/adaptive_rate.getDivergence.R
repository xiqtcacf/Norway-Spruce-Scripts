args <- commandArgs(T)

setwd(args[1])
#### get divergence input file for dfe-alpha
data <- read.table(args[2],header = T)
fold4 <- read.table(args[3],header = T)
cat("1",sum(data$all_sites), sum(data$fixed_diff),file="divergence1")
cat("0",sum(fold4$all_sites), sum(fold4$fixed_diff),file="divergence2")

#### get sfs input file for dfe-alpha
data_sfs <- read.table(args[4],header = F)
fold4_sfs <- read.table(args[5],header = F)

data_sfs$V1 <- NULL
data_sfs$V1 <- sum(data$all_sites) - rowSums(data_sfs)
data_sfs_new <- data_sfs[,c(ncol(data_sfs),1:(ncol(data_sfs)-1))]

fold4_sfs$V1 <- NULL
fold4_sfs$V1 <- sum(fold4$all_sites) - rowSums(fold4_sfs)
fold4_sfs_new <- fold4_sfs[,c(ncol(fold4_sfs),1:(ncol(fold4_sfs)-1))]

cat("1",file="sfs1")
cat("50",file="sfs2")
write.table(data_sfs_new, "sfs3",sep="\t", quote=FALSE, col.names = F, row.names = FALSE)
write.table(fold4_sfs_new, "sfs4",sep="\t", quote=FALSE, col.names = F, row.names = FALSE)
