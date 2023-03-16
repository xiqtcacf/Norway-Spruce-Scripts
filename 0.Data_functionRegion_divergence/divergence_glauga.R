library(dplyr)

#setwd("/domus/h1/xiwa/xi/processedData/analysis/3_ME_MajorRevison_natural_selection/divergence_glauga")
args <- commandArgs(T)

########
##### how to get the fixed differences 1/1 for glauga - (1/1 for glauga overlaping with abies)

abies_function <- read.table(args[1], header = T)
glauga_function <- read.table(args[2],header = T)

functions <- merge(abies_function, glauga_function, by= c("Chr","Pos")) ### merge two species according to scaffold and position
fixed_difference_function <- anti_join(glauga_function, functions, by = c("Chr", "Pos")) ###anti_join to identify the records from the original table that did not exist in our updated table.

write.table(fixed_difference_function, "fixed_difference_function.glauga.txt", sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)




