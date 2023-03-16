library(ggplot2)
library(reshape2)
library("gridExtra")
library("ggpubr")
library(car)
library(reshape)
library(dplyr)
library(grid)
library(gridExtra)
library(shape)

setwd('/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/dfe/')

#######
#### bootstraps 100 95% CI

Nes_0fold <- read.table('100bootstraps.0fold.4fold.NeS.txt',header = F)
omega_0fold <- read.table('100bootstraps.0fold.4fold.alpha_omega.out',header = F)
quantile(Nes_0fold$V3,probs=c(.025,.975))
quantile(Nes_0fold$V6,probs=c(.025,.975))
quantile(Nes_0fold$V9,probs=c(.025,.975))
quantile(Nes_0fold$V12,probs=c(.025,.975))
quantile(omega_0fold$V6,probs=c(.025,.975))
quantile(omega_0fold$V8,probs=c(.025,.975))

Nes_intron<- read.table('100bootstraps.intron.4fold.NeS.txt',header = F)
omega_intron <- read.table('100bootstraps.intron.4fold.alpha_omega.out',header = F)
quantile(Nes_intron$V3,probs=c(.025,.975))
quantile(Nes_intron$V6,probs=c(.025,.975))
quantile(Nes_intron$V9,probs=c(.025,.975))
quantile(Nes_intron$V12,probs=c(.025,.975))
quantile(omega_intron$V6,probs=c(.025,.975))
quantile(omega_intron$V8,probs=c(.025,.975))

Nes_conserved<- read.table('100bootstraps.conserved.4fold.NeS.txt',header = F)
omega_conserved <- read.table('100bootstraps.conserved.4fold.alpha_omega.out',header = F)
quantile(Nes_conserved$V3,probs=c(.025,.975))
quantile(Nes_conserved$V6,probs=c(.025,.975))
quantile(Nes_conserved$V9,probs=c(.025,.975))
quantile(Nes_conserved$V12,probs=c(.025,.975))
quantile(omega_conserved$V6,probs=c(.025,.975))
quantile(omega_conserved$V8,probs=c(.025,.975))

Nes_promoters<- read.table('100bootstraps.promoters.4fold.NeS.txt',header = F)
omega_promoters <- read.table('100bootstraps.promoters.4fold.alpha_omega.out',header = F)
quantile(Nes_promoters$V3,probs=c(.025,.975))
quantile(Nes_promoters$V6,probs=c(.025,.975))
quantile(Nes_promoters$V9,probs=c(.025,.975))
quantile(Nes_promoters$V12,probs=c(.025,.975))
quantile(omega_promoters$V6,probs=c(.025,.975))
quantile(omega_promoters$V8,probs=c(.025,.975))

Nes_intergenic <- read.table('100bootstraps.intergenic.4fold.NeS.txt',header = F)
omega_intergenic <- read.table('100bootstraps.intergenic.4fold.alpha_omega.out',header = F)
quantile(Nes_intergenic$V3,probs=c(.025,.975))
quantile(Nes_intergenic$V6,probs=c(.025,.975))
quantile(Nes_intergenic$V9,probs=c(.025,.975))
quantile(Nes_intergenic$V12,probs=c(.025,.975))
quantile(omega_intergenic$V6,probs=c(.025,.975))
quantile(omega_intergenic$V8,probs=c(.025,.975))


#### plot
## basic R
Nes <- readr::read_delim('plot.DFE.txt')
alpha <- read.table('plot.alpha.omega.txt',header = T)
pdf("/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/Figure1.1.dfe-alpha.pdf",
    width = 10, height = 10)
old.par <- par(mfrow=c(1, 1))
barplot(as.matrix(Nes), main="", ylab="", xlab="", beside=T, 
        col=c('dark red','dark green','yellow','pink','dark blue'),
        ylim=c(-0.05,1.08))
legend(19, 0.8, c("Zerofold","Intronic","Conserved","Promoters","Intergenic", "95% CI"), cex=1, 
       fill=c('dark red','dark green','yellow','pink','dark blue', 'black'))
arrows(x0=1.5,y0=0.623,y1=0.771,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=7.5,y0=0.076,y1=0.092,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=13.5,y0=0.084,y1=0.096,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=19.5,y0=0.043,y1=0.215,angle=90,code=3,length=0.05,col="black",lwd = 1.5)

arrows(x0=25.5,y0=-0.037,y1=0.227,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=26.5,y0=0.293,y1=0.517,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=27.5,y0=0.141,y1=0.338,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=28.5,y0=0.313,y1=0.491,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=29.5,y0=0.032,y1=0.256,angle=90,code=3,length=0.05,col="black",lwd = 1.5)

arrows(x0=31.5,y0=-0.024,y1=0.193,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=32.5,y0=0.414,y1=1.07,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=33.5,y0=0.164,y1=0.510,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=34.5,y0=0.456,y1=0.965,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=35.5,y0=0.033,y1=0.345,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
dev.off()

pdf("/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/Figure1.2.dfe-alpha.pdf",
    width = 5, height = 5)
old.par <- par(mfrow=c(1, 1))
plot(alpha$Number,alpha$a, ylim=c(-0.1,0.6), pch=16, type = "p", col = c('dark red','dark green','yellow','pink','dark blue'),cex=1.5)
arrows(x0=1,y0=-0.037,y1=0.227,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=2,y0=0.293,y1=0.517,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=3,y0=0.141,y1=0.338,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=4,y0=0.313,y1=0.491,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=5,y0=0.032,y1=0.256,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
dev.off()

pdf("/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/Figure1.3.dfe-alpha.pdf",
    width = 5, height = 5)
old.par <- par(mfrow=c(1, 1))
plot(alpha$Number,alpha$omega, ylim=c(-0.1,1),pch=16, type='p', col = c('dark red','dark green','yellow','pink','dark blue'),cex=1.5)
arrows(x0=1,y0=-0.024,y1=0.193,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=2,y0=0.414,y1=1.07,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=3,y0=0.164,y1=0.510,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=4,y0=0.456,y1=0.965,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
arrows(x0=5,y0=0.033,y1=0.345,angle=90,code=3,length=0.05,col="black",lwd = 1.5)
dev.off()
