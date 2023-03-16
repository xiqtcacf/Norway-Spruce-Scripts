library(ggplot2)
library(reshape2)
library("gridExtra")
library(dplyr)

setwd('/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/adaptive_rate_ME')

#### plot figure for adaptive rates, bins 10, recombination, gene density, 0fold VS 4fold, 0fold VS introns
recom_4fold <- read.table('plot.AdaptiveRate.recombination10k',header = T)
gene_4fold <- read.table('plot.AdaptiveRate.gene_density',header = T)

#######
### linear, non-linear model choose
library(stats)
### recombination rate
##Ka/Ks - recombination: linear model
model1 <- lm(Ka.Ks ~ recombination, data = recom_4fold)
summary(model1)
m.exp1 <- nls(Ka.Ks ~ I(a * exp(b * recombination)), data = recom_4fold, start = list(a = 1, b = 0), trace = T)
summary(m.exp1)
AIC(model1, m.exp1)
BIC(model1, m.exp1)
## Ka+ - recombination: expotential model
model11 <- lm(Ka. ~ recombination, data = recom_4fold)
m.exp11 <- nls(Ka. ~ I(a * exp(b * recombination)), data = recom_4fold, start = list(a = 1, b = 0), trace = T)
AIC(model11, m.exp11)
BIC(model11, m.exp11)
### gene density 
## Ka/Ks - gene density
model2 <- lm(Ka.Ks ~ gene_density,data = gene_4fold)
m.exp2 <- nls(Ka.Ks ~ I(a * exp(b * gene_density)), data = gene_4fold, start = list(a = 1, b = 0),trace = T)
AIC(model2, m.exp2)
BIC(model2, m.exp2)
## Ka+ - gene density
model22 <- lm(Ka. ~ gene_density,data = gene_4fold)
m.exp22 <- nls(Ka. ~ I(a * exp(b * gene_density)), data = gene_4fold, start = list(a = 1, b = 0),trace = T)
AIC(model22, m.exp22)
BIC(model22, m.exp22)

#### r co-efficiently
cor.test(recom_4fold$Ka.Ks, recom_4fold$recombination, 
         method="spearman")
cor.test(recom_4fold$Ka., recom_4fold$recombination, 
         method="spearman")
cor.test(gene_4fold$Ka.Ks, gene_4fold$gene_density, 
         method="spearman")
cor.test(gene_4fold$Ka., gene_4fold$gene_density, 
         method="spearman")

#### plot
pdf("/home/xi/Desktop/ME_Revisions_Spruce/2023.MERevision-Xi-Spruce-manuscript/Figure2.adaptive_rate.Ka_Ks.Ka+.recombination.gene.pdf",
    width = 10, height = 10)
### recombination rate
old.par <- par(mfrow=c(2, 2))
model = lm(Ka.Ks ~ recombination,
           data = recom_4fold)
summary(model)  
Anova(model, type="II")  
int =  model$coefficient["(Intercept)"]
slope =model$coefficient["recombination"]
plot(Ka.Ks ~ recombination,
     data = recom_4fold,
     pch=16,
     xlab = "Recombination rate", 
     ylab = "Ka/Ks",
     xlim=c(0,0.01),
     cex=2.5,
     cex.axis=1.5,
     cex.lab=1.5,
     col=c("grey"))
abline(int, slope, lty=1, lwd=1, col="blue")
legend("bottomright", 
       legend=c(as.expression(bquote("r=0.952***")),
                as.expression(bquote(~ r^2 ~"=0.860"))),
       #as.expression(bquote(~ r^2 ~"=0.026"))),
       #as.expression(bquote("P<2.2e-16***"))),
       cex=1.5,
       bty = "n")
#       title="Linear", text.font=1)
m.exp <- nls(Ka. ~ I(a * exp(b * recombination)), data = recom_4fold, start = list(a = 1, b = 0), trace = T)
summary(m.exp)
plot(Ka. ~ recombination,
     data = recom_4fold,
     pch=16,
     xlab = "Recombination rate", 
     ylab = "Ka+",
     xlim=c(0,0.01),
     cex=2.5,
     cex.axis=1.5,
     cex.lab=1.5,
     col=c("grey"))
lines(recom_4fold$recombination, predict(m.exp), lty = 1, col = "blue", lwd = 1)
RSS.p <- sum(residuals(m.exp)^2)  # Residual sum of squares
TSS <- sum((recom_4fold$Ka. - mean(recom_4fold$Ka.))^2)  # Total sum of squares
1 - (RSS.p/TSS)  # R-squared measure
legend("bottomright", 
       legend=c(as.expression(bquote("r=0.164")),
                as.expression(bquote(~ r^2 ~"=0.174"))),
       #as.expression(bquote(~ r^2 ~"=0.026"))),
       #as.expression(bquote("P<2.2e-16***"))),
       cex=1.5,
       bty = "n")
#      title="Curvilinear", text.font=1)

####gene density
m.exp <- nls(Ka.Ks ~ I(a * exp(b * gene_density)), data = gene_4fold, start = list(a = 1, b = 0), trace = T)
summary(m.exp)
plot(Ka.Ks ~ gene_density,
     data = gene_4fold,
     pch=16,
     xlab = "Gene density", 
     ylab = "Ka/Ks",
     xlim=c(0,1),
     cex=2.5,
     cex.axis=1.5,
     cex.lab=1.5,
     col=c("grey"))
lines(gene_4fold$gene_density, predict(m.exp), lty = 1, col = "blue", lwd = 1)
RSS.p <- sum(residuals(m.exp)^2)  # Residual sum of squares
TSS <- sum((gene_4fold$Ka.Ks - mean(gene_4fold$Ka.Ks))^2)  # Total sum of squares
1 - (RSS.p/TSS)  # R-squared measure
legend("topright", 
       legend=c(as.expression(bquote("r=-0.163")),
                as.expression(bquote(~ r^2 ~"=0.164"))),
       #as.expression(bquote(~ r^2 ~"=0.026"))),
       #as.expression(bquote("P<2.2e-16***"))),
       cex=1.5,
       bty = "n")
#      title="Curvilinear", text.font=1)
m.exp2 <- nls(Ka. ~ I(a * exp(b * gene_density)), data = gene_4fold, start = list(a = 1, b = 0), trace = T)
summary(m.exp2)
plot(Ka. ~ gene_density,
     data = gene_4fold,
     pch=16,
     xlab = "Gene density", 
     ylab = "Ka+",
     xlim=c(0,1),
     cex=2.5,
     cex.axis=1.5,
     cex.lab=1.5,
     col=c("grey"))
lines(gene_4fold$gene_density, predict(m.exp2), lty = 1, col = "blue", lwd = 1)
RSS.p <- sum(residuals(m.exp2)^2)  # Residual sum of squares
TSS <- sum((gene_4fold$Ka. - mean(gene_4fold$Ka.))^2)  # Total sum of squares
1 - (RSS.p/TSS)  # R-squared measure
legend("topright", 
       legend=c(as.expression(bquote("r=-0.345")),
                as.expression(bquote(~ r^2 ~"=0.099"))),
       #as.expression(bquote(~ r^2 ~"=0.026"))),
       #as.expression(bquote("P<2.2e-16***"))),
       cex=1.5,
       bty = "n")
#      title="Curvilinear", text.font=1)
dev.off()

