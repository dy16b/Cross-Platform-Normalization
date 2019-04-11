
### This code was written by Dr. Serin Zhang, as a part of her Statistics PhD
### dissertation project at Florida State University

# setwd("C:/Users/Disa Yu/Dropbox/Cross-platform-normalization/R_code")

library("genefilter")
library("MatchMixeR")

### create a "true" expression data and fit OLS/FLMER models 
set.seed(4321)
m <- 10000; n <- 100
## TrueMeans is the true mean expression for each gene.
TrueMeans <- rnorm(m, 7.5, .6)

## Generate the training data (matched data)
matchdat.X <- matrix(rnorm(m*n, 0, .5), m) + TrueMeans
beta0 <- rnorm(m, 2.5, 2.0); beta1 <- runif(m, .7, 1.0)
betas <- cbind(beta0=beta0, beta1=beta1)
matchdat.Y <- beta0 + matchdat.X * beta1 

## fit models on the first matched data : OLS & FLMER   
OLSmod <- MatchMixeR:::OLS(matchdat.X, matchdat.Y)
#FLMERmod <-ModIV(matchdat.X, matchdat.Y)
FLMERmod <-MM(matchdat.X, matchdat.Y)

### Initial Method comparison
## create research data(test data) "WITHOUT" DEGs - "smaller" matched data 
#n1 <- 5
n1 <- 30
sim.match.X <- matrix(rnorm(m*n1, 0, .3), m) + TrueMeans
sim.match.Y <- beta0 + sim.match.X * beta1 

#sim.match.X5 <- sim.match.X; sim.match.Y5 <- sim.match.Y
#sim.match.X30 <- sim.match.X; sim.match.Y30 <- sim.match.Y

## Apply fitted OLS&FLMER on the first matched data to the smaller matched data   
XtransOLS <- sim.match.X*OLSmod$betamat[, "Slope"] + OLSmod$betamat[, "Intercept"]
XtransFLMER <- sim.match.X*FLMERmod$betamat[, "Slope"] + FLMERmod$betamat[, "Intercept"]
## other methods normalization on smaller matched data
DWDmod <- dwd(sim.match.X, sim.match.Y)
EBmod <- eb(sim.match.X, sim.match.Y)
#GQmod <- gq(sim.match.X, sim.match.Y)
XPNmod <- xpn(sim.match.X, sim.match.Y)

## make plots : Coloumn-wise R2 comparison: DWD beats our methods 
corr.raw <- as.numeric(lapply(1:n1,function(i) cor(sim.match.X[,i], sim.match.Y[,i])))
corr.OLS <- as.numeric(lapply(1:n1,function(i) cor(XtransOLS[,i],sim.match.Y[,i])))
corr.FLMER <- as.numeric(lapply(1:n1,function(i) cor(XtransFLMER[,i],sim.match.Y[,i])))
corr.DWD <- as.numeric(lapply(1:n1,function(i) cor(DWDmod$x[,i],DWDmod$y[,i])))
corr.EB <- as.numeric(lapply(1:n1,function(i) cor(EBmod$x[,i],EBmod$y[,i])))
#corr.GQ <- as.numeric(lapply(1:n1,function(i) cor(GQmod$x[,i],GQmod$y[,i])))
corr.XPN <- as.numeric(lapply(1:n1,function(i) cor(XPNmod$x[,i],XPNmod$y[,i])))
r2.raw <- (corr.raw)^2
r2.OLS <- (corr.OLS)^2
r2.FLMER <- (corr.FLMER)^2
r2.DWD <- (corr.DWD)^2
r2.EB <- (corr.EB)^2
#r2.GQ <- (corr.GQ)^2
r2.XPN <- (corr.XPN)^2
r2.results <- cbind(MM=r2.FLMER,DWD=r2.DWD,XPN=r2.XPN,EB=r2.EB)
boxplot(r2.results,ylab="R(column wise)",ylim = c(0.999,1),cex.lab=1.3, cex.axis=1.3) 

#r2.sim.n5 <- r2.results
r2.sim.n30 <- r2.results

## RSS comparison: OLS vs FLMER
RSS.OLS <- mean((XtransOLS-sim.match.Y)^2)
RSS.FLMER <- mean((XtransFLMER-sim.match.Y)^2)
## mean sum of difference squres for other methods
Diff.raw <- mean((sim.match.X-sim.match.Y)^2)
Diff.EB <- mean((EBmod$x - EBmod$y)^2)
Diff.DWD <- mean(as.matrix(DWDmod$x - DWDmod$y)^2)
#Diff.GQ <- mean((GQmod$x - GQmod$y)^2)
Diff.XPN <- mean(as.matrix(XPNmod$x - XPNmod$y)^2)

#diff.n5 <- c(raw = Diff.raw, MM = RSS.FLMER,DWD = Diff.DWD,XPN = Diff.XPN,EB = Diff.EB)
diff.n30 <- c(raw = Diff.raw, MM = RSS.FLMER,DWD = Diff.DWD,XPN = Diff.XPN,EB = Diff.EB)

save(r2.sim.n30, diff.n30, file="sim_results.RData")

### DE analysis
##Generate GroupA with n3 platformX samples & n4 platformY samples
#GroupB with n5 platformX samples & n6 platformY samples
n2 <- 30; n3 <- 0 # GroupA&PlatformX ; GroupB&PlatformX  
n4 <- 0; n5 <- 30 # GroupA&PlatformY ; GroupB&PlatformY
TP <- matrix(0,30,6);  FP <- matrix(0,30,6)  #exclude GQ
colnames(TP)<-c("raw","OLS","MM","DWD","XPN","EB")
colnames(FP)<-c("raw","OLS","MM","DWD","XPN","EB")

## create research data with 1000 DEGs (By adding Effs to GroupB)
m1 <- 1000
Effs <- rep(0, m); Effs[1:m1] <- runif(m1, .2, 3)
for (i in 1:30)
{
  print(i)
  sim.A.X <- matrix(rnorm(m*n2, 0, .25), m) + TrueMeans 
  sim.B.X <- matrix(rnorm(m*n3, 0, .25), m) + TrueMeans + Effs
  sim.X <- cbind(sim.A.X,sim.B.X) 
  TrueExp.A.Y <- matrix(rnorm(m*n4, 0, .25), m) + TrueMeans
  TrueExp.B.Y <- matrix(rnorm(m*n5, 0, .25), m) + TrueMeans + Effs
  TrueExp.Y <- cbind(TrueExp.A.Y,TrueExp.B.Y)
  sim.Y <- beta0 + TrueExp.Y * beta1 
  
  ## Apply the trained betas to the research data
  Xtrans.OLS_DE <- sim.X*OLSmod$betamat[, "Slope"] + OLSmod$betamat[, "Intercept"]
  Xtrans.FLMER_DE <- sim.X*FLMERmod$betamat[, "Slope"] + FLMERmod$betamat[, "Intercept"]
  ## other methods normalization on research data with DEGs 
  DWDmod_DE <- dwd(sim.X,sim.Y)
  EBmod_DE <- eb(sim.X,sim.Y)
  #GQmod_DE <- gq(sim.X,sim.Y)
  XPNmod_DE <-xpn(sim.X,sim.Y)
  ## check stat. power and type I error
  simData <- cbind(sim.X, sim.Y)
  simDataOLS <- cbind(Xtrans.OLS_DE, sim.Y)
  simDataFLMER <- cbind(Xtrans.FLMER_DE, sim.Y)
  simDataDWD <- as.matrix(cbind(DWDmod_DE$x,DWDmod_DE$y))
  simDataEB <- as.matrix(cbind(EBmod_DE$x,EBmod_DE$y))
  #simDataGQ <- as.matrix(cbind(GQmod_DE$x,GQmod_DE$y))
  simDataXPN <- as.matrix(cbind(XPNmod_DE$x,XPNmod_DE$y))
  simFac <- as.factor(c(rep(0, n2), rep(1, n3),rep(0, n4), rep(1, n5)))
  rr.raw <- rowttests(simData, simFac)
  rr.OLS <- rowttests(simDataOLS, simFac)
  rr.FLMER <- rowttests(simDataFLMER, simFac)
  rr.DWD <- rowttests(simDataDWD, simFac)
  rr.EB <- rowttests(simDataEB, simFac)
  #rr.GQ <- rowttests(simDataGQ, simFac)
  rr.XPN <- rowttests(simDataXPN, simFac)
  padj.raw <- p.adjust(rr.raw[, "p.value"], "holm")
  padj.OLS <- p.adjust(rr.OLS[, "p.value"], "holm")
  padj.FLMER <- p.adjust(rr.FLMER[, "p.value"], "holm")
  padj.DWD <- p.adjust(rr.DWD[, "p.value"], "holm")
  padj.EB <- p.adjust(rr.EB[, "p.value"], "holm")
  #padj.GQ <- p.adjust(rr.GQ[, "p.value"], "holm")
  padj.XPN <- p.adjust(rr.XPN[, "p.value"], "holm")
  TP[i,1] <- sum(padj.raw[1:m1]<0.05); FP[i,1] <- sum(padj.raw[(m1+1):m]<0.05)
  TP[i,2] <- sum(padj.OLS[1:m1]<0.05); FP[i,2] <- sum(padj.OLS[(m1+1):m]<0.05)
  TP[i,3] <- sum(padj.FLMER[1:m1]<0.05); FP[i,3] <- sum(padj.FLMER[(m1+1):m]<0.05)
  TP[i,4] <- sum(padj.DWD[1:m1]<0.05); FP[i,4] <- sum(padj.DWD[(m1+1):m]<0.05)
  TP[i,6] <- sum(padj.EB[1:m1]<0.05); FP[i,6] <- sum(padj.EB[(m1+1):m]<0.05)
  #TP[i,7] <- sum(padj.GQ[1:m1]<0.05); FP[i,7] <- sum(padj.GQ[(m1+1):m]<0.05)
  TP[i,5] <- sum(padj.XPN[1:m1]<0.05); FP[i,5] <- sum(padj.XPN[(m1+1):m]<0.05)
}

DE <- matrix(1000,30,6) 
P <- TP/(TP+FP)
R <- TP/DE
F1 <- 2*(P*R)/(P+R) 


