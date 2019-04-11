
### This code was written by Dr. Serin Zhang, as a part of her Statistics PhD
### dissertation project at Florida State University

library("genefilter")
library("MatchMixeR")

# setwd("C:/Users/Disa Yu/Dropbox/Cross-platform-normalization")
# getwd()

load("./Real_data/matched_TCGA_mat.Rdata")
load("./Real_data/NormalvsTumor_mat.Rdata")

## na in array data -> 0 
matched_array_mat[is.na(matched_array_mat)] <- 0

## choose 100 samples(out of 523 matched samples) as training set
train.ID <- sample(1:523,100,replace=FALSE)
X <- matched_array_mat[,train.ID] ; Y <- matched_seq_mat[,train.ID]
train.sample <- colnames(X)
OLSmod <- MatchMixeR:::OLS(X, Y)
FLMERmod <-MM(X, Y)

### DE analysis 
## Get DEGs from Yuhang's list (SigGeneList_1.xls)
SigGeneList <- read.csv("./Real_data/SigGeneList_nocomments.csv")

SigGene <- SigGeneList[,c(1,3,7)]
# using cut off p-value < 0.001 ---> # of DEGs = 9555
SigGene_padj <- SigGene[SigGene$padj< .001,]
colnames(SigGene_padj)[1] <- "SigGeneList_padj"
SigG1 <- as.character(SigGene_padj$SigGeneList_padj)
glist_1 <- strsplit(SigG1, "|", fixed=TRUE) 
SigG1_ID <-t(data.frame(glist_1))
SigG1_ID <- data.frame(SigG1_ID)
colnames(SigG1_ID)[1] <- "geneID"
common.gene <- data.frame(row.names(matched_array_mat)); colnames(common.gene) <-"geneID" 
DEG1 <- merge(cbind(common.gene,c(1:16146)),SigG1_ID,by="geneID")
DEG_pa <- as.character(DEG1$geneID)

## remove train samples from the test sample pool  
normal_array_1 <- normal_array[ ,setdiff(colnames(normal_array), train.sample)]
tumor_array_1 <- tumor_array[,setdiff (colnames(tumor_array), train.sample)]
normal_seq_1 <- normal_seq[ ,setdiff(colnames(normal_seq), train.sample)]
tumor_seq_1 <- tumor_seq[,setdiff (colnames(tumor_seq), train.sample)]
s1 <- length(colnames(normal_array_1)); s2 <-length(colnames(tumor_array_1)) 
s3 <- length(colnames(normal_seq_1)); s4 <-length(colnames(tumor_seq_1)) 

TP <- matrix(0,30,8);  FP <- matrix(0,30,8); com <- matrix(0,30,4)
colnames(TP) <- c("raw","OLS","MM","DWD","XPN","EB","GQ","Fisher")
colnames(FP) <- c("raw","OLS","MM","DWD","XPN","EB","GQ","Fisher")
colnames(com) <- c("array.P", "seq.P", "com.P", "com.TP")

for (i in 1:5)
{ print(i)
  n1 <- 20; n2 <- 0 # GroupA&PlatformX ; GroupB&PlatformX  
  n3 <- 0; n4 <- 20 # GroupA&PlatformY ; GroupB&PlatformY
  normal_arrayID <-  sample(1:s1,n1,replace=FALSE)
  tumor_arrayID <-  sample(1:s2,n2,replace=FALSE)
  normal_seqID <-  sample(1:s3,n3,replace=FALSE)
  tumor_seqID <-  sample(1:s4,n4,replace=FALSE)
  test.array <- cbind(normal_array_1[,normal_arrayID],tumor_array_1[,tumor_arrayID])
  test.array[is.na(test.array)] = 0
  test.seq <- cbind(normal_seq_1[,normal_seqID],tumor_seq_1[,tumor_seqID])
  
  XtransOLS <- test.array*OLSmod$betamat[, "Slope"] + OLSmod$betamat[, "Intercept"]
  XtransFLMER <- test.array*FLMERmod$betamat[, "Slope"] + FLMERmod$betamat[, "Intercept"]
  test.array<- data.frame(test.array);test.seq<- data.frame(test.seq)
  DWDmod <- dwd(test.array, test.seq)
  EBmod <- eb(test.array, test.seq)
  GQmod <- gq(test.array,test.seq)  # Check: is gq from genefilter package?
  XPNmod <- xpn(test.array, test.seq)
  
  ## check stat. power and type I error
  DataArray <- as.matrix(test.array)
  DataSeq <- as.matrix(test.seq)
  DataRaw <- as.matrix(cbind(test.array, test.seq))
  DataOLS <- as.matrix(cbind(XtransOLS, test.seq))
  DataFLMER <- as.matrix(cbind(XtransFLMER, test.seq))
  DataDWD <- as.matrix(cbind(DWDmod$x,DWDmod$y))
  DataXPN <- as.matrix(cbind(XPNmod$x,XPNmod$y))
  DataEB <- as.matrix(cbind(EBmod$x,EBmod$y))
  DataGQ <- as.matrix(cbind(GQmod$x,GQmod$y))
  fac_1 <- as.factor(c(rep(0, n1), rep(1, n2),rep(0, n3), rep(1, n4)))
  fac_array <- as.factor(c(rep(0, n1), rep(1, n2)))
  fac_seq <- as.factor(c(rep(0, n3), rep(1, n4)))                 
  rr.raw <- rowttests(DataRaw, fac_1)
  rr.OLS <- rowttests(DataOLS, fac_1)
  rr.FLMER <- rowttests(DataFLMER, fac_1)
  rr.DWD <- rowttests(DataDWD, fac_1)
  rr.XPN <- rowttests(DataXPN, fac_1)
  rr.EB <- rowttests(DataEB, fac_1)
  rr.GQ <- rowttests(DataGQ, fac_1)
  rr.array <- rowttests(DataArray, fac_array)
  rr.seq <- rowttests(DataSeq, fac_seq)   
  
  p.array <- rr.array[rr.array$p.value < 0.05,]
  TP.array <- as.numeric(dim(p.array[row.names(p.array) %in% DEG_pa,])[1])
  FP.array <- as.numeric(dim(p.array)[1]) - TP.array  
  
  p.seq <- rr.seq[rr.seq$p.value < 0.05,]  
  TP.seq <- as.numeric(dim(p.seq[row.names(p.seq) %in% DEG_pa,])[1])
  FP.seq <- as.numeric(dim(p.seq)[1]) - TP.seq
  
  # common TP bet.array&seq results 
  tp.array <- p.array[row.names(p.array) %in% DEG_pa,]
  tp.seq <- p.seq[row.names(p.seq) %in% DEG_pa,]
  com.p <- as.numeric(dim(p.array[row.names(p.array) %in% row.names(p.seq),])[1])
  com.tp <- as.numeric(dim(tp.array[row.names(tp.array) %in% row.names(tp.seq),])[1])
  array.p <- as.numeric(dim(p.array)[1])
  seq.p <- as.numeric(dim(p.seq)[1])
  
  sumlogP <- -2*(log(rr.array$p.value) + log(rr.seq$p.value)) #sum of log 
  combineP <- 1-pchisq(sumlogP,4)
  rr.combine <- cbind(rr.array,combineP)
  p.fisher <- rr.combine[rr.combine$combineP < 0.05,]
  TP.fisher <- as.numeric(dim(p.fisher[row.names(p.fisher) %in% DEG_pa,])[1])
  FP.fisher<- as.numeric(dim(p.fisher)[1]) - TP.fisher
 
  
  p.raw <- rr.raw[rr.raw$p.value < 0.05,]
  TP.raw <- as.numeric(dim(p.raw[row.names(p.raw) %in% DEG_pa,])[1])
  FP.raw <- as.numeric(dim(p.raw)[1]) - TP.raw 
  
  
  p.OLS <- rr.OLS[rr.OLS$p.value < 0.05,]
  TP.OLS <- as.numeric(dim(p.OLS[row.names(p.OLS) %in% DEG_pa,])[1])
  FP.OLS <- as.numeric(dim(p.OLS)[1]) - TP.OLS 
  
  
  p.FLMER <- rr.FLMER[rr.FLMER$p.value < 0.05,]
  TP.FLMER <- as.numeric(dim(p.FLMER[row.names(p.FLMER) %in% DEG_pa,])[1])
  FP.FLMER <- as.numeric(dim(p.FLMER)[1]) - TP.FLMER 
  
  p.DWD <- rr.DWD[rr.DWD$p.value < 0.05,]
  TP.DWD <- as.numeric(dim(p.DWD[row.names(p.DWD) %in% DEG_pa,])[1])
  FP.DWD <- as.numeric(dim(p.DWD)[1]) - TP.DWD  
  
  p.XPN <- rr.XPN[rr.XPN$p.value < 0.05,]
  TP.XPN <- as.numeric(dim(p.XPN[row.names(p.XPN) %in% DEG_pa,])[1])
  FP.XPN <- as.numeric(dim(p.XPN)[1]) - TP.XPN
  
  p.EB <- rr.EB[rr.EB$p.value < 0.05,]
  TP.EB <- as.numeric(dim(p.EB[row.names(p.EB) %in% DEG_pa,])[1])
  FP.EB <- as.numeric(dim(p.EB)[1]) - TP.EB
  
  p.GQ <- rr.GQ[rr.GQ$p.value < 0.05,]
  TP.GQ <- as.numeric(dim(p.GQ[row.names(p.GQ) %in% DEG_pa,])[1])
  FP.GQ <- as.numeric(dim(p.GQ)[1]) - TP.GQ
  
  TP[i,1]<-TP.raw; FP[i,1]<-FP.raw
  TP[i,2]<-TP.OLS; FP[i,2]<-FP.OLS
  TP[i,3]<-TP.FLMER; FP[i,3]<-FP.FLMER
  TP[i,4]<-TP.DWD; FP[i,4]<-FP.DWD
  TP[i,5]<-TP.XPN; FP[i,5]<-FP.XPN
  TP[i,6]<-TP.EB; FP[i,6]<-FP.EB
  TP[i,7]<-TP.GQ; FP[i,7]<-FP.GQ
  TP[i,8]<-TP.fisher; FP[i,8]<-FP.fisher
  
  com [i,1]<- array.p;com [i,2]<-seq.p; com [i,3]<-com.p; com [i,4]<-com.tp
}

ex6 <- cbind(TP,FP,com) #repeat from line 62 with different n1-n4 numbers: ex1(10,10,10,10), ex2(15,15,5,5), ex3(15,5,5,15), ex4(20,0,5,15), ex5(25,5,0,10), ex6(20,0,0,20)
class(ex6) # "matrix"
dim(ex6) # 30 20

tp <- mean(ex6[,2:8])
fp <- mean(ex6[,10:16])
P <- tp/(tp + fp)
R <- tp/9555
F1 <- 2*(P*R)/(P+R)

# Getting tp, fp, and F1 values for MM in ex6
tp_MM = ex6[,3]
fp_MM = ex6[,11]

# Getting tp, fp, and F1 values for DWD in ex6
tp_DWD = ex6[,4]
fp_DWD = ex6[,12]

ratio6 = "ex6(20,0,0,20)"
# cat("iteration = ", iter <- iter + 1, "\n")

cat("The ratio is ", ratio6, "\n",
    "The value of tp is ", tp, "\n", 
    "The value of fp is ", fp, "\n",
    "The value of F1 is ", F1)
