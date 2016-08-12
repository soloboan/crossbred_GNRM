setwd("../")
X <- read.table("pureB12crossbred.tped")
X <- X[,-1:-4] - 1
X <- X[,seq(1,ncol(X),2)] + X[,seq(2,ncol(X),2)]
X <- t(X)
dim(X)

afrq <- read.table("pureB12crossbred.afrq")[,7]
afrq2 <- 2*afrq
het <- sqrt(2*afrq*(1-afrq))
W <- scale(x=X,center=afrq2,scale=het)
GoverallHET <- tcrossprod(W)/length(afrq)
hist(diag(GoverallHET),breaks=100,col="red")
hist(GoverallHET[lower.tri(GoverallHET)],breaks=100,col="green",xlim=c(-0.5,1.5))

GoverallHET[1:5,1:5]

##### breed specific allele frequency
fam <- read.table('pureB12crossbred.fam')
breeds<-as.vector(unique(fam[,1]))
for(i in 1:length(breeds)){
  famX <- as.numeric(rownames(fam[which(fam[,1]==breeds[i]),]))
  afrqX <- read.table(paste(breeds[i],".afrq",sep=""))[,7]
  afrq2X <- 2*afrqX
  hetX <- sqrt(2*afrqX*(1-afrqX))
  Wx <- scale(x=X[famX,],center=afrq2X,scale=hetX)
  if(i==1){Wxs <- Wx} else {Wxs <- rbind(Wxs,Wx)}
  rm(Wx)
}
Wxs <- data.matrix(Wxs)
GindHET <- tcrossprod(Wxs)/ncol(Wxs)
hist(diag(GindHET),breaks=100,col="blue")
hist(GindHET[lower.tri(GindHET)],breaks=100,col="purple",xlim=c(-0.5,1.5))

hist(diag(GoverallHET)-diag(GindHET),breaks=100)


##### crossbred with admixture
fam<- read.table("pureB12crossbred.tfam")
pureb <- c("pureB2","pureB1")
for(j in 1:length(pureb)){
  famX <- as.numeric(rownames(fam[which(fam[,1]==pureb[j]),]))
  afrqX <- read.table(paste(pureb[j],".afrq",sep=""))[,7]
  
  #### center and scale genotypes
  afrq2X <- 2*afrqX
  hetX <- sqrt(2*afrqX*(1-afrqX))
  Wx <- scale(x=X[famX,],center=afrq2X,scale=hetX)
  if(j==1){PBX <- Wx} else {PBX <- rbind(PBX,Wx)}
  rm(Wx)
  
  ### storing allele frequecies for crossbreeds
  afrqX <- data.frame(afrqX)
  colnames(afrqX) <- pureb[j]
  if(j==1){afrqbreeds <- afrqX} else {afrqbreeds <- cbind(afrqbreeds,afrqX)}
}

breeds <- as.vector(unique(fam[,1]))
CBbreeds <- breeds[!breeds %in% pureb]
for(i in 1:length(CBbreeds)){
  famX <- as.numeric(rownames(fam[which(fam[,1]==CBbreeds[i]),]))
  QX <- read.table(paste("pureB12crossbred.2.Q",sep=""))
  QX <- QX[famX,]
  PX <- 2*t(afrqbreeds)
  
  QPX <- data.matrix(QX) %*% data.matrix(PX)
  Wxy <-X[famX,]-QPX
  hetX <- sqrt(2*rowMeans(afrqbreeds)*(1-rowMeans(afrqbreeds)))
  Wx <- scale(x=Wxy,center=F,scale=hetX)
  if(i==1){CBX <- Wx} else {CBX <- rbind(CBX,Wx)}
  rm(Wxy,Wx)
}

CBPBX <- data.matrix(rbind(CBX,PBX))
GCBPBX <- tcrossprod(CBPBX)/ncol(CBPBX)

hist(diag(GCBPBX),breaks=100,col="brown")
hist(GCBPBX[lower.tri(GCBPBX)],breaks=100,col="cyan",xlim=c(-0.5,1.5))

hist(diag(GoverallHET)-diag(GindHET),breaks=100,xlim=c(-0.2,0.6))
hist(diag(GoverallHET)-diag(GCBPBX),breaks=100,xlim=c(-0.2,0.6))
hist(diag(GindHET)-diag(GCBPBX),breaks=100,xlim=c(-0.3,0.4))

hist(GindHET[lower.tri(GindHET)]-GCBPBX[lower.tri(GCBPBX)],
     breaks=100,xlim=c(-0.3,0.2))

#######################
eigen <- read.table("pureB12crossbred.eigenvec")
plot(eigen$V3,eigen$V4,pch=20,col=as.numeric(eigen$V1),xlab='PCA 1',
     ylab='PCA 2')

#cs=as.numeric(eigen$V1)
