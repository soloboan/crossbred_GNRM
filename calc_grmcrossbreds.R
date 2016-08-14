grmcrossbred <- function(genofile,genoformat,ana_type='yang',centeringMethod='all',scalingK,
                         animqrp,admixEst,missinggeno,outputformat='matrix',outputname){
  
  if(genoformat!="tped"){
    if(missing(genofile))
      stop(cat("\n  Need to specify the input genotype file"))
  } else if(genoformat=="tped"){
    genofile1 <- genofile[1]
    genofile2 <- genofile[2]
    if(missing(genofile1))
      stop(cat("\n  Need to specify the tped file"))
    if(missing(genofile2))
      stop(cat("\n  Need to specify the tfam file "))
  }
  if(missing(genoformat))
    stop(cat("\n  Specify genotype file format ... \n'1. ped'\n'2. tped'\n'3. genotypes'\n'4. Robj_ped'\n'5. Robj_genotypes'"))
  if(missing(ana_type))
    stop('\n  Specify the grm option ... \n"1. vanRaden"\n"2. yang"')
  if(missing(centeringMethod))
    stop('\n  Specify the centering (-2p) datatype for the grm  ... \n"1. all"\n"2. BSAF"\n"3. admixture"')
  if(missing(scalingK))
    stop('\n  Specify the centering (sum(2qp) | number of loci) datatype for the grm  ... \n"1. alldata"\n"2. purebreddata"\n')
  if(missing(outputformat))
    stop("\n  Specify the format of the output file ... \n'1. ASREML'\n'2. dense'\n'3. matrix' ")
  if(missing(outputname))
    stop("\n  Specify the name of the output file ")
  if(missing(missinggeno))
    stop("\n Specify mssing genotype code in the dataset")
  
  #############################################################################
  gF <- c("ped","tped","genotypes","Robj_ped","Robj_genotypes")
  oF <- c("ASREML","dense","matrix")
  opt <- c("vanRaden","yang")
  scalegrm <- c('alldata','purebreddata')
  centergrm <- c('all','BSAF','admixture')
  
  if(!genoformat %in% gF)
    stop(cat("\n Specify the genoformat correctly\n'1. ped'\n'2. tped'\n'3. genotypes'\n'4. Robj_ped'\n'5. Robj_genotypes'"))
  if (!outputformat %in% oF)
    stop(cat("\n Specify the outputformat correctly \n'1. ASREML'\n'2. dense'\n'3. matrix'"))
  if(!ana_type %in% opt)
    stop(cat('\n Specify the correct ana_type option \n one of the following ...',
             '\n"1. vanRaden"\n"2. yang",\n'))
  
  #####
  #checking if files specified exist
  if(genoformat=="ped" | genoformat=="genotypes"){
    if(file.exists(genofile)==F)
      stop(cat(paste(genofile," does not exit, re-specify the correct path or filename",sep="")))
  } else if(genoformat=="Robj_ped" | genoformat=="Robj_genotypes"){
    if(exists(genofile)==F)
      stop(cat(paste(genofile," does not exist as an R-object",sep="")))
  } else if(genoformat=="tped"){
    if(file.exists(genofile[1])==F)
      stop(cat(paste(genofile[1]," does not exit, re-specify the correct path or filename",sep="")))
    if(file.exists(genofile[2])==F)
      stop(cat(paste(genofile[2]," does not exit, re-specify the correct path or filename",sep="")))
  }
  
  ####
  # checking if file exist for computing vanRaden_SAF
  if(ana_type=="admixture"){
    if(file.exists(admixEst)==F)
      stop(cat(paste(admixEst," does not exit, re-specify the correct path or filename",sep="")))
  }
  
  #############################################################################
  if (genoformat=="ped"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile,header=F,nrows=2,na.strings=missinggeno))
    geno  <- read.table(genofile,header=F,colClasses=c(rep("character",6),rep("numeric",nc-6)),na.strings=missinggeno)
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID   <- as.vector(geno[,2])
    geno  <- geno[,-1:-6]
    geno  <- geno-1
    geno  <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]  
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="genotypes"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile,header=F,nrows=2,na.strings=missinggeno))
    geno  <- read.table(genofile,header=F,colClasses=c("character",rep("numeric",nc-1)),na.strings=missinggeno)
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID  <- as.vector(geno[,2])
    geno <- geno[,-1:-6]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="tped"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile[1],header=F,nrows=2,na.strings=missinggeno))
    geno  <- read.table(paste(genofile[1],sep=""),header=F,colClasses=c(rep("character",4),rep("numeric",nc-4)),na.strings=missinggeno)
    fam  <- read.table(paste(genofile[2],sep=""),header=F,colClasses=c(rep("character",6)))
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID   <- as.vector(fam[,2])
    geno  <- geno[,-1:-4]
    geno  <- geno-1
    geno  <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    geno <- t(geno)
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1)
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="Robj_ped"){
    geno <- get(genofile)
    cat('\n....... file processing for GRM computation started ..........\n\n')
    IID <- as.vector(geno[,1])
    geno <- geno[,-1]
    geno <- geno-1
    geno <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; z}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="Robj_genotypes"){
    cat('\n....... file processing for GRM computation started ..........\n\n')
    geno <- get(genofile)
    IID <- as.vector(geno[,1])
    geno <- geno[,-1]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
  }
  
  ########################################################################
  # Methods used in computing GRM (basd on vanRaden (2008) and yang et al. (2010))
  if(ana_type=="vanRaden"){
    if(centeringMethod=='all'){
      Z <- scale(x=geno,center=T,scale=F)
      K <- sum(apply(X=geno,FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
      G <- tcrossprod(Z)/K
    } else if(centeringMethod=='BSAF'){
      grpANIM <- as.vector(read.table(animqrp,colClasses=c('character'),stringsAsFactors=F)[,1])
      breeds<-as.vector(unique(grpANIM))
      for(i in 1:length(breeds)){
        grp <- which(grpANIM==breeds[i])
        Zgrp <- scale(x=geno[grp,],center=T,scale=F)
        if(i==1){Z <- Zgrp} else {Z <- rbind(Z,Zgrp)}
        rm(Zgrp)
      }
      Z <- data.matrix(Z)
      if(scalingK=='alldata'){
        K <- sum(apply(X=geno,FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
        G <- tcrossprod(Z)/K
      } else if (scalingK=='purebreddata'){
        grp <- which(grpANIM!='crossbred')
        K <- sum(apply(X=geno[grp,],FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
        G <- tcrossprod(Z)/K
      }
    } else if(centeringMethod=='admixture'){
      grpANIM <- as.vector(read.table(animqrp,colClasses=c('character'),stringsAsFactors=F)[,1])
      breeds<-as.vector(unique(grpANIM))
      ### determine purebreed afreq
      pureB.det <- aggregate(x=1:length(grpANIM),by=list(grpANIM),min)
      pureB.det <- pureB.det[order(pureB.det[,2]),1]
      pureB.det <- pureB.det[which(pureB.det!='crossbred')]
      for(i in 1:length(pureB.det)){
        grp <- which(grpANIM==pureB.det[i])
        afreq <- data.frame(colMeans(geno[grp,]))
        if(i==1){pureBfreq <-afreq} else {pureBfreq <- cbind(pureBfreq,afreq)}
        rm(afreq)
      }
      for(i in 1:length(breeds)){
        if(breeds[i]=='crossbred'){
          grp <- which(grpANIM==breeds[i])
          admix <- read.table(paste(admixEst,sep=""))
          admix <- admix[grp,]
          admix.mean <- data.matrix(admix) %*% data.matrix(t(pureBfreq))
          Zgrp <- geno[grp,]-admix.mean
        } else {
          grp <- which(grpANIM==breeds[i])
          Zgrp <- scale(x=geno[grp,],center=T,scale=F)
        }
        if(i==1){Z <- Zgrp} else {Z <- rbind(Z,Zgrp)}
        rm(Zgrp)
      }
      Z <- data.matrix(Z)
      if(scalingK=='alldata'){
        K <- sum(apply(X=geno,FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
        G <- tcrossprod(Z)/K
      } else if (scalingK=='purebreddata'){
        grp <- which(grpANIM!='crossbred')
        K <- sum(apply(X=geno[grp,],FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
        G <- tcrossprod(Z)/K
      }
    }
    cat('....... G computed according to vanRaden (2008) ..........\n')
  }
  
  if (ana_type=="yang"){
    if(centeringMethod=='all'){
      Z <- scale(x=geno,center=T,scale=T)
      K <- ncol(geno)
      G <- tcrossprod(Z)/K
    } else if(centeringMethod=='BSAF'){
      grpANIM <- as.vector(read.table(animqrp,colClasses=c('character'),stringsAsFactors=F)[,1])
      breeds<-as.vector(unique(grpANIM))
      for(i in 1:length(breeds)){
        grp <- which(grpANIM==breeds[i])
        checkafreq <- colMeans(geno[grp,])/2
        checkafreq <- length(checkafreq[which(checkafreq==0)])
        if(checkafreq>0)
          stop(paste(checkafreq,'markers have frequency (p) = 0 for ',breeds[i],'\n This wont work under Yang et al. (2010) \n suggestion ==> use vanraden (2008)'),'\n')
        Zgrp <- scale(x=geno[grp,],center=T,scale=T)
        if(i==1){Z <- Zgrp} else {Z <- rbind(Z,Zgrp)}
        rm(Zgrp)
      }
      Z <- data.matrix(Z)
      K <- ncol(geno)
      G <- tcrossprod(Z)/K
    } else if(centeringMethod=='admixture'){
      grpANIM <- as.vector(read.table(animqrp,colClasses=c('character'),stringsAsFactors=F)[,1])
      breeds<-as.vector(unique(grpANIM))
      ### determine purebreed afreq
      pureB.det <- aggregate(x=1:length(grpANIM),by=list(grpANIM),min)
      pureB.det <- pureB.det[order(pureB.det[,2]),1]
      pureB.det <- pureB.det[which(pureB.det!='crossbred')]
      for(i in 1:length(pureB.det)){
        grp <- which(grpANIM==pureB.det[i])
        afreq <- data.frame(colMeans(geno[grp,]))
        if(i==1){pureBfreq <-afreq} else {pureBfreq <- cbind(pureBfreq,afreq)}
        rm(afreq)
      }
      for(i in 1:length(breeds)){
        if(breeds[i]=='crossbred'){
          grp <- which(grpANIM==breeds[i])
          admix <- read.table(paste(admixEst,sep=""))
          admix <- admix[grp,]
          admix.mean <- data.matrix(admix) %*% data.matrix(t(pureBfreq))
          Zgrp <- geno[grp,]-admix.mean
          agreqCRB <- colMeans(geno[grp,])/2; agreqCRB <- sqrt(2*agreqCRB*(1-agreqCRB))
          checkafreq <- colMeans(geno[grp,])/2
          checkafreq <- length(checkafreq[which(checkafreq==0)])
          if(checkafreq>0)
            stop(paste(checkafreq,'markers have frequency (p) = 0 for ',breeds[i],'\n This wont work under Yang et al. (2010) \n suggestion ==> use vanraden (2008)'),'\n')
          Zgrp <- scale(Zgrp,center=F,scale=agreqCRB)
        } else {
          grp <- which(grpANIM==breeds[i])
          checkafreq <- sum(apply(X=geno[grp,],FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
          checkafreq <- length(checkafreq[which(checkafreq==0)])
          if(checkafreq>0)
            stop('Some alleles have frequency (p) = 0 for ',grp,'\n This wont work under Yang et al. (2010) \n suggestion ==> use vanraden (2008) \n')
          Zgrp <- scale(x=geno[grp,],center=T,scale=T)
        }
        if(i==1){Z <- Zgrp} else {Z <- rbind(Z,Zgrp)}
        rm(Zgrp)
      }
      Z <- data.matrix(Z)
      K <- ncol(geno)
      G <- tcrossprod(Z)/K
    }
    cat('....... G computed according to Yang et al. (2010) ..........\n')
  }
  
  #############
  if(outputformat=="ASREML"){
    cat('....... Output file preparation started ..........\n')
    require(reshape2)
    rowiseG <- G
    rowiseG[upper.tri(rowiseG)] <- NA
    colnames(rowiseG) <- IID; rownames(rowiseG) <- IID
    rowiseG <- melt(rowiseG)
    rowiseGnr <- G
    colnames(rowiseGnr) <- 1:length(IID); rownames(rowiseGnr) <- 1:length(IID)
    rowiseGnr <- melt(rowiseGnr)
    rowiseG <- cbind(rowiseGnr[,1:2],rowiseG)
    rm(rowiseGnr)
    rowiseG <- na.omit(rowiseG)
    rowiseG <- rowiseG[order(rowiseG[,1],rowiseG[,2]),]
    write.table(rowiseG,paste(outputname,".grm",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
    cat('+++++++++ output file columns ++++++++++\n')
    cat(' recodedID1 recodedID2 originalID1  originalID2 G_values \n\n')
  }  
  if(outputformat=="dense"){
    cat('....... Output file preparation started ..........\n')
    require(reshape2)
    rowiseG <- G
    colnames(rowiseG) <- IID; rownames(rowiseG) <- IID
    rowiseG <- melt(rowiseG)
    rowiseGnr <- G
    colnames(rowiseGnr) <- 1:length(IID); rownames(rowiseGnr) <- 1:length(IID)
    rowiseGnr <- melt(rowiseGnr)
    rowiseG <- cbind(rowiseGnr[,1:2],rowiseG)
    rm(rowiseGnr)
    write.table(rowiseG,paste(outputname,".grm",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
    cat('+++++++++ output file columns ++++++++++\n')
    cat(' recodedID1 recodedID2 originalID1  originalID2 G_values \n\n')
  } 
  if(outputformat=="matrix"){
    cat('....... Output file preparation started ..........\n')
    Gmat <- G
    colnames(Gmat) <- IID; rownames(Gmat) <- IID
    write.table(Gmat,paste(outputname,".grm",sep=""),col.names=T,row.names=F,quote=T,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
  } 
  
  plotname <-paste("Genomic_relationship_",outputname,".png",sep="")
  dpi=300;width.cm<-15;height.cm<-13;width.in<-width.cm/2.54;height.in<-height.cm/2.54
  png(file=plotname,width=width.in*dpi,height=height.in*dpi,pointsize=8,units="px",res=dpi)
  par(mfrow=c(1,2))
  par(mar=c(5,4,5,2))
  br<-length(IID)*0.5; if(br>100){br=100}
  hist(diag(G),breaks=br,col=sample(2:6,1),xlab="diagonal element of G-matrix",
       main=paste("grm for ",outputname,"\n","based on ",ana_type,sep=""))
  hist(G[lower.tri(G)],breaks=br,col=sample(2:6,1),
       xlab="off-diagonal element of G-matrix",main="")
  dev.off()
  colnames(G) <- IID; rownames(G) <- IID
  return(G)
}
