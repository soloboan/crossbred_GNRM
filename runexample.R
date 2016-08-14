genofile=c('../pureB12crossbred.tped','../pureB12crossbred.tfam')
genoformat='tped'
#ana_type='vanRaden'
ana_type=c('vanRaden','yang')
#centeringMethod='admixture'
centeringMethod=c('all','BSAF','admixture')
#scalingK='alldata'
scalingK=c('alldata','purebreddata')
animqrp='../pureB12crossbred.txt'
admixEst='../pureB12crossbred.2.Q'
#outputformat=c('ASREML','dense','matrix')
outputformat='ASREML'
outputname='results'
missinggeno='NA'

source('../calc_grmcrossbreds.R')
for(s in 1:length(ana_type)){
  for(m in 1:length(centeringMethod)){
    for(j in 1:length(scalingK)){
      if(ana_type[s]=='yang' & centeringMethod[m]!='all')
        next()
      datgrm <- grmcrossbred(genofile,genoformat,
                             ana_type=ana_type[s],
                             centeringMethod=centeringMethod[m],
                             scalingK=scalingK[j],
                             animqrp,admixEst,
                             missinggeno,
                             outputformat,
                             outputname=paste('PBCR_',ana_type[s],'-',centeringMethod[m],'-',scalingK[j],sep=''))
      assign(paste('PBCR_',ana_type[s],'-',centeringMethod[m],'-',scalingK[j],sep=''),datgrm)
    }
  }
}

var(`PBCR_vanRaden-all-alldata`[lower.tri(`PBCR_vanRaden-all-alldata`)])
var(`PBCR_vanRaden-all-purebreddata`[lower.tri(`PBCR_vanRaden-all-purebreddata`)])
var(`PBCR_vanRaden-BSAF-alldata`[lower.tri(`PBCR_vanRaden-BSAF-alldata`)])
var(`PBCR_vanRaden-BSAF-purebreddata`[lower.tri(`PBCR_vanRaden-BSAF-purebreddata`)])
var(`PBCR_vanRaden-admixture-alldata`[lower.tri(`PBCR_vanRaden-admixture-alldata`)])
var(`PBCR_vanRaden-admixture-purebreddata`[lower.tri(`PBCR_vanRaden-admixture-purebreddata`)])


diag(`PBCR_vanRaden-all-alldata`) - diag(`PBCR_vanRaden-all-purebreddata`)
diag(`PBCR_vanRaden-BSAF-alldata`) - diag(`PBCR_vanRaden-BSAF-purebreddata`)
diag(`PBCR_vanRaden-all-alldata`) - diag(`PBCR_vanRaden-all-purebreddata`)
diag(`PBCR_vanRaden-all-alldata`) - diag(`PBCR_vanRaden-all-purebreddata`)



`PBCR_vanRaden-all-alldata`[lower.tri(`PBCR_vanRaden-all-alldata`)]-`PBCR_vanRaden-all-purebreddata`[lower.tri(`PBCR_vanRaden-all-purebreddata`)]




vanR <- grmcrossbred(genofile,genoformat,
                  ana_type='vanRaden',centeringMethod='all',
                  scalingK,animqrp,admixEst,
                  missinggeno,
                  outputformat,outputname)

datbsaf <- grmcrossbred(genofile,genoformat,
                    ana_type,centeringMethod,
                    scalingK,animqrp,admixEst,
                    missinggeno,
                    outputformat,outputname)

hist(diag(dat)-diag(datbsaf))

var(dat[lower.tri(dat)])
is.na(dat[lower.tri(dat)])==T
var(datbsaf[lower.tri(datbsaf)])
