####MESGA
k=63
n=80
load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/1Simulation function/simulated.data80236Pat",n,".Rdata",sep=""))
AllGeneName = rownames(BigObsMat.list[[1]])
candi.geneList = NULL
for(m in 1:100){
  candi.geneTemp <- rownames(listof.pat.mut.list.test[[m]][[1]])
  candi.geneList[[m]] <- candi.geneTemp
}
load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/MESGA/MESGA.data",k,".Rdata",sep=""))
######
ResultME<-NULL
for (i in 1:100){
ResultM = ResultMSEGA[[i]]
GeneIndexMat1 = ResultM[c(1:3),1:10]
GeneIndexMat1<-matrix(match(GeneIndexMat1,candi.geneList[[i]]),nrow=3)
ResultME[[i]]<-GeneIndexMat1
}
##########De
ResultDe<-NULL
for (i in 1:100){
  matrixDe<-matrix(NA,nrow=3,ncol=10)
  aa<-read.csv(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/dendrix/Denrix",k,i-1,".csv",sep=""),header =F)
  matrixDe[1,]<-substr(aa[,1], start=7, stop=nchar(as.character(aa[1,1]))-1)
  matrixDe[2,]<-substr(aa[,2], start=3, stop=nchar(as.character(aa[1,2]))-1)
  matrixDe[3,]<-substr(aa[,3], start=3, stop=nchar(as.character(aa[1,3]))-3)
  matrixDe<-matrix(match(matrixDe,candi.geneList[[i]]),nrow=3)
  ResultDe[[i]]<-matrixDe
}

####Co
ResultCo<-NULL
for (i in 1:100){
  aa<-read.csv(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/comet/topgenesets/genesets",k,i-1,".csv",sep=""),header =F)
  matrixCo<-matrix(NA,nrow=3,ncol=10)
  matrixCo<-t(aa[c(1:10),c(1:3)])
  geneset<-matrixCo+1
  ResultCo[[i]]<-geneset
}

####re
ResultRe<-NULL
for (i in 1:100){
  matrixRe<-matrix(NA,nrow=3,ncol=10)
  bb<-read.csv(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/re_test/topgenesets/genesetstr_test",k,i-1,".csv",sep=""),header =F)
  matrixRe<-t(bb[c(1:10),c(1:3)])
  geneset<-matrixRe+1
  ResultRe[[i]]<-geneset
}

GeneIndexMat=NULL
load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/TG/Res",k,"With1Cat200Pat.C.Rdata",sep=""))
for(i in 1:100){
ResultMat = ResultLists[[i]]$results
GeneIndexMat[[i]] = ResultMat[1:3,1:10]
}

GeneforTime= NULL
for(i in 1:100){
  GeneforTime[[i]]<-cbind(ResultME[[i]],ResultDe[[i]],ResultCo[[i]],ResultRe[[i]],GeneIndexMat[[i]],c(1,2,3))
}  
  
filename=paste("ResTime",k,".Rdata",sep="")
#######################
save(GeneforTime,ResultME,ResultDe,ResultCo,ResultRe,GeneIndexMat,file=filename)




