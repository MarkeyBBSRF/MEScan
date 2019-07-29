setwd("/home/tzh232/TiMEx/")
library(TiMEx)
funcPv = function(tempGene,ObsMat){
  pvalue = testCliqueAsGroup(geneIdx = tempGene,mat=t(ObsMat))$pvalueLRT
  #print(pvalue)
  return(pvalue)
}


k<-61
n=60
load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/1Simulation function/simulated.data",k,"Pat",n,"WithTP53.Rdata",sep=""))
load(paste("ResTime",k,".Rdata",sep=""))
filename = paste("TiMExResults3Simu",k,"With200Pat.Rdata",sep="")
TiMExRes = NULL
rankingT = rep(0,100)
for(i in 1:100){
  
  #ResultMat = ResultLists[[i]]$results
  #GeneIndexMat = ResultMat[1:3,1:140] 
  GeneIndexMat = as.matrix(unique(GeneforTime[[i]],MARGIN = 2))  # 
  #############################
  pat.mut.list.test <- listof.pat.mut.list.test[[i]]
  candi.gene <- rownames(pat.mut.list.test[[1]])
  mut.gene <- mut.gene.mat[,i]
  Mut.Gind <- which(candi.gene %in% mut.gene)
  ObsMat = BigObsMat.list[[i]][candi.gene,]
  #############
  pvecTry = try(apply(GeneIndexMat,2,funcPv,ObsMat),TRUE) #get pvalue for first 50 set
  if(isTRUE(class(pvecTry)=="try-error")) { next } else { pvec=pvecTry} 
  res = rbind(GeneIndexMat,pvec)
  OrderRes = res[,order(pvec)]
  TiMExRes[[i]] <- OrderRes
  rankingTemp = which(colSums(Mut.Gind-OrderRes[1:3,])==0)   # whether all of three genes are in the first set.
  rankingT[i] = rankingTemp
  print(c(i,rankingT[i]))
}
####################
save(rankingT,TiMExRes,frac,file=filename)
