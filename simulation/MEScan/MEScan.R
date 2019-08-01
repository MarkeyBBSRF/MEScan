## test c implementation for one category
library(mescan)

### extract information
candi.geneList = NULL
ObsMatList     = NULL
RateMatList = NULL
for(k in 1:100){
  candi.geneTemp <- rownames(listof.pat.mut.list.test[[k]][[1]])
  candi.geneList[[k]] <- candi.geneTemp
  ObsMatList[[k]] <- BigObsMat.list[[k]][candi.geneTemp,]
  RateMatList[[k]] <- EtaMat[candi.geneTemp,]
}

TG.vecNew1 = rep(0,100)
Opt.setNew1 = matrix(rep(0,3*100),ncol=100)
rankingNew1 = rep(0,100)
lambdaVec = rep(0,100)
ResultLists = NULL
##run 100 times
for(i in 1:100){
  ObsMat <- ObsMatList[[i]]
  RateMat <- RateMatList[[i]]
  lambda = quantile(RateMat^(1/2), probs=0.05) 
  lambdaVec[i]= lambda
  mut.gene <- mut.gene.mat[,i]
  Mut.Gind = which(mut.gene %in% candi.geneList[[i]])
  ReTab1 = scan_genesets(c(1:20), 3, ObsMat,RateMat,lambda=lambda)
  ResultLists[[i]]<-ReTab1
  TG.vecNew1[i] = ReTab1$Larg_TG
  Opt.setNew1[,i] = ReTab1$Opt_set
  rankingNew1[i] = which(colSums(abs(ReTab1$results[1:3,]-Mut.Gind))==0)
  print(c(i,TG.vecNew1[i],lambda))
  print(rankingNew1[i])
}
######################
filename="Res51With1Cat200Pat.C.Rdata"
#######################
save(ResultLists,TG.vecNew1,Opt.setNew1,rankingNew1,lambdaVec,ResultLists,file=filename)
