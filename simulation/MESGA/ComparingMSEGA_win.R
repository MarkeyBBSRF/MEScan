source("MEGSA.R")

for (i in c(57,60)){
load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/1Simulation function/simulated.data",i,"Pat40WithTP53.Rdata",sep=""))
filename = paste("MESGA.data",i,".Rdata",sep="")
rankingM = rep(0,100)
ResultMSEGA = NULL
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  result<-rbind(id.all,temp.score)
  ResultMSEGA[[iter]]<-result[,order(as.numeric(result[4,]), decreasing = TRUE)][,c(1:50)]
  
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac, ResultMSEGA,file=filename)
}

i<-62

for (i in c(50,53)){
  load(paste("/Volumes/Transcend/MEsimulation/rerunsimulation/1Simulation function/simulated.data",i,"Pat20.Rdata",sep=""))
  filename = paste("MESGA.data",i,".Rdata",sep="")
  rankingM = rep(0,100)
  ResultMSEGA = NULL
  for(iter in 1:100){
    mut.gene = mut.gene.mat[,iter]
    pat.mut.list = listof.pat.mut.list.test[[iter]]
    tempmut = sapply(pat.mut.list,rowSums)
    tempmut[tempmut==0] <- FALSE
    tempmut[tempmut>0] <- TRUE
    mutationMat = t(tempmut)
    gname = colnames(mutationMat)
    id.all = combn(gname,3)
    ncombn = dim(id.all)[2]
    mutationList = NULL
    for(cb in 1:ncombn){
      mutationList[[cb]]<- mutationMat[,id.all[,cb]]
    }
    temp.score = rep(0,ncombn)
    for(j in 1:ncombn){
      temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
      if(j %% 50 == 1 ){print(j)}
    }
    sort.score = sort(temp.score,decreasing = T)
    result<-rbind(id.all,temp.score)
    ResultMSEGA[[iter]]<-result[,order(as.numeric(result[4,]), decreasing = TRUE)][,c(1:50)]
    
    pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
    rankingM[iter] = which(sort.score==temp.score[pos])[1]
    print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
  }
  save(rankingM,frac, ResultMSEGA,file=filename)
}

load("simulated.data51Pat40.Rdata")
filename = "MESGA.data51.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)


load("simulated.data52Pat60.Rdata")
filename = "MESGA.data52.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data53Pat20.Rdata")
filename = "MESGA.data53.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data54Pat40.Rdata")
filename = "MESGA.data54.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data55Pat60.Rdata")
filename = "MESGA.data55.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data56Pat20WithTP53.Rdata")
filename = "MESGA.data56.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data57Pat40WithTP53.Rdata")
filename = "MESGA.data57.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data58Pat60WithTP53.Rdata")
filename = "MESGA.data58.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data59Pat20WithTP53.Rdata")
filename = "MESGA.data59.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data60Pat40WithTP53.Rdata")
filename = "MESGA.data60.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data61Pat60WithTP53.Rdata")
filename = "MESGA.data61.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data80333Pat80.Rdata")
filename = "MESGA.data62.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data80236Pat80.Rdata")
filename = "MESGA.data63.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data8033300Pat80WithTP53.Rdata")
filename = "MESGA.data64.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

load("simulated.data8023600Pat80WithTP53.Rdata")
filename = "MESGA.data65.Rdata"
rankingM = rep(0,100)
for(iter in 1:100){
  mut.gene = mut.gene.mat[,iter]
  pat.mut.list = listof.pat.mut.list.test[[iter]]
  tempmut = sapply(pat.mut.list,rowSums)
  tempmut[tempmut==0] <- FALSE
  tempmut[tempmut>0] <- TRUE
  mutationMat = t(tempmut)
  gname = colnames(mutationMat)
  id.all = combn(gname,3)
  ncombn = dim(id.all)[2]
  mutationList = NULL
  for(cb in 1:ncombn){
    mutationList[[cb]]<- mutationMat[,id.all[,cb]]
  }
  temp.score = rep(0,ncombn)
  for(j in 1:ncombn){
    temp.score[j] = funEstimate(mutationMat = mutationList[[j]],tol = 1e-7)$S
    if(j %% 50 == 1 ){print(j)}
  }
  sort.score = sort(temp.score,decreasing = T)
  pos=which(apply(id.all,2,function(x) all(x %in% mut.gene)))
  rankingM[iter] = which(sort.score==temp.score[pos])[1]
  print(c(paste("iteration",iter,sep=" "),rankingM[iter]))
}
save(rankingM,frac,file=filename)

