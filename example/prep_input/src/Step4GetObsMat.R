
stopifnot(file.exists(mutable))
load(mutable)
all = rbind(nonsilent.mutation.table,silent.mutation.table)
gene.name = unique(nonsilent.mutation.table[,1])
Pat.name = unique(all[,8])

G = length(gene.name)
N = length(Pat.name)
BigObsMat = matrix(rep(0,G*N),ncol=N)
rownames(BigObsMat) <- gene.name
colnames(BigObsMat) <- Pat.name
for(iter in 1:dim(nonsilent.mutation.table)[1]){
  #print(iter)
  tempG <- nonsilent.mutation.table[iter,1]
  tempP <- nonsilent.mutation.table[iter,8]
  BigObsMat[tempG,tempP] <- 1
}

save(BigObsMat,file=file.path(path,"obsMat.Rdata"))
