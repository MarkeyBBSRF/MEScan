load(file.path(path,"obsMat.Rdata"))
#load(file.path(path,'Bij.Rdata'))
bij = bijAndQvec$bij
gene.idindex = which(sapply(bij,function(x) sum(is.na(x)))==0)
gname1 = gid[gene.idindex]
#########################
gene.id2 = which(rowSums(BigObsMat)>=2)
gname.from.mat = rownames(BigObsMat)
gname2 = gname.from.mat[gene.id2]


################## The genes we used for finding pattern
use.gene = intersect(gname1,gname2)
ifelse(length(use.gene) < 2000,warning("# of genes < 2000"),sprintf("# of genes = %d",length(use.gene)))
bij.use = bij[use.gene]
#########################
mut.mat = BigObsMat[use.gene,]
#########################
G0 = length(bij.use)
N0 = length(bij.use[[1]])
rate.mut = matrix(rep(0,G0*N0),ncol=N0)
for(g in 1:G0){
  rate.mut[g,] = 1-bij.use[[g]]
}
rownames(rate.mut) <- rownames(mut.mat)
colnames(rate.mut) <- colnames(mut.mat)

save(mut.mat,rate.mut,file = file.path(path,"MEScanInput.Rdata"))
