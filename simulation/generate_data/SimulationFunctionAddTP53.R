source("GenerateMutationTable.R")
## Simulations are different cause we add the nonsilent mutation with extrem case, type of highest mutation rate/lowest mutation rate
load("Bij.rdata") #load("bijAndQvec.Rdata")
load("pat.rate.list.Rdata")
load("ov.simulation.table.rdata")#load("final.table.Rdata")

#####################################
bij <- bijAndQvec$bij
rm(bijAndQvec)
G0 = length(bij)
N0 = length(bij[[1]])
Bm = matrix(rep(0, G0 * N0), ncol = N0)
for (g in 1:G0) {
  Bm[g, ] <- bij[[g]]
}
############################# find the names of genes
gene.rep.expr.file = "gene.rep.expr.RData"
gene.ptr <- load(gene.rep.expr.file)
gene.rep.expr <- get(gene.ptr)
rm(gene.ptr)
gid <- unlist(sapply(gene.rep.expr, function(x) x[[1]]))
########################################### 
sample.name = names(bij[[1]])
#################################################### 
rownames(Bm) <- gid
colnames(Bm) <- sample.name
############################################### 
use.id = which(rowSums(is.na(Bm)) == 0)
Bmat = Bm[use.id, ]
EtaMat = 1 - Bmat
mut.gene.name1 = rownames(EtaMat)


#################Creat simulation setting:
Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=20,frac=c(1/3,1/3,1/3),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=56,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=40,frac=c(1/3,1/3,1/3),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=57,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=60,frac=c(1/3,1/3,1/3),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=58,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=20,frac=c(1/2,1/3,1/6),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=59,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=40,frac=c(1/2,1/3,1/6),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=60,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=60,frac=c(1/2,1/3,1/6),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=61,
               EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

#####################
Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=80,frac=c(1/2,1/3,1/6),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=8023600,
                   EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)
###
Res = SimulGenTp53(seed.pat=12345,nPat=200,NumofaddPat=80,frac=c(1/3,1/3,1/3),NoiseGeneNum=17,extrem.ind.mut=5,trialNum=8033300,
                   EtaMat=EtaMat,pat.name=pat.name,final.table = final.table, mut.gene.name1=mut.gene.name1)

##############################
## seed.pat = set seed for picking patients
##nPat = # of patient picked for simulation
## frac = fraction of three mutually exclusive genes
## NoiseGeneNum = # of noise gene added
## extrem.ind.mut = minimun # of noise mutations that we want for noise genes
## trialNum = trials name 
# TP53 = "ENSG00000141510" was inside the function
SimulGenTp53 <- function(seed.pat,nPat,NumofaddPat=NumofaddPat,frac=frac,NoiseGeneNum=17,extrem.ind.mut=5,trialNum=101,
                     EtaMat=EtaMat,pat.name=pat.name,final.table=final.table,mut.gene.name1=mut.gene.name1){
  ##################################random pick 400 patients here!!!!!!!!!!!
  pat.name = colnames(EtaMat)
  set.seed(seed.pat)
  pat.pick = sample(pat.name, nPat)
  final.table = final.table[final.table$Tumor_Sample_Barcode %in% pat.pick,]  ## only look at these pick patients' mutations
  ############################### find all genes that have at least 2 mutations
  atleast2.pos = which(table(final.table[, 1]) >= 2)
  mut.gene.name2 <- names(table(final.table[, 1])[atleast2.pos])
  ###################### All the genes that have at least 2 mutations and no NA values of rate
  mut.gene.name = intersect(mut.gene.name1, mut.gene.name2)
  EtaMat  = EtaMat[mut.gene.name,pat.pick]     #update EtaMat for only picked patients 
  ############################################################
  ############################################################
  ############################################################
  NumofaddPat = NumofaddPat
  listof.pat.rate.list.test <- NULL
  listof.pat.mut.list.test <- NULL
  listof.candi.gene <- NULL
  BigObsMat.list <- NULL
  mut.gene.mat <- matrix(rep(NA,3*100),ncol=100)
  frac = frac
  TP53 = "ENSG00000141510"
  ########### random pick three genes for mutually exclusive pattern
  for (iter in 1:100) {
    set.seed(iter)
    mut.gene = sample(setdiff(mut.gene.name,TP53), 3)
    mut.gene.mat[,iter]<-mut.gene
    ############################################################################################# 
    mut.gene.pos = which(final.table[, 1] %in% mut.gene)  # pos of mut.gene in table
    exist.pat = unique(final.table[mut.gene.pos, 8])  # pat already have mutation
    NumofMut = table(final.table[mut.gene.pos, 1])  # How many mutations already have for each gene 
    mut.gene.reorder = names(NumofMut)
    #################################################
    set.seed(iter)
    add.pat = sample(setdiff(pat.pick, exist.pat), NumofaddPat)
    pos.mut = lapply(mut.gene.reorder,function(x)(which(final.table[,1]==x)))
    i1 = as.integer(frac*NumofaddPat)[1]
    i2 = as.integer(frac*NumofaddPat)[2]
    i3 = NumofaddPat -i1 - i2
    add.pos.table1=NULL
    add.pos.table2=NULL
    add.pos.table3=NULL
    final.tableG1 = final.table[pos.mut[[1]],]
    final.tableG2 = final.table[pos.mut[[2]],]
    final.tableG3 = final.table[pos.mut[[3]],]
    len1 = dim(final.tableG1)[1]; len2 = dim(final.tableG2)[1]; len3 = dim(final.tableG3)[1] 
    temp1 = do.call("rbind",rep(list(final.tableG1),as.integer(i1/len1)+1))
    add.pos.table1 = temp1[1:i1,]
    temp2 = do.call("rbind",rep(list(final.tableG2),as.integer(i2/len2)+1))
    add.pos.table2 = temp2[1:i2,]
    temp3 = do.call("rbind",rep(list(final.tableG3),as.integer(i3/len3)+1))
    add.pos.table3 = temp3[1:i3,]
    #############################################
    add.pos.table = rbind(add.pos.table1, add.pos.table2, add.pos.table3)
    #### 
    add.pos.table[, 8] <- add.pat
    simu.final.table = rbind(final.table, add.pos.table)  # simulated mutation table 
    ###################### Organize it to be the form used in TGscore function
    gen.name = rownames(pat.rate.list[[1]])
    pat.mut.list = Reorganize.pat.mut(final.table = simu.final.table, gen.name = gen.name, Pat.name = pat.pick)
    ### load data############## the data have: pat.mut.list simu.final.table mut.gene NumofMut add.pat atleast2.gene.id = which(
    ### rowSums(Reduce('+',pat.mut.list))>=2 ) get the genes have at least 2 mutations and not NA mutation rate
    pat.mut.list <- lapply(pat.mut.list, function(x) x[mut.gene.name, ])
    pat.rate.list = pat.rate.list[pat.pick]
    pat.rate.list <- lapply(pat.rate.list, function(x) x[mut.gene.name, ])
    ############################ 
    gen.name = rownames(pat.mut.list[[1]])
    ### Generate big observation table
    temp.mut = lapply(pat.mut.list, function(x) as.numeric(rowSums(x) > 0))
    BigObsMat = do.call(cbind, temp.mut)
    rownames(BigObsMat) <- gen.name
    colnames(BigObsMat) <- pat.pick
    BigObsMat.list[[iter]] <- BigObsMat
    #######################################
    mutual.ind = which(gen.name %in% mut.gene)
    ####################################### 
    extrem.ind = which(rowSums(BigObsMat) >=extrem.ind.mut & rowSums(BigObsMat)<=nPat*0.5)
    TP53.ind = as.numeric(which(rowSums(BigObsMat)>=nPat*0.5))
    ####################################### 
    set.seed(iter)
    add.ind = sample(setdiff(extrem.ind, mutual.ind), NoiseGeneNum)
    replace.pos = sample(1:length(add.ind),1)
    add.ind[replace.pos] <- TP53.ind
    ind = c(mutual.ind, add.ind)
    ###################################
    temp.pat.rate.list.test = lapply(pat.rate.list, function(x) x[ind, ])
    temp.pat.mut.list.test = lapply(pat.mut.list, function(x) x[ind, ])
    candi.geneTemp = rownames(temp.pat.mut.list.test[[1]])
    listof.candi.gene[[iter]] <- candi.geneTemp
    listof.pat.rate.list.test[[iter]] <- temp.pat.rate.list.test
    listof.pat.mut.list.test[[iter]] <- temp.pat.mut.list.test
  } 
  ########################
  filename = paste("simulated.data",trialNum,"Pat",NumofaddPat,"WithTP53.Rdata",sep="")
  save(seed.pat,NumofaddPat,frac,listof.pat.mut.list.test, listof.pat.rate.list.test, 
       mut.gene.mat,NumofMut, BigObsMat.list, EtaMat,listof.candi.gene,file = filename) 
  ########################
  listRes = list(NumofaddPat=NumofaddPat,frac=frac,EtaMat=EtaMat,mut.gene.mat=mut.gene.mat,
                 listof.pat.mut.list.test=listof.pat.mut.list.test,listof.pat.rate.list.test=listof.pat.rate.list.test,
                 BigObsMat.list=BigObsMat.list,listof.candi.gene=listof.candi.gene)
  #########################
  return(listRes)
}







