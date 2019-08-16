# Generate 2 files
# 1. Step2IntermediateData.Rdata is used for the simulation data generation.
# 2. Step2.rda

source("MADGiC.R")
library(abind)
mutable = inputTab
stopifnot(file.exists(mutable))
message(sprintf("Using mutation table %s",mutable))
path = outputDir
dir.create(path,recursive = TRUE)
message(sprintf("Saving output to %s",path))
load(mutable)

#nonsilent.mutation.table = nonsilent.mutation.table[1:5000,]
#silent.mutation.table = silent.mutation.table[1:2000,]
#######################################################
exome.file="exome_36.RData"
gene.rep.expr.file = "gene.rep.expr.RData"
gene.names.file="gene_names.txt"
prior.file="prior.RData"
alpha=0.2
beta=6
N=20
replication.file=NULL
expression.file=NULL
########################################################
### number of samples
S=length(unique(c(silent.mutation.table[,8], nonsilent.mutation.table[,8])))          

# get ensembl ids of all protein-coding genes
#all.gene.name <- as.vector(as.matrix(read.table( gene.names.file, header=FALSE)[,1]));
all.gene.name <- unique(c(nonsilent.mutation.table[,1],silent.mutation.table[,1]))
gene.ptr <- load(gene.rep.expr.file)
gene.rep.expr <- get(gene.ptr)
rm(gene.ptr)
gid <- unlist(sapply(gene.rep.expr, function(x) x[[1]]))

exome.ptr <- load(exome.file)
exome <- get(exome.ptr)
rm(exome.ptr)
exome.constants <- exome.SIFT <- exome.nonsil <- exome
for(i in 1:length(exome)){
  exome.constants[[i]] <- exome[[i]][,1:7]
  exome.SIFT[[i]] <- exome[[i]][,c(1,2,7:11)]
  exome.nonsil[[i]] <- exome[[i]][,c(1,2,7,12:15)]
}
rm(exome)

################################################################################
if(length(expression.file)>0){
  print(paste0("Loading expression file ", expression.file))
  expression <- read.table(expression.file, header=FALSE, stringsAsFactors=FALSE)
  colnames(expression) <- c("gid", "ex")
  matching <- sum(expression$gid %in% gid)
  print(paste0("Found ", matching, " genes with expression data"))
  
  # reorder to match gid object
  x <- match(gid, expression$gid)
  x <- x[!is.na(x)]
  expression <- expression[x,]
  
  # find tertiles of expression
  cutoffs <- quantile(expression$ex, c(0.333,0.666), na.rm=T)
  expression$cat <- 1
  expression$cat[expression$ex > cutoffs[1] & expression$ex < cutoffs[2]] <- 2
  expression$cat[expression$ex > cutoffs[2]] <- 3
  
  for (j in 1:length(exome.constants)){
    on.chr <- unlist(lapply(gene.rep.expr, function(x) x[[2]]==j))  # pull out all genes on chrom i
    gene.chr <- gene.rep.expr[on.chr]
    names <- unlist(lapply(gene.chr, function(x) x[[1]])) 
    x <- match(names, expression$gid)
    
    if(!is.null(gene.chr)) {
      cat <- expression$cat[x]
      gene.chr <- lapply(1:length(gene.chr), function(i) list(gene.chr[[i]][[1]], gene.chr[[i]][[2]], gene.chr[[i]][[3]], 
                                                              gene.chr[[i]][[4]], cat[i], gene.chr[[i]][[6]] ))
    }
    gene.rep.expr[on.chr] <- gene.chr
  }
}

if(length(replication.file)>0){
  print(paste0("Loading replication timing file ", replication.file))
  replication <- read.table(replication.file, header=FALSE, stringsAsFactors=FALSE)
  colnames(replication) <- c("gid", "repl")
  matching <- sum(replication$repl %in% gid)
  print(paste0("Found ", matching, " genes with replication timing data"))
  
  # reorder to match gid object
  x <- match(gid, replication$gid)
  x <- x[!is.na(x)]
  replication <- replication[x,]
  
  # find tertiles of replication
  cutoffs <- quantile(replication$repl, c(0.333,0.666), na.rm=T)
  replication$cat <- 1
  replication$cat[replication$repl > cutoffs[1] & replication$repl < cutoffs[2]] <- 2
  replication$cat[replication$repl > cutoffs[2]] <- 3
  
  for (j in 1:length(exome.constants)){
    on.chr <- unlist(lapply(gene.rep.expr, function(x) x[[2]]==j))  # pull out all genes on chrom i
    gene.chr <- gene.rep.expr[on.chr]
    names <- unlist(lapply(gene.chr, function(x) x[[1]])) 
    x <- match(names, expression$gid)
    
    if(!is.null(gene.chr)) {
      cat <- replication$cat[x]
      gene.chr <- lapply(1:length(gene.chr), function(i) list(gene.chr[[i]][[1]], gene.chr[[i]][[2]], gene.chr[[i]][[3]], 
        cat[i], gene.chr[[i]][[5]] , gene.chr[[i]][[6]] ))
}
    gene.rep.expr[on.chr] <- gene.chr
  }
}

#########################################################################
nosam= 1:length(gid)
for(i in 1:length(gid) )
  nosam[i]= sum(as.character(nonsilent.mutation.table[,1])==gid[i], na.rm=T)

nonsilent_passenger_gene=gid[nosam<=1]     ### names of genes having at most one nonsilent mutations
rm(nosam)

res=generate.sel.exome(nonsilent_passenger_gene,  gene.rep.expr,  exome.constants)   ### sequence of genes having at most one nonsilent mutatons
nonsil.type.const=preprocess.BM(res, gene.rep.expr)
rm(res)

all.gene.name=as.matrix(all.gene.name)         ### names of genes used for detecting silent mutations
sil.type.const=preprocess.BM(exome.constants, gene.rep.expr)  ### sequence of genes used for silent mutation detection - exome.constants

nonsilent_passenger_exclusive_gene= nonsilent_passenger_gene[ !(nonsilent_passenger_gene %in% all.gene.name )  ]
silent_exclusive_gene= all.gene.name[ !(all.gene.name %in% nonsilent_passenger_gene )  ]
both_gene= nonsilent_passenger_gene[ nonsilent_passenger_gene %in% all.gene.name   ]

rm(all.gene.name)

nonsil.mutab=nonsilent.mutation.table
sil.mutab= silent.mutation.table


##### calculate p_{7},p_{8}  ##########################################################################################################
x=nonsil.mutab[!is.na(match(nonsil.mutab[,1],nonsilent_passenger_gene)),4]
p_inframe=sum(x=="In_frame")/(sum(x=="In_frame")+ sum(x=="Frame_shift") )
p_frameshift=sum(x=="Frame_shift")/(sum(x=="In_frame")+ sum(x=="Frame_shift") )
if(p_inframe==0 | p_frameshift==0 | sum(x=="In_frame")+ sum(x=="Frame_shift") <=5 )
{
  p_inframe= 1/3
  p_frameshift= 2/3
}
rm(x)
rm(nonsilent.mutation.table)
#####  calculate p_{i}, i=1,2...,6 and selection bias r ##############################################################################

temp1 <- mut.type.converter(nonsil.mutab, exome.SIFT, nonsil.type.const[[1]], nonsil.type.const[[2]], gene.rep.expr)
nonsil.mut.type.sampl.sum <-temp1[[1]]
nonsil.mut.type <- temp1[[2]]
rm(temp1)

rm(silent.mutation.table)
temp2 <- mut.type.converter(sil.mutab, exome.SIFT, sil.type.const[[1]], sil.type.const[[2]], gene.rep.expr)
sil.mut.type.sampl.sum <- temp2[[1]]
sil.mut.type <- temp2[[2]]
rm(temp2)

# now *.mut.type.sampl.sum has a third dimension for the three expression categories
temp=  fit.background.model(nonsil.mutab, nonsil.mut.type.sampl.sum, sil.mut.type.sampl.sum, nonsil.type.const, sil.type.const, gene.rep.expr)
p=temp[[1]]
r=temp[[2]]  #selection bias
s=temp[[3]]  # rep timing
epsilon=temp[[4]] # expression
delta=temp[[5]]  # OR genes

rm(temp)

##### Calculating the likelihood of background mutations #############################################################################

exome.constants.mult=multiply.p.s.epsilon(exome.constants,p,s,epsilon, delta,gene.rep.expr)  # replace e_{k},f_{k},c_{k},d_{k} with (e_{k}*p_{t_{k}}, f_{k}*p_{v_{k}}, c_{k}*p_{t_{k}}, d_{k}*p_{v_{k}})*s_{k}
###save temp data
save(p,r,s,epsilon,delta,exome.constants,file = file.path(path,"Step2IntermediateData.Rdata"))
##################
rm(exome.constants)

res=generate.sel.exome(nonsilent_passenger_exclusive_gene, gene.rep.expr, exome.constants.mult)
nonsil.exclus=multiplyr(res,r)

sil.exclus=generate.sel.exome(silent_exclusive_gene,gene.rep.expr,exome.constants.mult)

res=generate.sel.exome(both_gene,gene.rep.expr,exome.constants.mult)
both=multiplyr(res,r)

a= mut.lik.change(nonsil.mut.type,nonsil.exclus,TRUE,FALSE,p_inframe,p_frameshift,p,r,s,epsilon,delta, gene.rep.expr)
b= mut.lik.change(nonsil.mut.type,both,TRUE,TRUE,p_inframe,p_frameshift,p,r,s,epsilon,delta,gene.rep.expr)
c= mut.lik.change(sil.mut.type,sil.exclus,FALSE,FALSE,p_inframe,p_frameshift,p,r,s,epsilon,delta,gene.rep.expr)
d= mut.lik.change(sil.mut.type,both,FALSE,TRUE,p_inframe,p_frameshift,p,r,s,epsilon,delta,gene.rep.expr)

B= rbind(a,b,c,d)
B=B[B[,2]!="0",]
muttable=B       ###  the change of the likelihood due to existing background mutations
rm(B)
rm(a)
rm(b)
rm(c)
rm(d)

A <- getA(nonsil.exclus, sil.exclus, both, p, r, s, epsilon, delta, gene.rep.expr)
A <- signif(A,digits=14)
rm(nonsil.exclus)
rm(sil.exclus)
rm(both)
uniqueA=unique(A)
tableA=as.vector(table(A))    ###  the likelihood of no background mutations
rm(A)

save.image(file = file.path(path,"Step2.rda"))