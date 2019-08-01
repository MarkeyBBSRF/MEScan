suppressPackageStartupMessages(library(mescan))
library(locfdr)

c_batch_run_tgscore <- function(gs_idx,size=3,n=10, mut_data,mut_rate,lambda,seed=1){
  if(n<2||n%%1!=0)
  {return("number of gene in the sets is not valid")}
  samples = matrix(0,nrow=n,ncol=size)
  for(i in 1:n){
    samples[i,] = sample(gs_idx,size)
  }

  samples = samples - 1
  rc = rc = dim(mut_data)
  mut_data <- mut_data
  mut_rate <- mut_rate
  tgscores=as.double(rep(0,nrow(samples)))
  tmp = .C("r_calculate_TG",
           as.integer(as.vector(mut_data)),
           as.integer(rc),
           as.vector(as.numeric(mut_rate)),
           as.integer(rc),
           as.integer(as.vector(as.matrix(samples))),
           as.integer(dim(samples)),
           as.double(lambda),tgscores=as.double(rep(0,nrow(samples)))
  )
  Tvec = tmp$tgscores
  return(Tvec)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=6){
  stop("Rscript get_threshold.R Data min max prefix nt df") 
}

library(doMC)


################
f = args[1]
load(f)

mut_data <- as.matrix(mut.mat)
rate_data <- as.matrix(rate.mut)

seeds = sample(1:1e5,1)
ns = as.numeric(args[5])
df = 0
if(is.null(args[6])){
  df = 15
}else{
  df = as.numeric(args[6])
}
genesetSize = as.numeric(args[2]):as.numeric(args[3])
# registerDoMC(args[5])
# res = foreach(size = genesetSize,.combine = rbind) %dopar% 
# {
fileConn<-file(args[4],'w')
cutVec = NULL
processed = NULL
#res = matrix(0,nrow=length(genesetSize),ncol = 2)
for (size in genesetSize){
  print(size)
  set.seed(seeds)
  tg_lambda = quantile(sqrt(rate_data), probs=0.05)
  TGc = c_batch_run_tgscore(1:nrow(mut_data),size=size,n=ns,mut_data,rate_data,tg_lambda,seed=seeds)
  #save(new_TGc, seeds, file= paste("random_TGvec_",size,".Rdata",sep=""))
  TGvec = TGc[TGc>0]
  fdr.res = try(locfdr(zz=TGvec,nulltype = 1,df=df,bre=500,mlests = c(median(TGvec),1)),
    silent = T) #pick the highest degree without error?
  if(class(fdr.res) == "try-error"){

    next
  }else{
    processed = c(processed,size)
  }
  # ########Fit CDF and get FDR 
  x = fdr.res$mat[,"x"]
  f0 = fdr.res$mat[,"f0"]
  f1 = fdr.res$mat[,"f"]
  fdr = fdr.res$mat[,"fdr"]
  p0 = fdr.res$fp0[3,3]; p1=1-p0
  fdrtheo = fdr.res$mat[,"fdrtheo"]
  upper = f1*fdr
  lower = f1
  len = length(x)
  FDR = rep(0,len)

  for(i in 1:499){
    FDR[i] = sum(upper[i:len])/sum(lower[i:len])
  }
  #############
  cut.id=which(FDR<0.05)[1]
  #sum(TGvec>=x[cut.id])/length(TGvec)
  Cutoff = x[cut.id]
  writeLines(sprintf("%d\t%.3f",size,Cutoff), fileConn,sep="\n")
  cutVec = c(cutVec,Cutoff)
}
print(cutVec)
if ( length(processed) < 2 ){
  print("Only processed less than 2 sizes")
  break
}
# lmR = lm(cutVec~processed)
# slmR = summary(lmR)$coefficients
# coef = slmR[2,1]
# ix = slmR[1,1]
# for(i in 2:10){
#   if (i %in% processed) next
#   Cutoff = ix+coef*i
#   print (Cutoff)
#   writeLines(sprintf("%d\t%.3f",i,Cutoff), fileConn,sep="\n")
# }

close(fileConn)




