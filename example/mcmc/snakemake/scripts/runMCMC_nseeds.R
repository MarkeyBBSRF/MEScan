suppressPackageStartupMessages(library(mescan))
anyNA <- function (x) 
{
    any(is.na(x))
      }
stopifnot_error <- function(err_message, ...)
{
  n <- length(ll <- list(...))
  for (i in 1:n)
    if (!(is.logical(r <- ll[[i]]) && !anyNA(r) && all(r))) {
      stop(err_message)
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=7){
  stop("Rscript runMCMC_nseeds data genesetSize iters tunning burning nseeds output") 
}


# is.number = grepl("^[[:digit:]]+L",args)
# args[is.number] = as.numeric(args[is.number])

print(args)
data = args[1]

geneset_size = as.numeric(args[2])
num_iters = as.numeric(args[3])
mcmc_tuning = as.numeric(args[4])
burning = as.numeric(args[5])
nseeds = as.numeric(args[6])
output_dir = args[7]
seeds = 1:nseeds-1
# MCMC parameter
re.dir <- output_dir


# Input matrix
load(data)
stopifnot_error(all(exists('mut.mat'),exists('rate.mut')),err_message="Check the data name!")

mut_data <- as.matrix(mut.mat)
rate_data <- as.matrix(rate.mut)

tg_lambda = quantile(sqrt(rate_data), probs=0.05)


if(num_iters %/% 1000000>=1){
    pretty_lab = paste0((num_iters %/% 1000000),'M')
}else if (num_iters %/% 1000>=10){
    pretty_lab = paste0((num_iters %/% 1000),'K')
}else{
    pretty_lab = num_iters
}


#---------------------SETUP ENV------------------#
wrapper <- function(x){
  seed = x
  genesets.f = paste0(re.dir,'/raw/',"gs",geneset_size,".",
    "s",seed,'.raw')
  genesets.stat = gsub("raw$","stat",genesets.f)
  dir.create(paste0(re.dir,'/raw/'),showWarnings=F,recursive=T)
  re <- mcmc_core(mut_data,rate_data,
                 num_iters = num_iters,burning = 0,mcmc_tuning = mcmc_tuning,
                 tg_lambda = tg_lambda,
                 geneset_size = geneset_size,
                 mcmc_output = genesets.f,
                 seed = -1,debug=0)
  rm(re)

}

library(snowfall)

bParallel=TRUE
sfInit( parallel=bParallel, cpus=length(seeds) )
stopifnot( sfCpus() == length(seeds) )
stopifnot( sfParallel() == bParallel )
sfExportAll()
sfLibrary(mescan)
tmp = sfLapply(seeds, wrapper)
sfStop()

