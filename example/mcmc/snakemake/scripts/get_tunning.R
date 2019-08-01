suppressPackageStartupMessages(library(mescan))
stopifnot_error <- function(err_message, ...)
{
  n <- length(ll <- list(...))
  for (i in 1:n)
    if (!(is.logical(r <- ll[[i]]) && all(r))) {
      stop(err_message)
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=5){
  stop("Rscript get_tunning.R data minSize maxSize iters tunning") 
}


is.number = grepl("^[[:digit:]]+L",args)
args[is.number] = as.numeric(args[is.number])

data = args[1]

minSize = args[2]
maxSize = args[3]
iters = args[4]
tunning = args[5]

stopifnot_error(file.exists(data),err_message="The ME input file does not exist!")
#tryCatch(stopifnot(file.exists(data),error=stop("The ME input file does not exist!")))


# Input matrix
load(data)

stopifnot_error(all(exists('mut.mat'),exists('rate.mut')),err_message="Check the data name!")

mut_data <- as.matrix(mut.mat)
rate_data <- as.matrix(rate.mut)
tg_lambda = quantile(sqrt(rate_data), probs=0.05)
sizes = minSize:maxSize
tmp = mcmc_tuning(mut_data,rate_data,
                 num_iters = iters,burning = 0,mcmc_tuning = tunning,
                 tg_lambda = tg_lambda,seed = -1, range=sizes)
