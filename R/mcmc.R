#' Perform Markov Chain Monte Carlo (MCMC) sampling using
#' Metropolis-Hasting algorithm.
#' @param mut_data Mutation table created using MADGIC algorithm.
#' @param mut_rate Mutation rate table create uisng MADGIC algorithm.
#' @param num_iters Number of MCMC samples.
#' @param burning Number of Burn-in samples.
#' @param mcmc_tuning MCMC tuning parameter.
#' @param tg_lambda lambda.
#' @param geneset_size Size of the gene sets to sample.
#' @param seed Seed for RGN.
#' @param previous_stats Recovery from the previous MCMC chain file.
#' @return Tab-delimited file contain gene sets and TG scores

mcmc_core <- function(mut_data, mut_rate, num_iters = 1000,
    burning = 0, mcmc_tuning, tg_lambda, mcmc_output,
    geneset_size, seed = -1, previous_stats = NULL,
    debug = -1) {
    genes <- -1
    if (!is.null(previous_stats)) {
        stats <- read.delim(previous_stats, header = F,
            stringsAsFactors = F, sep = "\t")
        stopifnot(geneset_size == stats[2, 1])
        seed <- stats[1, 1]
        genes <- as.numeric(strsplit(stats[3, 1], split = " ")[[1]])
        stopifnot(geneset_size == length(genes))
    }
    .C("mcmc", as.integer(as.vector(mut_data)), as.integer(dim(mut_data)),
        as.vector(as.numeric(mut_rate)), as.integer(dim(mut_rate)),
        as.integer(num_iters), as.integer(burning),
        as.numeric(mcmc_tuning), as.numeric(tg_lambda),
        tgscores = as.numeric(rep(0, num_iters)), gfile = as.character(mcmc_output),
        max_geneset = as.integer(rep(0, nrow(mut_data))),
        max_iter = as.integer(0), size = as.integer(geneset_size),
        seed = as.integer(seed), debug = as.integer(debug),
        genes = as.integer(genes))
}

#' Select tuning parameters for Markov Chain Monte Carlo (MCMC) sampling
#' based on the size of the gene sets
#' @param mut_data Mutation table created using MADGIC algorithm.
#' @param mut_rate Mutation rate table create uisng MADGIC algorithm.
#' @param num_iters Number of MCMC samples.
#' @param burning Number of Burn-in samples.
#' @param mcmc_tuning MCMC initial tuning parameter, .
#' @param tg_lambda lambda
#' @param seed Seed for RGN.
#' @return tuning parameter after the selection

mcmc_tuning <- function(mut_data, mut_rate, num_iters = 1000,
    burning = 0, mcmc_tuning, tg_lambda, seed = -1,
    range = c(2:3)) {
    .C("select_tuning", as.integer(as.vector(mut_data)),
        as.integer(dim(mut_data)), as.vector(as.numeric(mut_rate)),
        as.integer(dim(mut_rate)), as.integer(num_iters),
        as.numeric(mcmc_tuning), as.numeric(tg_lambda),
        seed = as.integer(seed), min = as.integer(min(range)),
        max = as.integer(max(range)))
}

