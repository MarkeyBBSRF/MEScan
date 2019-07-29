#' Return TG scores for the gene sets combinations from \code{gs_idx},
#' with size of \code{n}.
#' @param gs_idx A vector of genes index
#' @param n The size of the gene sets
#' @param mut_data Mutation table created using MADGIC algorithm.
#' @param mut_rate Mutation rate table create uisng MADGIC algorithm.
#' @param lambda lambda.
#' @return A list contains the gene set with maximun TG score, the maximum TG score and all the scanned gene sets
scan_genesets <- function(gs_idx, n = 2, mut_data, rate_data,
    lambda) {
    if (n < 2 || n > length(gs_idx) || n%%1 != 0) {
        return("number of gene in the sets is not valid")
    }

    comb <- combn(gs_idx, n)
    comb.t <- t(comb - 1)
    rc <- rc <- dim(mut_data)
    mut_data <- mut_data
    mut_rate <- rate_data

    tmp <- .C("r_calculate_TG", as.integer(as.vector(mut_data)),
        as.integer(rc), as.vector(as.numeric(mut_rate)),
        as.integer(rc), as.integer(as.vector(as.matrix(comb.t))),
        as.integer(dim(comb.t)), as.double(lambda),
        tgscores = as.double(rep(0, nrow(comb.t))))
    Tvec <- tmp$tgscores
    Larg_TG <- max(Tvec)  #maximum test statistic
    index <- which(Tvec == Larg_TG)  #index of maximum test statistic
    Opt_set <- comb.t[index, ] + 1  #best combination
    results <- rbind(comb, Tvec)  #combination and Test statistic together
    results <- results[, order(results[n + 1, ], decreasing = T)]
    return((list(max_tg = Larg_TG, mat_tg_geneset = Opt_set,
        all_tgs = results)))
}
