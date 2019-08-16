library(mescan)

##
# survey tg scores for geneset size n
##

load("MEScanInput.RData")
# calculate lambda
lambda = quantile(sqrt(mut_rate), probs=0.05)
# calculate tg scores to measure the mutual exclusiveness for all the size n genesets given m genes
set.seed(1)
# as an example, we randomly sample 5 genes from mut_tab
gene_index = sample(1:nrow(mut_tab),5)
# we then check the tg scores for geneset size = 3
size = 3
tgs = scan_genesets(gene_index,n=size,mut_data = mut_tab,rate_data = mut_rate,lambda = lambda)
# convert gene index to ID.
tg_mat = t(tgs$all_tgs)
tg_mat[,1:size] = rownames(mut_tab)[tg_mat[,1:size]]

