library(mescan)

#load("MEInput.Rdata")

# mut_tab = mut.mat
# mut_rate = rate.mut
# genes = sample(rownames(mut_tab),100)
# pat = sample(colnames(mut_tab),100)
# mut_tab = mut_tab[genes,pat]
# mut_rate = mut_rate[genes,pat]
# write.table(mut_tab,'mut.table.txt',sep='\t',row.names=T,col.names=T,quote=F)
# write.table(mut_rate,'mut.rate.txt',sep='\t',row.names=T,col.names=T,quote=F)

##
# survey tg scores for geneset size n
##

mut_tab = read.delim("mut.table.txt",row.names = 1)
mut_rate = read.delim("mut.rate.txt",row.names = 1)
mut_tab = as.matrix(mut_tab)
mut_rate = as.matrix(mut_rate)
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

