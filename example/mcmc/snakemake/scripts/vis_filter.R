load("/scratch/lji226/incoming/ChiWang/maybe_test/homosapien.ens2gs.rda")
suppressWarnings(library(gplots))
source("/scratch/lji226/incoming/ChiWang/maybe_test/scripts/preview_geneset.R")


cosmic = read.delim("/scratch/lji226/incoming/ChiWang/maybe_test/synapse/Census_allThu Oct 25 17_33_53 2018.tsv", stringsAsFactors=F)
head(cosmic$GeneSymbol)
args = commandArgs(trailingOnly=TRUE)
load(args[1])
mut.m=mut.mat
all_gs = NULL

inputFile=args[2]
con  <- file(inputFile, open = "r")
dataList <- list()
counter = 0
scores=NULL
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    #print(oneLine)
    tmp <- (strsplit(oneLine,' '))[[1]]
    len = length(tmp)
    #if((len-2) < 4) next
    score <- tmp[len]
    myVector <- tmp[1:(len-2)]
    dataList <- c(dataList,list(myVector))
    scores <- c(scores,as.numeric(score))
    counter = counter + 1

 } 
 close(con)

index = order(scores,decreasing=T)

cl=list()
for(a in index){
    cl = c(cl,dataList[a])
}

pdf(gsub(".txt$",'.vizCosmicNew.pdf',inputFile),h=4,w=20)
all_gs = NULL
glist = NULL
limit = ifelse(length(cl)<as.numeric(args[3]),length(cl),as.numeric(args[3]))
for(i in 1: limit){

    c1 = as.numeric(cl[[i]])+1
    tmp = rownames(mut.m[c1,])

    print(tmp)
    tmp1 = homo_dict[tmp]
    tmp1 = paste(c(tmp1,scores[index[i]]),collapse=' ')
    all_gs = c(all_gs,tmp1)
    glist= c(tmp,glist)
    preview_gs(mut.m[c1,],NULL,CosmicGene = cosmic$GeneSymbol,output = NA,map=homo_dict) 
}
dev.off()
write.table(all_gs,gsub(".txt$",'.GeneSymbol.txt',inputFile),row.names=F,col.names=F,quote=F,sep='')
# glist = unique(glist)

# hgsn = homo_dict[glist]
# out = cbind("ensembl" = glist, "genesymbol" = hgsn)
#write.table(out,'ov_map.csv',row.names=F,quote=F,sep=',')
