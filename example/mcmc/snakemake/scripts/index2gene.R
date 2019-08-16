
args = commandArgs(trailingOnly=TRUE)
MEInput=args[1]
load(MEInput)
all_gs = NULL

inputFile=args[2]
con  <- file(inputFile, open = "r")
dataList <- list()
scores=NULL
geneid=rownames(mut.mat)
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0 ) {
    tmp <- (strsplit(oneLine,' '))[[1]]
    len = length(tmp)
    score <- as.numeric(tmp[len])
    genes <- geneid[as.numeric(tmp[-c(len,len-1)])+1]
    genes = paste(genes,collapse=';')
    dataList <- c(dataList,list(c(genes,score)))

} 
close(con)

output = do.call(rbind,dataList)

write.table(output,gsub(".txt$",'.genes.txt',inputFile),row.names=F,col.names=F,quote=F,sep='\t')
