#### Example form of final.table input!!!!!!!!!!!!!!!!!!!!!!! The input need to be exact!!!!!!  Ensembl_gene_id Chrom Start_Position
#### Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode dCG type.vec type.reorder 1 ENSG00000148584 10
#### 52258002 SNP C C G TCGA-59-2352-01A-01W-0799-08 0 4 4 2 ENSG00000166535 12 8882236 SNP T T A TCGA-61-2002-01A-01W-0722-08 0 2 2 3
#### ENSG00000166535 12 8907744 SNP C C G TCGA-09-2045-01A-01W-0799-08 0 4 4 This function converge the above table to a list of 463 patients
#### and each is a matrix G*7 with indicator of having mutation or not.(1/0)
Reorganize.pat.mut <- function(final.table = final.table, gen.name = gen.name, Pat.name = Pat.name) {
    ### Some of the gene in final.table is not in the gene sets gen.name, remove them!!!
    all.mut.gene = unique(final.table[, 1])
    ind.notin = !(all.mut.gene %in% gen.name)
    gene.remove = all.mut.gene[ind.notin]
    ind.remove = which(final.table[, 1] %in% gene.remove)
    ####################### 
    if (length(ind.remove) == 0) {
        final.table.update = final.table
    } else {
        final.table.update = final.table[-ind.remove, ]
    }
    ########## Generate pat.mut.list!!!!!!!!!!!!!!!!!!!!!!
    pat.mut.list <- rep(list(0), length(Pat.name))
    names(pat.mut.list) <- Pat.name
    pat.row = final.table.update[, 8]
    NumofG = length(gen.name)
    ############################################ 
    for (pat in Pat.name) {
        ind.pos <- which(pat.row %in% pat)
        tempMat <- matrix(rep(0, NumofG * 7), ncol = 7)
        rownames(tempMat) <- gen.name
        ############################################################ 
        if (length(ind.pos) == 0) {
            pat.mut.list[[pat]] <- tempMat
        } else {
            ######################################################## 
            indOfMutation <- final.table.update[ind.pos, c(1, 11)]
            tempMat <- matrix(rep(0, NumofG * 7), ncol = 7)
            rownames(tempMat) <- gen.name
            for (i in 1:(dim(indOfMutation)[1])) {
                tempMat[indOfMutation$Ensembl_gene_id[i], indOfMutation$type.reorder[i]] <- 1
            }
            pat.mut.list[[pat]] <- tempMat
        }
        print(pat)
    }
    return(pat.mut.list)
} 
