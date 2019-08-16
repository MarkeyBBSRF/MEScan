source("../MADGiC.R")
library("abind")
library("biomaRt")
library(data.table)
proj = Proj

###############################
maf$Tumor_Sample_Barcode = sub("-(?:[^-]*-){3}[^-]*$","",maf$Tumor_Sample_Barcode,perl=T)
exome.file="../exome_36.RData"
gene.rep.expr.file = "../gene.rep.expr.RData"
gene.names.file="../gene_names.txt"
prior.file="../prior.RData"
alpha=0.2
beta=6
N=20
replication.file=NULL
expression.file=NULL
##############################
maf$Strand[maf$Strand=="-1"]="-"
maf$Strand[maf$Strand=="1"]="+"
maf.table <- maf 
for(i in 1:ncol(maf.table))
  maf.table[,i]=as.character(maf.table[,i])
maf.table=as.matrix(maf.table)

# remove entries with missing chromosome names
maf.table <- maf.table[!is.na(maf.table[,5]),]
maf.table[,5] = gsub("^chr","",maf.table[,5])

# ensembl version depends on NCBI build
ncbi.version <- as.character(unique(maf.table[,4]))
if(grepl("38",ncbi.version)) ncbi.version="38"
if(grepl("37",ncbi.version)) ncbi.version="37"
if(grepl("36",ncbi.version)) ncbi.version="36"

require(rtracklayer)

if (ncbi.version=="38") {
  is.change <- 1:nrow(maf.table)
  # convert 38 to 37
  chain <- import.chain("../hg38ToHg19.over.chain")
  gr <- GRanges(
    seqnames=paste("chr",as.character(maf.table[is.change,5]), sep=""),
    ranges=IRanges(start=as.numeric(as.character(maf.table[is.change,6])), end=as.numeric(as.character(maf.table[is.change,7]))),
    strand=maf.table[is.change,8])
  grnew <- liftOver(gr, chain)
  drop = which(elementNROWS(grnew)>1)
  cat("Remove those genes",maf.table[drop,1])
  keep = which(elementNROWS(grnew)<2)
  maf.table[keep,5] <- substring(as.character(seqnames(grnew)[keep]), first=4)
  maf.table[keep,6] <- as.numeric(start(grnew)[keep])
  maf.table[keep,7] <- as.numeric(end(grnew)[keep])
  ncbi.version <- "37"
}

if (ncbi.version=="37") {
  
  # convert 37 to 36
  chain <- import.chain("../hg19ToHg18.over.chain")
  maf.table = maf.table[!is.na(maf.table[,6]),]
  is.37 <- 1:nrow(maf.table)
  gr <- GRanges(
    seqnames=paste("chr",as.character(maf.table[is.37,5]), sep=""),
    ranges=IRanges(start=as.numeric(as.character(maf.table[is.37,6])), end=as.numeric(as.character(maf.table[is.37,7]))),
    strand=maf.table[is.37,8])
  grnew <- liftOver(gr, chain)
  drop = which(elementNROWS(grnew)>1)
  cat("Remove those genes",maf.table[drop,1])
  keep = which(elementNROWS(grnew)<2)
  maf.table[keep,5] <- substring(as.character(seqnames(grnew)[keep]), first=4)
  maf.table[keep,6] <- as.numeric(start(grnew)[keep])
  maf.table[keep,7] <- as.numeric(end(grnew)[keep])
  ncbi.version <- "36"
}
rm(chain)
rm(gr)
rm(grnew)

maf.table <- maf.table[!is.na(maf.table[,5]),]
nonsilent.maf.table <- maf.table[maf.table[,9] != "Silent",]   ### table of nonsilent mutations
silent.maf.table <- maf.table[maf.table[,9] == "Silent",]   ### table of nonsilent mutations
rm(maf.table)

if (ncbi.version == "36") { 
  ens.version <- "may2009.archive.ensembl.org"
  ensembl=useMart("ENSEMBL_MART_ENSEMBL",host=ens.version, dataset = "hsapiens_gene_ensembl") 
} else {
  stop(paste("build number ", ncbi.version, " not supported", sep=""))
}

# change coding of indels in nonsilent.maf.table
type.col <- which(colnames(nonsilent.maf.table)=="Variant_Type")
class.col <- which(colnames(nonsilent.maf.table)=="Variant_Classification")

which.indels <- which(nonsilent.maf.table[,type.col] %in% c("INS", "DEL"))
if (length(which.indels)>0) {
  indel.labels <- nonsilent.maf.table[which.indels, class.col]
  if (sum(indel.labels %in% c("Frame_Shift_Del", "Frame_Shift_Ins")) > 0) indel.labels[indel.labels %in% c("Frame_Shift_Del", "Frame_Shift_Ins")] <- "Frame_shift"
  if (sum(indel.labels %in% c("In_Frame_Del", "In_Frame_Ins")) > 0)  indel.labels[indel.labels %in% c("In_Frame_Del", "In_Frame_Ins")] <- "In_frame"
  nonsilent.maf.table[which.indels, type.col] <- indel.labels
}

rm(indel.labels)
Extra.type <- which(!(nonsilent.maf.table[,10] %in% c("SNP", "DNP", "TNP", "ONP", "Frame_shift", "In_frame")))
for (i in 1:length(Extra.type)){
  alleles <- c(nonsilent.maf.table[Extra.type[i],11], nonsilent.maf.table[Extra.type[i], 12],nonsilent.maf.table[Extra.type[i], 13]) 
  maxchar <- max(nchar(alleles))
  minchar <- min(nchar(alleles))
  if (("-" %in% alleles) | minchar != maxchar) {
    nonsilent.maf.table[Extra.type[i], 10] <- "Frame_shift"
    if (maxchar %% 3 == 0) nonsilent.maf.table[Extra.type[i], 10] <- "In_frame"
  } else{
    nonsilent.maf.table[Extra.type[i], 10] <- "SNP"
    if (maxchar > 1) nonsilent.maf.table[Extra.type[i], 10] <- "DNP"
    if (maxchar > 2) nonsilent.maf.table[Extra.type[i], 10] <- "TNP"
    if (maxchar > 3) nonsilent.maf.table[Extra.type[i], 10] <- "ONP"
  }
}
rm(Extra.type)
rm(alleles)
rm(maxchar)

# save only what is needed for computations
silent.mutation.table <- cbind(Ensembl_gene_id=
    convert.hgnc.to.ensembl(silent.maf.table,ensembl), 
    silent.maf.table[,c(5,6,10:13,16)])
nonsilent.mutation.table <- cbind(Ensembl_gene_id=
    convert.hgnc.to.ensembl(nonsilent.maf.table,ensembl), 
    nonsilent.maf.table[,c(5,6,10:13,16)])
silent.mutation.table <- silent.mutation.table[silent.mutation.table[,2] %in% c(1:24,"X","Y"),]
nonsilent.mutation.table <- nonsilent.mutation.table[nonsilent.mutation.table[,2] %in% c(1:24,"X","Y"),]
# remove NA emsembl_gene_id rows
silent.mutation.table <- silent.mutation.table[!is.na(silent.maf.table[,1]),]
nonsilent.mutation.table <- nonsilent.mutation.table[!is.na(nonsilent.mutation.table[,1]),]
###save mutation table data
save(silent.mutation.table, nonsilent.mutation.table,file= sprintf("mutationtable.%s.Rdata",proj))
