suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

parser <- OptionParser(usage = "%prog [maf file] [output directory]")
arguments <- parse_args(parser, positional_arguments = T);
if (length(arguments$args) != v2 ){
  print_help(parser)
  stop()
}

inputMAF <- arguments$args[1]
outputDir <- arguments$args[2]

inputTab = fread(inputMAF,header=T, stringsAsFactors=FALSE,data.table = F)

source("../Step1ForObtainingSilentAndNonsilentTable.R")
message("log: Running step 2")
source("src/Step2GetMutationRate.R")
message("log: Running step 3")
source("src/Step3GetBij.R")
message("log: Running step 4")
source("src/Step4GetObsMat.R")
message("log: Running step 5")
source("src/Step5GetMatrixOfRateAndMut.R")
