#!/usr/bin/env Rscript

# ----------------- deal with argument -------------------
# argument parsing for R, adapted from 
# https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
input=args[1] # Conf.csv
output=args[2] # MaxSpecies.R

# ----------------- do stuff -------------------
Conf <- read.csv(input,header=TRUE,row.names=1)
Conf.t <- t(Conf)
Conf0 <- Conf.t[rowSums(Conf.t) > 0,]
ConfP <- Conf0/rowSums(Conf0)
SpeciesMax <- cbind.data.frame("Species"=colnames(ConfP)[apply(ConfP,1,which.max)],"MaxValue"=apply(ConfP,1,max))
write.csv(SpeciesMax,output,quote=FALSE)
