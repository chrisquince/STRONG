#!/usr/bin/env Rscript

Conf <- read.csv("Conf.csv",header=TRUE,row.names=1)
Conf.t <- t(Conf)
Conf0 <- Conf.t[rowSums(Conf.t) > 0,]
ConfP <- Conf0/rowSums(Conf0)
SpeciesMax <- cbind.data.frame("Species"=colnames(ConfP)[apply(ConfP,1,which.max)],"MaxValue"=apply(ConfP,1,max))
write.csv(SpeciesMax,"SpeciesMax.csv",quote=FALSE)
