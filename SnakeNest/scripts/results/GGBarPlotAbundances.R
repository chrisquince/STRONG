#!/usr/bin/env Rscript
# R_barplot
# Takes a csv file of haplotype coverages and another of the variation in this and produces a normalised barplot with error bars
#
# file GGBarPlotcoverages
#
# inputs
# 	coverages 		- the haplotypes coverages
#   varAbunances	- the variation in these abunances
# 
# outputs
# 	coverageBarPlot - the barplot file produced

# Version    Author       Date      Affiliation
# 1.00       J K Summers  03/06/20  Quince Lab - Warwick Medical School - University of Warwick
#									Code adapted from a sample supplied by Chris Quince

library(reshape2)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

if (length(args) !=6){
  stop("Five arguments must be supplied (intensity file, intensity variance, normalisation, read length, assembly kmer length and file name for plot)", call = FALSE)
}

coverages <- args[1]
varcoverages <- args[2]
norm <- args[3]
Read_length <- as.numeric(args[4])
kmer_size <- as.numeric(args[5])
savePlot <- args[6]

# read data
norm <- read.table(norm, header = TRUE, row.names = 1)
gammacoverages <- read.csv(coverages, header = FALSE, row.names = 1)
varGammacoverages <- read.csv(varcoverages, header = TRUE, row.names = 1)

# Super hidden information : normalisation column is always in the same order as bayespath : bayespath order is derived from subgraphs, derived from assembly, derived from SAMPLE
# Normalisation order derive from collate.py script, which keep the same order as SAMPLE
# bayespath start at 0 R index start at 1, so we need to increment
colnames(gammacoverages)= colnames(norm)[unlist(gammacoverages[1,]+1)]
gammacoverages = gammacoverages[-1,,drop=FALSE]

# sometimes bayespath does not output coverage for all samples, normalisation need to be done only on these samples
norm = norm[,match(colnames(gammacoverages),colnames(norm),nomatch = 0),drop=FALSE]

# transform intensity in coverage:
# nb_reads_per_genome = (gamma * genome_length)
# Ck =  (gamma * genome_length) * read_length / genome_length

# transform kmer coverage to coverage
#Ck*Read_length/(Read_length-kmer_size+1) = C
kmer_intensity_to_nuc_coverage = (Read_length*Read_length/(Read_length-kmer_size+1))

# transforme intensity in normalised kmer coverage
gammaP <- kmer_intensity_to_nuc_coverage*t(gammacoverages)/t(norm)[,1]
# be carefull to apply scalling to variance not std
gammaVarP <- sqrt(kmer_intensity_to_nuc_coverage*t(varGammacoverages)/t(norm)[,1])

#gammaPS <- cbind.data.frame(gammaP)
gammaPSMelt <- melt(gammaP, id.var = 1)
colnames(gammaPSMelt) <- c('Samples', 'Strain', 'coverage')

#gammaVarPS <- cbind.data.frame(gammaVarP)
gammaVarPSMelt <- melt(gammaVarP, id.var = 1)
colnames(gammaVarPSMelt) <- c('Samples','Strain','Sdcoverage')

gammaVarPSMelt$Samples <- NULL
gammaVarPSMelt$Strain <- NULL

gammaB <- cbind.data.frame(gammaPSMelt, gammaVarPSMelt)
gammaB$Strain <- as.factor(gammaB$Strain)


coverageBarPlot <- ggplot(data = gammaB, aes(x = Samples, y = coverage, colour = Strain, fill = Strain, group = Strain))
coverageBarPlot <- coverageBarPlot + geom_bar(stat = "identity", position = position_dodge())
coverageBarPlot <- coverageBarPlot + geom_errorbar(aes(ymin = coverage - Sdcoverage, ymax = coverage + Sdcoverage), width = .2, position = position_dodge(.9), colour = "Black")
coverageBarPlot <- coverageBarPlot + theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(size = 18), legend.text = element_text(size = 16), plot.title = element_text(hjust = 0.5,size = 20))
coverageBarPlot = coverageBarPlot +  ggtitle(print(gsub("F_Intensity.csv", "", basename(coverages))))
coverageBarPlot = coverageBarPlot + ylab("coverage per genome")
ggsave(savePlot)

