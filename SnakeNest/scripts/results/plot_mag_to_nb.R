 #!/usr/bin/env Rscript
library(ggplot2)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
coverage = args[1]
haplo_nb = args[2]
output = args[3]

# get coverage
cov_tot = rowSums(read.table(coverage,header=TRUE,row.name=1))
# get haplotype nb 
hap_nb = read.table(haplo_nb,row.name=1)

# plot it 
data = data.frame("cov"=cov_tot[rownames(hap_nb)],"nb"=hap_nb[,])
png(output)
ggplot(data=data, aes(x=cov, y=nb))+geom_point()+ geom_smooth()+xlab("Total mag coverage over all samples")+ylab("Number of haplotypes resolved for each mag") + theme_bw() + ggtitle("Haplotype number versus total mag coverage")
dev.off()





