#!/usr/bin/env Rscript
library(ggplot2)
library(ggtree)
library(stringr)
library(stringi)
library(treeio)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
bac_tree = args[1]
ar_tree = args[2]
output = args[3]

# get coverage
cov_tot = rowSums(read.table(coverage,header=TRUE,row.name=1))
# get haplotype nb 
hap_nb = read.table(haplo_nb,row.name=1)

# plot it 
data = data.frame("cov"=cov_tot[rownames(hap_nb)],"nb"=hap_nb[,])
png(output)
ggplot(data=data, aes(x=cov, y=nb))+geom_point()+ stat_smooth(span=2)+xlab("Total mag coverage over all samples")+ylab("Number of haplotypes resolved for each mag") + theme_bw() + ggtitle("Haplotype number versus total mag coverage")
dev.off()



tree_original<-read.tree("gtdbtk.bac120.classify.tree")
p <- ggtree(tree_original) + geom_tiplab(align=TRUE) 


get_mag_metadata = function(file)
{
	df = read.table(file,header=TRUE,row.name=1,sep="\t")
	### extract tip taxonomy ###
	# remove stupid prefix 
	taxonomy = str_remove_all(df[,'classification'], ".__")
	# split into ranks 
	taxonomy = stri_list2matrix(strsplit(taxonomy, ";"), byrow=TRUE)

	# meaningfull colnames/rownames
	colnames(taxonomy) = c('Domain','Phylum','Class','Order','familly','genus','strain')
	rownames(taxonomy) = rownames(df)

	### get RED ###
	# replace N/A by 1
	red = as.numeric(str_replace(df$red_value,"N/A","1"))
	return(data.frame(cbind(taxonomy,RED = red)))
}
get_node_metadata = function(file,tree)
{
	df = read.table(file,row.name=1,sep="\t")
	### extract tip taxonomy ###
	# remove stupid prefix 
	taxonomy = str_remove_all(df[,1], ".__")
	# split into ranks 
	taxonomy = t(data.frame(strsplit(taxonomy, ";")))
	# meaningfull colnames/rownames
	colnames(taxonomy) = c('Domain','Phylum','Class','Order','familly','genus','strain')
	rownames(taxonomy) = rownames(df)
	# only return tips in the tree
	representatives = tree$tip.label
	taxonomy = taxonomy[which(rownames(taxonomy) %in% representatives),]
	return(data.frame(taxonomy))
}



df_res = get_mag_metadata("gtdbtk.bac120.summary.tsv")
gtdb_taxa = get_node_metadata("/mnt/gpfs/seb/Applications/miniconda3/envs/gtdbtk/share/gtdbtk-0.3.2/db/taxonomy/gtdb_taxonomy.tsv",tree_original)


phylums = unique(df_res$Phylum)
Phylum_trees <- vector(mode="list", length=length(phylums))

pdf("test.pdf")
for (i in 1:length(phylums))
{
	nodes_mags = rownames(df_res)[df_res$Phylum==phylums[i]]
	nodes_gtdb = rownames(gtdb_taxa)[gtdb_taxa$Phylum==phylums[i]]
	clade <- MRCA(tree_original,c(nodes_mags,nodes_gtdb))
	# define mapping of species to nodes
	classes = unique(gtdb_taxa[gtdb_taxa$Phylum==phylums[i],]$Class)
	classe_group = vector(mode="list", length=length(classes))
	names(classe_group)=classes
	for (name in classes)
	{
		classe_group[[name]]=rownames(gtdb_taxa)[which(gtdb_taxa$Class==name)]
	}
	# add colouring depending on mapping
	p_tree = tree_subset(tree_original,clade , levels_back = 0)
	p = ggtree(p_tree)
	p2 = groupOTU(p, classe_group, 'Class') + aes(color=Class) +  theme(legend.position="right")

	# cosmetic fix
	p2$data$Class = as.character(p2$data$Class)

	p2$data[p2$data$Class==0,]=phylums[i]



	p <- ggtree(tree) + geom_tiplab()
	p_clade = viewClade(p, clade)
}
dev.off()
}




p$data[match(rownames(df_res),p$data$label,nomatch=0),]

tree_original <- groupOTU(tree_original, rownames(df_res))


p = ggtree(tree_original, aes(colour = group)) 
p2 = ggtree(test, aes(colour = group)) 


library(treeio)
library(tidytree)
library(ape)
library(ggtree)





data(iris)
rn <- paste0(iris[,5], "_", 1:150)
rownames(iris) <- rn
d_iris <- dist(iris[,-5], method="man")

tree_iris <- ape::bionj(d_iris)
grp <- list(setosa     = rn[1:50],
            versicolor = rn[51:100],
            virginica  = rn[101:150])










