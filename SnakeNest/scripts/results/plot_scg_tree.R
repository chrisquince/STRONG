#!/usr/bin/env Rscript
rm(list = ls()) 
library(tidyr)
library(stringr)
library(ggtree)
library(ggplot2)
library(wesanderson)



# custom reworking of gheatmap from ggtree, mostly add number on heatmap
gheatmap2 = function (p, data, offset = 0, width = 1)
{
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
    if (width==0){width = 1/ncol(data)}
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
    df <- p$data
    df <- df[df$isTip, ]
    start <- max(df$x, na.rm = TRUE) + offset
    dd <- as.data.frame(data)
    i <- order(df$y)
    i <- i[!is.na(df$y[i])]
    lab <- df$label[i]
    dd <- dd[lab, , drop = FALSE]
    dd$y <- sort(df$y)
    dd$lab <- lab
    dd <- gather(dd, variable, value, -c(lab, y))
    i <- which(dd$value == "")
    if (length(i) > 0) {
        dd$value[i] <- NA
    }
    dd$variable <- factor(dd$variable, levels = colnames(data))
    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from = dd$variable, to = V2)
    mapping <- unique(mapping)
    dd$x <- V2
    dd$width <- width

    # draw the heatmap 
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), 
            width = width, color = "white",inherit.aes = FALSE)
    # add annotation on tiles 
    p2 <- p2 + geom_text(data = dd, aes(x,y,label=value, fontface="bold"),size=7/sqrt(ncol(data)),inherit.aes = FALSE)
    # add color 
    p2 <- p2 + scale_fill_gradientn(colours = pal)
    p2 <- p2 + theme(legend.position = "right")
    attr(p2, "mapping") <- mapping
    return(p2)
}


# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
newick <- args[1]
dist_matrix <- args[2]
savePlot <- args[3]

# -------- build a tree ------------
tree_original<-read.tree(newick)
p <- ggtree(tree_original) + geom_tiplab(align=TRUE) 
p = p + theme_tree2(axis.text.x = element_text(angle = 45, hjust = 1))

# add color if mag
tips = p$data[p$data$isTip,]$label
has_haplo = grepl("haplo",tips)
has_eval = grepl("eval",tips)
has_gtdb = grepl("gtdb",tips)
is_mag = list("Haplotype" = tips[which(has_haplo)], "Mag" = tips[which(!(has_haplo|has_eval|has_gtdb))])
if (sum(has_eval)!=0){is_mag[["Evaluation"]]=tips[has_eval]}
if (sum(has_gtdb)!=0){is_mag[["GTDB Reference"]]=tips[has_gtdb]}

p = groupOTU(p, is_mag, 'SCG') + aes(color=SCG)

# add heatmap using dist mat 
offset=1.5*max(p$data$x)
if (offset==0){offset=1}
data = read.table(dist_matrix,header=TRUE,row.name=1)
# reorder data so that the distmat is in the same order :

y_order = order(p$data$y[p$data$isTip])
label_order = p$data$label[p$data$isTip][y_order]
data = t(t(data)[,label_order])[,rev(label_order)]


r = gheatmap2(p, data,offset=offset,width=2)+ scale_x_ggtree()+ylim(0.5,max(p$data$y)+0.5) + labs(fill="Percent substitution")

# title please
mag = str_split(newick, "/")[[1]][2]
r = r + ggtitle(mag)

# -------- just draw the fig ------------
pdf(savePlot)
print(r)
dev.off()


