#!/usr/bin/env Rscript
rm(list = ls()) 
library(ggtree)
library(ggplot2)
library(tidyr)
library(wesanderson)



# custom reworking of gheatmap from ggtree, mostly add number on heatmap
gheatmap2 = function (p, data, offset = 0, width = 1)
{
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
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
    p2 <- p2 + geom_text(data = dd, aes(x,y,label=value, fontface="bold"),inherit.aes = FALSE)
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
is_mag = list("Haplotype" = tips[which(has_haplo)], "Mag" = tips[which(!has_haplo)])
p = groupOTU(p, is_mag, 'SCG') + aes(color=SCG)

# add heatmap using dist mat 
data = read.table(dist_matrix,header=TRUE,row.name=1)
r = gheatmap2(p, data, offset=1.5*max(p$data$x),width=2)+ scale_x_ggtree()+ylim(0.5,max(p$data$y)+0.5) + labs(fill="Percent substitution")

# -------- just draw the fig ------------
pdf(savePlot)
print(r)
dev.off()


