ggheatmap <- function(exprs) {
  require(ggdendro)
  require(cowplot)
  require(data.table)
  require(stringr)
  
  miRNAs <- exprs[,1]
  exprs <- t(scale(t(exprs[,-1])))
  
  dendro.col <- as.dendrogram(hclust(dist(exprs)))
  col.order <- order.dendrogram(dendro.col)
  
  dendro.row <- as.dendrogram(hclust(dist(t(exprs))))
  row.order <- order.dendrogram(dendro.row)
  
  exprs <- data.table(miRNAs, t(scale(t(exprs)))[col.order, row.order])
  exprs <- melt(exprs, id = 1, variable.name = 'sample.name')
  
  dendro.data.x <- dendro_data(dendro.row)
  dendro.data.y <- dendro_data(dendro.col)
  
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  
  heatmap.exprs <- ggplot(exprs, aes(x = sample.name, y = miRNAs)) + geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "green", mid = "black", high = "red", midpoint = 0) + xlab('Sample') + theme(axis.line = element_blank(), axis.text.y = element_text(size = rel(0.8)))
  dendrogram.x <- ggplot(segment(dendro.data.x)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    theme_none + theme(axis.title.x=element_blank())
  dendrogram.y <- ggplot(segment(dendro.data.y)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none
  return(list(heatmap.exprs, dendrogram.x, dendrogram.y))
}