Entrez2GOResults <-
  function(gene.list,
           output.prefix,
           watermark = TRUE,
           watermark.content = 'Suck my dick!') {
    # Get gene classification and counts table with barplot
    require(data.table)
    require(cowplot)
    require(clusterProfiler)
    require(stringr)
    ontologies <- c('BP', 'MF', 'CC')
    list.ggo <-
      lapply(ontologies, function(x)
        groupGO(
          gene = gene.list,
          OrgDb = 'org.Hs.eg.db',
          ont  = x,
          level = 3,
          readable = TRUE
        ))
    
    fwrite(list.ggo[[1]]@result,
           file = paste('ggo.', output.prefix, '.BP.csv', sep = ''))
    fwrite(list.ggo[[2]]@result,
           file = paste('ggo.', output.prefix, '.MF.csv', sep = ''))
    fwrite(list.ggo[[3]]@result,
           file = paste('ggo.', output.prefix, '.CC.csv', sep = ''))
    
    list.ego <- lapply(ontologies, function(x)
      enrichGO(
        gene = gene.list,
        OrgDb = 'org.Hs.eg.db',
        ont  = x,
        readable = TRUE
      ))
    
    fwrite(list.ego[[1]]@result,
           file = paste('ego.', output.prefix, '.BP.csv', sep = ''))
    fwrite(list.ego[[2]]@result,
           file = paste('ego.', output.prefix, '.MF.csv', sep = ''))
    fwrite(list.ego[[3]]@result,
           file = paste('ego.', output.prefix, '.CC.csv', sep = ''))
    
    list.barplot.ggo <- lapply(list.ggo, function(x)
      barplot(
        x,
        drop = TRUE,
        showCategory = 12,
        order = TRUE
      ) + scale_x_discrete(
        label = function(x)
          str_wrap(x, width = 50)
      ))
    
    list.barplot.ego <- lapply(list.ego, function(x)
      dotplot(x, showCategory = 12) + scale_y_discrete(
        label = function(x)
          str_wrap(x, width = 50)
      ))
    
    summary.barplot.ggo <-
      ggdraw() + draw_plot(list.barplot.ggo[[1]], .5, .5, .49, .5) + draw_plot(list.barplot.ggo[[2]], 0, 0, .49, .5) + draw_plot(list.barplot.ggo[[3]], .5, 0, .49, .5) + draw_plot_label(c("BP", "CC", "MF"), c(0.5, 0, 0.5), c(1, 0.5, 0.5), size = 15) + draw_text('GO classification', .25, .75, 50)
    if (watermark) {
      summary.barplot.ggo <-
        summary.barplot.ggo + draw_label(
          watermark.content,
          angle = 45,
          size = 80,
          alpha = .2
        )
    }
    save_plot(
      paste('figures/barplot.ggo.', output.prefix, '.png', sep = ''),
      summary.barplot.ggo,
      base_height = 10,
      base_width = 16
    )
    
    summary.barplot.ego <-
      ggdraw() + draw_plot(list.barplot.ego[[1]], .5, .5, .49, .5) + draw_plot(list.barplot.ego[[2]], 0, 0, .49, .5) + draw_plot(list.barplot.ego[[3]], .5, 0, .49, .5) + draw_plot_label(c("BP", "CC", "MF"), c(0.5, 0, 0.5), c(1, 0.5, 0.5), size = 15) + draw_text('GO enrichment', .25, .75, 50)
    if (watermark) {
      summary.barplot.ego <-
        summary.barplot.ego + draw_label(
          watermark.content,
          angle = 45,
          size = 80,
          alpha = .2
        )
    }
    save_plot(
      paste('figures/barplot.ego.', output.prefix, '.png', sep = ''),
      summary.barplot.ego,
      base_height = 10,
      base_width = 16
    )
  }