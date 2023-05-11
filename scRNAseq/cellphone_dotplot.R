library(ggplot2)


dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename,
                    pvalue = 0.05,
                    mean = 0,
                    width = 30,
                    height = 25,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair

  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
 
  pval = unlist(sel_pval)
  pval[pval==0] = 0.00000009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  # pr[pr==0] = 1
  plot.data = cbind(plot.data,pr)
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean','log2mean')
  plot.data$pair <- factor(plot.data$pair,levels = sort(as.character(unique(plot.data$pair))))
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  ## Filter p value and mean expression
  plot.data_filter <- plot.data[plot.data$pvalue < 0.05 & plot.data$mean > 0,]
  write.table(plot.data_filter,'./LR_filter.txt',row.names = F,col.names = T,quote = F,sep='\t')
  
  ggplot(plot.data_filter,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=log2mean)) +
    scale_color_gradientn('Log2 mean (Gene 1, Gene 2)', colors=my_palette) +
    theme_bw() +
    theme(
          axis.text=element_text(size=14, colour = "black",family = 'bold'),
          axis.text.x = element_text(angle = 90, hjust = 1, family = 'bold'),
          axis.text.y = element_text(size=12, colour = "black", family = 'bold'),
          axis.title=element_blank(),
          text = element_text(family='bold'),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}
