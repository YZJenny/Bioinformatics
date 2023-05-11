Output.Cytoscape = function(meta_file, pvalues_file, count_filename,
                         meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  
  #######   Network
  print(meta_file)
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  
  dim(all_intr)
  
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = sort(unique(meta[,2]))
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    count1=length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  # here we disable this file generation
  all_count = all_count[-1,]
  all_count = as.data.frame(all_count)
  all_count[,"count"] = as.numeric(all_count[,"count"])
  #all_count = all_count %>% arrange(-count)
  
  write.table(all_count, 'count_network.txt', sep='\t', quote=F, row.names = F)
  return(all_count)
}

