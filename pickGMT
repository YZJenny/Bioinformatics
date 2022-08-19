library(xlsx)
pick <- read.csv('/local/yanzijun/pickC2.csv',header = F)[,1]
geneset <- read.gmt('/mdshare/node9/yzj/publicData/GMT/c2.all.v6.2.symbols.gmt') 

for(i in 1:length(pick)){
  p <- pick[i]
  if(p %in% unique(geneset$term)){
    sub <- geneset[which(geneset$term == p),]
    line=paste0(p,paste0(sub$gene,collapse ='\t'),collapse ='\t')
    write.table(line,'/local/yanzijun/ZHY.txt',append = TRUE,row.names=F,col.names = F,quote = F)
  }
}
