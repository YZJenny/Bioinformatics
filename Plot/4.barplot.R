library(ggplot2)
size <- function(df){
  if(nrow(df) >= 1 & nrow(df) <= 10 ){
    s=11
  } else if(nrow(df) > 10  ) {
    s=7
  } 
  return(s)
} 
fdr=c(0.05,0.1)

david_wp='/Users/yzj/Desktop/0515/DATA/DAVID/'
term_wp <- '/Users/yzj/Desktop/0515/DATA/DAVID/term_list/'
out_wp <- '/Users/yzj/Desktop/0515/0.Data/3.BARPLOT/'
files=list.files(term_wp)
for(i in 1:length(files)){
  for(j in 1:length(fdr)){
    if(grepl(pattern = ".txt",x = files[i])){
      point_term=read.table(paste(term_wp,files[i],sep=""),header=F, sep='\t',stringsAsFactor=F)[,1]
      print(length(point_term))
      name=unlist(strsplit(files[i],'_'))
      if(length(name) == 5){
        david=name[1]
        fc_cutoff=name[3]
      }else if (length(name) == 6){
        david=paste(name[1],name[2],sep='_')
        fc_cutoff=name[4] 
      }
      david_res <- read.table(paste(david_wp,'fdr_',fdr[j],'/',david,'/up_',fc_cutoff,'.txt',sep=''),
                              header = T,sep='\t',stringsAsFactors = FALSE,quote = "")[,c(2,5)]
      df <- david_res[which(david_res$Term %in% point_term),]
      df$X.log10.pValue. <- -log(df$PValue,10)
      df$Term=factor(df$Term,levels = df$Term[order(df$X.log10.pValue.)])
      s=size(df)
      p <- ggplot(df,aes(Term,X.log10.pValue.))+
        #coord_cartesian(expand = F)+
        geom_bar(stat='identity',width = 0.9,colour='black',fill='blue')+
        coord_flip() + #horizontal
        labs(y='-log10(pValue)',x="")+
        theme(panel.grid=element_blank(),
              panel.background = element_blank(),
              axis.ticks.y= element_blank())+ #delete backgroud
        guides(fill=FALSE)+ #delete legend
        theme(axis.text.x = element_text(size = 8,face = 'bold'),
              axis.text.y = element_text(size = s,face = 'bold'), 
              axis.title.x = element_text(size = 8,face = 'bold'), 
              axis.title.y = element_text(size = 6,face = 'bold'))
       ggsave(filename=paste(out_wp,'fdr_',fdr[j],'/',files[i],'.pdf',sep=''),plot=p,width=8,height=6)
    }
  }
}
