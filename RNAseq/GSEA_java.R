height <- function(df){
  if(nrow(df) >= 1 & nrow(df) <= 20 ){
    height=7
  } else if(nrow(df) > 20 & nrow(df) <= 50 ) {
    height=10
  }
  return(height)
}

get_GSEA <- function(GMT,EXP,Cond1,Cond2,OUT_PATH){
  ## Construct EXP file
  col2 <- data.frame(Name=rownames(EXP),Description=rownames(EXP))
  EXP <- cbind(col2,EXP[,colnames(EXP)==Cond1],EXP[,colnames(EXP)==Cond2])
  colnames(EXP)[c(3,4)] <- c(Cond1,Cond2)
  write.table(EXP,paste(OUT_PATH,'EXP.txt',sep=''),col.names = T,row.names = F,quote = F,sep='\t')
  
  ## Construct CLS file
  L1 <- matrix(c(2,2,1),nrow = 1)
  L2 <- matrix(c('#',Cond1,Cond2),nrow = 1)
  L3 <-  matrix(c(rep(0,1),rep(1,1)),nrow = 1)
  CLS <- list(L1,L2,L3)
  for (i in 1:length(CLS)) {
    write.table(CLS[[i]], paste(OUT_PATH,Cond1,'_',Cond2,'.cls',sep=''),row.names = F,col.names = F,quote = F,sep=' ', append = T)
  }
  
  ###################
  ## Run GSEA
  ###################
  system(
    paste('/home/disk/fyh/tools/GSEA_Linux_4.0.2/gsea-cli GSEA -out ',OUT_PATH,' -rnd_seed 13579 -res ',OUT_PATH,'EXP.txt -cls ',OUT_PATH,Cond1,'_',Cond2,'.cls -gmx /home/yzj/publicData/GMT/',
          GMT,' -plot_top_x 2000 -permute gene_set -collapse false -set_min 2 -set_max 2000',sep='')
  )
  system(
    paste('mv ',OUT_PATH,'my_analysis* ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],sep=''))
  system(
    paste('mv ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond1,'_*.xls ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond1,'.xls',sep='')
  )
  system(
    paste('mv ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond2,'_*.xls ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond2,'.xls',sep='')
  )
  
  
  ###################
  ## BarPlot
  ###################
  if(file.exists(paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],sep=''))){
    A<-read.table(file=paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond1,'.xls',sep=''),sep="\t",head=T,row.names=1,fill=T)
    A<- A[A$FDR.q.val < 0.05,]
    B<-read.table(file=paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',Cond2,'.xls',sep=''),sep="\t",head=T,row.names=1,fill=T)
    B<- B[B$FDR.q.val < 0.05,]
    
    if(nrow(A) > 25){
      A=A[1:25,]
    }else{
      A=A
    }
    
    if(nrow(B) > 25){
      B=B[1:25,]
    }else{
      B=B
    }
    
    gsea_dot <- rbind(A,B)
    print(dim(gsea_dot))
    print(summary(gsea_dot$NES))
    
    gsea_dot=gsea_dot[order(gsea_dot$FDR.q.val,decreasing=F),]
    
    if(nrow(gsea_dot) > 0){
      ## mkdir
      if(!file.exists(paste(OUT_PATH,'Figure.',unlist(strsplit(GMT,'\\.'))[1],sep=''))){
        dir.create(file.path(OUT_PATH, paste('Figure.',unlist(strsplit(GMT,'\\.'))[1],sep='')),recursive = TRUE)
      }
      
      FIG_PATH=paste(OUT_PATH,'Figure.',unlist(strsplit(GMT,'\\.'))[1],'/',sep='')
      
      h=height(gsea_dot)
      p<- ggplot(gsea_dot, aes(x = NES, y = fct_reorder(GS.br..follow.link.to.MSigDB, NES),fill=NES))+
        geom_bar(stat='identity',colour='black',width = 0.78,position = position_dodge2(0.7))+
        theme_bw()+
        ylab("")+
        theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
        theme(axis.ticks = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_text(face = 'bold'),
              axis.title.x = element_text(face = 'bold'))+
        scale_fill_gradient2(low='blue',high='red',mid = 'white',midpoint = 0)
      
      pdf(paste(FIG_PATH,"DotPlot.pdf",sep=''),width = 20,height = h)
      print(p)
      dev.off()
    }else{
      print('no GSEA res!')
    } 
  }else{
    print('no Files!')
  }
}
