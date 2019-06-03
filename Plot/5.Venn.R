library(VennDiagram)
chip_wp='Desktop/0515/DATA/Mycchip/deg/'
deg_wp='Desktop/0515/DATA/DEG/'
output_wp='Desktop/0515/DATA/Mycchip/overlapDEG/'
outres_wp='Desktop/0515/7.Mycchip/1.VENN/'

fdr=c(0.05,0.1)
fc_cutoff=c(1,1.3)
cat=c('all','up','dn')
chip_files <- list.files(chip_wp)

for(i in 1:length(chip_files)){
  chip <- read.table(paste(chip_wp,chip_files[i],sep=''),header=F, sep='\t',stringsAsFactor=F)[,1]
  for(j in 1:length(fdr)){
    for(m in 1:length(fc_cutoff)){
      for(n in 1:length(cat)){
        deg <- read.table(paste(deg_wp,'fdr_',fdr[j],'/',cat[n],'_',fc_cutoff[m],'.txt',sep=''),header=F, sep='\t',stringsAsFactor=F)[,1]
        out=paste(unlist(strsplit(chip_files[i],'.txt'))[1],'_',cat[n],'_',fc_cutoff[m],sep='')
        venn.diagram(list(Mycchip=chip,RNAseq=deg),cat.pos=c(-20,  20),fill=c("blue2","yellow2"),paste(outres_wp,'fdr_',fdr[j],'/',out,".tiff",sep=''))
        ovlap=intersect(chip,deg)
        write.table(ovlap,paste(output_wp,'fdr_',fdr[j],'/',out,".txt",sep=''),sep="\t", quote=F, row.names=F, col.names=F)
        }
    }
  }
}
