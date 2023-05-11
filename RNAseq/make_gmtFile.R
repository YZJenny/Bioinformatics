library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(GO.db)

############
## 1. 提取KEGG
############
hsa_term <- clusterProfiler::download_KEGG("hsa")
names(hsa_term)
PATH2ID <- hsa_term$KEGGPATHID2EXTID
PATH2NAME <- hsa_term$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("ID", "ENTREZID", "DESCRPTION")

id2symbol=mapIds(org.Hs.eg.db,PATH_ID_NAME$ENTREZID,'SYMBOL',"ENTREZID")
symbol <- data.frame(ENTREZID=names(id2symbol),gene=id2symbol)

print(all(PATH_ID_NAME$ENTREZID==symbol$ENTREZID))
term2symbol <- cbind(PATH_ID_NAME,symbol)[,c(3,5)]

### 写出gmt
term2symbol_list <- tapply(term2symbol[,2],as.factor(term2symbol[,1]),function(x) x)
write.gmt <- function(geneSet,gmt_file){
  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    cat('\tNA\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')
  }
  sink()
}

write.gmt(geneSet=term2symbol_list,gmt_file='publicData/GMT/KEGG_clusterprofile.gmt')

############
## 2. 提取GO
############
keytypes(GO.db)
columns(GO.db)

k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ID2Gene <- select(org.Hs.eg.db,
       keys = k,
       columns = c("GO"),
       keytype="SYMBOL")

head(ID2Gene)

k <- keys(GO.db, keytype = "GOID")
ID2Des <- select(GO.db,
         keys = k,
         columns = c("TERM","ONTOLOGY"),
         keytype="GOID")
head(ID2Des)

ID2Gene_BP <- ID2Gene[ID2Gene$ONTOLOGY=='BP',]
ID2Des_BP <- ID2Des[ID2Des$ONTOLOGY=='BP',]
colnames(ID2Des_BP)[1] <- 'GO'

term2symbol <- merge(ID2Des_BP[,1:2],ID2Gene_BP[,1:2],by='GO')[,2:3]
dim(term2symbol)
head(term2symbol)



### 写出gmt
term2symbol_list <- tapply(term2symbol[,2],as.factor(term2symbol[,1]),function(x) x)
write.gmt <- function(geneSet,gmt_file){
  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    cat('\tNA\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')
  }
  sink()
}

write.gmt(geneSet=term2symbol_list,gmt_file='publicData/GMT/GO_clusterprofile.gmt')
