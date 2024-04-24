###### 预测某基因的转录因子

# 加载所需的包
library(biomaRt)
library(TFBSTools)
library(JASPAR2018)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd('/Users/zijunyan/Desktop/ZF/')

# 指定Ensembl数据库和感兴趣的基因
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes = c("WEE1") #输入要预测的基因，可以一次性输入多个基因，基因越多，运行时间越久。
# 把基因名转换为ensembl_id
gene_id<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
               filters = "external_gene_name", # 指定转化基因的格式
               values = genes, # 基因列表
               mart = ensembl)
gene_id
gene_id = gene_id$ensembl_gene_id

# 获取基因组信息
gene_info <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), 
                   filters = "ensembl_gene_id", 
                   values = gene_id, 
                   mart = ensembl)

# 获取基因组数据
genome <- BSgenome.Hsapiens.UCSC.hg38

###定义提取启动子序列的函数
# 定义提取启动子序列的函数
get_promoter_sequence <- function(chromosome, start, end, strand, upstream_distance = 2000) {
  if (strand == 1) {
    promoter_start <- max(start - upstream_distance, 1)
    promoter_end <- start - 1
  } else {
    promoter_start <- end + 1
    promoter_end <- min(end + upstream_distance, seqlengths(genome)[[chromosome]])
  }
  promoter_seq <- getSeq(genome, paste0("chr", chromosome), promoter_start, promoter_end)
  return(promoter_seq)
}

####4.构建JASPAR数据库
# 使用JASPAR数据库
db <- file.path(system.file("extdata", package="JASPAR2018"), 
                "JASPAR2018.sqlite")
opts <- list()
opts[["species"]] <- 9606
# opts[["type"]] <- "ChIP-seq"
opts[["all_versions"]] <- FALSE
# pfm <- getMatrixSet(db, opts)
pfm <- getMatrixSet(JASPAR2018, opts)
head(pfm)

###5.预测转录因子并保存成csv文件
# 打开一个csv文件进行写入
file_conn <- file("TFBS.csv", "w")
# 写入CSV文件头部
writeLines("Gene,TF,Predicted_TFBS", file_conn)

# 循环处理每个PWM对象
for (name in names(pfm)) {
  # 提取基因的上游序列
  for (i in 1:dim(gene_info)[1]) {
    promoter_seq <- get_promoter_sequence(gene_info[i, "chromosome_name"], 
                                          gene_info[i, "start_position"], 
                                          gene_info[i, "end_position"], 
                                          gene_info[i, "strand"])
    
    # 预测转录因子结合位点
    predicted_TFBS <- matchPWM(pfm[[name]]@profileMatrix, promoter_seq)
    
    # 获取基因名
    gene_name <- getBM(attributes = "external_gene_name", 
                       filters = "ensembl_gene_id",
                       values = gene_info[i, "ensembl_gene_id"], 
                       mart = ensembl)
    
    # 如果预测结果不为空且基因名存在
    if (!is.null(predicted_TFBS) && nrow(gene_name) > 0) {
      gene_name <- gene_name$external_gene_name[1]
      
      # 将预测结果写入CSV文件
      writeLines(paste(gene_name, pfm[[name]]@name, predicted_TFBS, sep = ","), file_conn)
    }
  }
}

# 关闭文件连接
close(file_conn)
