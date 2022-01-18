{
  ## cellphoneDB
  
  source('/home/yjingjing/project/hfz2020adjust/pla_adjust/6th_TPA/cellphoneDb/scr/SEURATEXT/cellphone/cellphone_util.R')
  source('/home/yjingjing/project/hfz2020adjust/pla_adjust/6th_TPA/cellphoneDb/scr/SEURATEXT/cellphone/cellphone_dotplot.R')
  source('/home/yjingjing/project/hfz2020adjust/pla_adjust/6th_TPA/cellphoneDb/scr/SEURATEXT/cellphone/cellphone_heat.R')
  library("stringr")
  
  # cellphoneDb_TBdeci_NM.R
  tag = c('NM','PE')
  # base_dir = '/home/disk/yangjing/10th_hfz/pla/cellphonedb'
  base_dir = '/home/yjingjing/cellPhonedb'
  # setwd("/home/disk/yangjing/10th_hfz/pla/cellphonedb")
  setwd('~/hfz/')
  deci_pbmc = readRDS('/home/yjingjing/project/hfz2020adjust/deci_adjust/9th_discard_deci1/deci_pbmc_disTB.20200903.RDS')
  
  for(i in tag){
    if(i=='NM'){
      # outDir <- paste(base_dir,'/',i,'_deciTB',sep='')
      # if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
      
      outDir <- '~/hfz/'
      deci <- subset(deci_pbmc,subset = tag==c("NM"))
      # pla_TB <- subset(pla_pbmc_disSC_left,,subset = tag==c("NM")& annoClust==c("VCT",'EVT','SCT') )
      # TB_deci <- merge(deci,y=c(pla_TB)) # decidua
    }else{
      # outDir <- file.path(base_dir,'/',i,'_deciTB',sep='')
      # if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
      outDir <- '~/hfz/'
      deci <- subset(deci_pbmc,subset = tag==c("PE"))
      # pla_TB <- subset(pla_pbmc_disSC_left,,subset = tag==c("PE") & annoClust==c("VCT",'EVT','SCT') )
      # TB_deci <- merge(deci,y=c(pla_TB)) # decidua
      
    }
    
    deci@meta.data$Cell <- row.names(deci@meta.data)
    head(deci@meta.data)
    deci_meta <- deci@meta.data[,c('Cell','new')]
    colnames(deci_meta) = c('Cell',"cell_type")
    head(deci_meta,row.names=FALSE)
    write.table(deci_meta,file=paste(outDir,'/',i,"_meta.txt",sep=''),row.names = FALSE,quote=FALSE,sep='\t')
    
    ### get count mtx 
    deci_count <- as.matrix(deci@assays$RNA@data)
    deci_count = as.data.frame(deci_count)
    deci_count = rownames_to_column(deci_count,var = "Gene")
    # count_mtx$Gene = toupper(count_mtx$Gene)
    write.table(deci_count,file=paste(outDir,'/',i,"_count.txt",sep=''),row.names = FALSE,quote=FALSE,sep='\t')
    
    # open terminal
    # setwd(outDir)
    # system(            
    #     source /home/yjingjing/cpdb-venv/bin/activate
    # )
    # cellphonedb method statistical_analysis test_meta.txt test_counts.txt
    
    # # cellphonedb method statistical_analysis /home/disk/yangjing/10th_hfz/pla/cellphonedb/plaNM_meta.txt /home/disk/yangjing/10th_hfz/pla/cellphonedb/plaNM_count.txt --counts-data gene_name
    # system(
    #     paste("cellphonedb method statistical_analysis  ", outDir,'/',i,"_count.txt  ", outDir,'/',i,"_count.txt --counts-data =gene_name",sep='')
    #  )
    # # cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt
    # cellphonedb plot heatmap_plot yourmeta.txt
    # getwd()
  }
}

{
  ## Cell Chat
  # devtools::install_github("sqjin/cellchat")
  # devtools::install_github("jokergoo/ComplexHeatmap")
  library(Seurat)
  library(NMF)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(svglite)
  options(stringsAsFactors = FALSE)
  
  ## I: Data input & processing and initialization of CellChat object
  # Load data
  # Make sure you are using the normalized data without negative values as input. 
  # data_source = 'pla'
  # sample_type = 'NM'
  # data_source_list = c('deci', 'pla')
  data_source_list = c('deci')
  sample_type_list = c('NM', 'PE')
  # sample_type_list = c('deci4-PE', 'deci6-PE', 'deci8-PE', 'deci9-PE')
  for (data_source in data_source_list){
    if (data_source == 'deci'){
      rds = '/home/disk/liuqiming/hfz_anno_data/deci_remove_unkown.rds'
    }else if(data_source == 'pla'){
      rds = '/home/disk/liuqiming/hfz_anno_data/pla.rds'
    }
    data = readRDS(rds)
    data@meta.data$celltype = as.character(data@meta.data$celltype)
    # summary(data)
    # str(data)
    for (sample_type in sample_type_list){
      print(sprintf('=============================%s-%s=============================', data_source, sample_type))
      print('preprocessing...')
      input_cells = rownames(data@meta.data)[which(data@meta.data$status == sample_type)]
      # input_cells = rownames(data@meta.data)[which(data@meta.data$orig.ident == sample_type)]
      data.input = data@assays$RNA@data[, input_cells]
      identity = data.frame(group = data@meta.data[input_cells, 'celltype'], row.names = input_cells)
      unique(identity$group)
      
      # Create a CellChat object
      cellchat = createCellChat(data = data.input)
      cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
      cellchat <- setIdent(cellchat, ident.use = "labels")
      levels(cellchat@idents)
      groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
      
      # Set the ligand-receptor interaction database
      CellChatDB <- CellChatDB.human
      # showDatabaseCategory(CellChatDB)
      dplyr::glimpse(CellChatDB$interaction)
      
      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
      cellchat@DB <- CellChatDB.use # set the used database in the object
      
      # Preprocessing the expression data for cell-cell communication analysis
      cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
      future::plan("multiprocess", workers = 4) # do parallel
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- projectData(cellchat, PPI.human)
      
      
      ## II: Inference of cell-cell communication network
      # Compute the communication probability and infer cellular communication network
      cellchat <- computeCommunProb(cellchat)
      # Infer the cell-cell communication at a signaling pathway level
      cellchat <- computeCommunProbPathway(cellchat)
      # Calculate the aggregated cell-cell communication network
      cellchat <- aggregateNet(cellchat)
      
      print('saving calculating results...')
      # save Legand_recepter result
      {
        save_path = paste0('/home/wzk/hfz/cellchat/', data_source, sample_type)
        if (!dir.exists(save_path)) dir.create(save_path)
        setwd(save_path)
        LR_list = names(cellchat@net$prob[1,1,])
        cell_types = sort(unique(identity[,1]))
        LR_prob_df = as.data.frame(matrix(nrow = length(LR_list), ncol = length(cell_types)^2), row.names = LR_list)
        LR_pval_df = as.data.frame(matrix(nrow = length(LR_list), ncol = length(cell_types)^2), row.names = LR_list)
        colname_list = c()
        for (i in cell_types){
          for (j in cell_types){
            colname_list = c(colname_list, paste0(i, '@', j))
          }
        }
        colnames(LR_prob_df) = colname_list
        colnames(LR_pval_df) = colname_list
        for (LR in LR_list){
          prob_df = cellchat@net$prob[, , LR]
          pval_df = cellchat@net$pval[, , LR]
          prob_vec = c()
          pval_vec = c()
          for (tag in cell_types){
            prob_vec = c(prob_vec, as.numeric(prob_df[tag,]))
            pval_vec = c(pval_vec, as.numeric(pval_df[tag,]))
          }
          LR_prob_df[LR, ] = prob_vec
          LR_pval_df[LR, ] = pval_vec
        }
        write.csv(LR_prob_df, file = paste0(data_source, sample_type, '_LR_prob.csv'))
        write.csv(LR_pval_df, file = paste0(data_source, sample_type, '_LR_pval.csv'))
      }
      
      # save pathway result
      {
        save_path = paste0('/home/wzk/hfz/cellchat/', data_source, sample_type)
        if (!dir.exists(save_path)) dir.create(save_path)
        setwd(save_path)
        pathway_list = names(cellchat@netP$prob[1,1,])
        cell_types = sort(unique(identity[,1]))
        pathway_prob_df = as.data.frame(matrix(nrow = length(pathway_list), ncol = length(cell_types)^2), row.names = pathway_list)
        # pathway_pval_df = as.data.frame(matrix(nrow = length(pathway_list), ncol = length(cell_types)^2), row.names = pathway_list)
        colname_list = c()
        for (i in cell_types){
          for (j in cell_types){
            colname_list = c(colname_list, paste0(i, '@', j))
          }
        }
        colnames(pathway_prob_df) = colname_list
        # colnames(pathway_pval_df) = colname_list
        for (pathway in pathway_list){
          prob_df = cellchat@netP$prob[, , pathway]
          # pval_df = cellchat@netP$pval[, , pathway]
          prob_vec = c()
          # pval_vec = c()
          for (tag in cell_types){
            prob_vec = c(prob_vec, as.numeric(prob_df[tag,]))
            # pval_vec = c(pval_vec, as.numeric(pval_df[tag,]))
          }
          pathway_prob_df[pathway, ] = prob_vec
          # pathway_pval_df[pathway, ] = pval_vec
        }
        write.csv(pathway_prob_df, file = paste0(data_source, sample_type, '_pathway_prob.csv'))
      }
      
      print('ploting pathway figures...')
      ## III: Visualization and systems analysis of cell-cell communication network
      dir_path = paste0('/home/wzk/hfz/cellchat/figures/', data_source, sample_type)
      if (!dir.exists(dir_path)) dir.create(dir_path)
      setwd(dir_path)
      cellchat@netP$pathways
      vertex.receiver = seq(1,(length(groupSize)+1)/2) # a numeric vector
      cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
      for (pathways.show in cellchat@netP$pathways){
        # pathways.show <- "TGFb"
        svglite(file =paste0(pathways.show, "_line.svg"))
        netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
        dev.off()
        svglite(file =paste0(pathways.show, "_circle.svg"))
        netVisual_aggregate(cellchat, signaling = c(pathways.show), layout = "circle", vertex.size = groupSize)
        dev.off()
        svglite(file =paste0(pathways.show, "_contribute.svg"))
        print(netAnalysis_contribution(cellchat, signaling = pathways.show))
        dev.off()
        svglite(file =paste0(pathways.show, "_signalingRole.svg"))
        # the slot 'netP' means the inferred intercellular communication network of signaling pathways
        netVisual_signalingRole(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10)
        dev.off()
      }
      
      ## Identify global communication patterns and major signals for specific cell groups
      # Identify outgoing communication pattern of secreting cells
      print('ploting pattern figures...')
      nPatterns = 5
      svglite(file =paste0('0_', data_source, sample_type, "_outgoing_patern_heatmap.svg"))
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, font.size = 8, width = 5, height = 12)
      dev.off()
      svglite(file =paste0('1_', data_source, sample_type, "_outgoing_patern_rivermap.svg"))
      print(netAnalysis_river(cellchat, pattern = "outgoing"))
      dev.off()
      svglite(file =paste0('2_', data_source, sample_type, "_outgoing_patern_dotmap.svg"))
      print(netAnalysis_dot(cellchat, pattern = "outgoing"))
      dev.off()
      
      svglite(file =paste0('0_', data_source, sample_type, "_incoming_patern_heatmap.svg"))
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, font.size = 8, width = 5, height = 12)
      dev.off()
      svglite(file =paste0('1_', data_source, sample_type, "_incoming_patern_rivermap.svg"))
      print(netAnalysis_river(cellchat, pattern = "incoming"))
      dev.off()
      svglite(file =paste0('2_', data_source, sample_type, "_incoming_patern_dotmap.svg"))
      print(netAnalysis_dot(cellchat, pattern = "incoming"))
      dev.off()
    }
  }
  #####
}

#####
{
  library(Seurat)
  library(NMF)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(svglite)
  options(stringsAsFactors = FALSE)
  
  data_source = 'deci-pla'
  deci_data = readRDS('/home/yjingjing/project/hfz2020adjust/deci_adjust/9th_discard_deci1/deci_pbmc_disTB.20200903.RDS')
  pla_data = readRDS('/home/yjingjing/project/hfz2020adjust/pla_adjust/9th_discard_pla4/pla_pbmc_disSC_left.20200903.RDS')
  # gene_list = intersect(rownames(deci_data@assays$RNA@data), rownames(pla_data@assays$RNA@data))
  for (pos in 1:length(deci_data@meta.data$annoClust)){
    deci_data@meta.data$new_tag[pos] = paste0('d', deci_data@meta.data$annoClust[pos])
  }
  for (pos in 1:length(pla_data@meta.data$annoClust)){
    pla_data@meta.data$new_tag[pos] = paste0('p', pla_data@meta.data$annoClust[pos])
  }
  meta_clonames = colnames(deci_data@meta.data)
  meta_clonames[7] = "RNA_snn_res"
  colnames(deci_data@meta.data) = meta_clonames
  colnames(pla_data@meta.data) = meta_clonames
  
  meta_data_df = rbind(deci_data@meta.data, pla_data@meta.data)
  sample_type_list = c('PE')
  for (sample_type in sample_type_list){
    print(sprintf('=============================%s-%s=============================', data_source, sample_type))
    print('preprocessing...')
    input_cells_1 = rownames(deci_data@meta.data)[which(deci_data@meta.data$tag == sample_type)]
    input_cells_2 = rownames(pla_data@meta.data)[which(pla_data@meta.data$tag == sample_type)]
    data.input = cbind(deci_data@assays$RNA@data[, input_cells_1], pla_data@assays$RNA@data[, input_cells_2])
    identity = data.frame(group = meta_data_df[c(input_cells_1, input_cells_2), 'new_tag'], row.names = c(input_cells_1, input_cells_2))
    unique(identity$group)
    
    # Create a CellChat object
    cellchat = createCellChat(data = data.input)
    cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
    cellchat <- setIdent(cellchat, ident.use = "labels")
    levels(cellchat@idents)
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    # Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human
    showDatabaseCategory(CellChatDB)
    dplyr::glimpse(CellChatDB$interaction)
    
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
    cellchat@DB <- CellChatDB.use # set the used database in the object
    
    # Preprocessing the expression data for cell-cell communication analysis
    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multiprocess", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    
    
    ## II: Inference of cell-cell communication network
    # Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat)
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    # Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    print('saving calculating results...')
    # save Legand_recepter result
    {
      save_path = paste0('/home/wzk/hfz/cellchat/', data_source, sample_type)
      if (!dir.exists(save_path)) dir.create(save_path)
      setwd(save_path)
      LR_list = names(cellchat@net$prob[1,1,])
      cell_types = sort(unique(identity[,1]))
      LR_prob_df = as.data.frame(matrix(nrow = length(LR_list), ncol = length(cell_types)^2), row.names = LR_list)
      LR_pval_df = as.data.frame(matrix(nrow = length(LR_list), ncol = length(cell_types)^2), row.names = LR_list)
      colname_list = c()
      for (i in cell_types){
        for (j in cell_types){
          colname_list = c(colname_list, paste0(i, '@', j))
        }
      }
      colnames(LR_prob_df) = colname_list
      colnames(LR_pval_df) = colname_list
      for (LR in LR_list){
        prob_df = cellchat@net$prob[, , LR]
        pval_df = cellchat@net$pval[, , LR]
        prob_vec = c()
        pval_vec = c()
        for (tag in cell_types){
          prob_vec = c(prob_vec, as.numeric(prob_df[tag,]))
          pval_vec = c(pval_vec, as.numeric(pval_df[tag,]))
        }
        LR_prob_df[LR, ] = prob_vec
        LR_pval_df[LR, ] = pval_vec
      }
      write.csv(LR_prob_df, file = paste0(data_source, sample_type, '_LR_prob.csv'))
      write.csv(LR_pval_df, file = paste0(data_source, sample_type, '_LR_pval.csv'))
    }
    
    # save pathway result
    {
      save_path = paste0('/home/wzk/hfz/cellchat/', data_source, sample_type)
      if (!dir.exists(save_path)) dir.create(save_path)
      setwd(save_path)
      pathway_list = names(cellchat@netP$prob[1,1,])
      cell_types = sort(unique(identity[,1]))
      pathway_prob_df = as.data.frame(matrix(nrow = length(pathway_list), ncol = length(cell_types)^2), row.names = pathway_list)
      # pathway_pval_df = as.data.frame(matrix(nrow = length(pathway_list), ncol = length(cell_types)^2), row.names = pathway_list)
      colname_list = c()
      for (i in cell_types){
        for (j in cell_types){
          colname_list = c(colname_list, paste0(i, '@', j))
        }
      }
      colnames(pathway_prob_df) = colname_list
      # colnames(pathway_pval_df) = colname_list
      for (pathway in pathway_list){
        prob_df = cellchat@netP$prob[, , pathway]
        # pval_df = cellchat@netP$pval[, , pathway]
        prob_vec = c()
        # pval_vec = c()
        for (tag in cell_types){
          prob_vec = c(prob_vec, as.numeric(prob_df[tag,]))
          # pval_vec = c(pval_vec, as.numeric(pval_df[tag,]))
        }
        pathway_prob_df[pathway, ] = prob_vec
        # pathway_pval_df[pathway, ] = pval_vec
      }
      write.csv(pathway_prob_df, file = paste0(data_source, sample_type, '_pathway_prob.csv'))
    }
    
    print('ploting pathway figures...')
    ## III: Visualization and systems analysis of cell-cell communication network
    dir_path = paste0('/home/wzk/hfz/cellchat/figures/', data_source, sample_type)
    if (!dir.exists(dir_path)) dir.create(dir_path)
    setwd(dir_path)
    cellchat@netP$pathways
    vertex.receiver = seq(1,(length(groupSize)+1)/2) # a numeric vector
    cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
    for (pathways.show in cellchat@netP$pathways){
      # pathways.show <- "TGFb"
      svglite(file =paste0(pathways.show, "_line.svg"))
      netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
      dev.off()
      svglite(file =paste0(pathways.show, "_circle.svg"))
      netVisual_aggregate(cellchat, signaling = c(pathways.show), layout = "circle", vertex.size = groupSize)
      dev.off()
      svglite(file =paste0(pathways.show, "_contribute.svg"))
      print(netAnalysis_contribution(cellchat, signaling = pathways.show))
      dev.off()
      svglite(file =paste0(pathways.show, "_signalingRole.svg"))
      # the slot 'netP' means the inferred intercellular communication network of signaling pathways
      netVisual_signalingRole(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10)
      dev.off()
    }
    
    ## Identify global communication patterns and major signals for specific cell groups
    # Identify outgoing communication pattern of secreting cells
    print('ploting pattern figures...')
    nPatterns = 5
    svglite(file =paste0('0_', data_source, sample_type, "_outgoing_patern_heatmap.svg"))
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, font.size = 8, width = 5, height = 12)
    dev.off()
    svglite(file =paste0('1_', data_source, sample_type, "_outgoing_patern_rivermap.svg"))
    print(netAnalysis_river(cellchat, pattern = "outgoing"))
    dev.off()
    svglite(file =paste0('2_', data_source, sample_type, "_outgoing_patern_dotmap.svg"))
    print(netAnalysis_dot(cellchat, pattern = "outgoing"))
    dev.off()
    
    svglite(file =paste0('0_', data_source, sample_type, "_incoming_patern_heatmap.svg"))
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, font.size = 8, width = 5, height = 12)
    dev.off()
    svglite(file =paste0('1_', data_source, sample_type, "_incoming_patern_rivermap.svg"))
    print(netAnalysis_river(cellchat, pattern = "incoming"))
    dev.off()
    svglite(file =paste0('2_', data_source, sample_type, "_incoming_patern_dotmap.svg"))
    print(netAnalysis_dot(cellchat, pattern = "incoming"))
    dev.off()
  }
}

{
  ## analyse cellchat res
  library(pheatmap)
  
  data_source_list = c('pla', 'deci')
  for (data_source in data_source_list){
    if (data_source == 'pla'){
      cc_list1 = c('VCTp.VCTp', 'VCTp.VCT', 'VCTp.SCT', 'VCTp.EVTp', 'VCTp.EVT', 'VCTp.FB', 'VCTp.VSM',
                   'VCT.VCTp', 'VCT.VCT', 'VCT.SCT', 'VCT.EVTp', 'VCT.EVT', 'VCT.FB', 'VCT.VSM',
                   'SCT.VCTp', 'SCT.VCT', 'SCT.SCT', 'SCT.EVTp', 'SCT.EVT', 'SCT.FB', 'SCT.VSM',
                   'EVTp.VCTp', 'EVTp.VCT', 'EVTp.SCT', 'EVTp.EVTp', 'EVTp.EVT', 'EVTp.FB', 'EVTp.VSM',
                   'EVT.VCTp', 'EVT.VCT', 'EVT.SCT', 'EVT.EVTp', 'EVT.EVT', 'EVT.FB', 'EVT.VSM', 
                   'FB.VCTp', 'FB.VCT', 'FB.SCT', 'FB.EVTp', 'FB.EVT', 'FB.FB', 'FB.VSM',
                   'VSM.VCTp', 'VSM.VCT', 'VSM.SCT', 'VSM.EVTp', 'VSM.EVT', 'VSM.FB', 'VSM.VSM')
      
      cc_list2 = c('VCTp.NK', 'VCTp.T', 'VCTp.B', 'VCTp.MC', 'VCTp.DC',
                   'VCT.NK', 'VCT.T', 'VCT.B', 'VCT.MC', 'VCT.DC',
                   'SCT.NK', 'SCT.T', 'SCT.B', 'SCT.MC', 'SCT.DC',
                   'EVT.NK', 'EVT.T', 'EVT.B', 'EVT.MC', 'EVT.DC',
                   'EVTp.NK', 'EVTp.T', 'EVTp.B', 'EVTp.MC', 'EVTp.DC',
                   'FB.NK', 'FB.T', 'FB.B', 'FB.MC', 'FB.DC',
                   'VSM.NK', 'VSM.T', 'VSM.B', 'VSM.MC', 'VSM.DC')
      
      cc_list3 = c('NK.VCTp', 'NK.VCT', 'NK.SCT', 'NK.EVTp', 'NK.EVT', 'NK.FB', 'NK.VSM',
                   'T.VCTp', 'T.VCT', 'T.SCT', 'T.EVTp', 'T.EVT', 'T.FB', 'T.VSM',
                   'B.VCTp', 'B.VCT', 'B.SCT', 'B.EVTp', 'B.EVT', 'B.FB', 'B.VSM',
                   'MC.VCTp', 'MC.VCT', 'MC.SCT', 'MC.EVTp', 'MC.EVT', 'MC.FB', 'MC.VSM',
                   'DC.VCTp', 'DC.VCT', 'DC.SCT', 'DC.EVTp', 'DC.EVT', 'DC.FB', 'DC.VSM')
      
      tag_list = c('PE', 'pla4-PE', 'pla5-PE', 'pla6-PE', 'pla7-PE')
    }else if(data_source == 'deci'){
      cc_list1 = c('FB1.FB1', 'FB1.FB2', 'FB1.SMC', 'FB1.DSC1', 'FB1.DSC2', 'FB1.Epi', 'FB1.LEC',
                   'FB2.FB1', 'FB2.FB2', 'FB2.SMC', 'FB2.DSC1', 'FB2.DSC2', 'FB2.Epi', 'FB2.LEC',
                   'SMC.FB1', 'SMC.FB2', 'SMC.SMC', 'SMC.DSC1', 'SMC.DSC2', 'SMC.Epi', 'SMC.LEC',
                   'DSC1.FB1', 'DSC1.FB2', 'DSC1.SMC', 'DSC1.DSC1', 'DSC1.DSC2', 'DSC1.Epi', 'DSC1.LEC',
                   'DSC2.FB1', 'DSC2.FB2', 'DSC2.SMC', 'DSC2.DSC1', 'DSC2.DSC2', 'DSC2.Epi', 'DSC2.LEC',
                   'Epi.FB1', 'Epi.FB2', 'Epi.SMC', 'Epi.DSC1', 'Epi.DSC2', 'Epi.Epi', 'Epi.LEC',
                   'LEC.FB1', 'LEC.FB2', 'LEC.SMC', 'LEC.DSC1', 'LEC.DSC2', 'LEC.Epi', 'LEC.LEC')
      
      cc_list2 = c('FB1.NK1', 'FB1.NK2', 'FB1.T1', 'FB1.T2', 'FB1.B', 'FB1.Gran', 'FB1.APC',
                   'FB2.NK1', 'FB2.NK2', 'FB2.T1', 'FB2.T2', 'FB2.B', 'FB2.Gran', 'FB2.APC',
                   'SMC.NK1', 'SMC.NK2', 'SMC.T1', 'SMC.T2', 'SMC.B', 'SMC.Gran', 'SMC.APC',
                   'DSC1.NK1', 'DSC1.NK2', 'DSC1.T1', 'DSC1.T2', 'DSC1.B', 'DSC1.Gran', 'DSC1.APC',
                   'DSC2.NK1', 'DSC2.NK2', 'DSC2.T1', 'DSC2.T2', 'DSC2.B', 'DSC2.Gran', 'DSC2.APC',
                   'Epi.NK1', 'Epi.NK2', 'Epi.T1', 'Epi.T2', 'Epi.B', 'Epi.Gran', 'Epi.APC',
                   'LEC.NK1', 'LEC.NK2', 'LEC.T1', 'LEC.T2', 'LEC.B', 'LEC.Gran', 'LEC.APC')
      # cc_list2 = c('FB1.NK1', 'FB2.NK2', 'SMC.T1', 'DSC1.T2', 'DSC2.B', 'Epi.Gran', 'LEC.APC',
      #              'FB1.NK1', 'FB2.NK1', 'SMC.NK1', 'DSC1.NK1', 'DSC2.NK1', 'Epi.NK1', 'LEC.NK1',)
      
      cc_list3 = c('NK1.FB1', 'NK1.FB2', 'NK1.SMC', 'NK1.DSC1', 'NK1.DSC2', 'NK1.Epi', 'NK1.LEC',
                   'NK2.FB1', 'NK2.FB2', 'NK2.SMC', 'NK2.DSC1', 'NK2.DSC2', 'NK2.Epi', 'NK2.LEC',
                   'T1.FB1', 'T1.FB2', 'T1.SMC', 'T1.DSC1', 'T1.DSC2', 'T1.Epi', 'T1.LEC',
                   'T2.FB1', 'T2.FB2', 'T2.SMC', 'T2.DSC1', 'T2.DSC2', 'T2.Epi', 'T2.LEC',
                   'B.FB1', 'B.FB2', 'B.SMC', 'B.DSC1', 'B.DSC2', 'B.Epi', 'B.LEC',
                   'Gran.FB1', 'Gran.FB2', 'Gran.SMC', 'Gran.DSC1', 'Gran.DSC2', 'Gran.Epi', 'Gran.LEC',
                   'APC.FB1', 'APC.FB2', 'APC.SMC', 'APC.DSC1', 'APC.DSC2', 'APC.Epi', 'APC.LEC')
      tag_list = c('PE', 'deci4-PE', 'deci6-PE', 'deci8-PE', 'deci9-PE')
    }
    for (tag in tag_list){
      print(sprintf('=======================================%s:%s-NM=======================================', data_source, tag))
      nmtag = 'NM'
      nm_wp = paste0('~/hfz/cellchat/', data_source, nmtag, '/')
      pe_wp = paste0('~/hfz/cellchat/', data_source, tag, '/')
      nm_lr_prob_df = read.csv(paste0(nm_wp, paste0(data_source, nmtag, '_LR_prob.csv')), sep = ',', header = T, row.names = 1)
      nm_lr_pval_df = read.csv(paste0(nm_wp, paste0(data_source, nmtag, '_LR_pval.csv')), sep = ',', header = T, row.names = 1)
      pe_lr_prob_df = read.csv(paste0(pe_wp, paste0(data_source, tag, '_LR_prob.csv')), sep = ',', header = T, row.names = 1)
      pe_lr_pval_df = read.csv(paste0(pe_wp, paste0(data_source, tag, '_LR_pval.csv')), sep = ',', header = T, row.names = 1)
      nm_pw_prob_df = read.csv(paste0(nm_wp, paste0(data_source, nmtag, '_pathway_prob.csv')), sep = ',', header = T, row.names = 1)
      pe_pw_prob_df = read.csv(paste0(pe_wp, paste0(data_source, tag, '_pathway_prob.csv')), sep = ',', header = T, row.names = 1)
      figure_wp = '~/hfz/cellchat/figures/pheatmap/'
      # nm_lr_prob_df = log10(nm_lr_prob_df + 1e-10)
      # nm_col = c()
      # for (i in rownames(nm_lr_prob_df)){
      #     if (mean(as.numeric(nm_lr_prob_df[i,])) > -10){
      #         nm_col = c(nm_col, i)
      #     }
      # }
      # nm_lr_prob_df = nm_lr_prob_df[nm_col,]
      # pe_lr_prob_df = log10(pe_lr_prob_df + 1e-10)
      # pe_col = c()
      # for (i in rownames(pe_lr_prob_df)){
      #     if (mean(as.numeric(pe_lr_prob_df[i,])) > -10){
      #         pe_col = c(pe_col, i)
      #     }
      # }
      # pe_lr_prob_df = pe_lr_prob_df[pe_col,]
      # i_lr_list = intersect(rownames(nm_lr_prob_df), rownames(pe_lr_prob_df))
      # nm_lr_list = rownames(nm_lr_prob_df)[which(!(rownames(nm_lr_prob_df) %in% i_lr_list))]
      # pe_lr_list = rownames(pe_lr_prob_df)[which(!(rownames(pe_lr_prob_df) %in% i_lr_list))]
      
      if (0){
        # plot cell-cell communication heatmap by pathways
        nm_pw_prob_df = log10(nm_pw_prob_df + 1e-10)
        pe_pw_prob_df = log10(pe_pw_prob_df + 1e-10)
        
        
        
        i_pw_list = intersect(rownames(nm_pw_prob_df), rownames(pe_pw_prob_df))
        nm_pw_list = rownames(nm_pw_prob_df)[which(!(rownames(nm_pw_prob_df) %in% i_pw_list))]
        pe_pw_list = rownames(pe_pw_prob_df)[which(!(rownames(pe_pw_prob_df) %in% i_pw_list))]
        
        for (col_name in c(cc_list1, cc_list2, cc_list3)){
          if (!(col_name %in% colnames(nm_pw_prob_df))){
            nm_pw_prob_df$a = NA
            colnames(nm_pw_prob_df)[length(colnames(nm_pw_prob_df))] = col_name
          }
          if (!(col_name %in% colnames(pe_pw_prob_df))){
            pe_pw_prob_df$a = NA
            colnames(pe_pw_prob_df)[length(colnames(pe_pw_prob_df))] = col_name
          }
        }
        
        cor_test_df_1 = as.data.frame(matrix(nrow = length(i_pw_list), ncol = 2), row.names = i_pw_list)
        colnames(cor_test_df_1) = c('stat', 'pval')
        for (row_pos in i_pw_list){
          tmp_res = cor.test(as.numeric(nm_pw_prob_df[row_pos, cc_list1]), as.numeric(pe_pw_prob_df[row_pos, cc_list1]))
          cor_test_df_1[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
          cor_test_df_1[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
        }
        cor_test_df_1 = cor_test_df_1[order(cor_test_df_1$stat), ]
        pw_list_1 = rownames(cor_test_df_1)[which(cor_test_df_1$stat < 0.5)]
        pw_list_1_nm = c(nm_pw_list, pw_list_1)
        pw_list_1_pe = c(pe_pw_list, pw_list_1)
        
        cor_test_df_2 = as.data.frame(matrix(nrow = length(i_pw_list), ncol = 2), row.names = i_pw_list)
        colnames(cor_test_df_2) = c('stat', 'pval')
        for (row_pos in i_pw_list){
          tmp_res = cor.test(as.numeric(nm_pw_prob_df[row_pos, cc_list2]), as.numeric(pe_pw_prob_df[row_pos, cc_list2]))
          cor_test_df_2[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
          cor_test_df_2[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
        }
        cor_test_df_2 = cor_test_df_2[order(cor_test_df_2$stat), ]
        pw_list_2 = rownames(cor_test_df_2)[which(cor_test_df_2$stat < 0.5)]
        pw_list_2_nm = c(nm_pw_list, pw_list_2)
        pw_list_2_pe = c(pe_pw_list, pw_list_2)
        
        cor_test_df_3 = as.data.frame(matrix(nrow = length(i_pw_list), ncol = 2), row.names = i_pw_list)
        colnames(cor_test_df_3) = c('stat', 'pval')
        for (row_pos in i_pw_list){
          tmp_res = cor.test(as.numeric(nm_pw_prob_df[row_pos, cc_list3]), as.numeric(pe_pw_prob_df[row_pos, cc_list3]))
          cor_test_df_3[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
          cor_test_df_3[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
        }
        cor_test_df_3 = cor_test_df_3[order(cor_test_df_3$stat), ]
        pw_list_3 = rownames(cor_test_df_3)[which(cor_test_df_3$stat < 0.5)]
        pw_list_3_nm = c(nm_pw_list, pw_list_3)
        pw_list_3_pe = c(pe_pw_list, pw_list_3)
        
        pdf(file = paste0(figure_wp, data_source, tag, '_cc_commu.pdf'), width = 10, height = 5)
        print(pheatmap(nm_pw_prob_df[pw_list_1_nm,cc_list1],cluster_rows = F, cluster_cols = F, main = 'Normal: non-immune cell communications',
                       gaps_col = seq(0,length(cc_list1)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black'))
        print(pheatmap(pe_pw_prob_df[pw_list_1_pe,cc_list1],cluster_rows = F, cluster_cols = F, main = 'PE: non-immune cell communications',
                       gaps_col = seq(0,length(cc_list1)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black'))
        print(pheatmap(nm_pw_prob_df[pw_list_2_nm,cc_list2],cluster_rows = F, cluster_cols = F, main = 'Normal: TBs to immune cell communications',
                       gaps_col = seq(0,length(cc_list2)-1, by = 5)[-1], clustering_method = 'average', na_col = 'black'))
        print(pheatmap(pe_pw_prob_df[pw_list_2_pe,cc_list2],cluster_rows = F, cluster_cols = F, main = 'PE: TBs to immune cell communications',
                       gaps_col = seq(0,length(cc_list2)-1, by = 5)[-1], clustering_method = 'average', na_col = 'black'))
        print(pheatmap(nm_pw_prob_df[pw_list_3_nm,cc_list3],cluster_rows = F, cluster_cols = F, main = 'Normal: immune cell to TBs communications',
                       gaps_col = seq(0,length(cc_list3)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black'))
        print(pheatmap(pe_pw_prob_df[pw_list_3_pe,cc_list3],cluster_rows = F, cluster_cols = F, main = 'PE: immune cell to TBs communications',
                       gaps_col = seq(0,length(cc_list3)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black'))
        dev.off()
      }
      if (1){
        # plot cell-cell communication heatmap by L-R pairs
        min_prob = 1
        for (col_name in colnames(nm_lr_prob_df)){
          tmp_min= min(nm_lr_prob_df[which(nm_lr_prob_df[, col_name] > 0), col_name])
          if (min_prob > tmp_min){
            min_prob = tmp_min
          } 
        }
        for (col_name in colnames(pe_lr_prob_df)){
          tmp_min= min(pe_lr_prob_df[which(pe_lr_prob_df[, col_name] > 0), col_name])
          if (min_prob > tmp_min){
            min_prob = tmp_min
          } 
        }
        min_prob = min_prob * 1e-1
        
        nm_lr_prob_df = log10(nm_lr_prob_df + min_prob)
        pe_lr_prob_df = log10(pe_lr_prob_df + min_prob)
        
        i_lr_list = intersect(rownames(nm_lr_prob_df), rownames(pe_lr_prob_df))
        nm_lr_list = rownames(nm_lr_prob_df)[which(!(rownames(nm_lr_prob_df) %in% i_lr_list))]
        pe_lr_list = rownames(pe_lr_prob_df)[which(!(rownames(pe_lr_prob_df) %in% i_lr_list))]
        
        for (col_name in c(cc_list1, cc_list2, cc_list3)){
          if (!(col_name %in% colnames(nm_lr_prob_df))){
            nm_lr_prob_df$a = NA
            colnames(nm_lr_prob_df)[length(colnames(nm_lr_prob_df))] = col_name
          }
          if (!(col_name %in% colnames(pe_lr_prob_df))){
            pe_lr_prob_df$a = NA
            colnames(pe_lr_prob_df)[length(colnames(pe_lr_prob_df))] = col_name
          }
          if (!(col_name %in% colnames(nm_lr_pval_df))){
            nm_lr_pval_df$a = NA
            colnames(nm_lr_pval_df)[length(colnames(nm_lr_pval_df))] = col_name
          }
          if (!(col_name %in% colnames(pe_lr_pval_df))){
            pe_lr_pval_df$a = NA
            colnames(pe_lr_pval_df)[length(colnames(pe_lr_pval_df))] = col_name
          }
        }
        {
          cor_test_df_1 = as.data.frame(matrix(nrow = length(i_lr_list), ncol = 2), row.names = i_lr_list)
          colnames(cor_test_df_1) = c('stat', 'pval')
          for (row_pos in i_lr_list){
            tmp_res = cor.test(as.numeric(nm_lr_prob_df[row_pos, cc_list1]), as.numeric(pe_lr_prob_df[row_pos, cc_list1]))
            cor_test_df_1[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
            cor_test_df_1[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
          }
          cor_test_df_1 = cor_test_df_1[order(cor_test_df_1$stat), ]
          lr_list_1 = rownames(cor_test_df_1)[which(cor_test_df_1$stat < 0.5)]
          lr_list_1_nm = c()
          for (lr in nm_lr_list){
            if (min(nm_lr_pval_df[lr, cc_list1][which(nm_lr_pval_df[lr, cc_list1] >= 0)]) <= 0.05){
              lr_list_1_nm = c(lr_list_1_nm, lr)
            }
          }
          lr_list_1_pe = c()
          for (lr in pe_lr_list){
            if (min(pe_lr_pval_df[lr, cc_list1][which(pe_lr_pval_df[lr, cc_list1] >= 0)]) <= 0.05){
              lr_list_1_pe = c(lr_list_1_pe, lr)
            }
          }
          lr_list_1 = c(lr_list_1_nm, lr_list_1_pe, lr_list_1)
          
          cor_test_df_2 = as.data.frame(matrix(nrow = length(i_lr_list), ncol = 2), row.names = i_lr_list)
          colnames(cor_test_df_2) = c('stat', 'pval')
          for (row_pos in i_lr_list){
            tmp_res = cor.test(as.numeric(nm_lr_prob_df[row_pos, cc_list2]), as.numeric(pe_lr_prob_df[row_pos, cc_list2]))
            cor_test_df_2[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
            cor_test_df_2[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
          }
          cor_test_df_2 = cor_test_df_2[order(cor_test_df_2$stat), ]
          lr_list_2 = rownames(cor_test_df_2)[which(cor_test_df_2$stat < 0.5)]
          lr_list_2_nm = c()
          for (lr in nm_lr_list){
            if (min(nm_lr_pval_df[lr, cc_list2][which(nm_lr_pval_df[lr, cc_list2] >= 0)]) <= 0.05){
              lr_list_2_nm = c(lr_list_2_nm, lr)
            }
          }
          lr_list_2_pe = c()
          for (lr in pe_lr_list){
            if (min(pe_lr_pval_df[lr, cc_list2][which(pe_lr_pval_df[lr, cc_list2] >= 0)]) <= 0.05){
              lr_list_2_pe = c(lr_list_2_pe, lr)
            }
          }
          lr_list_2 = c(lr_list_2_nm, lr_list_2_pe, lr_list_2)
          
          cor_test_df_3 = as.data.frame(matrix(nrow = length(i_lr_list), ncol = 2), row.names = i_lr_list)
          colnames(cor_test_df_3) = c('stat', 'pval')
          for (row_pos in i_lr_list){
            tmp_res = cor.test(as.numeric(nm_lr_prob_df[row_pos, cc_list3]), as.numeric(pe_lr_prob_df[row_pos, cc_list3]))
            cor_test_df_3[row_pos, 'stat'] = as.numeric(tmp_res$estimate)
            cor_test_df_3[row_pos, 'pval'] = as.numeric(tmp_res$p.value)
          }
          cor_test_df_3 = cor_test_df_3[order(cor_test_df_3$stat), ]
          lr_list_3 = rownames(cor_test_df_3)[which(cor_test_df_3$stat < 0.5)]
          lr_list_3_nm = c()
          for (lr in nm_lr_list){
            if (min(nm_lr_pval_df[lr, cc_list3][which(nm_lr_pval_df[lr, cc_list3] >= 0)]) <= 0.05){
              lr_list_3_nm = c(lr_list_3_nm, lr)
            }
          }
          lr_list_3_pe = c()
          for (lr in pe_lr_list){
            if (min(pe_lr_pval_df[lr, cc_list3][which(pe_lr_pval_df[lr, cc_list3] >= 0)]) <= 0.05){
              lr_list_3_pe = c(lr_list_3_pe, lr)
            }
          }
          lr_list_3 = c(lr_list_3_nm, lr_list_3_pe, lr_list_3)
        }
        
        {
          # getting differential expression FC
          # cc_list1
          lr_prob_df_1 = as.data.frame(matrix(nrow = length(lr_list_1), ncol = length(cc_list1)), row.names = lr_list_1)
          colnames(lr_prob_df_1) = cc_list1
          for (lr in lr_list_1){
            if (lr %in% nm_lr_list){
              lr_prob_df_1[lr, ] = (log10(min_prob) - nm_lr_prob_df[lr, cc_list1])/sqrt(abs(max(nm_lr_prob_df[lr, cc_list1][which(!(is.na(nm_lr_prob_df[lr, cc_list1])))])))
            }else if (lr %in% pe_lr_list){
              lr_prob_df_1[lr, ] = (pe_lr_prob_df[lr, cc_list1] - log10(min_prob))/sqrt(abs(max(pe_lr_prob_df[lr, cc_list1][which(!(is.na(pe_lr_prob_df[lr, cc_list1])))])))
            }else{
              lr_prob_df_1[lr, ] = (pe_lr_prob_df[lr, cc_list1] - nm_lr_prob_df[lr, cc_list1])/sqrt(abs(max(as.numeric(pe_lr_prob_df[lr, cc_list1][which(!(is.na(pe_lr_prob_df[lr, cc_list1])))]),
                                                                                                            as.numeric(nm_lr_prob_df[lr, cc_list1][which(!(is.na(nm_lr_prob_df[lr, cc_list1])))]))))
            }
          }
          # cc_list2
          lr_prob_df_2 = as.data.frame(matrix(nrow = length(lr_list_2), ncol = length(cc_list2)), row.names = lr_list_2)
          colnames(lr_prob_df_2) = cc_list2
          for (lr in lr_list_2){
            if (lr %in% nm_lr_list){
              lr_prob_df_2[lr, ] = (log10(min_prob) - nm_lr_prob_df[lr, cc_list2])/sqrt(abs(max(nm_lr_prob_df[lr, cc_list2][which(!(is.na(nm_lr_prob_df[lr, cc_list2])))])))
            }else if (lr %in% pe_lr_list){
              lr_prob_df_2[lr, ] = (pe_lr_prob_df[lr, cc_list2] - log10(min_prob))/sqrt(abs(max(pe_lr_prob_df[lr, cc_list2][which(!(is.na(pe_lr_prob_df[lr, cc_list2])))])))
            }else{
              lr_prob_df_2[lr, ] = (pe_lr_prob_df[lr, cc_list2] - nm_lr_prob_df[lr, cc_list2])/sqrt(abs(max(as.numeric(pe_lr_prob_df[lr, cc_list2][which(!(is.na(pe_lr_prob_df[lr, cc_list2])))]),
                                                                                                            as.numeric(nm_lr_prob_df[lr, cc_list2][which(!(is.na(nm_lr_prob_df[lr, cc_list2])))]))))
            }
          }
          # cc_list3
          lr_prob_df_3 = as.data.frame(matrix(nrow = length(lr_list_3), ncol = length(cc_list3)), row.names = lr_list_3)
          colnames(lr_prob_df_3) = cc_list3
          for (lr in lr_list_3){
            if (lr %in% nm_lr_list){
              lr_prob_df_3[lr, ] = (log10(min_prob) - nm_lr_prob_df[lr, cc_list3])/sqrt(abs(max(nm_lr_prob_df[lr, cc_list3][which(!(is.na(nm_lr_prob_df[lr, cc_list3])))])))
            }else if (lr %in% pe_lr_list){
              lr_prob_df_3[lr, ] = (pe_lr_prob_df[lr, cc_list3] - log10(min_prob))/sqrt(abs(max(pe_lr_prob_df[lr, cc_list3][which(!(is.na(pe_lr_prob_df[lr, cc_list3])))])))
            }else{
              lr_prob_df_3[lr, ] = (pe_lr_prob_df[lr, cc_list3] - nm_lr_prob_df[lr, cc_list3])/sqrt(abs(max(as.numeric(pe_lr_prob_df[lr, cc_list3][which(!(is.na(pe_lr_prob_df[lr, cc_list3])))]),
                                                                                                            as.numeric(nm_lr_prob_df[lr, cc_list3][which(!(is.na(nm_lr_prob_df[lr, cc_list3])))]))))
            }
          }
        }
        
        if (data_source == 'pla'){
          # pla and deci having different columns
          pdf(file = paste0(figure_wp, data_source, tag, '_cc_commu_lr.pdf'), width = 10, height = 8)
          {
            max1 = round(max(as.numeric(as.matrix(lr_prob_df_1))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_1)))))])+1)
            min1 = round(min(as.numeric(as.matrix(lr_prob_df_1))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_1)))))]))
            bk1 <- c(seq(min1,-0.01,by=0.01),seq(0,max1,by=0.01))
            print(pheatmap(lr_prob_df_1,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: non-immune cell communications',
                           gaps_col = seq(0,length(cc_list1)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min1,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max1,by=0.01)))),
                           legend_breaks=seq(min1,max1,2), breaks = bk1))
          }
          {
            max2 = round(max(as.numeric(as.matrix(lr_prob_df_2))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_2)))))])+1)
            min2 = round(min(as.numeric(as.matrix(lr_prob_df_2))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_2)))))]))
            bk2 <- c(seq(min2,-0.01,by=0.01),seq(0,max2,by=0.01))
            print(pheatmap(lr_prob_df_2,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: non-immune cell to immune cell communications',
                           gaps_col = seq(0,length(cc_list2)-1, by = 5)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min2,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max2,by=0.01)))),
                           legend_breaks=seq(min2,max2,2), breaks = bk2))
          }
          {
            max3 = round(max(as.numeric(as.matrix(lr_prob_df_3))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_3)))))])+1)
            min3 = round(min(as.numeric(as.matrix(lr_prob_df_3))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_3)))))]))
            bk3 <- c(seq(min3,-0.01,by=0.01),seq(0,max3,by=0.01))
            print(pheatmap(lr_prob_df_3,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: immune cell to non-immune cell communications',
                           gaps_col = seq(0,length(cc_list3)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min3,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max3,by=0.01)))),
                           legend_breaks=seq(min3,max3,2), breaks = bk3))
          }
          dev.off()
        }else{
          pdf(file = paste0(figure_wp, data_source, tag, '_cc_commu_lr.pdf'), width = 10, height = 9)
          {
            max1 = round(max(as.numeric(as.matrix(lr_prob_df_1))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_1)))))])+1)
            min1 = round(min(as.numeric(as.matrix(lr_prob_df_1))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_1)))))]))
            bk1 <- c(seq(min1,-0.01,by=0.01),seq(0,max1,by=0.01))
            print(pheatmap(lr_prob_df_1,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: non-immune cell communications',
                           gaps_col = seq(0,length(cc_list1)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min1,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max1,by=0.01)))),
                           legend_breaks=seq(min1,max1,2), breaks = bk1))
          }
          {
            max2 = round(max(as.numeric(as.matrix(lr_prob_df_2))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_2)))))])+1)
            min2 = round(min(as.numeric(as.matrix(lr_prob_df_2))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_2)))))]))
            bk2 <- c(seq(min2,-0.01,by=0.01),seq(0,max2,by=0.01))
            print(pheatmap(lr_prob_df_2,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: non-immune cell to immune cell communications',
                           gaps_col = seq(0,length(cc_list2)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min2,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max2,by=0.01)))),
                           legend_breaks=seq(min2,max2,2), breaks = bk2))
          }
          {
            max3 = round(max(as.numeric(as.matrix(lr_prob_df_3))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_3)))))])+1)
            min3 = round(min(as.numeric(as.matrix(lr_prob_df_3))[which(!(is.na(as.numeric(as.matrix(lr_prob_df_3)))))]))
            bk3 <- c(seq(min3,-0.01,by=0.01),seq(0,max3,by=0.01))
            print(pheatmap(lr_prob_df_3,cluster_rows = F, cluster_cols = F, main = 'PE-Normal: immune cell to non-immune cell communications',
                           gaps_col = seq(0,length(cc_list3)-1, by = 7)[-1], clustering_method = 'average', na_col = 'black',
                           color = c(colorRampPalette(colors = c("royalblue","snow"))(length(seq(min3,-0.01,by=0.01))),
                                     colorRampPalette(colors = c("snow","red"))(length(seq(0,max3,by=0.01)))),
                           legend_breaks=seq(min3,max3,2), breaks = bk3))
          }
          dev.off()
          
        }
      }
    }
  }
}
















