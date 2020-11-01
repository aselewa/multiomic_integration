library(Seurat)
library(Signac)
library(ggplot2)
library(future)

get_integration_anchors <- function(RNA, ATAC, k.anchor = 5, k.filter = 200, k.score = 30){
  
  genes.use <- VariableFeatures(RNA)
  refdata <- GetAssayData(RNA, assay = "RNA", slot = "data")[genes.use, ]

  plan(strategy = "multicore", workers = 10)
  options(future.globals.maxSize = 12 * 1024^3)  
  transfer.anchors <- FindTransferAnchors(reference = RNA, 
                                          query = ATAC, 
                                          features = genes.use, 
                                          reference.assay = "RNA", 
                                          query.assay = "ACTIVITY", 
                                          reduction = "cca", 
                                          k.anchor = k.anchor,
                                          k.filter = k.filter, 
                                          k.score = k.score, 
                                          normalization.method = "LogNormalize")
  
  return(transfer.anchors)
  
}

transfer_data <- function(ATAC, data, anchors){
  
  imputed.labs <- TransferData(anchorset = anchors, refdata = data, weight.reduction = ATAC[["lsi"]], dims = 1:10)

  rna.ident <- factor(imputed.labs$predicted.id, levels = sort(unique(data)))
  names(rna.ident) <- rownames(imputed.labs)
  ATAC$RNA_IDENT <- rna.ident
  ATAC$PREDICTED_SCORE_MAX <- imputed.labs$prediction.score.max
  
  return(ATAC)
  
}

run_combed <- function(){
  p1 <- DimPlot(heart_atac, group.by = 'rna.ident', cols = colz) + ggtitle("scATAC-seq")
  p2 <- DimPlot(heart_rna, group.by = 'seurat_clusters', cols = colz) + ggtitle("scRNA-seq")
  CombinePlots(list(p1=p1,p2=p2))
  
  imputed.rna <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = heart_atac[["lsi"]], dims = 1:10)
  heart_atac[["RNA"]] <- imputed.rna
  #mse <- sum((Matrix::as.matrix(imputed.rna@data) - Matrix::as.matrix(obj@assays$RNA@data))^2) / (dim(imputed.rna)[1] * dim(imputed.rna)[2])
  heart_rna$tech <- "10X_RNA"
  heart_atac$tech <- "10X_ATAC"
  coembed <- merge(x = heart_rna, y = heart_atac)
  
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  
  DimPlot(coembed, group.by = "tech")
}