setwd('/project2/gca/aselewa/integration_benchmarks/')
source('R/seurat_integration.R')

s.rna <- readRDS('data/Seurat/rna_seurat.rds')
s.atac <- readRDS('data/Seurat/atac_seurat.rds')

anchor.files <- list.files(path = 'data/anchors/', pattern = 'anchors_kfilter_*', full.names = T)
anchor.list <- lapply(anchor.files, FUN = readRDS)
n <- sapply(strsplit(anchor.files, split='_'), function(x){x[3]})
n <- gsub(pattern = '.rds', replacement = '', x = n)
names(anchor.list) <- as.numeric(n)

res <- seq(0.05, 1.5, 0.1)
lab.acc.list <- list()
nclusts.list <- list()
for(name in n){
  lab.acc <- rep(0, length(res))
  nclusts <- rep(0, length(res))
  for(i in 1:length(res)){
    s.rna <- FindClusters(s.rna, resolution = res[i], verbose = F)
    s.atac$clusters <- Idents(s.rna)
    s.atac <- suppressMessages(transfer_data(ATAC = s.atac, data = Idents(s.rna), anchors = anchor.list[[name]]))
    lab.acc[i] <- sum(s.atac$RNA_IDENT==s.atac$clusters)/length(s.atac$RNA_IDENT)
    nclusts[i] <- length(levels(s.rna))
  }
  lab.acc.list[[name]] <- lab.acc
  nclusts.list[[name]] <- nclusts
}

saveRDS(lab.acc.list, file = 'data/transfer_label_benchmarks/label_transfer_accuracy_kfilter_2.rds')
saveRDS(nclusts.list, file = 'data/transfer_label_benchmarks/label_transfer_nclusts_kfilter_2.rds')
saveRDS(res, file='data/transfer_label_benchmarks/resolutions.rds')