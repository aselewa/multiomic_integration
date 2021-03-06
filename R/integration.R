library(Seurat)
library(Signac)

setwd('/project2/gca/aselewa/integration_benchmarks/')
source('R/seurat_integration.R')

comb <- Matrix::readMM(file = 'data/filtered_feature_bc_matrix/matrix.mtx')
features <- readr::read_tsv(file = 'data/filtered_feature_bc_matrix/features.tsv', col_names = F)
barcodes <- readr::read_tsv(file = 'data/filtered_feature_bc_matrix/barcodes.tsv', col_names = F)

rna <- comb[features$X3 == "Gene Expression", ]
atac <- comb[features$X3 == "Peaks", ]
rna.feats <- features[features$X3 == "Gene Expression",]
atac.feats <- features[features$X3 == "Peaks",]

rownames(rna) <- rna.feats$X2
colnames(rna) <- barcodes$X1
rownames(atac) <- atac.feats$X2
colnames(atac) <- barcodes$X1

s.rna <- CreateSeuratObject(counts = rna, project = "integration_RNA", min.cells = 10)
s.rna[["percent.mt"]] <- PercentageFeatureSet(s.rna, pattern = "^MT-")
VlnPlot(s.rna, features = c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.1)

s.rna <- subset(s.rna, subset = nFeature_RNA > 600 & nFeature_RNA < 2200)
s.rna <- NormalizeData(s.rna)
s.rna <- FindVariableFeatures(s.rna)
s.rna <- ScaleData(s.rna)
s.rna <- RunPCA(s.rna, verbose = F)
s.rna <- FindNeighbors(s.rna, reduction = "pca", dims = 1:10)
s.rna <- RunUMAP(s.rna, reduction = "pca", dims = 1:10)
s.rna <- FindClusters(s.rna, reduction = "pca", resolution = 0.15)
DimPlot(s.rna, reduction = 'umap', label=T)

atac <- atac[,Cells(s.rna)]
s.atac <- CreateChromatinAssay(counts = atac, sep = c(':','-'), fragments = 'data/pbmc_unsorted_10k_atac_fragments.tsv.gz')
s.atac <- CreateSeuratObject(counts = s.atac, assay = "peaks", project = "integration_ATAC")
s.atac <- RunTFIDF(s.atac)
s.atac <- FindTopFeatures(s.atac, min.cutoff = 'q0')
s.atac <- RunSVD(s.atac)
s.atac <- FindNeighbors(object = s.atac, reduction = 'lsi', dims = 2:10)
s.atac <- RunUMAP(object = s.atac, reduction = 'lsi', dims = 2:10)
s.atac <- FindClusters(object = s.atac, verbose = FALSE, algorithm = 3, resolution = 0.1)
DimPlot(s.atac)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(s.atac) <- annotations
gene.activities <- GeneActivity(s.atac)

s.atac[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
s.atac <- NormalizeData(object = s.atac, assay = 'ACTIVITY',normalization.method = 'LogNormalize', scale.factor = median(s.atac$nCount_ACTIVITY))

s.rna <- RenameIdents(s.rna, '0' = 'CD14+ Monocytes','1'='Naive CD4+ T', '2'='Memory CD4+ T','3'='NK','4'='CD8+ T','5'='B','6'='NK')
s.rna$clusters <- factor(as.character(s.rna@active.ident), levels = sort(unique(as.character(s.rna@active.ident))))
s.atac$clusters <- factor(as.character(s.atac@active.ident), levels = sort(unique(as.character(s.atac@active.ident))))

# Get anchors for several values of k anchor and k filter

k.anchors <- c(1, 5, 10, 15, 20, 50)
for(k.anchor in k.filters){
  anchors <- get_integration_anchors(RNA = s.rna, ATAC = s.atac, k.anchor = k.anchor)
  saveRDS(anchors, paste0('data/anchors/anchors_kfilter_',k.anchor,'.rds'))
}

k.filters <- c(10, 100, 200, 500, Inf)
for(k.filter in k.filters){
  anchors <- get_integration_anchors(RNA = s.rna, ATAC = s.atac, k.filter = ifelse(is.infinite(k.filter), NA, k.filter))
  saveRDS(anchors, paste0('data/anchors/anchors_kfilter_',k.filter,'.rds'))
}

saveRDS(s.rna, 'data/rna_seurat.rds')
saveRDS(s.atac, 'data/atac_seurat.rds')

## Aggregate into pseudo cells

s.rna <- FindClusters(s.rna, resolution = 50) # fine-scale clustering

aggregate.counts <- list()
type <- c()
og.counts <- as.matrix(s.rna@assays$RNA@counts)
for(i in unique(s.rna$seurat_clusters)){
  sub.counts <- og.counts[,s.rna$seurat_clusters == i]
  aggregate.counts[[i]] <- rowSums(sub.counts)
  type <- c(type, as.character(s.rna$clusters[s.rna$seurat_clusters == i][1]))
}
aggregate.counts <- as.data.frame(aggregate.counts)
colnames(aggregate.counts) <- 1:ncol(aggregate.counts)

srna_agg <- CreateSeuratObject(counts = aggregate.counts, project = "aggregate_RNA")
srna_agg$clusters <- type
srna_agg <- NormalizeData(srna_agg)
srna_agg <- FindVariableFeatures(srna_agg)
srna_agg <- ScaleData(srna_agg)
srna_agg <- RunPCA(srna_agg, verbose = F)
srna_agg <- FindNeighbors(srna_agg, reduction = "pca", dims = 1:10)
srna_agg <- RunUMAP(srna_agg, reduction = "pca", dims = 1:10)

aggregate.counts <- list()
og.counts <- as.matrix(s.atac@assays$peaks@counts)
for(i in unique(s.rna$seurat_clusters)){
  sub.counts <- og.counts[ , s.rna$seurat_clusters == i]
  aggregate.counts[[i]] <- rowSums(sub.counts)
}
aggregate.counts <- as.data.frame(aggregate.counts)
colnames(aggregate.counts) <- 1:ncol(aggregate.counts)

satac_agg <- CreateChromatinAssay(counts = aggregate.counts, sep = c(':','-'))
satac_agg <- CreateSeuratObject(counts = satac_agg, assay = "peaks", project = "aggregate_peaks")
satac_agg$clusters <- type
satac_agg <- RunTFIDF(satac_agg)
satac_agg <- FindTopFeatures(satac_agg, min.cutoff = 'q0')
satac_agg <- RunSVD(satac_agg)
satac_agg <- FindNeighbors(object = satac_agg, reduction = 'lsi', dims = 2:10)
satac_agg <- RunUMAP(object = satac_agg, reduction = 'lsi', dims = 2:10)

saveRDS(srna_agg, 'data/agg_rna_seurat.rds')
saveRDS(satac_agg, 'data/agg_atac_seurat.rds')
