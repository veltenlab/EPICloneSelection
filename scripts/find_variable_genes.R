library(Seurat)
nestorowa <- read.table(gzfile("/users/lvelten/lvelten/Data/Nestorowa2015/GSE81682_HTSeq_counts.txt.gz"), header=T, row.names=1)
nesto.seurat <- CreateSeuratObject(nestorowa, min.cells = 5)
nesto.seurat[["percent.mt"]] <- PercentageFeatureSet(nesto.seurat, pattern = "Mt-")
nesto.seurat <- NormalizeData(nesto.seurat)
nesto.seurat <- FindVariableFeatures(nesto.seurat)
var.genes <- nesto.seurat@assays$RNA@var.features
write.csv(var.genes, '/users/lvelten/project/Methylome/analysis/selection_pipeline/misc/variable_genes_Nestorowa2015.csv')