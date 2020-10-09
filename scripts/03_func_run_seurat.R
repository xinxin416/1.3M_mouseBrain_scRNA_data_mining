library(ggplot2)
library(tidyverse)
library(Seurat)
library(reticulate)

### 基本流程包括标准化、特征选择、归一化，可以通过SCTransform实现
run_seurat <- function(rawcounts,ndim){
	obj_seurat <- CreateSeuratObject(counts = rawcounts)
    #obj_seurat <- ScaleData(obj_seurat)
    obj_seurat <- SCTransform(obj_seurat,variable.features.n = 2000, verbose = FALSE)  ### stored data in SCT assay
    obj_seurat <- RunPCA(obj_seurat,feature = VariableFeatures(obj_seurat),verbose = FALSE, npcs = ndim)
    obj_seurat <- FindNeighbors(obj_seurat, dims = 1:ndim)
    obj_seurat <- FindClusters(obj_seurat,resolution = 1)
    obj_seurat <- RunUMAP(obj_seurat, dims = 1:ndim,umap.method = "umap-learn",metric = "correlation")
    obj_seurat <- RunTSNE(obj_seurat, verbose = FALSE,dims = 1:ndim,check_duplicates = FALSE)
    return(obj_seurat)
}