## 不同细胞数样本进行Seurat分群分析

#! /public/home/fanlj/software/miniconda3/envs/mlf_gatk/bin/Rscript
library(ggplot2)
library(tidyverse)
library(Seurat)
library(reticulate)
library(future)
library(future.apply)
options(future.globals.maxSize = 15 * 1024^3)
library(SingleR)
library(scater)
np <- import("numpy")


wdir = "/home/chenhy/mouse_brain/sub_20w_rep/"
setwd(wdir)
source("./scripts/03_func_run_seurat.R")  ### func: run_seurat(rawcounts,ndim)

file_name <- list.files("./RData/sub20w_sample_1234/")
dir.create(file.path(paste0(wdir,"/Results/cluster_result_PC30/"),"sub20w_sample_1234"),showWarnings = FALSE)
dir.create(file.path(paste0(wdir,"/Results/singler_anno/"),"sub20w_sample_1234"),showWarnings = FALSE)

### 载入参考数据集
load("./RData/Campbell_ref_sce.RData")  ## sce_ref_Campbell
load("./RData/Chen_ref_sce.RData")  ## sce_ref_Chen

p_um = list()
p_ts = list()

j = 1
sink("./scripts/log_cluster_singler_20w.txt",append = TRUE)
for(file in file_name){

    cat("run cluster in file:")
    print(file)
    name <- str_split(file,'_')[[1]][3] 
    load(paste0("./RData/sub20w_sample_1234/",file))  ## rawcounts, sub_data

    ### 调用函数，聚类,singleR注释
    num_Seurat <- run_seurat(sub_data,30)  ### return obj_seurat done clustering

    ### 聚类结果可视化、注释
    p_um[[j]] <- DimPlot(num_Seurat, reduction = "umap" , label = T) + ggtitle(paste0("umap",format(as.numeric(name),big.mark = ',')," cells"))
    p_ts[[j]] <- DimPlot(num_Seurat, reduction = "tsne" , label = T) + ggtitle(paste0("tsne",format(as.numeric(name),big.mark = ',')," cells"))
    path3 <- paste0("./Results/cluster_result_PC30/sub20w_sample_1234/sub20w_cluster_",name,".pdf")
    pdf(path3)
    print(p_um[[j]])
    print(p_ts[[j]])
    dev.off()

    #### singler 注释,用两个数据集Campbell，Chen，两种类型，cell_type1，cell_ontology_class
    test_sce <- as.SingleCellExperiment(num_Seurat)

    pred_Campbell_type.brain <- SingleR(test = test_sce, ref = sce_ref_Campbell, labels = sce_ref_Campbell$cell_type1, 
        method = "cluster", clusters = colData(test_sce)$seurat_clusters)

    pred_Chen_go.brain <- SingleR(test = test_sce, ref = sce_ref_Chen, labels = sce_ref_Chen$cell_ontology_class, 
        method = "cluster", clusters = colData(test_sce)$seurat_clusters)

    pred_Chen_type.brain <- SingleR(test = test_sce, ref = sce_ref_Chen, labels = sce_ref_Chen$cell_type1, 
        method = "cluster", clusters = colData(test_sce)$seurat_clusters)

    #### 统计注释类群数
    cat("#number of cell types:","\n")
    cat("name","\tCampbell","\tChen\n")
    cat(name,"\t",length(unique(pred_Campbell_type.brain$labels)),"\t",length(unique(pred_Chen_type.brain$labels)),"\n")

    cat("#number of cell ontology class:","\n")
    cat("name","\tChen\n")
    cat(name,"\t",length(unique(pred_Chen_go.brain$labels)),"\n")

    #### 注释结果可视化
    num_Seurat[["Chen_go_labels"]] <- 
        pred_Chen_go.brain$labels[match(num_Seurat[[]][["seurat_clusters"]], rownames(pred_Chen_go.brain))]

    num_Seurat[["Campbell_ctype_labels"]] <- 
        pred_Campbell_type.brain$labels[match(num_Seurat[[]][["seurat_clusters"]], rownames(pred_Campbell_type.brain))]

    num_Seurat[["Chen_ctype_labels"]] <- 
        pred_Chen_type.brain$labels[match(num_Seurat[[]][["seurat_clusters"]], rownames(pred_Chen_type.brain))]


    path1 = paste0("./Results/singler_anno/sub20w_sample_1234/sub20w_",name,"_cells_anno.pdf")
    pdf(path1)
    print(DimPlot(num_Seurat, reduction = "umap" , label = T,group.by = "Chen_go_labels") + ggtitle(paste0("Chen GO ",format(as.numeric(name),big.mark = ',')," cells")))
    print(DimPlot(num_Seurat, reduction = "umap" , label = T,group.by = "Campbell_ctype_labels") + ggtitle(paste0("Campbell cell type ",format(as.numeric(name),big.mark = ',')," cells")))
    print(DimPlot(num_Seurat, reduction = "umap" , label = T,group.by = "Chen_ctype_labels") + ggtitle(paste0("Chen cell type ",format(as.numeric(name),big.mark = ',')," cells")))
    dev.off()

    path1 = paste0("./Results/singler_anno/sub20w_sample_1234/sub20w_",name,"_cells_predScore.pdf")
    pdf(path1,width = 14, height = 8)
    plotScoreHeatmap(pred_Chen_go.brain,main = "Chen GO")
    plotScoreHeatmap(pred_Campbell_type.brain,main = "Campbell cell type")
    plotScoreHeatmap(pred_Chen_type.brain,main = "Chen cell type")
    dev.off()

    j = j+1

}
sink()





