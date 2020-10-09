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
source("./scripts/03_func_run_seurat.R")  ### func: run_seurat(rawcounts,ndim)，return obj_seurat

seed = "3456"
file_name <- list.files("./RData/sub20w_sample_3456/")
dir.create(file.path(paste0(wdir,"/Results/diff_markers/"),"sub20w_sample_3456"),showWarnings = FALSE)

### 载入参考数据集
#load("./RData/Campbell_ref_sce.RData")  ## sce_ref_Campbell
load("./RData/Chen_ref_sce.RData")  ## sce_ref_Chen

p_um = list()
p_ts = list()
num_markers = list()
singler_anno_chen = list()

j = 1
sink("./scripts/log_diffGenes_20w.txt",append = TRUE)
cat("run diff markers in sample seed:")
print(seed)
for(file in file_name){

    cat("run cluster in file:")
    print(file)
    name <- str_split(file,'_')[[1]][3] 
    load(paste0("./RData/sub20w_sample_3456/",file))  ## rawcounts, sub_data
    ### 调用函数，运行Seurat 基本流程
    num_Seurat <- run_seurat(sub_data,30)  ### return obj_seurat done clustering

    ### 差异分析，确定各个簇marker gene,保存输出矩阵存为list，最后统一保存到该抽样结果文件中。
    cat("info of num_Seurat is:")
    print(num_Seurat)  ## 激活SCT assay
    #### 一个坑，参数中有点号，赋值不能有空格
    num_markers[[name]]<-FindAllMarkers(num_Seurat,features=VariableFeatures(object = num_Seurat),logfc.threshold=0.25,test.use="roc",only.pos=TRUE,return.thresh=0.05)
    cat(">significant markers from sample ",name," cells is:")
    print(dim(num_markers[[name]][1]))
    print(head(num_markers[[name]],2))

    ### singleR 注释，只用Chen 数据集
    sce_test <- SingleCellExperiment(list(counts=sub_data))
    sce_test <- logNormCounts(sce_test)
    cat("info of sce object is:")
    print(sce_test)

    cat("cell info comparing from seurat and sce object:")
    # print(head(num_Seurat@meta.data),6)
    # print(head(colData(sce_test),6))

    pred_Chen_type.brain <- SingleR(test = sce_test, ref = sce_ref_Chen, labels = sce_ref_Chen$cell_type1, 
        method = "cluster", clusters = num_Seurat@meta.data$seurat_clusters)
    singler_anno_chen[[name]] <- pred_Chen_type.brain

    #### 统计注释类群数
    cat("#number of cell types:","\n")
    cat("name","\tChen\n")
    cat(name,"\t",length(unique(pred_Chen_type.brain$labels)),"\n")

    #### 注释结果可视化
    num_Seurat[["Chen_ctype_labels"]] <- 
        pred_Chen_type.brain$labels[match(num_Seurat[[]][["seurat_clusters"]], rownames(pred_Chen_type.brain))]


    ### 聚类结果可视化、注释
    p_um[[j]] <- DimPlot(num_Seurat, reduction = "umap" , label = T) + ggtitle(paste0("umap",format(as.numeric(name),big.mark = ',')," cells"))
    p_ts[[j]] <- DimPlot(num_Seurat, reduction = "tsne" , label = T) + ggtitle(paste0("tsne",format(as.numeric(name),big.mark = ',')," cells"))
    path3 <- paste0("./Results/diff_markers/sub20w_sample_3456/sub20w_new_cluster_",name,".pdf")
    pdf(path3)
    print(p_um[[j]])
    print(p_ts[[j]])
    dev.off()

    ### singleR注释可视化
    path1 = paste0("./Results/diff_markers/sub20w_sample_3456/sub20w_new_",name,"_cells_anno.pdf")
    pdf(path1)
    print(DimPlot(num_Seurat, reduction = "umap" , label = T,group.by = "Chen_ctype_labels") + ggtitle(paste0("Chen cell type ",format(as.numeric(name),big.mark = ',')," cells")))
    plotScoreHeatmap(pred_Chen_type.brain,main = "Chen cell type")
    dev.off()


    j = j+1

}
save(num_markers,file = "./Results/diff_markers/sub20w_sample_3456_diff_markers.RData") 
save(singler_anno_chen,file = "./Results/diff_markers/sub20w_sample_3456_singler_anno_chen.RData")
sink()





