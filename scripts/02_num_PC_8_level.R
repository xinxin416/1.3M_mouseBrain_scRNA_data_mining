# 每个数量细胞PC分析
##! /public/home/fanlj/software/miniconda3/envs/mlf_gatk/bin/Rscript
library(ggplot2)
library(stringr)
library(Seurat)
library(reticulate)
library(future)
library(future.apply)
options(future.globals.maxSize = 15 * 1024^3)

wdir = "/public4/chy/206/chenhy/mouse_brain/sub_20w_rep/"
setwd(wdir)

seed = "3456"
file_name <- list.files("./RData/sub20w_sample_3456/")

p_EL = list()
p_JS = list()

#### 综合考虑，当细胞数达到10000时，显著PC开始>100，为降低计算时间，只计算从1000到10000的四个水平 c(1,2,4,7)
for(file in file_name){ 
	
	j = 1
	print(file)
	name <- str_split(file,'_')[[1]][3] 
	load(paste0("./RData/sub20w_sample_3456/",file))

	num_PC <- CreateSeuratObject(counts = sub_data)
	# print(num_PC[1:4,1:4])
	# all.genes <- rownames(num_PC)
	num_PC <- NormalizeData(num_PC, normalization.method = "LogNormalize", scale.factor = 10000)
	num_PC <- FindVariableFeatures(num_PC, selection.method = "vst",nfeatures = 2000)
	num_PC <- ScaleData(num_PC)
	#num_PC <- SCTransform(num_PC,variable.features.n = 2000, verbose = FALSE) 
	num_PC <- RunPCA(num_PC, features = VariableFeatures(num_PC),verbose = FALSE, npcs = 100)
	num_PC <- JackStraw(num_PC, num.replicate = 100,dims = 100)
	num_PC <- ScoreJackStraw(num_PC, dims = 1:100)
	p_JS[[j]] <- JackStrawPlot(num_PC, dims = 1:100) + ggtitle(paste0(format(as.numeric(name),big.mark = ',')," cells"))
	p_EL[[j]] <- ElbowPlot(num_PC, ndims = 100,reduction = "pca") + ggtitle(paste0(format(as.numeric(name),big.mark = ',')," cells"))

	### 绘图可视化
	dir.create(file.path(paste0(wdir,"Results/PC_results/"),"sub20w_sample_3456"),showWarnings = FALSE)
	path2 <- paste0("./Results/PC_results/sub20w_sample_3456/sub20w_new_",name,"_PC_JS.pdf")
	pdf(path2,width = 14, height = 8)
	print(p_JS[[j]])
	dev.off()

	path2 <- paste0("./Results/PC_results/sub20w_sample_3456/sub20w_new_",name,"_PC_EL.pdf")
	pdf(path2,width = 14, height = 8)
	print(p_EL[[j]])
	dev.off()

	j = j+1 
}


