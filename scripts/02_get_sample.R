library(ggplot2)
library(Seurat)
library(reticulate)
library(future)
library(future.apply)
options(future.globals.maxSize = 15 * 1024^3)


setwd("/public4/chy/206/chenhy/mouse_brain/RData")
load("mm1M_sub20w.RData")  ## obj_integrated(merge seurat object), hvg_features
norm.counts <- obj_integrated[["RNA"]]@counts
c_num <- c(1000,3000,5000,10000,20000,30000,40000,50000)

dir.create(file.path("/public4/chy/206/chenhy/mouse_brain/RData","mm1M_sub20w_sample_1234"),showWarnings = FALSE)
setwd("/public4/chy/206/chenhy/mouse_brain/RData/mm1M_sub20w_sample_1234")
for (i in 1:8){
	set.seed(1234)
	f.name <- paste("obj_integrated",c_num[i],"cell.RData",sep = "_")
	cols = sample(1:ncol(norm.counts), c_num[i], replace = FALSE)
	sub_data <- norm.counts[,cols]
	print(dim(sub_data))
	save(sub_data,file =f.name )
}

save(norm.counts,file = "obj_integrated_all_cell.RData")
