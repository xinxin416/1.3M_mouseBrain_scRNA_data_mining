

### 读取hdf5 文件，根据细胞列名进行拆分文件，形成10x 文件格式
### python 进行拆分，gzip 进行解压
### 拆分步骤
### 1. hdf5 文件中data 数据集通过scipy 将其压缩成稀疏矩阵，而后取子集存储为mtx 文件
### 2. 读取features、barcodes 文件，将barcodes 文件进行取子集
### 3. 依次创建拆分文件夹，每个目录下存放该三个文件，其中features（27998 行），barcodes(20409 行)，mtx 对应。
### 
### 通过for 循环拆分，根据barcodes 拆分得到64个文件。
### 创建Seurat 对象，做基本过滤，并将其合并，最终得到2个数据集，而后进行抽样。

### python HDF5 文件处理
import h5py
import numpy as np 
import pandas as pd 
import scipy.sparse as sp_sparse
from scipy import io,sparse 
from pathlib import Path
import os

os.chdir("../rawdata_scRNA")
f = h5py.File("1M_neurons_filtered_gene_bc_matrices_h5.h5",'r')

[print(f[item]) for item in f.keys()]

list(f['mm10'].keys())
indices = f['mm10']['indices']  ## 0下标的行对应data 元素
indptr = f['mm10']['indptr']  ## 0下标的列对应data 元素
shape = f['mm10']['shape']
mtx = f['mm10']['data'][:]  ###Nonzero UMI counts in column-major order

matrix = sp_sparse.csc_matrix((mtx, indices, indptr), shape=shape)

### 读取features, barcodes文件
fg = pd.read_table('features.tsv.gz',sep = '\t',names = ["Ensembl","geneID"])
fg.drop_duplicates(subset=['geneID']).shape  ## 无重复
fg.index = fg.index + 1 ## 设置索引从1开始
 
fb = pd.read_table('barcodes.tsv.gz')  ## 很多列
fc = fb.loc[:,'b']  ## 重命名失败
fc.index = fc.index +1 


### 创建目录，取子集
nrow = 20409
#left,right = 0, 20409
#left,right = 40818, 61227
for i in range(3,64):   ### range， 1：n-1
	path = "/home/chenhy/mouse_brain/rawdata_scRNA/x" + str(i) 
	print(path) 
	Path(path).mkdir(parents=True, exist_ok=True) 
	os.chdir(path)
	mtx_sub1 = matrix[:,left:right]
	io.mmwrite("matrix.mtx", mtx_sub1)
	fg.to_csv("features.tsv",sep = "\t",header =False,index = False)
	cells_sub = fc[left:right]  ### loc 函数切片符左闭右闭
	cells_sub.to_csv("barcodes.tsv",sep = "\t",header = False,index = False)
	left,right = right,right+nrow 
	print(left)
	print(right)
	os.chdir("../")

### 最后一个文件
left, right = 1285767,1306127
path = "/home/chenhy/mouse_brain/rawdata_scRNA/x" + str(64)
print(path)
Path(path).mkdir(parents=True, exist_ok=True)
os.chdir("./x64")
mtx_sub1 = matrix[:,left:right]
io.mmwrite("matrix.mtx", mtx_sub1)
fg.to_csv("features.tsv",sep = "\t",header =False,index = False)
cells_sub = fc[left:right]
cells_sub.to_csv("barcodes.tsv",sep = "\t",header = False,index = False)
print(left)
print(right)
os.chdir("../")


#### 压缩所有文件

gzip x*/*

#### R 环境，创建Seurat 对象
library(Seurat)

### 合并得到四个个rawcounts 文件，保存为RData（rawcounts)，dgmatrix格式。
setwd("../rawdata_scRNA")
for (j in 1:4){
	k = c(1,13,25,49)  ## 合并文件命名序列起始下标
	n = c(12,24,48,64)  ## 合并文件命名序列终止下标
	path = paste0("x",k[j])
	rawdata <- Read10X("x1")
	obj <- CreateSeuratObject(counts = rawdata,project = path,min.cells = 3) 
	obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
	obj <- subset(obj, subset = nFeature_RNA > 1000 & percent.mt<10)
	print(obj)
	mm1M = obj 
	### 根据设置的k，n 进行两次合并，注意y 文件下标需要 k+1
	for (i in (k[j]+1):n[j]){  ### 数学运算需要加括号，否则出现错误
		path = paste0("x",i)
		rawdata = Read10X(path)
		obj <- CreateSeuratObject(counts = rawdata,project = path,min.cells = 3) 
		obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
		obj <- subset(obj, subset = nFeature_RNA > 1000 & percent.mt<10)
		print(obj)
		### merge 原理是什么？类似按基因合并？
		mm.combined <- merge(x = mm1M, y = obj, project = "mouseBrain_1M")
		mm1M = mm.combined
	}
	print(mm1M)
	print(table(mm1M@meta.data$orig.ident))
	### 提取counts data
	rawcounts = mm1M[['RNA']]@counts

	fname = paste("../RData/mm1M_sub",k[j],n[j],"10X.RData",sep = "_")
	save(rawcounts,file = fname)
}


### 合并得到2个文件，（细胞数比例约为3:1）,按比例从这两个文件中抽取对应细胞数，100000（75000：25000）,200000（150000:30000）
library(Seurat)
library(tidyverse)
setwd("../RData")
for (i in 1:2){	
	sink("log_sample.txt")
	rate = c(0.75,0.25)
	ncells = c(100000,200000)
	fnames = c("10w","20w")
	prjnames = c("sub_1_48","sub_49_64")

	load("../RData/mm1M_sub_1_48_10X.RData")
	rawdt1 = rawcounts; rm(rawcounts)
	cat("dims of rawcounts 1 is")
	print(dim(rawdt1))
	load("../RData/mm1M_sub_49_64_10X.RData")
	rawdt2 = rawcounts; rm(rawcounts)
	cat("dims of rawcounts 2 is")
	print(dim(rawdt2))

	sub1 = rawdt1[,sample(ncol(rawdt1),(rate[1]*ncells[i]),replace = FALSE)]
	obj1 <- CreateSeuratobj1ect(counts = sub1,project = prjnames[1],min.cells = 3) 
	obj1[["percent.mt"]] <- PercentageFeatureSet(obj1, pattern = "^mt-")
	obj1 <- subset(obj1, subset = nFeature_RNA > 1000 & percent.mt<10)
	cat("sampling result info form data 1 is")
	print(obj1)

	sub2 = rawdt2[,sample(ncol(rawdt2),(rate[2]*ncells[i]),replace = FALSE)]
	obj2 <- CreateSeuratobj2ect(counts = sub2,project = prjnames[2],min.cells = 3) 
	obj2[["percent.mt"]] <- PercentageFeatureSet(obj2, pattern = "^mt-")
	obj2 <- subset(obj2, subset = nFeature_RNA > 1000 & percent.mt<10)
	cat("sampling result info form data 1 is")
	print(obj2)
	### merge 原理是什么？
	obj_integrated <- merge(x = obj1, y = obj2, project = fnames[i])

	fname = paste0("../RData/mm1M_sub",fnames[i],".RData")
	save(obj_integrated,file = fname)
	### 前三个数据集抽取70万细胞，再和第4个数据集合并，得到100万细胞
	rm(obj1,obj2,obj_integrated)
}
sink()

### 提取counts data
sink()
rawcounts = mm1M[['RNA']]@counts
fname = "../RData/mm1M_rawcounts_filtered.RData"
save(rawcounts,file = fname)

### 抽样,
load("../RData/mm1M_rawcounts_filtered.RData")
sink()
for (i in c(50000,100000,200000)){
	set.seed(1234)
	print(i)
	sub_counts <- rawcounts[,sample(ncol(rawcounts),i,replace = FALSE)]
	sub_path = paste0("mm1M_sub_",i/10000,'w')
	dir.create(file.path("../RData/",sub_path),showWarnings = FALSE)
	setwd(sub_path)
	c_num <- c(1000,3000,5000,10000,20000,30000,40000,50000)
	for (i in 1:8){
		set.seed(1234)
		f.name <- paste("obj_integrated",c_num[i],"cell.RData",sep = "_")
		cols = sample(1:ncol(sub_counts), c_num[i], replace = FALSE)
		sub_data <- sub_counts[,cols]
		print(dim(sub_data))
		save(sub_data,file =f.name )
	}
	setwd("../")
}
sink()

