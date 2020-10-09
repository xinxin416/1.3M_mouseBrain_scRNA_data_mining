library(SingleR)
library(scater)
library(reticulate)
np <- import("numpy")
library(SingleR)
library(scater)
library(Seurat)
library(tidyverse)

setwd("/home/chenhy/mouse_brain/")
wdir = "/home/chenhy/mouse_brain/"

setwd("./RData/Chen")
counts <- np$load("expr_rawcounts.npy") %>% t()
#counts <- t(counts)
genes <- read.table("features.tsv",sep = "\t",header = TRUE);names(genes) = "geneID"
cells <- read.table("barcodes.tsv",sep = "\t",header = TRUE);names(cells) = "barcodes"
cols <- read.table("cell_ontology_labels.tsv",sep = "\t",header = TRUE)
cat("dimention of counts is: ")
print(dim(counts))
colnames(counts) = cells$barcodes
rownames(counts) = genes$geneID

sce_ref_Chen <- SingleCellExperiment(list(counts=counts),
    colData=cols
)
sce_ref_Chen <- logNormCounts(sce_ref_Chen)  ### 数据较大，计算时间更长
sce_ref_Chen
save(sce_ref_Chen, file = "ref_sce.RData")  ### 占用内存更小
