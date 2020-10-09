import h5py
import numpy as np 
import pandas as pd 
import scipy.sparse as sp_sparse
from scipy import io,sparse 
from pathlib import Path
import os
from pathlib import Path
import logging
import sys

log = open("myprog.log", "a")
sys.stdout = log

f2 = h5py.File("../rawdata_scRNA/ACA_datas/Chen.h5",'r')
os.chdir("../RData/")
path = "Chen"
Path(path).mkdir(parents=True, exist_ok=True) 

[print(f2[item]) for item in f2.keys()]

os.chdir(path)
### 保存barcodes，features 信息
cells = f2['obs_names'][...]
tmp = pd.DataFrame(cells)
tmp.to_csv("barcodes.tsv",sep = "\t",header = True,index = False)

genes = f2['var_names'][...]
tmp = pd.DataFrame(genes)
tmp.to_csv("features.tsv",sep = "\t",header = True,index = False)

### 保存注释label 信息
print(list(f2['/obs'].keys()))
cell_ontology_label = f2['/obs']['cell_ontology_class'][...]
cell_type1 = f2['/obs']['cell_type1'][...]
cols = pd.DataFrame({"cell_ontology_class":cell_ontology_label,"cell_type1":cell_type1})
cols.to_csv("cell_ontology_labels.tsv",sep = "\t",header = True,index = False)
### sed -i "s/b//g; s/'//g" */*.tsv 

### 保存matrix 信息
indices = f2['/exprs']['indices']  ## 0下标的行对应data 元素
# (6682428,)
indptr = f2['/exprs']['indptr']  ## 0下标的列对应data 元素
# (3661,)
shape = f2['/exprs']['shape']
# array([ 3660, 23797])
mtx = f2['/exprs']['data'][:]  ###Nonzero UMI counts in column-major order
## rawcounts
matrix = sp_sparse.csr_matrix((mtx, indices, indptr), shape=shape).toarray()
#np.savetxt('expr_rawcounts.csv', matrix, delimiter=",")  ### 耗时长，文件很大
np.save("expr_rawcounts.npy",matrix) ## 速度更快，占用内存更小




