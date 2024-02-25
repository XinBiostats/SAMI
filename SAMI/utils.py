import os
import scanpy as sc
import pandas as pd
import numpy as np

def list_files(directory, substring, extension):
    files = sorted([file for file in os.listdir(directory) if substring in file and file.endswith(extension)])
    return files

def normalization(adata):
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)
    return adata

def add_preva_abund(adata):
    abundance=np.sum(adata.X!=0,axis=1)/adata.shape[1]
    prevalence=np.sum(adata.X!=0,axis=0)/adata.shape[0]
    adata.obs['abundance']=abundance
    adata.var['prevalence']=prevalence
    return adata

def adata_filter(adata,abundance=0.1,prevalence=0.1):
    adata_temp=add_preva_abund(adata)
    adata_filter = adata_temp[adata_temp.obs['abundance']>=abundance,adata_temp.var['prevalence']>=prevalence]
    return adata_filter

def cluster_assign(adata,cluster):
    adata.obs=pd.merge(adata.obs,cluster,left_index=True,right_index=True)
    return adata

def adj_pvalue(sort):
    n=sort.shape[0]
    sort['rank']=sort.index+1
    sort['adj_pvalue']=sort['pvalue']*n/sort['rank']
    return sort