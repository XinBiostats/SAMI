import os
import scanpy as sc
import pandas as pd
import numpy as np

def list_files(directory, substring, extension):
    """
    Lists files in a given directory that contain a specific substring and have a certain file extension.

    Parameters
    ----------
    directory : str
        Path to the directory where files are located.
    substring : str
        Substring that must be present in the filenames.
    extension : str
        File extension to filter files (e.g., '.csv', '.h5ad').

    Returns
    -------
    list
        A sorted list of filenames that match the specified criteria.
    """
    files = sorted([file for file in os.listdir(directory) if substring in file and file.endswith(extension)])
    return files

def normalization(adata):
    """
    Normalizes an AnnData object by total counts and applies log transformation.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing gene/metabolite expression data.

    Returns
    -------
    AnnData
        Normalized AnnData object with log-transformed values (log base 2).
    """
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)
    return adata

def add_preva_abund(adata):
    """
    Computes abundance and prevalence for features in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.

    Returns
    -------
    AnnData
        Modified AnnData object with two new attributes:
        - 'abundance' (obs-level): Fraction of nonzero values per sample.
        - 'prevalence' (var-level): Fraction of nonzero values per feature.
    """
    abundance=np.sum(adata.X!=0,axis=1)/adata.shape[1]
    prevalence=np.sum(adata.X!=0,axis=0)/adata.shape[0]
    adata.obs['abundance']=abundance
    adata.var['prevalence']=prevalence
    return adata

def adata_filter(adata,abundance=0.1,prevalence=0.1):
    """
    Filters an AnnData object based on abundance and prevalence thresholds.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    abundance : float, optional
        Minimum fraction of features that must be nonzero in a sample. Default is 0.1.
    prevalence : float, optional
        Minimum fraction of samples in which a feature must be detected. Default is 0.1.

    Returns
    -------
    AnnData
        Filtered AnnData object containing only features and samples that meet the thresholds.
    """
    adata_temp=add_preva_abund(adata)
    adata_filter = adata_temp[adata_temp.obs['abundance']>=abundance,adata_temp.var['prevalence']>=prevalence]
    return adata_filter

def cluster_assign(adata,cluster):
    """
    Assigns cluster labels to an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    cluster : pd.DataFrame
        DataFrame containing cluster assignments with index matching `adata.obs`.

    Returns
    -------
    AnnData
        AnnData object with updated cluster labels in `adata.obs`.
    """
    adata.obs=pd.merge(adata.obs,cluster,left_index=True,right_index=True)
    return adata

def adj_pvalue(sort):
    """
    Computes adjusted p-values using the Benjamini-Hochberg method for multiple testing correction.

    Parameters
    ----------
    sort : pd.DataFrame
        DataFrame containing a 'pvalue' column.

    Returns
    -------
    pd.DataFrame
        DataFrame with two additional columns:
        - 'rank': Rank of the p-values in ascending order.
        - 'adj_pvalue': Adjusted p-values using the Benjamini-Hochberg correction.
    """
    n=sort.shape[0]
    sort['rank']=sort.index+1
    sort['adj_pvalue']=sort['pvalue']*n/sort['rank']
    return sort

def adata_concat(*adata_list):
    """
    Concatenates multiple AnnData objects based on their common features.

    Parameters
    ----------
    *adata_list : tuple of AnnData
        One or more AnnData objects to be concatenated.

    Returns
    -------
    AnnData
        Concatenated AnnData object, including only common features across all input datasets.
    """

    common_vars = adata_list[0].var_names
    for adata in adata_list[1:]:
        common_vars = common_vars.intersection(adata.var_names)
    
    adata_sub_list = [adata[:, list(common_vars)] for adata in adata_list]
    
    adata_combined = sc.concat(adata_sub_list, join='inner', label='sample')

    adata_combined.var = adata_sub_list[0].var
    
    return adata_combined
