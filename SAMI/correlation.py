import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from scipy.stats import gaussian_kde,pearsonr
import seaborn as sns
from utils import *

def calculate_corr(adata1,adata2):
    """
    Computes the correlation between two spatial omics datasets based on shared pixels.

    This function identifies common pixels and features between two datasets,
    filters out shared features, and calculates the correlation between remaining features.

    Parameters
    ----------
    adata1 : AnnData
        First AnnData object.
    adata2 : AnnData
        Second AnnData object.

    Returns
    -------
    pd.DataFrame
        Dataframe containing correlation values between features in `adata1` and `adata2`.

    Saves
    -----
    CSV file containing correlation values:
        `../results/correlation/correlation.csv`
    """
    common_pixel = set(adata1.obs_names).intersection(set(adata2.obs_names))
    common_feat = set(adata1.var_names).intersection(set(adata2.var_names))
    data1 = adata1[pd.Index(common_pixel), ~adata1.var_names.isin(common_feat)].copy().to_df()
    data2 = adata2[pd.Index(common_pixel), ~adata2.var_names.isin(common_feat)].copy().to_df()
    data = pd.merge(data1,data2,left_index=True,right_index=True)
    corr_df = data.corr()
    corr_df = corr_df.iloc[:data1.shape[1],data1.shape[1]:]
    corr = corr_df.stack().reset_index().rename(columns={'level_0':'omic1','level_1':'omic2',0:'corr'})
    file_path = '../results/correlation/'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    corr.to_csv(os.path.join(file_path,'correlation.csv'), index=False)
    return corr

def corr_plot(adata1,adata2,mz1,mz2,show=False):
    """
    Plots the correlation between two selected features from different omics datasets.

    This function generates a scatter plot with a kernel density estimation (KDE)
    to visualize the correlation between the selected features.

    Parameters
    ----------
    adata1 : AnnData
        First AnnData object.
    adata2 : AnnData
        Second AnnData object.
    mz1 : str
        Feature name from `adata1`.
    mz2 : str
        Feature name from `adata2`.
    show : bool, optional
        If True, displays the plot instead of saving. Default is False.

    Saves
    -----
    PNG file containing the correlation plot:
        `../results/correlation/{mz2}_vs_{mz1}.png`
    """
    data1 = pd.DataFrame(adata1.X, columns=adata1.var_names, index=adata1.obs_names)
    data2 = pd.DataFrame(adata2.X, columns=adata2.var_names, index=adata2.obs_names)

    plt.figure(figsize=(12,10))
    plt.subplots_adjust(left=0.24,bottom=0.15)
    x = data1[mz1].values
    y = data2[mz2].values
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    plt.scatter(x,y,c=z,s=1)
    sns.regplot(x=x,y=y, scatter = False, ci = 95, fit_reg = True, color = 'darkorange') 
    plt.xlabel(mz1,fontsize=45,fontweight='bold')
    plt.ylabel(mz2,fontsize=45,fontweight='bold')
    plt.locator_params(axis='y',nbins=2)
    plt.locator_params(axis='x',nbins=2)
    plt.xticks(fontsize=35,fontweight='bold')
    plt.yticks(fontsize=35,fontweight='bold')
    corr, p_value = pearsonr(x, y)
    plt.text(max(x),min(y),f'R={corr:.2f}',ha='right',va='bottom',fontsize=45,fontweight='bold')
    
    file_path = '../results/correlation/'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(os.path.join(file_path,f'{mz2}_vs_{mz1}.png'))
    if show == False:
        plt.close()
    else:
        plt.show()

def create_network(adata,abundance,prevalence):
    """
    Constructs a correlation matrix for intra- and inter-omics interactions.

    This function filters the input dataset based on abundance and prevalence thresholds,
    then computes the correlation matrix between the remaining features.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    abundance : float
        Minimum abundance threshold for feature selection.
    prevalence : float
        Minimum prevalence threshold for feature selection.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the correlation matrix.
    """
    data = adata_filter(adata,abundance=abundance,prevalence=prevalence).to_df()
    corr_matrix = data.corr() ## calculate correlation intra and inter omics
    corr_df = pd.DataFrame(data=corr_matrix, index=data.columns, columns=data.columns)
    return corr_df
    
    
def corr_network(adata, corr_df):
    """
    Constructs and visualizes a correlation network for spatial omics data.

    This function creates a network graph based on the correlation matrix,
    setting a threshold of 0.5 for edge inclusion. Nodes are colored based
    on their omics type.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing omics information.
    corr_df : pd.DataFrame
        Correlation matrix dataframe.

    Saves
    -----
    PNG file of the correlation network:
        `../results/correlation/triple_network.png`
    """
    # create a graph object
    fig, ax = plt.subplots(figsize=(10,10))
    G = nx.Graph()

    for i in range(corr_df.shape[0]):
        G.add_node(corr_df.columns[i],omics=adata.var.loc[corr_df.columns[i]][0])
        for j in range(i+1, corr_df.shape[1]):
            if corr_df.iloc[i,j] > 0.5: #set a correlation threshold here
                G.add_edge(corr_df.columns[i],corr_df.columns[j],weight=corr_df.iloc[i, j])

    node_degree = dict(G.degree())
    node_size = [(v+1)*3 for v in node_degree.values()]
    #label_dict = {node: node if G.degree(node) > 100 else '' for node in G.nodes()}
    node_color = ['green' if adata.var.loc[node][0]=='metabolomics' else 'blue' if adata.var.loc[node][0]=='lipidomics' else 'red' for node in G.nodes()]

    edge_color=[]
    for u,v in G.edges:
        if adata.var.loc[u][0]==adata.var.loc[v][0]:
            edge_color.append('lightgray')
        else:
            edge_color.append('yellow')

    pos = nx.spring_layout(G,k=1,pos=None)
    # for node in pos:
    #     if adata.var.loc[node][0]=='metabolomics':
    #         pos[node][0]=-1*abs(pos[node][0])
    #         pos[node][1]=abs(pos[node][1])
    #     elif adata.var.loc[node][0]=='lipidomics':
    #         pos[node][0]=abs(pos[node][0])
    #         pos[node][1]=abs(pos[node][1])
    #     elif adata.var.loc[node][0]=='glycomics':
    #         pos[node][0]=pos[node][0]
    #         pos[node][1]=-1*abs(pos[node][1])

    # plot the graph
    #plt.title('Triple Domains Correlation Network',fontsize=15,fontweight='bold')
    nx.draw(G, pos,node_size=node_size,node_color=node_color,edge_color=edge_color,width=0.1)

    node_color_dict = {'metabolomics':'green','lipidomics':'blue','glycomics':'red'}
    edge_color_dict = {'intradomain':'lightgray','interdomain':'yellow'}
    node_legend = [plt.Line2D([], [], marker='o', markersize=10, color=color, linestyle='None',label=domain) for domain, color in node_color_dict.items()]
    edge_legend = [plt.Line2D([], [], color=color, linestyle='-',label=type) for type,color in edge_color_dict.items()]
    handles = node_legend+edge_legend
    plt.legend(handles=handles)
    
    file_path = '../results/correlation/'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(file_path+'triple_network.png')
    plt.show()