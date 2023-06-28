import os
import scanpy as sc
import matplotlib.pyplot as plt
import string
import pandas as pd
import numpy as np

def cluster(adataname,res=1):
    adata=sc.read(os.path.join('../datasets/',adataname))
    sc.pp.pca(adata)
    sc.pp.neighbors(adata,n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=res)
    
    res=str(round(res,1))
    name, ext = os.path.splitext(adataname)
    modelname = name+'_'+res+'.h5ad'
    file_path = '../results/clustering'
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    adata.write(os.path.join(file_path,modelname))
    
    
def plot_umap_cluster(modelname,size=50):
    adata = sc.read(os.path.join('../results/clustering',modelname))
    name,ext = os.path.splitext(modelname)
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=adata.obs[color].unique().tolist()
    n_cluster=str(len(colorlist))
    
    fig, axs = plt.subplots(1, 2, figsize=(30, 10))

    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist],
        ax=axs[0],
        show=False,
    )

    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color=color,
        title=name,
        size=size,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        ax=axs[1],
        show=False
    )
    
    name,ext = os.path.splitext(modelname)
    
    plt.savefig(os.path.join('../results/clustering',name+'_umap_'+n_cluster+'.png'))
    
def plot_cluster(modelname,size=50):
    adata = sc.read(os.path.join('../results/clustering',modelname))
    name,ext = os.path.splitext(modelname)
    color='leiden'
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=adata.obs[color].unique().tolist()
    n_cluster=str(len(colorlist))

    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color=color,
        title=name,
        size=size,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        show=False
    )
    
    name,ext = os.path.splitext(modelname)
    
    plt.savefig(os.path.join('../results/clustering',name+'_'+n_cluster+'.png'))
    
def plot_select_cluster(modelname,cluster,merge_cluster=False,size=50):
    adata = sc.read(os.path.join('../results/clustering',modelname))
    if merge_cluster == False:
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs['leiden'].unique().tolist()
        palette=[v if k in list(map(str,cluster)) else '#A0A0A0' for k, v in clusters_colors.items()]
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color='leiden',
            legend_loc='none',
            #title='brain2_pool'+'_'+str(1.4)+'_cluster'+str(cluster),
            title='',
            size=size,
            palette=palette,
            show=False
        )
    else:
        cluster_merge = sorted(adata.obs['merge'].unique().tolist())
        n_cluster = len(adata.obs['leiden'].unique().tolist())
        n_cluster_merge = len(adata.obs['merge'].unique().tolist())
        indices = [index for index, value in enumerate(cluster_merge) if value in cluster]
        palette = sc.pl.palettes.default_102[n_cluster+3:]
        palette=["#A0A0A0"]*n_cluster_merge
        for i in indices:
            palette[i]=sc.pl.palettes.default_102[n_cluster+3:][i]
        
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color='merge',
            legend_loc='none',
            #title='brain2_pool'+'_'+str(1.4)+'_cluster'+str(cluster),
            title='',
            size=size,
            palette=palette,
            show=False
        )
    plt.savefig(os.path.join('../results/clustering/','selected_cluster_'+str(cluster).replace("'","")+'.png')) 
    
    
def merge_cluster(modelname,merge_dict):
    adata = sc.read(os.path.join('../results/clustering/',modelname))
    adata.obs['merge']= adata.obs['leiden'].astype(int).map({value: key for key, values in merge_dict.items() for value in values}).fillna('other')
    adata.write(os.path.join('../results/clustering',modelname))
    
    
def cluster_mapping(test,ref,testadata,markerfile):
    markers_file = pd.ExcelFile(os.path.join('../results/markers/',markerfile))
    sheets_names = markers_file.sheet_names

    markers = pd.DataFrame()
    for sheet in sheets_names:
        sections, omics = sheet.split('_')
        temp = pd.read_excel(markers_file,sheet)
        temp[['omics']] = sheet.split('_')[1]
        temp[['sections']] = sections
        markers = pd.concat([markers,temp], ignore_index=True)

    cluster_num1 = len(np.unique(markers.loc[markers['sections']==test,'cluster']))
    cluster_num2 = len(np.unique(markers.loc[markers['sections']==ref,'cluster']))

    prop = np.zeros((cluster_num2,cluster_num1))
    intersect = np.empty((cluster_num2,cluster_num1),dtype=object)
    intersect.fill(set())
    for i in range(cluster_num2):
        set1 = set(markers.loc[(markers['cluster']==i)&(markers['sections']==ref),'feature'])
        for j in range(cluster_num1):
            set2 = set(markers.loc[(markers['cluster']==j)&(markers['sections']==test),'feature'])
            intersection = set1.intersection(set2)
            intersection_size = len(intersection)
            union_size = len(set1.union(set2))
            prop[i][j]=intersection_size
            intersect[i][j]=intersection

    highest_prop = np.argmax(prop, axis=0)
    mapping = []
    for index1, index2 in enumerate(highest_prop):
        mapping.append((index1,index2))

    # cluster1 = []
    # cluster2 = []
    # intersection = []

    # for x, y in mapping:
    #     cluster1.append(x)
    #     cluster2.append(y)
    #     intersection.append(intersect[y, x])

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))

    colormap = mapping
    for i, tup in enumerate(colormap):
        key = str(tup[1])
        if key in clusters_colors:
            value = clusters_colors[key]
            colormap[i]=(tup[0],value)

    new_palette = [tup[1] for tup in colormap]

    adata = sc.read(os.path.join('../results/clustering/',testadata))

    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color='leiden',
        title='',
        size=50,
        palette=new_palette,
        #legend_loc='none',
        show=False
        )

    plt.savefig(os.path.join('../results/clustering/',test+'_mapping_'+ref+'.png'))