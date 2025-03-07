import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def cluster_mapping(reference,query,reference_data,query_data,size=50,show=False):
    """
    Maps clusters between a reference dataset and a query dataset using marker genes.

    This function calculates the Jaccard Index between clusters from two datasets based on
    shared marker genes. It assigns the most similar reference cluster to each query cluster.
    The function then visualizes the spatial clustering results for both datasets, using a
    color-matching scheme to indicate the mapping.

    Parameters
    ----------
    reference : str
        Name of the reference dataset (used for reading marker genes and clustering results).
    query : str
        Name of the query dataset (used for reading marker genes and clustering results).
    reference_data : str
        Filename of the reference dataset `.h5ad` file.
    query_data : str
        Filename of the query dataset `.h5ad` file.
    size : int, optional
        Size of the points in the spatial plot. Default is 50.
    show : bool, optional
        If True, displays the plot instead of saving. Default is False.

    Saves
    -----
    1. Cluster mapping results as a CSV file:
        `../results/clustermapping/{query}_map_{reference}.csv`
        - Columns: `query_cluster`, `reference_cluster`, `Jaccard_Index`, `Overlap_Markers`
    2. Cluster mapping visualization as a PNG file:
        `../results/clustermapping/{query}_map_{reference}.png`
        - Left: Reference spatial clustering.
        - Right: Query dataset colored based on reference mapping.

    """
    markers1 = pd.read_csv(os.path.join('../results/markers/',f'{reference}_marker.csv')) #ref
    markers2 = pd.read_csv(os.path.join('../results/markers/',f'{query}_marker.csv')) #query
    cluster_num1 = len(markers1['cluster'].unique().tolist())
    cluster_num2 = len(markers2['cluster'].unique().tolist())
    prop = np.zeros((cluster_num1,cluster_num2))
    intersect = np.empty((cluster_num1,cluster_num2),dtype=object)
    intersect.fill(set())
    for i in range(cluster_num1):
        set1 = set(markers1.loc[(markers1['cluster']==i),'feature'])
        for j in range(cluster_num2):
            set2 = set(markers2.loc[(markers2['cluster']==j),'feature'])
            intersection = set1.intersection(set2)
            intersection_size = len(intersection)
            union_size = len(set1.union(set2))
            prop[i][j]=intersection_size/union_size
            intersect[i][j]=intersection
    highest_prop = np.argmax(prop, axis=0)

    mapping = []
    mapping_df = pd.DataFrame()
    for que, ref in enumerate(highest_prop):
        mapping_df = mapping_df.append({query:que,reference:ref,'Jaccard_Index':prop[ref,que],'Overlap_Markers':intersect[ref,que]},ignore_index=True)
        mapping.append((que,ref))
        
        
    adata1 = sc.read(os.path.join('../results/clustering/',reference_data))
    adata2 = sc.read(os.path.join('../results/clustering/',query_data))

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist = adata1.obs['leiden'].unique().tolist()

    colormap = mapping.copy()
    for i, tup in enumerate(colormap):
        key = str(tup[1])
        if key in clusters_colors:
            value = clusters_colors[key]
            colormap[i]=(tup[0],value)

    new_palette = [tup[1] for tup in colormap]
    
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    
    sc.pl.spatial(
        adata1,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color='leiden',
        title=reference,
        size=size,
        palette=[v for k, v in clusters_colors.items() if k in colorlist],
        ax=axs[0],
        show=False
        )
    
    sc.pl.spatial(
        adata2,
        img_key=None,
        library_id=None,
        spot_size=1, #important to add to avoid error
        color='leiden',
        title=f'Map {query} to {reference}',
        size=size,
        palette=new_palette,
        ax=axs[1],
        show=False
        )
    
    file_path = '../results/clustermapping/'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    plt.savefig(os.path.join(file_path,f'{query}_map_{reference}.png'),dpi=400)
    
    if show==False:
        plt.close()
    else:
        plt.show()
        
    mapping_df.to_csv(os.path.join(file_path,f'{query}_map_{reference}.csv'),index=False)
