import os
import scanpy as sc
import matplotlib.pyplot as plt
import string
import pandas as pd
import numpy as np

class Clusters:
    """
    Performs clustering analysis on spatial omics data using Scanpy.

    Attributes
    ----------
    data_path : str
        Path to the directory containing input .h5ad files.
    result_path : str
        Path to the directory where results will be saved.
    region : str
        Name of the region being analyzed.
    modality : str
        Type of omics data (e.g., pool, metabolomics, lipidomics, glycomics).
    resolution : float
        Resolution parameter for Leiden clustering.
    """
    def __init__(self, data_path, result_path, region, modality, resolution):
        self.data_path = data_path
        self.result_path = result_path
        self.region = region
        self.modality = modality
        self.resolution = round(resolution,1)

    def clustering(self):
        """
        Performs PCA, Harmony integration (if multiple regions exist), and Leiden clustering.
        Saves the clustered dataset as an .h5ad file.
        """
        adata=sc.read(os.path.join(self.data_path,f'{self.region}_{self.modality}.h5ad'))
        sc.pp.pca(adata)
        if len(np.unique(adata.obs['region']))>1:
            sc.external.pp.harmony_integrate(adata,key='region',max_iter_harmony=20)
            all_pcs = adata.obsm['X_pca_harmony'].shape[1]
            n_pcs = min(20, all_pcs)
            sc.pp.neighbors(adata,use_rep = 'X_pca_harmony', n_pcs=n_pcs)
        else:
            all_pcs = adata.obsm['X_pca'].shape[1]
            n_pcs = min(20, all_pcs)
            sc.pp.neighbors(adata,n_pcs=n_pcs)
            
        sc.tl.umap(adata)
        sc.tl.leiden(adata,resolution=self.resolution)

        file_path = os.path.join(self.result_path,f'clustering/{self.modality}')
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        adata.write(os.path.join(file_path,f'{self.region}_{self.resolution}.h5ad'))
    
    
    def plot_umap_cluster(self,size=50,show=False):
        """
        Plots the UMAP clustering and spatial clustering.

        Parameters
        ----------
        size : int, optional
            Size of the points in the spatial plot. Default is 50.
        show : bool, optional
            If True, displays the plot instead of saving. Default is False.
        """
        adata = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}.h5ad'))
        color='leiden'
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs[color].cat.categories.tolist()
        n_cluster=str(len(colorlist))

        fig, axs = plt.subplots(1, 2, figsize=(30, 10))

        sc.pl.umap(
            adata, 
            legend_loc='on data',
            color=[color], 
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            ax=axs[0],
            show=False,
        )

        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color=color,
            title=f'{self.region}_{self.modality}_{self.resolution}',
            size=size,
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            ax=axs[1],
            show=False
        )
        
            
        samples = np.unique(adata.obs['region']).tolist()
        for sample in samples:
            adata_sub = adata[adata.obs['region']==sample]
            x_label = np.min(adata_sub.obsm['spatial'][:,0])
            y_label = np.min(adata_sub.obsm['spatial'][:,1])
            axs[1].text(x_label,y_label, sample,fontsize=10, color='black')
    
        
        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}_umap_{n_cluster}.png'),dpi=400)
        
        if show==False:
            plt.close()
        else:
            plt.show()
    
    
    def plot_cluster(self,size=50,show=False):
        """
        Plots the spatial clustering of the dataset only.

        This function generates a spatial plot of the Leiden clustering results, where each cluster is color-coded.

        Parameters
        ----------
        size : int, optional
            Size of the points in the spatial plot. Default is 50.
        show : bool, optional
            If True, displays the plot instead of saving it. Default is False.

        Saves
        -----
        A spatial clustering plot as a PNG file:
            {result_path}/clustering/{modality}/{region}_{resolution}_{num_clusters}.png
        """
        adata = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}.h5ad'))
        color='leiden'
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs[color].cat.categories.tolist()
        n_cluster=str(len(colorlist))

        fig, axs = plt.subplots(1, 1, figsize=(15, 10))
        
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color=color,
            title=f'{self.region}_{self.modality}_{self.resolution}',
            size=size,
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            ax=axs,
            show=False
        )
        
        samples = np.unique(adata.obs['region']).tolist()
        for sample in samples:
            adata_sub = adata[adata.obs['region']==sample]
            x_label = np.min(adata_sub.obsm['spatial'][:,0])
            y_label = np.min(adata_sub.obsm['spatial'][:,1])
            axs.text(x_label,y_label, sample,fontsize=10, color='black')

        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}_{n_cluster}.png'),dpi=400)
        if show==False:
            plt.close()
        else:
            plt.show()
    
    def plot_select_cluster(self,cluster,size=50,show=False):
        """
        Highlights selected clusters in the spatial plot while graying out others.

        This function allows visualization of specific clusters by maintaining their original colors while displaying
        all other clusters in gray. It also annotates sample region names on the spatial plot.

        Parameters
        ----------
        cluster : int or list of int
            The cluster(s) to highlight in the plot. Other clusters will be grayed out.
        size : int, optional
            Size of the points in the spatial plot. Default is 50.
        show : bool, optional
            If True, displays the plot instead of saving it. Default is False.

        Saves
        -----
        A spatial clustering plot with selected clusters highlighted:
            {result_path}/clustering/{modality}/{region}_{resolution}_{selected_clusters}.png
        """
        adata = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}.h5ad'))
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        
        if isinstance(cluster, int):
            cluster = [cluster]
            
        colorlist=adata.obs['leiden'].cat.categories.tolist()
        palette=[v if k in map(str,cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
        
        fig, axs = plt.subplots(1, 1, figsize=(15, 10))
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color='leiden',
            legend_loc='none',
            title=f'{self.region}_{self.modality}_{self.resolution}_{cluster}',
            size=size,
            palette=palette,
            ax=axs,
            show=False
        )
        
        samples = np.unique(adata.obs['region']).tolist()
        for sample in samples:
            adata_sub = adata[adata.obs['region']==sample]
            x_label = np.min(adata_sub.obsm['spatial'][:,0])
            y_label = np.min(adata_sub.obsm['spatial'][:,1])
            axs.text(x_label,y_label, sample,fontsize=10, color='black')
            
        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region}_{self.resolution}_{cluster}.png'),dpi=400)
        if show==False:
            plt.close()
        else:
            plt.show()
    
class Cluster_Integration:
    """
    Integrates and visualizes clustered spatial omics datasets.

    This class performs integration of two datasets, applies Leiden clustering, 
    and generates various visualizations, including UMAP and spatial plots.

    Attributes
    ----------
    data_path : str
        Path to the directory containing input .h5ad files.
    result_path : str
        Path to the directory where results will be saved.
    adata : str
        Filename of the first dataset to integrate.
    adata_ref : str
        Filename of the reference dataset for integration.
    modality : str
        Type of omics data (e.g., metabolomics, lipidomics, glycomics).
    """
    def __init__(self,data_path, result_path, adata,adata_ref,modality):
        self.data_path = data_path
        self.result_path = result_path
        self.adata = adata
        self.adata_ref = adata_ref
        self.modality = modality
        self.region1, self.resolution1 =self.splitname(adata)
        self.region2, self.resolution2 =self.splitname(adata_ref)
        
    def integrate(self):
        """
        Integrates two spatial omics datasets using the ingest function in Scanpy.

        This function reads two datasets, assigns sample labels, 
        and applies Scanpy's `ingest()` function to integrate the datasets 
        based on Leiden clustering.

        Saves
        -----
        Integrated datasets as `.h5ad` files:
            {result_path}/clustering/{modality}/{region1}_integrated.h5ad
            {result_path}/clustering/{modality}/{region2}_integrated.h5ad
        """        
        adata1 = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.adata}'))
        adata2 = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.adata_ref}'))
        var_names = adata1.var_names.intersection(adata2.var_names)
        adata1.obs['sample'] = self.region1
        adata2.obs['sample'] = self.region2
        sc.tl.ingest(adata1,adata2,obs='leiden')
        
        adata1.write(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region1}_integrated.h5ad'))
        adata2.write(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region2}_integrated.h5ad'))
        
    def splitname(self,adata):
        """
        Extracts the region name and resolution from the .h5ad filename.

        Parameters
        ----------
        adata : str
            Filename of the .h5ad dataset.

        Returns
        -------
        tuple (str, str)
            Region name and resolution.
        """
        name,_ = os.path.splitext(adata)
        region, resolution = name.split('_')
        return region, resolution
    
    def plot_umap_cluster(self,region,size=50,show=False):
        """
        Plots the UMAP clustering and spatial clustering of an integrated dataset.

        Parameters
        ----------
        region : str
            Region name to visualize.
        size : int, optional
            Size of the points in the spatial plot. Default is 50.
        show : bool, optional
            If True, displays the plot instead of saving. Default is False.

        Saves
        -----
        UMAP and spatial clustering plot:
            {result_path}/clustering/{modality}/{region}_umap_integrated.png
        """        
        adata = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{region}_integrated.h5ad'))
        color='leiden'
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs[color].cat.categories.tolist()
        n_cluster=str(len(colorlist))

        fig, axs = plt.subplots(1, 2, figsize=(30, 10))

        sc.pl.umap(
            adata, 
            legend_loc='on data',
            color=[color], 
            title=region,
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            ax=axs[0],
            show=False,
        )

        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color=color,
            title=region,
            size=size,
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            ax=axs[1],
            show=False
        )
        
        samples = np.unique(adata.obs['region']).tolist()
        for sample in samples:
            adata_sub = adata[adata.obs['region']==sample]
            x_label = np.min(adata_sub.obsm['spatial'][:,0])
            y_label = np.min(adata_sub.obsm['spatial'][:,1])
            axs[1].text(x_label,y_label, sample,fontsize=10, color='black')

        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{region}_umap_integrated.png'),dpi=400)

        if show==False:
            plt.close()
        else:
            plt.show()
            
    def plot_overlap_umap(self,show=False):
        """
        Plots an overlap of two integrated datasets in UMAP space.

        This function concatenates two datasets and visualizes their overlap.

        Parameters
        ----------
        show : bool, optional
            If True, displays the plot instead of saving. Default is False.

        Saves
        -----
        Overlapping UMAP plot:
            {result_path}/clustering/{modality}/{region1}_{region2}_umap_overlap.png
        """
        adata1 = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region1}_integrated.h5ad'))
        adata2 = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region2}_integrated.h5ad'))
        adata_concat = adata2.concatenate(adata1, batch_categories=[self.region2,self.region1])
        fig,ax = plt.subplots(figsize=(10,8))
        plt.subplots_adjust(left=0.1,right=0.7)
        sc.pl.umap(adata_concat,color='batch',title=f'{self.region1} vs {self.region2}',size=5,legend_fontsize=20,ax=ax,show=False)
        
        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{self.region1}_{self.region2}_umap_overlap.png'),dpi=200)
        if show==False:
            plt.close()
        else:
            plt.show()
            
            
    def plot_select_cluster(self,region,cluster,size=50,show=False):
        adata = sc.read(os.path.join(self.result_path,f'clustering/{self.modality}/{region}_integrated.h5ad'))
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        
        if isinstance(cluster, int):
            cluster = [cluster]
            
        colorlist=adata.obs['leiden'].cat.categories.tolist()
        palette=[v if k in map(str,cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
        
        fig, axs = plt.subplots(1, 1, figsize=(15, 10))
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color='leiden',
            legend_loc='none',
            title=f'{region}_{self.modality}_integrated_{cluster}',
            size=size,
            palette=palette,
            ax=axs,
            show=False
        )
        
        samples = np.unique(adata.obs['region']).tolist()
        for sample in samples:
            adata_sub = adata[adata.obs['region']==sample]
            x_label = np.min(adata_sub.obsm['spatial'][:,0])
            y_label = np.min(adata_sub.obsm['spatial'][:,1])
            axs.text(x_label,y_label, sample,fontsize=10, color='black')
            
        plt.savefig(os.path.join(self.result_path,f'clustering/{self.modality}/{region}_{cluster}.png'),dpi=400) 
        if show==False:
            plt.close()
        else:
            plt.show()
