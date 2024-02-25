import os
import scanpy as sc
import matplotlib.pyplot as plt
import string
import pandas as pd
import numpy as np

class Clusters:
    def __init__(self, region, modality, resolution):
        self.region = region
        self.modality = modality
        self.resolution = round(resolution,1)

    def clustering(self):
        adata=sc.read(os.path.join('../datasets/',f'{self.region}_{self.modality}.h5ad'))
        sc.pp.pca(adata)
        sc.pp.neighbors(adata,n_pcs=20)
        sc.tl.umap(adata)
        sc.tl.leiden(adata,resolution=self.resolution)

        file_path = '../results/clustering'
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        adata.write(os.path.join(file_path,f'{self.region}_{self.modality}_{self.resolution}.h5ad'))
    
    
    def plot_umap_cluster(self,size=50,show=False):
        adata = sc.read(os.path.join('../results/clustering',f'{self.region}_{self.modality}_{self.resolution}.h5ad'))
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
        
        plt.savefig(os.path.join('../results/clustering',f'{self.region}_{self.modality}_{self.resolution}_umap_{n_cluster}.png'),dpi=400)
        
        if show==False:
            plt.close()
        else:
            plt.show()
    
    
    def plot_cluster(self,size=50,show=False):
        adata = sc.read(os.path.join('../results/clustering',f'{self.region}_{self.modality}_{self.resolution}.h5ad'))
        color='leiden'
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs[color].cat.categories.tolist()
        n_cluster=str(len(colorlist))

        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color=color,
            title=f'{self.region}_{self.modality}_{self.resolution}',
            size=size,
            palette=[v for k, v in clusters_colors.items() if k in colorlist],
            show=False
        )

        plt.savefig(os.path.join('../results/clustering',f'{self.region}_{self.modality}_{self.resolution}_{n_cluster}.png'),dpi=400)
        if show==False:
            plt.close()
        else:
            plt.show()
    
    def plot_select_cluster(self,cluster,size=50,show=False):
        adata = sc.read(os.path.join('../results/clustering',f'{self.region}_{self.modality}_{self.resolution}.h5ad'))
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs['leiden'].cat.categories.tolist()
        palette=[v if k==str(cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
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
            show=False
        )
        plt.savefig(os.path.join('../results/clustering/',f'{self.region}_{self.modality}_{self.resolution}_{cluster}.png'),dpi=400) 
        if show==False:
            plt.close()
        else:
            plt.show()
    
class Cluster_Integration:
    def __init__(self,adata,adata_ref):
        self.adata = adata
        self.adata_ref = adata_ref
        self.region1, self.modality1, self.resolution1 =self.splitname(adata)
        self.region2, self.modality2, self.resolution2 =self.splitname(adata_ref)
        
    def integrate(self):
        adata1 = sc.read(f'../results/clustering/{self.adata}')
        adata2 = sc.read(f'../results/clustering/{self.adata_ref}')
        var_names = adata1.var_names.intersection(adata2.var_names)
        adata1.obs['sample'] = self.region1
        adata2.obs['sample'] = self.region2
        sc.tl.ingest(adata1,adata2,obs='leiden')
        
        adata1.write(f'../results/clustering/{self.region1}_integrated.h5ad')
        adata2.write(f'../results/clustering/{self.region2}_integrated.h5ad')
        
    def splitname(self,adata):
        name,_ = os.path.splitext(adata)
        region, modality, resolution = name.split('_')
        return region, modality, resolution
    
    def plot_umap_cluster(self,region,size=50,show=False):
        adata = sc.read(os.path.join('../results/clustering',f'{region}_integrated.h5ad'))
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

        plt.savefig(os.path.join('../results/clustering',f'{region}_umap_integrated.png'),dpi=400)

        if show==False:
            plt.close()
        else:
            plt.show()
            
    def plot_overlap_umap(self,show=False):
        adata1 = sc.read(os.path.join('../results/clustering',f'{self.region1}_integrated.h5ad'))
        adata2 = sc.read(os.path.join('../results/clustering',f'{self.region2}_integrated.h5ad'))
        adata_concat = adata2.concatenate(adata1, batch_categories=[self.region2,self.region1])
        fig,ax = plt.subplots(figsize=(10,8))
        plt.subplots_adjust(left=0.1,right=0.7)
        sc.pl.umap(adata_concat,color='batch',title=f'{self.region1} vs {self.region2}',size=5,legend_fontsize=20,ax=ax,show=False)
        
        plt.savefig(f'../results/clustering/{self.region1}_{self.region2}_umap_overlap.png',dpi=200)
        if show==False:
            plt.close()
        else:
            plt.show()
            
            
    def plot_select_cluster(self,region,cluster,size=50,show=False):
        adata = sc.read(os.path.join('../results/clustering',f'{region}_integrated.h5ad'))
        clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
        colorlist=adata.obs['leiden'].cat.categories.tolist()
        palette=[v if k==str(cluster) else '#A0A0A0' for k, v in clusters_colors.items()]
        sc.pl.spatial(
            adata,
            img_key=None,
            library_id=None,
            spot_size=1, #important to add to avoid error
            color='leiden',
            legend_loc='none',
            title=f'{region}_integrated_{cluster}',
            size=size,
            palette=palette,
            show=False
        )
        plt.savefig(os.path.join('../results/clustering/',f'{region}_{cluster}.png'),dpi=400) 
        if show==False:
            plt.close()
        else:
            plt.show()