import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
os.environ['R_HOME'] = "/blue/li.chen1/xin.ma/conda/envs/maldi/lib/R"
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

class Pathway:
    def __init__(self,region,modality):
        if not os.path.exists('../results/pathway'):
            os.makedirs('../results/pathway')
            
        self.region = region
        if modality=='metabolomics' or modality=='glycomics':
            print('Metabolomics and Glycomics are combined for pathway enrichment analysis.')
            self.modality = 'metabolomics_glycomics'
        else:
            self.modality = 'lipidomics'
            
    def findpathway(self):
        r_source = ro.r['source']
        r_source('../SAMI/MetabAnalystR.R')  
        
        ann_dict = self.anndict()
        markers = pd.read_csv(os.path.join('../results/markers/',f'{self.region}_marker.csv'))
        clusters = markers['cluster'].unique().tolist()
        
        if self.modality=='metabolomics_glycomics':
            markers = markers.loc[markers['feature'].isin(ann_dict['metabolomics']+ann_dict['glycomics'])]
        else:
            markers = markers.loc[markers['feature'].isin(ann_dict['lipidomics'])]
        
        logical_arg = ro.BoolVector([True,False])
        clustersR = ro.vectors.IntVector(clusters)
        regionR = ro.StrVector([self.region])
        with localconverter(ro.default_converter + pandas2ri.converter):
            markersR = ro.conversion.py2rpy(markers)
            
        PathwayEnrichment = ro.r['PathwayEnrichment']
        
        if self.modality=='metabolomics_glycomics':
            PathwayEnrichment(region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[1],merge=logical_arg[1])
        else:
            PathwayEnrichment(region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[0],merge=logical_arg[1])
        
    def anndict(self):
        anndict = {}
        for file in os.listdir('../annotation'):
            if re.match(r'^\w+annotation\.csv',file):
                modality = file.split('_')[0]
                data_temp = pd.read_csv(os.path.join('../annotation', file),delimiter=',')
                feat_cols = data_temp['Name'].tolist()
                anndict[modality] = feat_cols
        return anndict
    

    def plot_dot(self,cluster,scale,height,top=None,show=False):
        ora_data = pd.read_csv(f'../results/pathway/ora_{self.modality}.csv')
        ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
        
        if top is None:
            data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
        else:
            data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

        if data.shape[0]!=0:

            fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(10,10))
            fig.subplots_adjust(left=0.45,wspace=0.35)
            #### p-value  #####   
            pathway = data['pathway'].tolist()
            log_trans = -np.log10(data['Raw.p'])

            #axis
            y_pos = np.arange(len(pathway))*5
            x_pos = np.arange(0,np.ceil(max(log_trans))*1.1,5 if max(log_trans)>5 else 1)
            ax[0].set_xticks(x_pos)
            ax[0].tick_params(axis='x', labelsize=15)
            ax[0].set_xlim(-np.ceil(max(log_trans))*0.05,np.ceil(max(log_trans))*1.05)
            ax[0].set_yticks(y_pos,pathway,fontsize=18,fontweight='bold')
            ax[0].invert_yaxis()

            size = data['enrichment_ratio']*scale
            x = log_trans.tolist()
            color = data['Raw.p'].tolist()

            scatter = ax[0].scatter(x,y_pos,s=size,c=color,cmap='PuRd_r',zorder=2)
            ax[0].grid(alpha=0.5,zorder=1)
            ax[0].set_xlabel('-log10(p_value)',fontsize=18,fontweight='bold')

            cbar = fig.colorbar(scatter,shrink=0.5)
            cbar.ax.tick_params(labelsize=15)
            cbar.ax.set_title('p-value',fontsize=15,fontweight='bold')
            loc = cbar.ax.get_position()
            new_loc =[loc.x0+0.02,loc.y0+0.15,loc.width,loc.height]
            cbar.ax.set_position(new_loc)

            fig.set_size_inches(15, 8)

            #### FDR  #####   
            pathway = data['pathway'].tolist()
            log_trans = -np.log10(data['FDR'])

            #axis
            y_pos = np.arange(len(pathway))
            x_pos = np.arange(0,max(0.001,np.ceil(max(log_trans))*1.1),5 if max(log_trans)>5 else 1)
            ax[1].set_xticks(x_pos)
            ax[1].tick_params(axis='x', labelsize=15)
            ax[1].set_xlim(-max(np.ceil(max(log_trans))*0.05,0.05),max(1,np.ceil(max(log_trans))*1.1))
            ax[1].set_yticks(y_pos,'',fontsize=10)
            ax[1].invert_yaxis()

            size = data['enrichment_ratio']*scale
            x = log_trans.tolist()
            color = data['FDR'].tolist()

            scatter = ax[1].scatter(x,y_pos,s=size,c=color,cmap='PuRd_r',zorder=2)
            ax[1].grid(alpha=0.5,zorder=1)
            ax[1].set_xlabel('-log10(FDR)',fontsize=18,fontweight='bold')

            handles, labels = scatter.legend_elements(prop="sizes", alpha=1,num=4,func=lambda x: x/scale)
            legend = ax[1].legend(handles, labels, bbox_to_anchor=(1.8,0.43), ncol=1, title="Enrichment\nRatio",fontsize=15,labelspacing=1,borderpad=0.5)
            legend.get_title().set_fontsize(15)
            legend.get_title().set_weight('bold')

            cbar = fig.colorbar(scatter,shrink=0.5)
            scatter.set_clim(vmax=1)
            cbar.ax.tick_params(labelsize=15)
            cbar.ax.set_title('FDR',fontsize=15,fontweight='bold')
            loc = cbar.ax.get_position()
            new_loc =[loc.x0+0.01,loc.y0+0.15,loc.width,loc.height]
            cbar.ax.set_position(new_loc)
            fig.set_size_inches(20,height)

            plt.savefig(os.path.join('../results/pathway/',f'{self.region}_{self.modality}_{cluster}.png'))
            if show == False:
                plt.close()
            else:
                plt.show()

    
    
    def pathway_network(self,cluster,top=None,show=False):
        ora_data = pd.read_csv(f'../results/pathway/ora_{self.modality}.csv')
        ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
        
        if self.modality=='metabolomics_glycomics':
            libname = 'smpdb_pathway'
        else:
            libname = 'sub_class'
        
        if top is None:
            data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
        else:
            data = ora_data.loc[(ora_data['modality']==self.modality)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

        if data.shape[0]!=0:

            lib = pd.read_csv(f'../lib/{libname}.csv',index_col=0)
            lib['member']=lib['member'].str.split(';')
            lib = lib.explode('member').reset_index(drop=True)

            pathway_list = data['pathway'].tolist()
            lib_temp = lib.loc[lib['name'].isin(pathway_list)]
            pathway_num = len(pathway_list)

            prop = np.zeros((pathway_num,pathway_num))
            for i in range(pathway_num):
                set1 = set(lib.loc[(lib['name']==pathway_list[i]),'member'])
                for j in range(pathway_num):
                    set2 = set(lib.loc[(lib['name']==pathway_list[j]),'member'])
                    intersection_size = len(set1.intersection(set2))
                    union_size = len(set1.union(set2))
                    prop[i][j]=intersection_size/union_size
            prop_df=pd.DataFrame(data=prop,index=pathway_list,columns=pathway_list)

            #network
            fig, ax = plt.subplots(figsize=(50,40))
            plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
            G = nx.Graph()

            for i in range(prop_df.shape[0]):
                G.add_node(prop_df.columns[i])
                for j in range(i+1, prop_df.shape[1]):
                    if prop_df.iloc[i,j] > 0.25: # set a correlation threshold here
                        G.add_edge(prop_df.columns[i],prop_df.columns[j],weight=prop_df.iloc[i,j]*50)

            pos = nx.circular_layout(G)
            
            ymax = max([coord[1] for coord in pos.values()])
            ymin = min([coord[1] for coord in pos.values()])
            crit = (ymax-ymin)*0.05
            
            for pathway1, coord1 in pos.items():
                if (ymax>coord1[1]>ymax-crit) | (ymin<coord1[1]<ymin+crit):
                    coord1[0] *= 1.5

            plt.xlim(min([coord[0] for coord in pos.values()])*1.5,max([coord[0] for coord in pos.values()])*1.5)

            label_pathways = data['pathway'].tolist()
            node_size = []
            node_color = []
            labels = {}
            for node in G.nodes():
                node_size.append(data.loc[data['pathway']==node,'enrichment_ratio']*800)
                node_color.append(data.loc[data['pathway']==node,'Raw.p'])
                if node in label_pathways:
                    labels[node] = node

            label_pos = {}
            for key,cord in pos.items():
                if key in label_pathways:
                    label_pos[key] = (cord[0],cord[1]-0.05)
                    
            weights = [G[u][v]['weight'] for u,v in G.edges()]

            nx.draw_networkx_labels(G, label_pos,labels, font_size=45,font_weight='bold', verticalalignment='top')
            nx.draw(G, pos, node_size=node_size, node_color=node_color,cmap='PuRd_r',edge_color='lightblue',width=weights,edgecolors='black',linewidths=5)

            plt.savefig(os.path.join('../results/pathway/',f'network_{self.region}_{self.modality}_{cluster}.png'))
            if show == False:
                plt.close()
            else:
                plt.show()