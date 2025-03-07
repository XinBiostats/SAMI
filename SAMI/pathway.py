import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import networkx as nx
os.environ['R_HOME'] = "/blue/li.chen1/xin.ma/conda/envs/maldi/lib/R"
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

class Pathway:
    """
    A class for performing pathway enrichment analysis and visualization for metabolomics, lipidomics, and glycomics data.

    Attributes
    ----------
    data_path : str
        Path to the directory containing input data.
    result_path : str
        Path to the directory where results will be saved.
    region : str
        Name of the region being analyzed.
    omics : str
        Type of omics data (e.g., metabolomics, lipidomics, glycomics).
    modality : str
        Specifies whether the analysis is performed on pooled data or a single omics dataset.
    """
    def __init__(self,data_path, result_path,region,omics,pool=True):
        self.data_path = data_path
        self.result_path = result_path
        if not os.path.exists(os.path.join(self.result_path,'pathway')):
            os.makedirs(os.path.join(self.result_path,'pathway'))
            
        self.region = region
        if pool:
            self.modality = 'pool'
        else:
            self.modality = omics
        
        if omics=='metabolomics' or omics=='glycomics':
            print('Metabolomics and Glycomics are combined for pathway enrichment analysis.')
            self.omics = 'metabolomics_glycomics'
        else:
            self.omics = 'lipidomics'
            
    def findpathway(self,direction_threshold):
        """
        Performs pathway enrichment analysis and assigns regulation direction based on fold change.

        Parameters
        ----------
        direction_threshold : float
            Threshold for determining up- or down-regulated pathways.

        Notes
        -----
        - Calls an R script (`MetabAnalystR.R`) for pathway enrichment.
        - Identifies pathways as up-regulated, down-regulated, or neutral based on the average log2 fold change of pathway compounds.
        - Saves the enriched pathway results as a CSV file.
        """
        r_source = ro.r['source']
        r_source('../SAMI/MetabAnalystR.R')  
        
        markers = pd.read_csv(os.path.join(self.result_path,f'markers/{self.modality}/{self.region}_marker.csv'))
        clusters = markers['cluster'].unique().tolist()
        
        if self.omics=='metabolomics_glycomics':
            markers = markers.loc[markers['omics'].isin(['metabolomics','glycomics'])]
            
            result_pathR = ro.StrVector([self.result_path])
            logical_arg = ro.BoolVector([True,False])
            modalityR = ro.StrVector([self.modality])
            clustersR = ro.vectors.IntVector(clusters)
            regionR = ro.StrVector([self.region])
            with localconverter(ro.default_converter + pandas2ri.converter):
                markersR = ro.conversion.py2rpy(markers)

            PathwayEnrichment = ro.r['PathwayEnrichment']
            
            PathwayEnrichment(result_path=result_pathR,modality=modalityR,region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[1])
            
            #### 2025-03-05 add direction ###
            pathways = pd.read_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'))
            measured_compounds = pd.read_csv(os.path.join(self.result_path,f'markers/{self.modality}/{self.region}_pvalue.csv'))
            
            pathways['pathway_avg_log2FC'] = None
            pathways['regulated'] = None
            for idx, row in pathways.iterrows():
                cluster = row['cluster']
                all_compounds = row['all_compounds'].split(';')
                select_measured_compounds = measured_compounds[(measured_compounds['cluster']==cluster)&(measured_compounds['feature'].isin(all_compounds))]
                pathway_fc = round(np.log2(select_measured_compounds['mean_a'].sum()/select_measured_compounds['mean_b'].sum()),4)
                pathways.at[idx, 'pathway_avg_log2FC'] = pathway_fc
                if pathway_fc > direction_threshold:
                    pathways.at[idx, 'regulated'] = 'up'
                elif pathway_fc < -direction_threshold:
                    pathways.at[idx, 'regulated'] = 'down'
                else:
                    pathways.at[idx, 'regulated'] = 'neutral'
                
            pathways.to_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'),index=False)
                    
        else:
            markers = markers.loc[markers['omics']=='lipidomics']
        
            result_pathR = ro.StrVector([self.result_path])
            logical_arg = ro.BoolVector([True,False])
            modalityR = ro.StrVector([self.modality])
            clustersR = ro.vectors.IntVector(clusters)
            regionR = ro.StrVector([self.region])
            with localconverter(ro.default_converter + pandas2ri.converter):
                markersR = ro.conversion.py2rpy(markers)

            PathwayEnrichment = ro.r['PathwayEnrichment']
        
            PathwayEnrichment(result_path=result_pathR,modality=modalityR,region=regionR,markers=markersR,clusters=clustersR,lipid=logical_arg[0])
            
            #### 2025-03-05 add direction ###
            pathways = pd.read_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'))
            measured_compounds = pd.read_csv(os.path.join(self.result_path,f'markers/{self.modality}/{self.region}_pvalue.csv'))
            
            pathways['pathway_avg_log2FC'] = None
            pathways['regulated'] = None
            for idx, row in pathways.iterrows():
                cluster = row['cluster']
                all_compounds = row['all_compounds'].split(';')
                select_measured_compounds = measured_compounds[(measured_compounds['cluster']==cluster)&(measured_compounds['feature'].isin(all_compounds))]
                pathway_fc = round(np.log2(select_measured_compounds['mean_a'].sum()/select_measured_compounds['mean_b'].sum()),4)
                pathways.at[idx, 'pathway_avg_log2FC'] = pathway_fc
                if pathway_fc > direction_threshold:
                    pathways.at[idx, 'regulated'] = 'up'
                elif pathway_fc < -direction_threshold:
                    pathways.at[idx, 'regulated'] = 'down'
                else:
                    pathways.at[idx, 'regulated'] = 'neutral'
                
            pathways.to_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'),index=False)
        
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
        """
        Generates a dot plot to visualize enriched pathways based on statistical significance (p-value and FDR).

        Parameters
        ----------
        cluster : int
            Cluster ID for which the pathway enrichment analysis is performed.
        scale : float
            Scaling factor for the dot sizes based on the enrichment ratio.
        height : float
            Height of the output figure.
        top : int, optional
            Number of top pathways to display. If None, all pathways are shown. Default is None.
        show : bool, optional
            If True, displays the plot instead of saving it. Default is False.

        Notes
        -----
        - The function plots two dot plots side by side:
          - Left: -log10(p-value) vs. pathways.
          - Right: -log10(FDR) vs. pathways.
        - Dot size represents the enrichment ratio.
        - Dot color represents p-value and FDR values.
        - Saves the figure as a PNG file in the pathway results directory.

        Output
        ------
        - A PNG file containing the dot plot is saved in the pathway results folder.
        """
        ora_data = pd.read_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'))
        ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
        
        if top is None:
            data = ora_data.loc[(ora_data['omics']==self.omics)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
        else:
            data = ora_data.loc[(ora_data['omics']==self.omics)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

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

            plt.savefig(os.path.join(self.result_path,f'pathway/{self.modality}/{self.region}_{self.omics}_{cluster}_dotplot.png'))
            if show == False:
                plt.close()
            else:
                plt.show()
                
    def plot_bar_for_direction(self,cluster,top=None,show=False):
        """
        Plots bar charts of significant up- and down-regulated pathways.

        Parameters
        ----------
        cluster : int
            Cluster ID to analyze.
        top : int, optional
            Number of top pathways to display. If None, all significant pathways are shown. Default is None.
        show : bool, optional
            If True, displays the plot instead of saving. Default is False.

        Notes
        -----
        - Uses different colors for up-regulated and down-regulated pathways.
        - Saves the plot as a PNG file.
        """
        ora_data = pd.read_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'))
        if top is None:
            up_data = ora_data[(ora_data['regulated'] == 'up')&(ora_data['cluster']==cluster)]
            down_data = ora_data[(ora_data['regulated'] == 'down')&(ora_data['cluster']==cluster)]
        else:
            up_data = ora_data[(ora_data['regulated'] == 'up')&(ora_data['cluster']==cluster)].head(top)
            down_data = ora_data[(ora_data['regulated'] == 'down')&(ora_data['cluster']==cluster)].head(top)

        fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)


        if up_data.shape[0]!=0:
            ax[0].barh(up_data['pathway'], -np.log10(up_data['Raw.p']), color='#e64b35')
            ax[0].text(1.05, 0.5, 'Up-Regulated', fontsize=15, fontweight='bold', transform=ax[0].transAxes, ha='right', va='center', rotation=270)
            #ax[0].set_ylabel('Pathway',fontsize=15,fontweight='bold')
            ax[0].invert_yaxis()
            ax[0].set_yticks(range(len(up_data['pathway'])))
            ax[0].set_yticklabels(up_data['pathway'], fontsize=15, weight='bold')
            ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
            for label in ax[0].get_xticklabels():
                label.set_fontsize(12)
                label.set_fontweight('bold')

        else:
            fig.delaxes(ax[0])

        if down_data.shape[0]!=0:
            ax[1].barh(down_data['pathway'], -np.log10(down_data['Raw.p']), color='#4dbbd5')
            ax[1].text(1.05, 0.5, 'Down-Regulated', fontsize=15, fontweight='bold', transform=ax[1].transAxes, ha='right', va='center', rotation=270)
            ax[1].set_xlabel('-log10(p-value)',fontsize=15,fontweight='bold')
            #ax[1].set_ylabel('Pathway',fontsize=15,fontweight='bold')
            ax[1].invert_yaxis()
            ax[1].set_yticks(range(len(down_data['pathway'])))
            ax[1].set_yticklabels(down_data['pathway'], fontsize=15, weight='bold')
            ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
            for label in ax[1].get_xticklabels():
                label.set_fontsize(12)
                label.set_fontweight('bold')
        else:
            fig.delaxes(ax[1])
            
        fig.text(0.02, 0.5, 'Pathway', fontsize=20, fontweight='bold', ha='center', va='center', rotation=90)
            
        plt.tight_layout(rect=[0.05, 0, 1, 1])
        plt.savefig(os.path.join(self.result_path,f'pathway/{self.modality}/{self.region}_{self.omics}_{cluster}_barplot.png'))
        if show == False:
            plt.close()
        else:
            plt.show()  
    
    def pathway_network(self,cluster,top=None,show=False):
        """
        Generates a pathway network visualization based on Jaccard similarity between pathways.

        Parameters
        ----------
        cluster : int
            Cluster ID for which the pathway network is constructed.
        top : int, optional
            Number of top pathways to include in the network. If None, all pathways are used. Default is None.
        show : bool, optional
            If True, displays the plot instead of saving it. Default is False.

        Notes
        -----
        - The function constructs a network where nodes represent pathways.
        - An edge between two pathways is added if they share common metabolites.
        - Edge weights are determined by the Jaccard similarity between pathway members.
        - Pathway enrichment ratio is used to determine node size.
        - Node color represents the p-value of the pathway.

        Output
        ------
        - A pathway network visualization is saved as a PNG file in the pathway results directory.
        - If `show=True`, the plot is displayed interactively.
        """
        ora_data = pd.read_csv(os.path.join(self.result_path,f'pathway/{self.modality}/ora_{self.region}_{self.omics}.csv'))
        ora_data['enrichment_ratio'] = ora_data['hits']/ora_data['expected']
        
        if self.omics=='metabolomics_glycomics':
            libname = 'smpdb_pathway'
        else:
            libname = 'sub_class'
        
        if top is None:
            data = ora_data.loc[(ora_data['omics']==self.omics)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)]
        else:
            data = ora_data.loc[(ora_data['omics']==self.omics)&(ora_data['cluster']==cluster)&(ora_data['region']==self.region)].head(top)

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

            plt.savefig(os.path.join(self.result_path,f'pathway/{self.modality}/network_{self.region}_{self.omics}_{cluster}.png'))
            if show == False:
                plt.close()
            else:
                plt.show()