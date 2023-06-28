import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx


def plot_dot(ora_data,cluster,scale):
    
    sample=ora_data['sample'].tolist()[0]
    
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(10,10))
    fig.subplots_adjust(left=0.45,wspace=0.3)
    #### p-value  #####   
    data = ora_data.loc[(ora_data['sample']==sample)&(ora_data['cluster']==cluster)]
    pathway = data['pathway'].tolist()
    log_trans = -np.log10(data['Raw p'])

    #axis
    y_pos = np.arange(len(pathway))
    x_pos = np.arange(0,np.ceil(max(log_trans))*1.1,5 if max(log_trans)>5 else 1)
    ax[0].set_xticks(x_pos)
    ax[0].tick_params(axis='x', labelsize=15)
    ax[0].set_xlim(-np.ceil(max(log_trans))*0.05,np.ceil(max(log_trans))*1.05)
    ax[0].set_yticks(y_pos,pathway,fontsize=18,fontweight='bold')
    ax[0].invert_yaxis()

    size = data['enrichment_ratio']*scale
    x = log_trans.tolist()
    color = data['Raw p'].tolist()

    scatter = ax[0].scatter(x,y_pos,s=size,c=color,cmap='PuRd_r',zorder=2)
    ax[0].grid(alpha=0.5,zorder=1)
    ax[0].set_xlabel('-log10(p_value)',fontsize=18,fontweight='bold')

    #handles, labels = scatter.legend_elements(prop="sizes", alpha=1,num=4,func=lambda x: x/scale)
    #legend = ax[0].legend(handles, labels, bbox_to_anchor=(1.4,0.4), ncol=1, title="Enrichment\nRatio",fontsize=10,labelspacing=2,borderpad=1.2)
    #legend.get_title().set_fontsize(10)

    cbar = fig.colorbar(scatter,shrink=0.5)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_title('p-value',fontsize=15,fontweight='bold')
    loc = cbar.ax.get_position()
    new_loc =[loc.x0+0.02,loc.y0+0.15,loc.width,loc.height]
    cbar.ax.set_position(new_loc)
    
    fig.set_size_inches(15, 8)
    
    #### FDR  #####   
    data = ora_data.loc[(ora_data['sample']==sample)&(ora_data['cluster']==cluster)]
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
    fig.set_size_inches(15, 8)
    
#     if 'lipid' in sample:
#         fig.suptitle('Pathway Enrichment for Lipidomics (Cluster '+str(cluster)+')',fontsize=20,fontweight='bold')
#     else:
#         fig.suptitle('Pathway Enrichment for Metabolomics & Glycomics (Cluster '+str(cluster)+')',fontsize=20,fontweight='bold')
    
    
    plt.savefig(os.path.join('../results/pathway/',sample+'_'+str(cluster)+'.png'))
    plt.show()
    
    
    
def pathway_network(ora_data,lib_name,cluster):
    sample=ora_data['sample'].tolist()[0]
    
    lib = pd.read_csv('../lib/'+lib_name+'.csv',index_col=0)
    lib['member']=lib['member'].str.split(';')
    lib = lib.explode('member').reset_index(drop=True)

    ora_data_temp = ora_data.loc[(ora_data['cluster']==cluster)&(ora_data['sample']==sample)]
    pathway_list = ora_data_temp['pathway'].tolist()
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
    fig, ax = plt.subplots(figsize=(10,10))
    plt.subplots_adjust(left=0.001,right=0.999)
    G = nx.Graph()

    for i in range(prop_df.shape[0]):
        G.add_node(prop_df.columns[i])
        for j in range(i+1, prop_df.shape[1]):
            if prop_df.iloc[i,j] > 0.25: # set a correlation threshold here
                G.add_edge(prop_df.columns[i],prop_df.columns[j],weight=prop_df.iloc[i, j])

    node_size = [v*10 for v in ora_data_temp['enrichment_ratio'].tolist()]
    node_color = ora_data_temp['Raw p'].tolist()
    pos = nx.spring_layout(G,k=3)

    plt.xlim(min([coord[0] for coord in pos.values()])*2,max([coord[0] for coord in pos.values()])*2)
    label_pos = {}
    for key,cord in pos.items():
        label_pos[key] = (cord[0],cord[1]-0.03)
    nx.draw_networkx_labels(G, label_pos, font_size=15,font_weight='bold', verticalalignment='top')
    nx.draw(G, pos,node_size=node_size, node_color=node_color,cmap='PuRd_r',edge_color='lightgray')
    
    #plt.title("Pathway Network for Metabolomics & Glycomics (Cluster "+str(cluster)+")",fontsize=20,fontweight='bold')
    plt.savefig(os.path.join('../results/pathway/',sample+'_network_'+lib_name+'_'+str(cluster)+'.png'))
    plt.show()