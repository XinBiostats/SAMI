import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import ranksums
from utils import *
import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import seaborn as sns

def findmarkers(adata,samplename,adj_pval_cutoff,top,annotation=None,merge_cluster=False):
    if merge_cluster==False:
        cluster_param='leiden'
    else:
        cluster_param='merge'
    
    pvalue_result=pd.DataFrame()
    marker_result=pd.DataFrame()
    cluster_list=sorted(adata.obs[cluster_param].unique().tolist())
    for i in range(len(cluster_list)):
        pvalue_temp=pd.DataFrame(columns=['feature','pvalue','avg_log2FC','abs_avg_log2FC','pct1','pct2','adj_pvalue','cluster','rank'])
        for feature in adata.var.index:
            a=adata[adata.obs[cluster_param]==str(cluster_list[i])][:,feature].X.copy()
            b=adata[adata.obs[cluster_param]!=str(cluster_list[i])][:,feature].X.copy()
            pct1=sum(a!=0)/len(a)
            pct2=sum(b!=0)/len(b)
            if (pct1>0.1)&(pct2>0.1):
                stat,pvalue=ranksums(a,b)
                pvalue=pvalue[0]
                mean_a=np.mean(a)
                mean_b=np.mean(b)
                avg_log2FC=np.log2(mean_a/mean_b)
                abs_avg_log2FC=abs(avg_log2FC)
                pvalue_feature=pd.DataFrame({'feature':feature,'pvalue':[pvalue],'avg_log2FC':avg_log2FC,'abs_avg_log2FC':abs_avg_log2FC,'pct1':pct1,'pct2':pct2,'cluster':cluster_list[i]})
                pvalue_temp=pd.concat([pvalue_temp,pvalue_feature],ignore_index=True)
            else:
                continue
        
        #results of all the features
        pvalue_result=pd.concat([pvalue_result,pvalue_temp],ignore_index=True)
               
        #results of markers
        result_temp=pvalue_temp.sort_values(by=['pvalue'],ascending=True).reset_index(drop=True)
        result_temp=adj_pvalue(result_temp).sort_values(by=['adj_pvalue'],ascending=True).loc[result_temp['adj_pvalue']<=adj_pval_cutoff,:]
        result_temp=result_temp.sort_values(by=['abs_avg_log2FC'],ascending=False).reset_index(drop=True).head(top)
        marker_result=pd.concat([marker_result,result_temp],ignore_index=True)
    
    #annotation
    if annotation is not None:
        pvalue_result = pd.merge(pvalue_result,annotation,left_on='feature',right_on='m.z. export',how='left')
        pvalue_result['Name'] = pvalue_result.apply(lambda row: str(row['feature']) if pd.isna(row['Name']) else row['Name'], axis=1)
        pvalue_result=pvalue_result.drop(columns=['m.z. export'])
        
        marker_result = pd.merge(marker_result,annotation,left_on='feature',right_on='m.z. export',how='left')
        marker_result['Name'] = marker_result.apply(lambda row: str(row['feature']) if pd.isna(row['Name']) else row['Name'], axis=1)
        marker_result=marker_result.drop(columns=['m.z. export'])
    else:
        pvalue_result['Name'] = pvalue_result['feature']
        marker_result['Name'] = marker_result['feature']
        
    pvalue_result=pvalue_result.drop(columns=['rank','pct1','pct2','adj_pvalue'])
    marker_result=marker_result.drop(columns=['rank'])
    
    
    file_path = '../results/markers/'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    if merge_cluster==False:
        if not os.path.exists(os.path.join(file_path,'pvalues.xlsx')):
            with pd.ExcelWriter(os.path.join(file_path,'pvalues.xlsx'), mode="w", engine="openpyxl") as writer:
                pvalue_result.to_excel(writer, samplename, index=False)
        else:
            with pd.ExcelWriter(os.path.join(file_path,'pvalues.xlsx'), mode="a", engine="openpyxl",if_sheet_exists="replace") as writer:
                pvalue_result.to_excel(writer, samplename, index=False)

        if not os.path.exists(os.path.join(file_path,'markers.xlsx')):
            with pd.ExcelWriter(os.path.join(file_path,'markers.xlsx'), mode="w", engine="openpyxl") as writer:
                marker_result.to_excel(writer, samplename, index=False)
        else:
            with pd.ExcelWriter(os.path.join(file_path,'markers.xlsx'), mode="a", engine="openpyxl",if_sheet_exists="replace") as writer:
                marker_result.to_excel(writer, samplename, index=False)
    else:
        if not os.path.exists(os.path.join(file_path,'pvalues_merge.xlsx')):
            with pd.ExcelWriter(os.path.join(file_path,'pvalues_merge.xlsx'), mode="w", engine="openpyxl") as writer:
                pvalue_result.to_excel(writer, samplename, index=False)
        else:
            with pd.ExcelWriter(os.path.join(file_path,'pvalues_merge.xlsx'), mode="a", engine="openpyxl",if_sheet_exists="replace") as writer:
                pvalue_result.to_excel(writer, samplename, index=False)

        if not os.path.exists(os.path.join(file_path,'markers_merge.xlsx')):
            with pd.ExcelWriter(os.path.join(file_path,'markers_merge.xlsx'), mode="w", engine="openpyxl") as writer:
                marker_result.to_excel(writer, samplename, index=False)
        else:
            with pd.ExcelWriter(os.path.join(file_path,'markers_merge.xlsx'), mode="a", engine="openpyxl",if_sheet_exists="replace") as writer:
                marker_result.to_excel(writer, samplename, index=False)
        
        
        
def circular_tree(markerfile,section,clusters,top=5):
    markers_file = pd.ExcelFile(os.path.join('../results/markers/',markerfile))
    sheets_names = markers_file.sheet_names

    markers = pd.DataFrame()
    for sheet in sheets_names:
        sections, omics = sheet.split('_')
        temp = pd.read_excel(markers_file,sheet)
        temp[['omics']] = sheet.split('_')[1]
        temp[['sections']] = sections
        markers = pd.concat([markers,temp], ignore_index=True)

    top = markers.groupby(['sections','omics','cluster'],as_index=False).head(top)
    top['c|s'] = top['cluster'].astype(str)+'|'+top['omics']
    top['c|s|f'] = top['cluster'].astype(str)+'|'+top['omics']+'|'+top['Name']
    top['center'] = 'center'
    top = top.loc[(top['sections']==section)&(top['cluster'].isin(clusters))]
    data1 = top.rename(columns={'c|s|f':'start','c|s':'target'}).astype(str)
    data2 = top.rename(columns={'c|s':'target','c|s|f':'start'}).astype(str)
    data3 = top.rename(columns={'c|s':'start','cluster':'target'}).drop_duplicates(subset=['start','target'],keep='first').astype(str)
    data4 = top.rename(columns={'c|s|f':'target','c|s':'start'}).drop_duplicates(subset=['start','target'],keep='first').astype(str)
    data5 = top.rename(columns={'cluster':'start','center':'target'}).drop_duplicates(subset=['start','target'],keep='first').astype(str)
    data6 = top.rename(columns={'center':'start','cluster':'target'}).drop_duplicates(subset=['start','target'],keep='first').astype(str)
    data = pd.concat([data1, data2, data3, data4, data5, data6], axis=0)

    data['label'] = data.apply(lambda x: '' if x['start']=='center' else x['feature'] if x['start'].count('|') == 2 else x['omics'] if x['start'].count('|') == 1 else x['start'] , axis=1)
    data['size'] = data.apply(lambda x: x['abs_avg_log2FC'] if x['start'].count('|') == 2 else 10, axis=1)
    color_list = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    if 'merge' not in markerfile:
        data['color'] = data.apply(lambda x: str(50) if 'meta' in x['start'] else str(69) if 'lipid' in x['start'] else str(75) if 'gly' in x['start'] else str(80) if 'center' in x['start'] else str(int(x['start'])) , axis=1)
    else:
        cluster_merge = markers['cluster'].unique().tolist()
        data['color'] = data.apply(lambda x: str(50) if 'meta' in x['start'] else str(69) if 'lipid' in x['start'] else str(75) if 'gly' in x['start'] else str(80) if 'center' in x['start'] else str(int(cluster_merge.index(x['start']))+44+3) , axis=1)
    data['color'] = data['color'].map(color_list)

    plt.figure(figsize=(15, 15))
    G = nx.Graph()

    nodes = pd.unique(data[['start', 'target']].values.ravel('K'))
    G.add_nodes_from(nodes)

    edges = data[['start', 'target']].values
    G.add_edges_from(edges)

    size_dict = dict(zip(data['start'], data['size']))
    nx.set_node_attributes(G, size_dict, 'size')

    color_dict = dict(zip(data['start'], data['color']))
    nx.set_node_attributes(G, color_dict, 'color')


    label_dict = dict(zip(data['start'], data['label']))
    nx.set_node_attributes(G, label_dict, 'label')

    pos = nx.nx_agraph.graphviz_layout(G, prog="twopi", args="")

    center = pos['center']
    label_pos = {}
    for key,cord in pos.items():
        label_pos[key] = ((cord[0]-center[0])*1.05+center[0],(cord[1]-center[1])*1.03+center[1])

    labels = nx.draw_networkx_labels(G, label_pos, label_dict, font_size=15, font_weight='bold')
    theta = {k: np.arctan2(v[1]-center[1], v[0]-center[0]) * 180/np.pi for k, v in label_pos.items() }
    for key,t in labels.items():
        if key.count('|')==2:
            if 90 < theta[key] or theta[key] < -90 :
                angle = 180 + theta[key]
                t.set_ha('right')
            else:
                angle = theta[key]
                t.set_ha('left')
            t.set_va('center')
            t.set_rotation(angle)
            t.set_rotation_mode('anchor')
        elif key.count('|')==1:
            if 90 < theta[key] or theta[key] < -90 :
                angle = 180 + theta[key]
                t.set_ha('center')
            else:
                angle = theta[key]
                t.set_ha('center')
            t.set_va('center')
            t.set_rotation(angle)
            t.set_rotation_mode('anchor')
        else:
            if 90 < theta[key] or theta[key] < -90 :
                angle = 180 + theta[key]
                t.set_ha('right')
            else:
                angle = theta[key]
                t.set_ha('left')
            t.set_va('center')
            t.set_rotation(angle)
            t.set_rotation_mode('anchor')

    nx.draw(G, pos, node_color=[d['color'] for n, d in G.nodes(data=True)],node_size=[float(d['size'])*10 for n, d in G.nodes(data=True)],edge_color='#D3D3D3')
    x_values, y_values = zip(*pos.values())
    x_max = max(x_values)
    x_min = min(x_values)
    y_max = max(y_values)
    y_min = min(y_values)
    x_margin = (x_max - x_min) * 0.3
    y_margin = (y_max - y_min) * 0.3
    plt.xlim(x_min - x_margin, x_max + x_margin)
    plt.ylim(y_min - y_margin, y_max + y_margin)
    plt.savefig(os.path.join('../results/markers/','circular_tree_'+str(clusters).replace("'","")+'.png'))
    plt.show()
    
    
    
def volcano_plot(pvaluefile,section,omics,cluster):
    pvalues_file = pd.ExcelFile(os.path.join('../results/markers/',pvaluefile))
    sheets_names = pvalues_file.sheet_names

    pvalues = pd.DataFrame()
    for sheet in sheets_names:
        sections, omics = sheet.split('_')
        temp = pd.read_excel(pvalues_file,sheet)
        temp[['omics']] = sheet.split('_')[1]
        temp[['sections']] = sections
        pvalues = pd.concat([pvalues,temp], ignore_index=True)
        
    pval = pvalues.loc[pvalues['pvalue']>0]
    pval = pval.loc[(pval['sections']==section)&(pval['omics']==omics)&(pval['cluster']==cluster)&(pval['pvalue']>=np.percentile(pval['pvalue'],5))&(pval['abs_avg_log2FC']<=np.percentile(pval['abs_avg_log2FC'],95))]
    
    x = pval['avg_log2FC']
    y = -np.log10(pval['pvalue'])
    
    colors = np.where(y <= -np.log10(0.05),'darkgray',np.where(x >= 0, 'red', 'blue'))
    
    fig, ax = plt.subplots(figsize=(10,10))
    fig.subplots_adjust(left=0.2,bottom=0.2)
    plt.scatter(x,y,c=colors,s=10)
    
    crit1=pval[(pval['avg_log2FC']<=0.5)&(pval['avg_log2FC']>=0)].sort_values('pvalue').head(2).index
    crit2=pval[(pval['avg_log2FC']>=-0.5)&(pval['avg_log2FC']<=0)].sort_values('pvalue').head(2).index
    ann_df = pval.loc[pval.index.isin(crit1)|pval.index.isin(crit2)]
    ann_dict = {}
    for index, row in ann_df.iterrows():
        x = row['avg_log2FC']
        y = -np.log10(row['pvalue'])
        annotation = row['feature']
        coord = (x, y)
        ann_dict[coord] = annotation

    texts=[]
    for key, values in ann_dict.items():
        texts.append(plt.text(key[0], key[1], values, size=23, weight='bold')) 
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1))
        
    x_min = np.min(pval['avg_log2FC'])
    x_max = np.max(pval['avg_log2FC'])
    x_abs_max = max(abs(x_min), abs(x_max))*1.1
    plt.xlim(-x_abs_max, x_abs_max)
    
    #plt.axvline(x=0.5, color='gray', linestyle='--',alpha=0.5)
    #plt.axvline(x=-0.5, color='gray', linestyle='--',alpha=0.5)
    plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--',alpha=0.5)
    
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=2))
    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=1))
    plt.xticks(fontsize=40,fontweight='bold')
    plt.yticks(fontsize=40,fontweight='bold')
    
    #plt.title('Volcano Plot for '+omics.capitalize()+' Cluster '+str(cluster),fontsize=23,fontweight='bold')
    plt.xlabel('log2(FC)',fontsize=40,fontweight='bold')
    plt.ylabel('-log10(pvalue)',fontsize=40,fontweight='bold')
    plt.savefig(os.path.join('../results/markers/','vol_'+section+'_'+omics+'_'+str(cluster).replace("'","")+'.png'))
    plt.show()

    

def density_plot(adata1,adata2,cluster,marker):
    adata=sc.read(os.path.join('../results/clustering/',adata1))
    cluster_old=adata.obs[['leiden']]
    adata=sc.read(os.path.join('../datasets/',adata2))
    adata = cluster_assign(adata,cluster_old)
    adata= adata[:,marker]
    data = pd.concat([adata.to_df(),adata.obs],axis=1)
    data['cluster']=np.where(data['leiden']==str(cluster),data['leiden'],'other')

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors['other']='#A0A0A0'
    
    fig,ax = plt.subplots(figsize=(6,2))
    sns.kdeplot(data=data,x=marker,hue='cluster',palette=clusters_colors,fill=True,common_norm=False,alpha=0.5,linewidth=0)
    plt.setp(ax.get_legend().get_texts(), fontsize='10') 
    plt.setp(ax.get_legend().get_title(), fontsize='10',fontweight='bold')

    ax.set_xlabel(marker,fontsize=15,fontweight='bold')
    ax.set_ylabel('Density',fontsize=15,fontweight='bold')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(3))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
    section = adata1.split('_')[0]
    plt.savefig(os.path.join('../results/markers/','density_'+section+'_'+str(cluster)+'_'+marker+'.png'))
    
    
