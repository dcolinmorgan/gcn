import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,glob,sys,importlib,pickle,tqdm
from joblib import Parallel,delayed
import multiprocessing
import graph_tool.all as gt
from graph_tool import *
import networkx as nx
from graph_tool.all import *
from pylab import *
# from varname import nameof

# patt='all'
# sys.path.insert(1, './run/gcn/')
# import gcn_func
# importlib.reload(sys.modules['gcn_func'])
# from gcn_func import bip, load_list_of_dicts, meas, time_bar,proc_dat,rev_tbar,group_time_plot,time_order_net,build_gcn#plotRidge,LayeredNetworkGraph,plot_sankey
# from gcn_func import dgtz,shan_entropy,get_SI,get_red,calc_MI,calc_CMI,calc_II,puc_cal
# ##BUILD BI-NETS
# relgene=pd.read_csv('all_arg_subset_genefamilies-cpm.tsv',sep='\t')
# relgene['gene']=relgene['# Gene Family'].str.split('|').str[0]
# relgene=relgene[relgene['gene']!='UniRef90_unknown']
# relgene=relgene[relgene['gene']!='UNMAPPED']
# relgene.index=relgene['# Gene Family']
# del relgene['gene'], relgene['# Gene Family']
# relgene['gen']=relgene.index.str.split('|').str[1].str.split('.').str[0].tolist()
# relgene['spec']=relgene.index.str.split('.').str[1]#.str.split('.').str[0].tolist()
# relgene['spec'].replace('_',' ')
# relgene.index=relgene.index.str.split('|').str[0]
# relgene=relgene.dropna()
# cc=relgene.groupby(['# Gene Family','spec']).sum()
# cc=cc.reset_index()
# cc=cc.rename(columns={'# Gene Family':'gene'})
# fnet=[]
# Bnet=[]
# min_deg=3
# # def buildddd(i,net,cc,min_deg=5):
# for i,net in enumerate(relgene.columns[:-2]):
#     dd=cc[['spec','gene',net]]
#     # dd=dd[dd[net]!=0]
#     dd=dd[dd[net]>1]
# #     # ee=nx.from_pandas_edgelist(dd,source='spec',target='gene',edge_attr=net)
# #     # remove = [node for node,degree in dict(ee.degree()).items() if degree <min_deg]
# #     # ee.remove_nodes_from(remove)
# #     # fnet.append(ee)

#     B = nx.Graph()
#     dd=dd[dd['spec'].str.contains('s__')]
#     dd=dd[dd['gene'].str.contains('UniRef')]
#     B.add_nodes_from(dd['spec'], bipartite=0)
#     B.add_nodes_from(dd['gene'], bipartite=1)
#     B.add_edges_from(tuple(dd[['spec','gene']].itertuples(index=False, name=None)))
#     remove = [node for node,degree in dict(B.degree()).items() if degree <min_deg]
#     B.remove_nodes_from(remove)
#     Bnet.append(B)
    
# # # with open('data/gcn/NX_Ad'+str(min_deg)+'_more_ARG.pkl', 'wb') as f:
# # #     pickle.dump(fnet, f)
# with open('data/gcn/B2X_Emore_ARG.pkl', 'wb') as f:
#     pickle.dump(Bnet, f)

    
# ##NESTEDNESS ANALYSIS
# print(sys.argv[1])
# cc=int(sys.argv[1])
# # rand=(sys.argv[2])
# # deg_rand=float(sys.argv[3])
# # print(cc)
# # cc=int(cc)
# dd=cc+10
# sys.path.insert(1, './run/gcn/')
# import gcn_func
# importlib.reload(sys.modules['gcn_func'])
# from gcn_func import bip, load_list_of_dicts, meas, calc_bipart, bipart,time_bar,proc_dat,rev_tbar,group_time_plot,time_order_net,build_gcn,shuffle_net,structural_analysis
# relgene=pd.read_csv('all_arg_subset_genefamilies-cpm.tsv',sep='\t',nrows=1)
# Ngraphs = load_list_of_dicts('data/gcn/NX_Emore_ARG.pkl')
# from networkx import bipartite
# # relgene=pd.read_csv('all_arg_subset_genefamilies-cpm.tsv',sep='\t',nrows=1)
# Ngraphs = load_list_of_dicts('data/gcn/NX_Emore_ARG.pkl')

# ARG_meta=pd.read_excel('run/gcn/ARG_treatment_infor_modified.xlsx',index_col=0)
# ARG_meta2=pd.read_excel('run/gcn/patients_Tx_batch3_for_DM.xlsx',index_col=None,skiprows=1,names=['id','group'])
# relgene.columns=relgene.columns.str.replace("-00", "-00ST")
# relgene.columns=relgene.columns.str.replace("-00STST", "-00ST")
# relgene.columns=relgene.columns.str.split('-').str[0]+'-'+relgene.columns.str.split('-').str[1]
# ARG_meta['id']=ARG_meta['id'].str.replace('-00ST','')
# META=pd.concat([pd.DataFrame(ARG_meta[['id','group']]),ARG_meta2],ignore_index=True)
# Bgraphs = load_list_of_dicts('data/gcn/BXX_Emore_ARG.pkl')
# B2graphs = load_list_of_dicts('data/gcn/B2X_Emore_ARG.pkl')

# PS=pd.DataFrame(columns=['G','C','AC','RAC','D','SB'])

# for ii,i in enumerate(relgene.columns[cc:dd]):
    # N,Q,I,R=Parallel(n_jobs=10)(structural_analysis(ii,i,graphs,ARG_meta,rand=rand,deg_rand=deg_rand))
# # for ii,i in enumerate(relgene.columns[10:]):
#     # structural_analysis(ii,i,graphs,ARG_meta)# for ii,i in enumerate(relgene.columns[cc:dd]))


# # # structural_analysis(ii,i,graphs,ARG_meta) for ii,i in enumerate(relgene.columns[3:])


# # # pd.DataFrame([N,Q,I,R]).to_csv('Emore_ARG_nested_an.txt',sep='\t')

# # tissue=sys.argv[1]
# # omic=sys.argv[2]

# # ##PIDC analysis
# # data=pd.read_csv('data/Pipeline_consolidate_220301/'+tissue+'_gene_sum_'+omic+'.txt',sep='\t',header=None,index_col=0)
# # genes=data.index
# # data=data[0:10]
# # puc_genes=[]
# # # print(data.shape)
# num_cores = multiprocessing.cpu_count()
# # # Parallel(n_jobs=num_cores)(delayed(puc_cal)(data,i,genes,str(tissue),str(omic)) for i in np.arange(1,len(data)))

# # for i in np.arange(1,len(data)):
# #     puc_cal(data,i,genes,str(tissue),str(omic))
# # # with open('data/Pipeline_consolidate_220301/gcn/'+tissue+'_'+omic+'cdf_ALL.pkl', 'wb') as f:
# # #     pickle.dump(cdf, f)
# # # with open('data/Pipeline_consolidate_220301/gcn/'+tissue+'_'+omic+'pucN_ALL.pkl', 'wb') as f:
# # #     pickle.dump(pucN, f)

    

# Parallel(n_jobs=10)(delayed(calc_bipart)(B2graphs,META,ii,i) for ii,i in enumerate(relgene.columns[1:]))


# files=glob.glob('data/gcn/comp_net/*post_nest_clust_BIP2.txt')
# for i,file in enumerate(files):
#     jj=os.path.basename(file).split('_post')[0]
#     if i==0:
#         LC=pd.read_csv(file,sep='\t',names=[jj+'tC',jj+'_tLC',jj+'_tNR',jj+'_tDC',jj+'_tBC',jj+'_tCC',jj+'_tSBP',jj+'bC',jj+'_bLC',jj+'_bNR',jj+'_bDC',jj+'_bBC',jj+'_bCC',jj+'_bSBP'])
#     else:
#         LC=pd.merge(LC,pd.read_csv(file,sep='\t',names=[jj+'tC',jj+'_tLC',jj+'_tNR',jj+'_tDC',jj+'_tBC',jj+'_tCC',jj+'_tSBP',jj+'bC',jj+'_bLC',jj+'_bNR',jj+'_bDC',jj+'_bBC',jj+'_bCC',jj+'_bSBP']),right_index=True,left_index=True,how='outer')
# LC.columns=LC.columns.str.replace("_00", "_00ST")
# LC.columns=LC.columns.str.replace("_00STST", "_00ST")
# LC.columns=LC.columns.str.replace("_03ST", "_02ST")

# LCC=LC.melt()
# LCC['trx']=LCC.variable.str.split('_').str[1]
# LCC['time']=LCC.variable.str.split('_').str[2]
# LCC['metric']=LCC.variable.str.split('_').str[3]
# LCC=LCC[LCC.value!=0]
# LCC.to_csv('data/gcn/comp_net/LCC_clust_BIP2.txt',sep='\t',header=True,index=True)


# %reset -f
# %config Completer.use_jedi = True
# %matplotlib widget
from scipy.stats import rankdata
from sklearn.preprocessing import normalize
import sklearn.utils as sku
import plotly.graph_objects as go
import plotly.express as px
import chart_studio.plotly as py
# import chart_studio
# chart_studio.tools.set_credentials_file(username='dcolinmorgan', api_key='9FS3nO6nWYFq5zT6BRHD')
import plotly
import scipy
import warnings
warnings.filterwarnings('ignore')

##network stuff
import algorithmx
import networkx as nx
import leidenalg as la
import igraph as ig
import community as community_louvain
import networkx.algorithms.community as nx_comm
from networkx.algorithms import bipartite
import graphistry
graphistry.register(api=3, username='dcolinmorgan', password='f5UwthGEF@F@xnP')

import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle,tqdm
from itertools import combinations,chain#,scipy,coolbox,pybedtools,
# from scipy.stats import linregress
# from scipy.ndimage import gaussian_filter
from tqdm import tqdm
from IPython.display import Image
import pandas as pd
import seaborn as sns
from scipy import stats
import networkx as nx
from pathlib import Path
# import pyvis
# from pyvis.network import Network
import networkit as nk
# from statannot import add_stat_annotation
from statannotations.Annotator import Annotator
# import biosppy
# from sklearn import metrics

os.chdir('/home/dcmorgan')
os.getcwd()

patt='all'
sys.path.insert(1, './run/gcn/')
import gcn_func
importlib.reload(sys.modules['gcn_func'])
from gcn_func import pdize_net, plot_comm, load_list_of_dicts

sys.path.insert(1, './nestedness_analysis/')
import nestedness_metrics_other_functions
from nestedness_metrics_other_functions import from_edges_to_matrix
# importlib.reload(sys.modules['EO_functions_bipartite'])
# import extremal_bi
sys.path.insert(1, './nestedness_analysis/network-generator/')
import generate_synthetic_networks
from netgen import NetworkGenerator

from networkx import bipartite
from functools import reduce
def bipart(G,method,nodes):
    return pd.DataFrame.from_dict(eval('bipartite.'+method)(G,nodes),orient='index',columns=[str(method)])

# for j in (OTHER_00ST,OTHER_01ST,OTHER_02ST,CLA_00ST,CLA_01ST,CLA_02ST,LEVO_00ST,LEVO_01ST,LEVO_02ST):
OTHER_00ST = load_list_of_dicts('data/gcn/OTHER_00ST.pkl')
OTHER_01ST = load_list_of_dicts('data/gcn/OTHER_01ST.pkl')
OTHER_02ST = load_list_of_dicts('data/gcn/OTHER_02ST.pkl')

CLA_00ST = load_list_of_dicts('data/gcn/CLA_00ST.pkl')
CLA_01ST = load_list_of_dicts('data/gcn/CLA_01ST.pkl')
CLA_02ST = load_list_of_dicts('data/gcn/CLA_02ST.pkl')

LEVO_00ST = load_list_of_dicts('data/gcn/LEVO_00ST.pkl')
LEVO_01ST = load_list_of_dicts('data/gcn/LEVO_01ST.pkl')
LEVO_02ST = load_list_of_dicts('data/gcn/LEVO_02ST.pkl')


O_0b,O_0c,O_0d=pdize_net(OTHER_00ST)
O_1b,O_1c,O_1d=pdize_net(OTHER_01ST)
O_2b,O_2c,O_2d=pdize_net(OTHER_02ST)

C_0b,C_0c,C_0d=pdize_net(CLA_00ST)
C_1b,C_1c,C_1d=pdize_net(CLA_01ST)
C_2b,C_2c,C_2d=pdize_net(CLA_02ST)

L_0b,L_0c,L_0d=pdize_net(LEVO_00ST)
L_1b,L_1c,L_1d=pdize_net(LEVO_01ST)
L_2b,L_2c,L_2d=pdize_net(LEVO_02ST)

comm_met=['multilevel','label_propagation','spinglass','infomap','leiden']#,'walktrap','edge_betweenness','leading_eigenvector','label_propagation']#,'optimal_modularity'] #10

def ready_net(C):#C in ['C','L','O']:
    LL=eval(str(C)+'_0b').merge(eval(str(C)+'_1b'),right_on=['source','target'],left_on=['source','target'],how='outer').merge(eval(str(C)+'_2b'),right_on=['source','target'],left_on=['source','target'],how='outer').fillna(0)
    LL.columns=['00ST','01ST','02ST']
    # modes = pd.concat([LL['mode'] for df in (df1, df2, df3, df4)], ignore_index=True).unique()
    LL.reset_index(inplace=True)
    LL['source']=LL['source'].str.split('_').str[2]
    
    return LL
CLA=ready_net('C')
LEVO=ready_net('L')
OTHER=ready_net('O')
species=pd.unique(np.concatenate([CLA['source'],LEVO['source'],OTHER['source']]))
colors = sns.color_palette('hls', len(species))
palette = {mode: color for mode, color in zip(species, colors)}

Parallel(n_jobs=5)(delayed(plot_comm)('CLA','LEVO','OTHER',palette,j) for j in comm_met)



# OTHER=(O_0b).merge(O_1b,right_on=['source','target'],left_on=['source','target'],how='outer').merge(O_2b,right_on=['source','target'],left_on=['source','target'],how='outer').fillna(0)
# OTHER.columns=['00ST','01ST','02ST']

# CLA=(C_0b).merge(C_1b,right_on=['source','target'],left_on=['source','target'],how='outer').merge(C_2b,right_on=['source','target'],left_on=['source','target'],how='outer').fillna(0)
# CLA.columns=['00ST','01ST','02ST']

# LEVO=(L_0b).merge(L_1b,right_on=['source','target'],left_on=['source','target'],how='outer').merge(L_2b,right_on=['source','target'],left_on=['source','target'],how='outer').fillna(0)
# LEVO.columns=['00ST','01ST','02ST']

# for i,LL in enumerate([OTHER,CLA,LEVO]):
#     jj=pd.DataFrame(LL['00ST'].reset_index())
#     jj.columns=['start','end','value']
#     jj=jj[jj['value']>0]
#     kk=pd.DataFrame(np.unique(jj[['start','end']].melt()['value']))
#     kk.columns=['name']
#     kk['id']=kk['name']
#     kk['types']=np.concatenate([np.ones(len(np.unique(jj.end))),np.zeros(len(np.unique(jj.start)))]).astype('int')
#     G_0=igraph_from_pandas(edges_table=jj, vertices_table=kk, source_cl='start', target_cl='end', vertex_attrs=list(kk.columns), vertex_id_cl='name', directed=True)

#     jj=pd.DataFrame(LL['01ST'].reset_index())
#     jj.columns=['start','end','value']
#     jj=jj[jj['value']>0]
#     kk=pd.DataFrame(np.unique(jj[['start','end']].melt()['value']))
#     kk.columns=['name']
#     kk['id']=kk['name']
#     kk['types']=np.concatenate([np.ones(len(np.unique(jj.end))),np.zeros(len(np.unique(jj.start)))]).astype('int')
#     G_1=igraph_from_pandas(edges_table=jj, vertices_table=kk, source_cl='start', target_cl='end', vertex_attrs=list(kk.columns), vertex_id_cl='name', directed=True)


#     jj=pd.DataFrame(LL['02ST'].reset_index())
#     jj.columns=['start','end','value']
#     jj=jj[jj['value']>0]
#     kk=pd.DataFrame(np.unique(jj[['start','end']].melt()['value']))
#     kk.columns=['name']
#     kk['id']=kk['name']
#     kk['types']=np.concatenate([np.ones(len(np.unique(jj.end))),np.zeros(len(np.unique(jj.start)))]).astype('int')
#     G_2=igraph_from_pandas(edges_table=jj, vertices_table=kk, source_cl='start', target_cl='end', vertex_attrs=list(kk.columns), vertex_id_cl='name', directed=True)

#     print([G_0.is_bipartite(),G_1.is_bipartite(),G_2.is_bipartite()])

#     optimiser = la.Optimiser()
#     G_coupling = ig.Graph.Formula('1 -- 2 -- 3');
#     G_coupling.es['weight'] = 0.1; # Interslice coupling strength
#     G_coupling.vs['slice'] = [G_0, G_1, G_2]
    
#     # G_coupling.vs['types'] = [G_1.vs['types'], G_2.vs['types'], G_3.vs['types']]
#     node_size=0
#     gamma=0

#     # membership, improvement = la.find_partition_temporal([G_1, G_2, G_3],la.CPMVertexPartition,interslice_weight=0.1,resolution_parameter=gamma)
#     layers, interslice_layer, G_full = la.time_slices_to_layers([G_0, G_1, G_2],interslice_weight=0.1)
#     partitionsA = [la.CPMVertexPartition(H,weights='weight',resolution_parameter=gamma)for H in layers]
#     partitionsA.append(la.CPMVertexPartition(interslice_layer, resolution_parameter=0, weights='weight'))
#     # diff = optimiser.optimise_partition_multiplex(partitions + [interslice_partition]);

#     # membership, improvement = la.find_partition_multiplex([G_1, G_2, G_3],la.ModularityVertexPartition)
#     layers, interslice_layer, G_full = la.slices_to_layers(G_coupling)
#     partitionsB = [la.CPMVertexPartition(H, node_sizes='node_size',weights='weight',resolution_parameter=gamma)for H in layers]
#     partitionsB.append(la.CPMVertexPartition(interslice_layer, resolution_parameter=0,node_sizes='node_size', weights='weight'))
# # diff = optimiser.optimise_partition_multiplex(partitions + [interslice_partition]);
    
#     if i==0:
#         j='OTHER'
#     elif i==1:
#         j='CLA'
#     elif i==2:
#         j='LEVO'
    
#     g = G_0.to_graph_tool()
#     in_hist = vertex_hist(g, "in")

#     y = in_hist[0]
#     err = sqrt(in_hist[0])
#     err[err >= y] = y[err >= y] - 1e-2

#     figure(figsize=(6,4))
#     errorbar(in_hist[1][:-1], in_hist[0], fmt="o", yerr=err,
#             label="00ST_in_deg")

#     g = G_1.to_graph_tool()
#     in_hist = vertex_hist(g, "in")

#     y = in_hist[0]
#     err = sqrt(in_hist[0])
#     err[err >= y] = y[err >= y] - 1e-2

#     # figure(figsize=(6,4))
#     errorbar(in_hist[1][:-1], in_hist[0], fmt="o", yerr=err,
#             label="01ST_in_deg")

#     g = G_2.to_graph_tool()
#     in_hist = vertex_hist(g, "in")

#     y = in_hist[0]
#     err = sqrt(in_hist[0])
#     err[err >= y] = y[err >= y] - 1e-2

#     # figure(figsize=(6,4))
#     errorbar(in_hist[1][:-1], in_hist[0], fmt="o", yerr=err,
#             label="02ST_in_deg")
#     plt.legend(loc='best')#, bbox_to_anchor=(1.03, 1))
#     gca().set_yscale("log")
#     gca().set_xscale("log")
#     gca().set_ylim(1e-1, 1e4)
#     gca().set_xlim(0.8, 1e2)
#     subplots_adjust(left=0.2, bottom=0.2)
#     xlabel("$k_{in}$")
#     ylabel("$NP(k_{in})$")

#     tight_layout()
#     savefig('run/gcn/img/deg/in_degree_'+j+'.png')

#     for i in np.arange(0,len(partitionsA)-1):
#         plot_time_net(partitionsB,i,(j))
    