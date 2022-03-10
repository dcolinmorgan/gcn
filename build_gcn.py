import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,glob,sys,importlib,pickle,tqdm
from joblib import Parallel
import networkx as nx
patt='all'
sys.path.insert(1, './run/gcn/')
import gcn_func
importlib.reload(sys.modules['gcn_func'])
from gcn_func import bip, load_list_of_dicts, meas, time_bar,proc_dat,rev_tbar,group_time_plot,time_order_net,build_gcn#plotRidge,LayeredNetworkGraph,plot_sankey

###BUILD BI-NETS
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
# min_deg=5
# def buildddd(i,net,cc,min_deg=5):
# for i,net in enumerate(relgene.columns[:-2]):
#     dd=cc[['spec','gene',net]]
#     dd=dd[dd[net]!=0]
#     ee=nx.from_pandas_edgelist(dd,source='spec',target='gene',edge_attr=net)
#     remove = [node for node,degree in dict(ee.degree()).items() if degree <min_deg]
#     ee.remove_nodes_from(remove)
#     fnet.append(ee)

#     B = nx.Graph()
#     B.add_nodes_from(dd['spec'], bipartite=0)
#     B.add_nodes_from(dd['gene'], bipartite=1)
#     B.add_edges_from(tuple(dd[['spec','gene']].itertuples(index=False, name=None)))
#     remove = [node for node,degree in dict(B.degree()).items() if degree <min_deg]
#     B.remove_nodes_from(remove)
#     Bnet.append(B)
    
# with open('data/gcn/NX_Emore_ARG.pkl', 'wb') as f:
#     pickle.dump(fnet, f)
# with open('data/gcn/BX_Emore_ARG.pkl', 'wb') as f:
#     pickle.dump(Bnet, f)

    
##NESTEDNESS ANALYSIS
print(sys.argv[1])
cc=int(sys.argv[1])
rand=(sys.argv[2])
deg_rand=float(sys.argv[3])
# print(cc)
# cc=int(cc)
dd=cc+10
from gcn_func import bip, load_list_of_dicts, meas, plotRidge,structural_analysis

relgene=pd.read_csv('all_arg_subset_genefamilies-cpm.tsv',sep='\t',nrows=1)
graphs = load_list_of_dicts('data/gcn/NX_Emore_ARG.pkl')
ARG_meta0=pd.read_excel('run/gcn/ARG_treatment_infor_modified.xlsx',index_col=0)
ARG_meta2=pd.read_excel('run/gcn/patients_Tx_batch3_for_DM.xlsx',index_col=None,skiprows=1,names=['id','group'])
relgene.columns=relgene.columns.str.replace("-00", "-00ST")
relgene.columns=relgene.columns.str.replace("-00STST", "-00ST")
relgene.columns=relgene.columns.str.split('-').str[0]+'-'+relgene.columns.str.split('-').str[1]
ARG_meta=pd.concat([pd.DataFrame(ARG_meta0[['id','group']]),ARG_meta2],ignore_index=True)
ARG_meta['id']=ARG_meta['id'].str.replace('-00ST','')

for ii,i in enumerate(relgene.columns[cc:dd]):
    # N,Q,I,R=Parallel(n_jobs=10) (
    structural_analysis(ii,i,graphs,ARG_meta,rand=rand,deg_rand=deg_rand)
# for ii,i in enumerate(relgene.columns[10:]):
    # structural_analysis(ii,i,graphs,ARG_meta)# for ii,i in enumerate(relgene.columns[cc:dd]))


# structural_analysis(ii,i,graphs,ARG_meta) for ii,i in enumerate(relgene.columns[3:])


# pd.DataFrame([N,Q,I,R]).to_csv('Emore_ARG_nested_an.txt',sep='\t')