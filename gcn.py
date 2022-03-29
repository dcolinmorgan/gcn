import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
from tqdm import tqdm
import pandas as pd
import seaborn as sns
# from scipy import stats
import networkx as nx
# import networkit as nk
from networkx.algorithms import bipartite
from joblib import Parallel
from gcn_func import bip, load_list_of_dicts, meas,structural_analysis
# import argparse
from gcn_func import dgtz, shan_entropy, get_red, calc_MI, calc_CMI, calc_II, puc_cal


# parser = argparse.ArgumentParser(description='Process some integers.')
# parser.add_argument('dr', type=int,
#                  help='subsample fraction used to shuffle network')
# args = parser.parse_args()
# patt='all'

# if patt=='all':
#     relgene=pd.read_csv('~/data/gcn/ht_genefamilies-cpm.tsv',sep='\t')
# # elif patt=='other':
# #     relgene=pd.read_csv('~/ht_subset_genefamilies-cpm.tsv',sep='\t')
# elif patt=='50': ##humann2, others are humann3
#     # relgene=pd.read_csv('/groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/50_genefamilies.tsv',sep='\t')
#     relgene=pd.read_csv('~/ht_50_genefamilies-cpm.tsv',sep='\t')
#     patt=str(patt)
    
# if not os.path.isfile('data/gcn/relgene_'+patt+'.txt'):
#     relgene['gene']=relgene['# Gene Family'].str.split('|').str[0]
#     relgene=relgene[relgene['gene']!='UniRef90_unknown']
#     relgene=relgene[relgene['gene']!='UNMAPPED']
#     relgene.index=relgene['# Gene Family']
#     del relgene['gene'], relgene['# Gene Family']
#     # relgene=relgene/relgene.sum(axis=0)
#     relgene['gen']=relgene.index.str.split('|').str[1].str.split('.').str[0].tolist()
#     relgene['spec']=relgene.index.str.split('.').str[1]#.str.split('.').str[0].tolist()
#     relgene['spec'].replace('_',' ')
#     relgene.index=relgene.index.str.split('|').str[0]
#     relgene=relgene.dropna()
#     relgene.rename(columns={
#     'R243-LSF-DNA_Abundance-RPKs':'R0243-LSF-DNA_Abundance-RPKs',
#      'R247-LSY-DNA_Abundance-RPKs':'R0247-LSY-DNA_Abundance-RPKs',
#      'R253-HHF-DNA_Abundance-RPKs':'R0253-HHF-DNA_Abundance-RPKs',
#      'R256-TSW-DNA_Abundance-RPKs':'R0256-TSW-DNA_Abundance-RPKs'},inplace=True)
#     # relgene.columns[relgene.columns.str.contains('RPKs')].to_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
#     relgene[1:2].to_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
# else:
#     relgene=pd.read_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
    
# if not os.path.isfile('data/gcn/cc_'+patt+'.txt'):
#     cc=relgene.groupby(['# Gene Family','spec']).sum()
#     # dd=relgene.groupby(['# Gene Family','gen']).sum()
#     cc=cc.reset_index()
#     # dd=dd.reset_index()
#     cc=cc.rename(columns={'# Gene Family':'gene'})#,'spec':0,'gene':1})
#     cc.to_csv('data/gcn/cc_'+patt+'.txt',sep='\t')
# else:
#     cc=pd.read_csv('data/gcn/cc_'+patt+'.txt',sep='\t')

# primary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Primary Data')
# uni_bact=primary[['Age','Hypertension Category by 24h BP w/o considering antihypertensive med']]
# uni_bact=uni_bact.rename(columns={"Hypertension Category by 24h BP w/o considering antihypertensive med": "HT"})

# ff=[]
# C=[]
# dd=relgene.columns[relgene.columns.str.contains('RPKs')]
# # for i,net in enumerate(dd):
#     # ff,C=bip(cc,net,ff,C,patt)

# # cc.columns.str.contain()
# if not os.path.isfile('data/gcn/fBX_'+patt+'_HT.pkl'):
#     ff=Parallel(n_jobs=48) (bip(cc,net,ff,C,patt) for i,net in enumerate(dd))
#     C=[tup[1] for tup in ff]
#     ff=[tup[0] for tup in ff]
#     # xx=[tup[2] for tup in ff]
#     with open('data/gcn/NX_'+str(patt)+'_HT.pkl', 'ab+') as f:
#         pickle.dump(ff, f)
#     with open('data/gcn/BX_'+str(patt)+'_HT.pkl', 'ab+') as f:
#         pickle.dump(C, f)
#     # with open('data/gcn/fBX_'+str(patt)+'_HT.pkl', 'ab+') as f:
#     #     pickle.dump(xx, f)

# # graphs = load_list_of_dicts('data/gcn/BX_'+patt+'_hypert.pkl')
# # relgene=pd.read_csv('data/gcn/relgene'+patt+'.txt',sep='\t',nrows=1)

# # def plott(measur,uni_bact,relgene,graphs):
#     # measur='nx.betweenness_centrality'
# # df=meas(measur,uni_bact,relgene,graphs)
# # plotVio(df,'NoHT',measur,patt)
# dd=[#'nx.harmonic_centrality',
#     # 'nx.communicability',
#     'nx.betweenness_centrality','nx.degree_centrality','nx.degree','nx.closeness_centrality','nx.node_redundancy']
# # df=Parallel(n_jobs=7) (meas(measur,uni_bact,relgene,graphs,patt) for measur in dd)
# # C=[tup[1] for tup in ff]
# # df=[tup[0] for tup in df]
    
# # measur='nx.degree_centrality'
# # df=meas(measur,uni_bact,relgene,graphs)
# # plotVio(df,'NoHT',measur,patt)
# # plotVio(df,'HT',measur,patt)

# # measur='nx.closeness_centrality'
# # df=meas(measur,uni_bact,relgene,graphs)
# # plotVio(df,'NoHT',measur,patt)
# # plotVio(df,'HT',measur,patt)

# # measur='nx.node_redundancy'
# # df=meas(measur,uni_bact,relgene,graphs)
# # plotVio(df,'NoHT',measur,patt)
# # plotVio(df,'HT',measur,patt)


# patt='all'
# sys.path.insert(1, './run/gcn/')
# import gcn_func
# importlib.reload(sys.modules['gcn_func'])
# from gcn_func import bip, load_list_of_dicts, meas, time_bar,proc_dat,rev_tbar,group_time_plot,time_order_net,build_gcn,shuffle_net,structural_analysis,makeSYNCSAnet,group4SYNCSA
# JEFF=[]
# DR=sys.argv[1]

# relgene=pd.read_csv('all_arg_subset_genefamilies-cpm.tsv',sep='\t',nrows=1)
# graphs = load_list_of_dicts('data/gcn/NX_Emore_ARG.pkl')

# ARG_meta=pd.read_excel('run/gcn/ARG_treatment_infor_modified.xlsx',index_col=0)
# ARG_meta2=pd.read_excel('run/gcn/patients_Tx_batch3_for_DM.xlsx',index_col=None,skiprows=1,names=['id','group'])
# relgene.columns=relgene.columns.str.replace("-00", "-00ST")
# relgene.columns=relgene.columns.str.replace("-00STST", "-00ST")
# relgene.columns=relgene.columns.str.split('-').str[0]+'-'+relgene.columns.str.split('-').str[1]
# ARG_meta['id']=ARG_meta['id'].str.replace('-00ST','')
# META=pd.concat([pd.DataFrame(ARG_meta[['id','group']]),ARG_meta2],ignore_index=True)


# from joblib import Parallel
# dd=Parallel(n_jobs=20) (makeSYNCSAnet(ii,i,graphs,JEFF,deg_rand=DR) for ii,i in enumerate(relgene.columns[1:-1]))
# dd=pd.DataFrame()
# ff=pd.DataFrame()
# for ii,i in enumerate(relgene.columns[1:-1]):

# dd=makeSYNCSAnet(relgene,graphs,JEFF,META,deg_rand=DR)
# dd.to_csv('~/SYNCSA_eval/'+str(DR)+'_dd_.csv')

# dd=pd.read_csv('SYNCSA_eval/'+str(DR)+'_dd_.csv',sep=',',index_col=0)

# names=pd.unique(dd.columns.str.split('_').str[1]+'_'+dd.columns.str.split('_').str[2])[1:]

# # Parallel(n_jobs=10) (
# for i in (names):
#     group4SYNCSA(i,dd,DR)

# sys.path.insert(1, './run/gcn/')
# import gcn_func
# importlib.reload(sys.modules['gcn_func'])
# from gcn_func import research_orthologs

# norm_Px=pd.read_csv('data/Pipeline_consolidate_220301/Norm_Px.txt',sep='\t')
# PS=pd.DataFrame(columns=['mus','homo','spec'])
# for i in norm_Px['Protein']:
#     P,S=research_orthologs(str(i),'Homo sapiens')
#     cc = pd.DataFrame({'mus': i,
#                            'homo': P,
#                            'spec': S}, 
#                            index = [0])
#     PS=PS.append(cc)
#     cc.to_csv('mus_homo_int.txt',sep='\t',mode='a')
# PS.to_csv('mus_homo.txt',sep='\t')


# selection=MS['uniprot_id']
# jeff_ind=mg.querymany(selection, scopes='uniprot', species='human',as_dataframe=True)
# jeff_ind.to_csv('jeff_indtmp.txt',sep='\t')

import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle,tqdm,math
from itertools import combinations
from tqdm import tqdm
from IPython.display import Image
import pandas as pd
import seaborn as sns
from scipy import stats
import networkx as nx
from pathlib import Path
from joblib import Parallel
from statannotations.Annotator import Annotator
from scipy.stats import zscore
import mygene
mg = mygene.MyGeneInfo()

os.chdir('/home/dcmorgan')
os.getcwd()

# from https://hmdb.ca/
# from https://hmdb.ca/
# norm_MS=pd.read_csv('data/Pipeline_consolidate_220301/raw/metabolites.master-raw.csv')
# norm_MS=pd.read_csv('data/Pipeline_consolidate_220301/output/data/normalised.eigenMS_median.csv')
norm_MS=pd.read_csv('data/Pipeline_consolidate_220301/output/data/normalised.eigenMS_median.csv')

MS_map=pd.read_csv('data/Pipeline_consolidate_220301/raw/metabolites.annotation.master-NetID.csv')
MS_map2=pd.read_csv('data/Pipeline_consolidate_220301/test.txt',index_col=0)
prot=MS_map2.merge(MS_map,left_on='accession',right_on='HMDB')
len(np.unique(prot['HMDB']))
MS=norm_MS.merge(prot[['ColID','uniprot_id']])
# MS_map=pd.read_csv('data/Pipeline_consolidate_220301/raw/metabolites.annotation.master-NetID.csv')
# MS_map2=pd.read_csv('data/Pipeline_consolidate_220301/test.txt',index_col=0)
# prot=MS_map2.merge(MS_map,left_on='accession',right_on='HMDB')
# len(np.unique(prot['HMDB']))
# MS=norm_MS.merge(prot[['ColID','uniprot_id']])
del norm_MS,prot,MS_map,MS_map2

jeff_ind=pd.read_csv('jeff_indtmp.txt',sep='\t')
jeff_ind=jeff_ind.reset_index()[['query','symbol']]
del MS['mz'],MS['rt'],MS['ColID'],MS['negpos']
MS=MS.reset_index().merge(jeff_ind,left_on='uniprot_id',right_on='query')
np.max(MS,axis=1)

Lx=MS[MS['exp']=='Lx']
SP=Lx.groupby('symbol').mean().filter(regex='SP')
SP=SP.loc[~(SP==0).all(axis=1)]
SP.to_csv('data/Pipeline_consolidate_220301/SP_gene_sum_Lx.txt',sep='\t',header=False)

UP=Lx.groupby('symbol').mean().filter(regex='UP')
UP=UP.loc[~(UP==0).all(axis=1)]
UP.to_csv('data/Pipeline_consolidate_220301/UP_gene_sum_Lx.txt',sep='\t',header=False)

LF=Lx.groupby('symbol').mean().filter(regex='LF')
LF=LF.loc[~(LF==0).all(axis=1)]
LF.to_csv('data/Pipeline_consolidate_220301/LF_gene_sum_Lx.txt',sep='\t',header=False)


Mx=MS[MS['exp']=='Mx']
SP=Mx.groupby('symbol').mean().filter(regex='SP')
SP=SP.loc[~(SP==0).all(axis=1)]
SP.to_csv('data/Pipeline_consolidate_220301/SP_gene_sum_Mx.txt',sep='\t',header=False)

UP=Mx.groupby('symbol').mean().filter(regex='UP')
UP=UP.loc[~(UP==0).all(axis=1)]
UP.to_csv('data/Pipeline_consolidate_220301/UP_gene_sum_Mx.txt',sep='\t',header=False)

LF=Mx.groupby('symbol').mean().filter(regex='LF')
LF=LF.loc[~(LF==0).all(axis=1)]
LF.to_csv('data/Pipeline_consolidate_220301/LF_gene_sum_Mx.txt',sep='\t',header=False)
