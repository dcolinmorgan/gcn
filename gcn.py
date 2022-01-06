import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
# from tqdm import tqdm
import pandas as pd
import seaborn as sns
# from scipy import stats
import networkx as nx
# import networkit as nk
from networkx.algorithms import bipartite
from joblib import Parallel
from gcn_func import bip, load_list_of_dicts, meas, plotRidge

patt='all'

if patt=='all':
    relgene=pd.read_csv('~/data/gcn/ht_genefamilies-cpm.tsv',sep='\t')
# elif patt=='other':
#     relgene=pd.read_csv('~/ht_subset_genefamilies-cpm.tsv',sep='\t')
elif patt=='50': ##humann2, others are humann3
    # relgene=pd.read_csv('/groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/50_genefamilies.tsv',sep='\t')
    relgene=pd.read_csv('~/ht_50_genefamilies-cpm.tsv',sep='\t')
    patt=str(patt)
    
if not os.path.isfile('data/gcn/relgene_'+patt+'.txt'):
    relgene['gene']=relgene['# Gene Family'].str.split('|').str[0]
    relgene=relgene[relgene['gene']!='UniRef90_unknown']
    relgene=relgene[relgene['gene']!='UNMAPPED']
    relgene.index=relgene['# Gene Family']
    del relgene['gene'], relgene['# Gene Family']
    # relgene=relgene/relgene.sum(axis=0)
    relgene['gen']=relgene.index.str.split('|').str[1].str.split('.').str[0].tolist()
    relgene['spec']=relgene.index.str.split('.').str[1]#.str.split('.').str[0].tolist()
    relgene['spec'].replace('_',' ')
    relgene.index=relgene.index.str.split('|').str[0]
    relgene=relgene.dropna()
    relgene.rename(columns={
    'R243-LSF-DNA_Abundance-RPKs':'R0243-LSF-DNA_Abundance-RPKs',
     'R247-LSY-DNA_Abundance-RPKs':'R0247-LSY-DNA_Abundance-RPKs',
     'R253-HHF-DNA_Abundance-RPKs':'R0253-HHF-DNA_Abundance-RPKs',
     'R256-TSW-DNA_Abundance-RPKs':'R0256-TSW-DNA_Abundance-RPKs'},inplace=True)
    # relgene.columns[relgene.columns.str.contains('RPKs')].to_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
    relgene[1:2].to_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
else:
    relgene=pd.read_csv('data/gcn/relgene_'+patt+'.txt',sep='\t')
    
if not os.path.isfile('data/gcn/cc_'+patt+'.txt'):
    cc=relgene.groupby(['# Gene Family','spec']).sum()
    # dd=relgene.groupby(['# Gene Family','gen']).sum()
    cc=cc.reset_index()
    # dd=dd.reset_index()
    cc=cc.rename(columns={'# Gene Family':'gene'})#,'spec':0,'gene':1})
    cc.to_csv('data/gcn/cc_'+patt+'.txt',sep='\t')
else:
    cc=pd.read_csv('data/gcn/cc_'+patt+'.txt',sep='\t')

primary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Primary Data')
uni_bact=primary[['Age','Hypertension Category by 24h BP w/o considering antihypertensive med']]
uni_bact=uni_bact.rename(columns={"Hypertension Category by 24h BP w/o considering antihypertensive med": "HT"})

ff=[]
C=[]
dd=relgene.columns[relgene.columns.str.contains('RPKs')]
# for i,net in enumerate(dd):
    # ff,C=bip(cc,net,ff,C,patt)

# cc.columns.str.contain()
if not os.path.isfile('data/gcn/BX_'+patt+'_hypert.pkl'):
    ff=Parallel(n_jobs=48) (bip(cc,net,ff,C,patt) for i,net in enumerate(dd))
    C=[tup[1] for tup in ff]
    ff=[tup[0] for tup in ff]
    with open('data/gcn/NX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
        pickle.dump(ff, f)
    with open('data/gcn/BX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
        pickle.dump(C, f)

graphs = load_list_of_dicts('data/gcn/BX_'+patt+'_hypert.pkl')
# relgene=pd.read_csv('data/gcn/relgene'+patt+'.txt',sep='\t',nrows=1)

# def plott(measur,uni_bact,relgene,graphs):
    # measur='nx.betweenness_centrality'
# df=meas(measur,uni_bact,relgene,graphs)
# plotVio(df,'NoHT',measur,patt)
dd=[['nx.harmonic_centrality','nx.communicability','nx.betweenness_centrality','nx.degree_centrality','nx.closeness_centrality','nx.node_redundancy']]
df=Parallel(n_jobs=7) (meas(measur,uni_bact,relgene,graphs,patt) for measur in dd)
# C=[tup[1] for tup in ff]
df=[tup[0] for tup in df]
    
# measur='nx.degree_centrality'
# df=meas(measur,uni_bact,relgene,graphs)
# plotVio(df,'NoHT',measur,patt)
# plotVio(df,'HT',measur,patt)

# measur='nx.closeness_centrality'
# df=meas(measur,uni_bact,relgene,graphs)
# plotVio(df,'NoHT',measur,patt)
# plotVio(df,'HT',measur,patt)

# measur='nx.node_redundancy'
# df=meas(measur,uni_bact,relgene,graphs)
# plotVio(df,'NoHT',measur,patt)
# plotVio(df,'HT',measur,patt)