import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
# from tqdm import tqdm
import pandas as pd
import seaborn as sns
# from scipy import stats
from scipy.stats import rankdata
import networkx as nx
# import networkit as nk
from networkx.algorithms import bipartite
from joblib import Parallel
from gcn_func import bip, load_list_of_dicts, meas, plotRidge, time_bar, rev_tbar, proc_dat

primary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Primary Data')
diet=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Diet Data')
blood_stool=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='blood and stool biomarkers')
secondary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Secondary Data')
MRI=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='MRI scores')
uni_bact=primary[['Age','Hypertension Category by 24h BP w/o considering antihypertensive med']]
uni_bact=uni_bact.rename(columns={"Hypertension Category by 24h BP w/o considering antihypertensive med": "HT"})
# uni_bact.to_csv('data/gcn/uni_bact.txt',sep='\t')
patt='all'
relgene=pd.read_csv('data/gcn/relgene_all.txt',sep='\t')
graphs = load_list_of_dicts('data/gcn/BX_'+patt+'_HT.pkl')
HTXX=uni_bact[uni_bact.index.isin(relgene.columns[1:-2].str.split('-').str[0])]
HTXX['index']=np.arange(len(HTXX))


# for i,net in tqdm.tqdm(enumerate(BX_graphs)):
for i,net in (enumerate(HTXX[HTXX['HT']!=5]['index'].values)):
  cc=nx.convert_matrix.to_pandas_edgelist(graphs[i])
  # cc['weight']=np.random.randn(len(cc))
  rrr=str(HTXX[HTXX['index']==i]['Age'].item())+'_'+str(HTXX[HTXX['index']==i]['HT'].item())#+'_'+str(HTXX[HTXX['index']==i]['sex'].item())
  cc.rename(columns={cc.columns[2]:rrr},inplace=True)
  if i==0:
    dd=cc
  else:
    dd=dd.merge(cc,on=['source','target'],how='outer')
# dd.dropna(how='any')
   
    
    
dd.set_index(['source', 'target'], inplace=True) #>> run only first time editing dd
# dd = dd/np.max(dd,axis=0)
# dd=dd/np.sum(dd,axis=0)
### dd=np.argsort(dd)



# time_bar(noHT,5,'rank','all')
# time_bar(noHT,5,'value','all')
# time_bar(HT1,5,'rank','all')
# time_bar(HT1,5,'value','all')
# time_bar(HT2,5,'rank_diff','all')
# time_bar(HT2,5,'diff','all')

noHT=dd.filter(regex='_0_').dropna(how='all')
noHT=proc_dat(noHT)
rev_tbar(noHT,10,'noHT')

HT1=dd.filter(regex='_1_').dropna(how='all')
HT1=proc_dat(HT1)
rev_tbar(HT1,10,'HT1')

HT2=dd.filter(regex='_2_').dropna(how='all')
HT2=proc_dat(HT2)
rev_tbar(HT2,10,'HT2')
    
# Parallel(n_jobs=3) (rev_tbar(data,10) for data in ([noHT,HT1,HT2]))