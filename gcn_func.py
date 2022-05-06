import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
# from tqdm import tqdm
from scipy.stats import rankdata
from functools import reduce
import pandas as pd
import networkx as nx
import seaborn as sns
from joblib import delayed, wrap_non_picklable_objects
from pathlib import Path
import plotly
from networkx import bipartite
from numba import jit
from joblib import Parallel
import sklearn.utils as sku
import plotly.graph_objects as go
import plotly.express as px
# j=sys.argv[1]
from urllib import request
import xml.etree.ElementTree as ET
import urllib
sys.path.insert(1, './nestedness_analysis/')
import nestedness_metrics_other_functions
from nestedness_metrics_other_functions import from_edges_to_matrix
# importlib.reload(sys.modules['EO_functions_bipartite'])
import extremal_bi

@delayed
@wrap_non_picklable_objects

def bip(cc,net,ff,C,patt):

    # print(net)
    # dd=cc[['spec','gene',net]]
    dd=pd.read_csv('data/gcn/cc_'+patt+'.txt',index_col=False,sep='\t',usecols=['spec','gene',net])
    # try:
    dd=dd[dd[net]!=0]
    # except:
    #     pass
    # ee=nx.from_pandas_edgelist(dd,source='spec',target='gene')
    # remove = [node for node,degree in dict(ee.degree()).items() if degree <5]
    # ee.remove_nodes_from(remove)
    # ff.append(ee)
    
    B = nx.Graph()
    B.add_nodes_from(dd['spec'], bipartite=0)
    B.add_nodes_from(dd['gene'], bipartite=1)
    B.add_weighted_edges_from(tuple(dd[['spec','gene',net]].itertuples(index=False, name=None)))
    remove = [node for node,degree in dict(B.degree()).items() if degree <5]
    B.remove_nodes_from(remove)
    # C.append(B)
    
    xx=nx.from_pandas_edgelist(dd,source='spec',target='gene',edge_attr=net)
    remove = [node for node,degree in dict(xx.degree()).items() if degree <5]
    xx.remove_nodes_from(remove)

    # with open('data/gcn/NX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
    #     pickle.dump(ff, f)
    # with open('data/gcn/BX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
    #     pickle.dump(C, f)
    return xx,B

def load_list_of_dicts(filename, create_using=nx.Graph):
    with open(filename, 'rb') as f:
        list_of_dicts = pickle.load(f)   
    graphs = [create_using(graph) for graph in list_of_dicts]
    return graphs


# @delayed
# @wrap_non_picklable_objects
def meas(measur,uni_bact,relgene,graphs,patt):
    HTXX=uni_bact[uni_bact.index.isin(relgene.columns[1:-2].str.split('-').str[0])]
    HTXX['index']=np.arange(len(HTXX))
    # measur=eval(measur)
    S = [eval(measur)(graphs[i]) for i in HTXX[HTXX['HT']==0]['index'].values]
    T = [eval(measur)(graphs[i]) for i in HTXX[HTXX['HT']!=0]['index'].values]
    
    if measur!='nx.degree':
        non=pd.DataFrame(S).melt()
        yes=pd.DataFrame(T).melt()
    elif measur=='nx.degree':
        non=pd.DataFrame(S.pop())
        non=non.rename(columns={0:'variable',1:'value'})
        yes=pd.DataFrame(T.pop())
        yes=yes.rename(columns={0:'variable',1:'value'})
    non['type']='NoHT'
    non.dropna(inplace=True)
    non=non[non.value!=0]
    non=non[~non['variable'].str.contains('UniRef90')]
    non.value=non.value/np.sum(non.value)
    yes['type']='HT'
    yes.dropna(inplace=True)
    yes=yes[yes.value!=0]
    yes=yes[~yes['variable'].str.contains('UniRef90')]
    yes.value=yes.value/np.sum(yes.value)
    df=non.append(yes)
    # df=df.dropna()
    df['gen']=df.variable.str.split('_').str[2]
    df.to_csv("data/gcn/"+patt+"_"+str(measur)+".txt",sep='\t')

    plt.figure(figsize=(10,30))
    sns.set_theme(style="whitegrid")

    sns.violinplot(data=df, y="gen", x="value",hue="type",
                   split=True, inner="quart", linewidth=1,
                   orient="h")
    sns.despine(left=True)
    plt.savefig("data/gcn/"+patt+"_"+str(measur)+"_violin.png",dpi=300,bbox_inches = "tight")
    return df



def time_bar(data,XX,rank='rank',species='all'):
    if rank=='rank':
        data['rank']=rankdata(data.value,method='min')
    elif rank=='rank_diff' or rank=='diff':
        data['vx']=rankdata(data.value_x,method='min')
        data['vy']=rankdata(data.value_y,method='min')
        data['rank_diff']=data['vx'].astype('int')-data['vy'].astype('int')
        data['diff']=data['value_x']-data['value_y']
    # elif rank=='value':
        # rank=data.value
        
    if species!='all':
        data=data[data['species']==species]
    # clust = ll.groupby(['species','target','time'], as_index=False)['diff'].sum()
    
    df = data[['species','target','time',rank]]#.sort_values(['time'], ascending=[True]).groupby(['species','time']).max(5)
    jeff=pd.DataFrame(df.groupby(['species','time'])[rank].nlargest(XX))
    jeff.reset_index(inplace=True)


    for cc in np.unique(jeff.species):
        jeff2=jeff[jeff['species']==cc]
        if species!='all':
            jeff2=df.loc[jeff2['level_2']]
        else:
            jeff2=df.iloc[jeff2['level_2']]
        plt.figure(figsize=(15,5))
        ax = sns.histplot(jeff2, x='time', hue='target', weights=rank,
                     multiple='stack', palette='icefire', shrink=0.6,bins=len(pd.unique(jeff2.time))+5)
        ax.set_ylabel(str(rank)+'_HT')
        ax.set_title(cc)
        # Fix the legend so it's not on top of the bars.
        # legend = ax.get_legend()
        plt.legend([],[], frameon=False)
        Path("data/gcn/img/"+cc).mkdir(parents=True, exist_ok=True)
        plt.savefig("data/gcn/img/"+cc+"/"+str(data)+"_"+cc+"_"+str(rank)+".png",dpi=300,bbox_inches = "tight")

def proc_dat(noHT):
    # noHT=jj.filter(regex=str(focus)).dropna(how='all')
    noHT.columns=noHT.columns.str.split('_').str[0]
    noHT=noHT.groupby(by=noHT.columns, axis=1).mean()
    noHT=noHT.dropna(how='any')
    noHT.reset_index(inplace=True)
    jj=noHT.melt(['source','target'])
    jj.rename(columns={'variable':'time'},inplace=True)
    jj['t']=jj['time']
    # jj['time']=jj['time'].astype('int')+2000
    # jj['time'] = pd.to_datetime(jj['time'], format='%Y')
    # jj=jj[jj['value']>5]
    jj['species']=jj['source'].str.split('_').str[2]
    jj=jj.dropna(how='any')
    return jj
    
# @delayed
# @wrap_non_picklable_objects
def rev_tbar(jj,XX,gg,species='all'):
    
    data=jj[['species','target','time','t','value']]
    # df=data.copy()
    # data.reset_index(inplace=True)
    data['sum']=pd.DataFrame(data.groupby(['species','t','target'])['value'].transform('sum'))
    # jeff.reset_index(inplace=True)
    del data['value']
    data.drop_duplicates(inplace=True)
    data.reset_index(inplace=True)
    del data['index'],data['time']
    jeff=pd.DataFrame(data.groupby(['species','t'])['sum'].nlargest(XX))
    jeff.reset_index(inplace=True)

    jeffA=data.iloc[jeff['level_2']]
    tim_len=len(np.unique(jeffA['t']))
    if species!='all':
        jeff=jeff[jeff['species']==species]
    JJ=pd.DataFrame()
    rr=[]
    for q,ee in enumerate((np.unique(jeff.species))):
        jeff2=jeffA[jeffA['species']==ee]#.explode('target')
        dd=pd.DataFrame(jeff2['target'].to_numpy().reshape(int(len(jeff2)/tim_len),tim_len,order='F'))

        if len(dd.melt())==(tim_len*XX):
            JJ=JJ.append(dd)
            rr=np.append(rr, ee)
    jeffA=jeffA.sort_values(['species', 't'], ascending=[True, True])      
    labels,levels=pd.factorize(sku.shuffle(JJ.melt()['value']))
    cc=pd.DataFrame(np.array(labels).reshape((XX)*len(rr),tim_len,order='F'))

    for i in np.arange(0,len(cc),XX+1):
        for col in cc:
            cc.iloc[i:i+XX,col] = cc.iloc[i:i+XX,col].sort_values(ignore_index=True)
            cc.loc[i+XX]=0

    plt.figure(figsize=(10,30))
    ax=sns.heatmap(cc,cmap='rocket_r',annot=True, fmt="d",cbar=False,xticklabels=False,
                   yticklabels=False).set(ylabel='         -         '.join(rr))

    # plt.show()
    data.to_csv('data/gcn/'+str(gg)+'.csv',sep='\t')
    # Path("data/gcn/img/"+cc).mkdir(parents=True, exist_ok=True)
    plt.savefig("data/gcn/img/full_"+str(gg)+"_10.png",dpi=300,bbox_inches = "tight")


def group_time_plot(noHT,steps,XX,spec_spec):
        
    noHT.columns=noHT.columns.str.split('_').str[0]
    noHT.columns=pd.qcut((noHT.columns).astype('int'), steps, labels=False)
    noHT=noHT.groupby(by=noHT.columns, axis=1).mean()
    noHT=noHT.dropna(how='all')


    noHT.reset_index(inplace=True)
    jj=noHT.melt(['source','target'])
    jj.rename(columns={'variable':'time'},inplace=True)
    jj['t']=jj['time']
    # jj['time']=jj['time'].astype('int')+2000
    # jj['time'] = pd.to_datetime(jj['time'], format='%Y')
    # jj=jj[jj['value']>5]
    jj['species']=jj['source'].str.split('_').str[2]
    jj=jj.dropna(how='any')

    jj['rank']=rankdata(jj.value,method='min')

    XX=50 #10
    # df = noHT[['species','target','time','rank']]
    del jj['value'], jj['t'], jj['source']
    if spec_spec=='1':
        
        jeff=pd.DataFrame(jj.groupby(['species','time'])['rank_diff'].nlargest(XX))
        jeff=jeff.dropna(how='any')

        jeff.reset_index(inplace=True)
        jeff2=jj.loc[jeff['level_2']]
    else:
        jeff=pd.DataFrame(jj.groupby(['time'])['rank'].nlargest(XX))
        jeff=jeff.dropna(how='any')

        jeff.reset_index(inplace=True)
        jeff2=jj.loc[jeff['level_1']]

    plt.figure(figsize=(15,5))
    ax = sns.histplot(jeff2, x='time', hue='target', weights='rank',
                 multiple='stack', palette='icefire', shrink=0.6,bins=len(pd.unique(jeff2['time']))+5)
    ax.set_ylabel('rank_noHT')
    # ax.set_title(cc)
    # Fix the legend so it's not on top of the bars.
    # legend = ax.get_legend()
    plt.legend([],[], frameon=False)

def time_order_net(control,case,thresh=10**-6,group='source',groups=6,rounder=1,math='mean'):
    def preproc(data):
        data.columns=data.columns.str.split('_').str[0]
        data.columns=pd.qcut((data.columns).astype('int'), groups, labels=False)
        noHTm=data.groupby(by=data.columns, axis=1).mean()
        noHTm=noHTm.dropna(how='all')
        noHTm.reset_index(inplace=True)
        noHTv=data.groupby(by=data.columns, axis=1).var()
        noHTv=noHTv.dropna(how='all')
        noHTv.reset_index(inplace=True)
        return noHTm,noHTv
    noHTm,noHTv=preproc(control)
    HTm,HTv=preproc(case)
    
    if math=='mean':
        BB=noHTm[noHTm[0]>thresh].dropna().groupby(group).mean()-HTm[HTm[0]>thresh].dropna().groupby(group).mean()
    elif math=='median':
        BB=noHTm[noHTm[0]>thresh].dropna().groupby(group).median()-HTm[HTm[0]>thresh].dropna().groupby(group).median()
    BB=np.round(BB,rounder)
    aa='(BB[0]>='
    bb='(BB[0]<='
    for i in np.arange(groups)[1:]:
        cc='BB['+str(i)+'])&(BB['+str(i)+']>='
        aa=aa+str(cc)
        dd='BB['+str(i)+'])&(BB['+str(i)+']<='
        bb=bb+str(dd)
        
    grow=BB[eval(bb[:-9])]
    die=BB[eval(aa[:-9])]
    
    def proc_run(BBgrow,grow):
        if len(BBgrow)>0:
            BBgrow[groups]=BBgrow[0]-BBgrow[groups-1]
            BBgrow=BBgrow[BBgrow[groups]!=0]
            BBgrow.sort_values(by=groups,inplace=True)
            del BBgrow[groups]
            BBgrow.to_csv('data/gcn/comp_net/'+str(group)+'_'+str(thresh)+'_'+str(math)+'_'+str(groups)+'_'+grow+'.txt',sep='\t')
        else:
            BBgrow=0
        return BBgrow
    BBgrow=proc_run(grow,'grow')
    BBdie=proc_run(die,'die')
    return BBgrow,BBdie,noHTm,HTm


    
    
def build_gcn(i,net,cc,min_deg=5):

#     relgene=pd.read_csv(path,sep='\t')
#     # relgene=pd.read_csv('50_genefamilies-cpm.tsv')
#     # relgene=pd.read_csv('hmp_subset_genefamilies-cpm.tsv',sep='\t',nrows=100)

#     relgene['gene']=relgene['# Gene Family'].str.split('|').str[0]
#     relgene=relgene[relgene['gene']!='UniRef90_unknown']
#     relgene=relgene[relgene['gene']!='UNMAPPED']
#     relgene.index=relgene['# Gene Family']
#     del relgene['gene'], relgene['# Gene Family']
#     # relgene=relgene/relgene.sum(axis=0)

#     # relgene=relgene/relgene.sum(axis=0)
#     relgene['gen']=relgene.index.str.split('|').str[1].str.split('.').str[0].tolist()
#     relgene['spec']=relgene.index.str.split('.').str[1]#.str.split('.').str[0].tolist()
#     relgene['spec'].replace('_',' ')
#     relgene.index=relgene.index.str.split('|').str[0]
#     relgene=relgene.dropna()
#     cc=relgene.groupby(['# Gene Family','spec']).sum()
#     cc=cc.reset_index()
#     cc=cc.rename(columns={'# Gene Family':'gene'})

#     ff=[]
#     C=[]
    # for i,net in enumerate(relgene.columns[1:-2]):
        # pd.read_csv()
    dd=cc[['spec','gene',net]]
    dd=dd[dd[net]!=0]
    ee=nx.from_pandas_edgelist(dd,source='spec',target='gene',edge_attr=net)
    remove = [node for node,degree in dict(ee.degree()).items() if degree <min_deg]
    ee.remove_nodes_from(remove)
    # ff.append(ee)

    B = nx.Graph()
    B.add_nodes_from(dd['spec'], bipartite=0)
    B.add_nodes_from(dd['gene'], bipartite=1)
    B.add_edges_from(tuple(dd[['spec','gene']].itertuples(index=False, name=None)))
    remove = [node for node,degree in dict(B.degree()).items() if degree <min_deg]
    B.remove_nodes_from(remove)
    # C.append(B)
    return ee,B

    # with open('data/gcn/NX_Emore_'+name+'.pkl', 'wb') as f:
    #     pickle.dump(ff, f)
    # with open('data/gcn/BX_Emore_'+name+'.pkl', 'wb') as f:
    #     pickle.dump(C, f)
    

def buildSYNCSA(dd): 
    names=pd.unique(dd.columns.str.split('_').str[1]+'_'+dd.columns.str.split('_').str[2])[1:]
    for i in names:
        # ff.columns = ff.columns.str.strip('_x')
        # ff.columns = ff.columns.str.strip('_y')
        # i=i.split('_')[1]+'_'+i.split('_')[2]
        ff=dd.loc[:,dd.columns.str.contains(i)]
        ff[['source','target']]=dd[['source','target']]
        # ff=ff[ff['source'].str.contains('s__')]
        # ff=ff[ff['target'].str.contains('UniRef')]
        ff.source, ff.target = np.where(ff.source.str.contains('UniRef'), [ff.target, ff.source], [ff.source, ff.target])

        ff.groupby('source').sum().transpose().to_csv('comm_'+i+'.csv')
        ff.reset_index(inplace=True)
        ff.set_index(['source', 'target'], inplace=True)
        del ff['index']
        ff.columns=(ff.columns.str.split('_').str[1]+'_'+ff.columns.str.split('_').str[2])
        gg=ff.groupby(by=ff.columns, axis=1).sum()
        traits=gg[[i]].reset_index().pivot('source','target',i).dropna(how='all',axis=1).replace(np.nan,0)
        traits.to_csv('trait_'+i+'.csv')
    

def buildNestedNess():
    C=pd.DataFrame(columns=['N','Q','I','type'])
    D=[]
    files=glob.glob('*.npz')
    for i,j in enumerate(files):
        d=np.load(j)
        C.loc[i]=[float(d['N']),float(d['Q']),float(d['I']),j.split('_')[1]+'_'+j.split('_')[2]+'_'+j.split('_')[3].split('.')[0]]
        
def structural_analysis(ii,i,graphs,ARG_meta,rand,deg_rand):
    # aa= rrr[['from','to','value']].values
    ccc=nx.convert_matrix.to_pandas_edgelist(graphs[ii])
    ee=nx.convert_matrix.to_pandas_edgelist(graphs[ii])
    # cc['weight']=np.random.randn(len(cc))
    pww=i
    j=(i.split('-')[1])
    i=(i.split('-')[0])
    rrr=str(ARG_meta[ARG_meta['id']==i].index.item())+'_'+str(ARG_meta[ARG_meta['id']==i]['group'].item())+'_'+str(j)
    ccc.rename(columns={ccc.columns[2]:rrr},inplace=True)
    
#     a,b=pd.factorize(ccc['source']) 
#     c,d=pd.factorize(ccc['target'])
#     rrr=pd.DataFrame()
    
#     rrr['from']=a
#     rrr['to']=c
#     rrr['value']=1
    sss=str(ARG_meta[ARG_meta['id']==i]['group'].item())+'_'+str(j)
    Path('nest/'+sss).mkdir(parents=True, exist_ok=True)
    # rrr[['from','to','value']].to_csv('~/nest/'+sss+'/'+str(ccc.columns[2])+'.csv',sep=' ',index=False,header=False)
    # ccc.to_csv('~/nest/'+sss+'/'+str(ccc.columns[2])+'.txt',sep='\t',index=False,header=False)
    # np.transpose(pd.DataFrame([a,b])).to_csv('~/nest/'+sss+'/source_'+str(pww)+'.txt',sep='\t',index=False,header=False)
    # np.transpose(/'pd.DataFrame([c,d])).to_csv('~/nest/'+sss+'/target_'+str(pww)+'.txt',sep='\t',index=False,header=False)
    # aa= ccc[['from','to','value']].values
    aa=ccc.values
    
    if rand==True: ## to randomize
        aa=pd.DataFrame(aa)
        ddd=aa.sample(frac=np.float(deg_rand), replace=False, random_state=1) ##degree randomized
        rrr=aa[~aa.isin(ddd)].dropna(how='all')
        ddd.reset_index(inplace=True)
        del ddd['index']
        ssst=shuffle(ddd)
        aa=pd.concat([rrr,ssst])
        aa=np.array(aa).astype(int)
        ccc.values=aa
    # nodes_cols = int(max(aa[j,1] for j in range(aa.shape[0]))+1)
    # nodes_rows= int(max(aa[j,0] for j in range(aa.shape[0]))+1)
    # matrix=np.zeros((nodes_rows,nodes_cols),dtype='int')
    # matrix=np.zeros((len(np.unique(ccc['source'])),len(np.unique(ccc['target']))),dtype='int')
    # for j in range(aa.shape[0]):
    #     matrix[aa[j,0],aa[j,1]] = 1
    # M=matrix
    M=ccc.pivot('source',columns='target',values=ccc.columns[2]).fillna(0)#.round().astype('int').to_numpy())
    # BB=ccc.pivot('source',columns='target',values=ccc.columns[2]).fillna(0).round().astype('int')#.to_numpy())
    
    M=((M>0)*1).to_numpy()
    # M.to_numpy()
    with open('~/nest/'+sss+'/net_'+str(pww)+'.txt', "ab") as f:
        np.savetxt(f,M,delimiter = '\t')
    
    cols_degr=M.sum(axis=0)
    row_degr=M.sum(axis=1)
    R,C=M.shape #rows and cols
    #Nestednes
    # In-block nestedness with B=1
    Cn_=[np.repeat(1, R),np.repeat(1, C)]
    max_blockN=max(max(Cn_[0]),max(Cn_[1]))+1
    lambdasN=extremal_bi.call_lambda_i(M,cols_degr,row_degr,Cn_[1],Cn_[0],max_blockN,True)
    N=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasN[0],lambdasN[1],True)

    #Modularity Extremal
    C_=extremal_bi.recursive_step(M,cols_degr,row_degr,.7,3,False)
    max_blockQ=max(max(C_[0]),max(C_[1]))+1
    lambdasQ=extremal_bi.call_lambda_i(M,cols_degr,row_degr,C_[1],C_[0],max_blockQ,False)
    Q=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasQ[0],lambdasQ[1],False)

    # Inblock nestedness extremal
    Ci_=extremal_bi.recursive_step(M,cols_degr,row_degr,.7,3,True)
    max_blockI=max(max(Ci_[0]),max(Ci_[1]))+1
    lambdasI=extremal_bi.call_lambda_i(M,cols_degr,row_degr,Ci_[1],Ci_[0],max_blockI,True)
    I=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasI[0],lambdasI[1],True)
    zzz=[str(N),str(Q),str(I),sss]
    print(zzz)
    # np.savetxt('ARG_nest_test.txt', zzz, delimiter = '\t', fmt='%s')
    dfq=pd.DataFrame({'rows': pd.Series(C_[0]), 'cols': pd.Series(C_[1])})
    dfi=pd.DataFrame({'rows': pd.Series(Ci_[0]), 'cols': pd.Series(Ci_[1])})
    # dfn=pd.DataFrame({'rows': pd.Series(Cn_[0]), 'cols': pd.Series(Cn_[1])})
    dfq.to_csv('~/nest/'+sss+'/modularity_'+str(pww)+'.txt',mode='a',sep='\t',index=False,header=False)
    dfi.to_csv('~/nest/'+sss+'/in-block_'+str(pww)+'.txt',mode='a',sep='\t',index=False,header=False)
    # dfn.to_csv('~/nest/'+sss+'/nested_'+str(pww)+'.txt',mode='a',sep='\t',index=False,header=False)
    with open("~/nest/randC_ARG_nest.txt", "ab") as f:
        np.savetxt(f,np.column_stack([str(N),str(Q),str(I),sss,pww]),delimiter = '\t', fmt='%s')
    return N,Q,I,sss,pww


def shuffle_net(df, n=1, axis=0):     
    df = df.copy()
    for _ in range(n):
        df.apply(np.random.shuffle, axis=axis)
    return df
#https://stackoverflow.com/questions/15772009/shuffling-permutating-a-dataframe-in-pandas

def create_data(path,rand):    
    C=pd.read_csv(path,header=0)#,'type','init'],sep='\t')
    #### for randomized samples
    C=C[C['name']!='name']
    C['R']=C['R'].str.replace("False", "0")
    # pd.unique(C['name'])
    C=C[C['R']==rand]
    del C['R']

    # ####
    C['type']=C['name'].str.split('_').str[1]+'_'+C['name'].str.split('_').str[2]
    C['type']=C['type'].str.replace("_00", "_00ST")
    # # C=C[~C['type'].str.contains("03")]
    C['type']=C['type'].str.replace("_03ST", "_02ST")
    C['type']=C['type'].str.replace("_00STST", "_00ST")
    C['sample']=C['name'].str.split('_').str[0]
    C=C[C['N']!='0.0']
    # C=C[~C['name'].duplicated(keep='last')]
    C=C[~C[['type','sample']].duplicated(keep='last')]
    del C['name']
    C.reset_index(inplace=True)
    del C['index']
    D=C.pivot(index='sample', columns='type', values=['N','I','Q'])
    D=D.astype('float')
    return D

def form_tests(data,var,level):
    E0=data[var].reset_index()
    E0=E0[['sample',level+'_00ST',level+'_01ST',level+'_02ST']]
    E0['var']=var
    return E0
def merge_form(data,level):
    E0=form_tests(data,'N',level)
    E1=form_tests(data,'I',level)
    E2=form_tests(data,'Q',level)
    E=E0.append(E1)
    E=E.append(E2)
    return E
def output_data(D):
    E=merge_form(D,'CLA')
    G=merge_form(D,'LEVO')
    F=merge_form(D,'OTHER')
    H=E.merge(G,on=['sample','var'])
    H=H.merge(F,on=['sample','var'])
    # H.set_index(['var','sample'],inplace=True)
    # del H['var_x'],H['var_y']#,H0['type']
    return H


def makeSYNCSAnet(relgene,graphs,JEFF,META,deg_rand):
# for i,net in tqdm.tqdm(enumerate(BX_graphs)):
# for ii,i in tqdm():
    for ii,i in (enumerate(relgene.columns[1:])):
        ccc=nx.convert_matrix.to_pandas_edgelist(graphs[ii])
        ee=nx.convert_matrix.to_pandas_edgelist(graphs[ii])
        # cc['weight']=np.random.randn(len(cc))
        pww=i
        j=(i.split('-')[1])
        i=(i.split('-')[0])
        try:
            rrr=str(META[META['id']==i].index.item())+'_'+str(META[META['id']==i]['group'].item())+'_'+str(j)
            ccc.rename(columns={ccc.columns[2]:rrr},inplace=True)
            ddd=ccc[ccc['source'].str.contains('UniRef')]
            ddd[['source','target']] = ddd[['target','source']]
            ccc=ccc[~ccc['source'].str.contains('UniRef')].append(ddd)
            if deg_rand!=0:
                aa=pd.DataFrame(ccc)
                pcc=aa.sample(frac=np.float(deg_rand), replace=False, random_state=1) ##degree randomized
                pol=aa[~aa.isin(pcc)].dropna(how='all')
                pcc.reset_index(inplace=True)
                del pcc['index']
                lll=shuffle_net(pcc)
                ccc=pd.concat([pol,lll])
                del aa,pol,pcc,lll

            # a,b=pd.factorize(ccc['source']) 
            # c,d=pd.factorize(ccc['target'])
            # rrr=pd.DataFrame()
            # rrr['from']=a
            # rrr['to']=c
            # rrr['value']=1
            # sss=str(META[META['id']==i]['group'].item())+'_'+str(j)
            # Path('~/nest/'+sss).mkdir(parents=True, exist_ok=True)
            # rrr[['from','to','value']].to_csv('~/nest/'+sss+'/'+str(ccc.columns[2])+'.csv',sep=' ',index=False,header=False)
            # ee.rename(columns={ee.columns[2]:sss},inplace=True)
            print(ii)

            if ii==0:
                dd=ccc
                # ff=ee
            else:
                dd=dd.merge(ccc,on=['source','target'],how='outer')
                # ff=ff.merge(ee,on=['source','target'],how='outer')
            del ddd,rrr,ee,ccc
        except:
            print('no match for '+str(i))
    return dd
    

# names=pd.unique(dd.columns.str.split('_').str[1]+'_'+dd.columns.str.split('_').str[2])[1:]
# for i in tqdm(names):

def group4SYNCSA(i,dd,DR):
    # ff.columns = ff.columns.str.strip('_x')
    # ff.columns = ff.columns.str.strip('_y')
    # i=i.split('_')[1]+'_'+i.split('_')[2]
    ff=dd.loc[:,dd.columns.str.contains(i)]
    ff[['source','target']]=dd[['source','target']]
    ff=ff[ff['source'].str.contains('s__')]
    ff=ff[ff['target'].str.contains('UniRef')]
    comm=ff.groupby('source').sum().transpose()
    comm.to_csv('~/SYNCSA_eval/'+str(DR)+'_rand_comm_'+i+'.csv')
    ff.reset_index(inplace=True)
    ff.set_index(['source', 'target'], inplace=True)
    del ff['index']
    ff.columns=(ff.columns.str.split('_').str[1]+'_'+ff.columns.str.split('_').str[2])
    gg=ff.groupby(by=ff.columns, axis=1).sum()
    # traits=gg[[i]].reset_index().pivot('source','target',i).dropna(how='all',axis=1).replace(np.nan,0)
    traits=gg[[i]].reset_index().groupby(['source','target']).mean().reset_index().pivot('source','target',i).dropna(how='all',axis=1).replace(np.nan,0)

    traits.to_csv('~/SYNCSA_eval/'+str(DR)+'_rand_trait_'+i+'.csv')


def research_orthologs(uniID, species):
    '''
    research orthologs of a specific UniProtID in the inparanoid database available on the web
    save the info in the neo4j database
    '''
    # import urllib
    
    url="http://inparanoid.sbc.su.se/cgi-bin/gene_search.cgi?idtype=all;all_or_selection=all;scorelimit=0.05;rettype=xml;id=" + uniID
    # request = request(url)
    # for n in range(10):

    response = request.urlopen(url)
    txml = response.read()
    root = ET.fromstring(txml)
    for cluster in root.iter('cluster'):
        # print(cluster)
        c = cluster.attrib['nr']			#Cluster ID
        for p in cluster.iter('protein'):
            pid = p.attrib['prot_id']
            spec = p.attrib['speclong']
            if spec in species:
                P=pid
                S=spec
    return P,S
def dgtz(A):
    bis=np.round(np.sqrt(len(A))).astype(int)
    # bisA=np.abs(np.round((np.round(np.min(A),-1)-np.round(np.max(A),-1))/bis))
    # binA=np.arange(np.round(np.min(A),-1),np.round(np.max(A),-1),bisA)
    bisA=np.abs(np.round(np.min(A))-(np.max(A))/bis)
    binA=np.arange((np.min(A)),(np.max(A)),bisA)
    tmpA=np.digitize(A,bins=binA)
    return tmpA,bis

def shan_entropy(c):
    c_normalized = c / (np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log2(c_normalized))  
    return H

def get_SI(c_X,c_Z,c_XZ):
    return np.sum(np.nan_to_num(np.divide(c_XZ,c_Z)*(np.log10(np.divide(c_XZ,(np.matmul(c_X,c_Z))))), posinf=0, neginf=0),axis=0)

def get_red(SI0,SI1,c_Z):
    minSI=np.min([SI0,SI1])
    return np.sum((c_Z)*minSI)
            
def calc_MI(H_X,H_Y,H_XY):
    return H_X + H_Y - H_XY
    
def calc_CMI(H_XZ,H_YZ,H_XYZ,H_Z):
    return H_XZ+H_YZ+H_XYZ-H_Z
def calc_II(CMI_XY,MI_XY):
    return CMI_XY-MI_XY

# cdf=[]
# pucN=[]
# data=pd.read_csv('data/Pipeline_consolidate_220301/UP_gene_sum_Lx.txt',sep='\t',header=None,index_col=0)
# data=data[1:12]
# for i in np.arange(1,len(data)):
# @jit(nopython=True)
def puc_cal(data,i,genes,cdf,puc_genes):
    if cdf!=[]:
        tiss=str(cdf)
        cdf=[]
    if puc_genes!=[]:
        omic=str(puc_genes)
        puc_genes=[]
    for j in (np.arange(1,len(data))):#,position=1, desc="j", leave=False, colour='magenta'):
        pucN=[];
        for k in (np.arange(1,len(data))):#,position=2, desc="k", leave=False, colour='cyan'):
            # try:
                # time.sleep(0.5)
            puc=0
            if (i!=j)&(i!=k)&(j!=k):
                print(genes[i],genes[j],genes[k])
                A=np.array(data.iloc[i]).tolist()
                B=np.array(data.iloc[j]).tolist()
                C=np.array(data.iloc[k]).tolist()
                # print([A,B,C])
                # print(A)
        # A=np.array(data[1]).tolist()
        # B=np.array(data[2]).tolist()
        # C=np.array(data[3]).tolist()

                tmpA,ba=dgtz(A);
                tmpB,bb=dgtz(B);
                tmpC,cc=dgtz(C)
                A=tmpA/np.sum(tmpA);B=tmpB/np.sum(tmpB);C=tmpC/np.sum(tmpC);
                b=np.round(np.sqrt(len(A))).astype(int)

                c_XY = np.histogramdd([A,B],bins=(b,b))[0]
                c_XZ = np.histogramdd([A,C],bins=(b,b))[0]
                c_YZ = np.histogramdd([B,C],bins=(b,b))[0]
                c_XYZ = np.histogramdd([A,B,C],bins=(b,b,b))[0]
                c_Z = np.histogramdd([C],bins=(b))[0]
                c_Y = np.histogramdd([B],bins=(b))[0]
                c_X = np.histogramdd([A],bins=(b))[0]

                # H=[]
                # for ii,i in enumerate([c_X,c_Y,c_Z,c_XY,c_XZ,c_YZ,c_XYZ]):
                    # i.split('_')[1]
                H_X= shan_entropy(c_X)
                H_Y = shan_entropy(c_Y)
                H_Z = shan_entropy(c_Z)
                H_XY = shan_entropy(c_XY)
                H_XZ = shan_entropy(c_XZ)
                H_YZ = shan_entropy(c_YZ)
                H_XYZ = shan_entropy(c_XYZ)


                MIxy=calc_MI(H_X,H_Y,H_XY)
                MIxz=calc_MI(H_X,H_Z,H_XZ)
                MIyz=calc_MI(H_Z,H_Y,H_YZ)

                CMIx=calc_CMI(H_XZ,H_XY,H_XYZ,H_X)
                CMIy=calc_CMI(H_XY,H_YZ,H_XYZ,H_Y)
                CMIz=calc_CMI(H_XZ,H_YZ,H_XYZ,H_Z)

                [calc_II(CMIz,MIxy),
                calc_II(CMIy,MIxz),
                calc_II(CMIx,MIyz)]

                SI_xy=get_SI(c_X,c_Y,c_XY)
                SI_xz=get_SI(c_X,c_Z,c_XZ)
                SI_yz=get_SI(c_Y,c_Z,c_YZ)
                
                rZ=get_red(SI_xz,SI_yz,c_Z)
                rY=get_red(SI_xy,SI_yz,c_Y)
                rX=get_red(SI_xy,SI_xz,c_X)
    #                 # [rX[0],rY[0],rZ[0]]

    #                 Szxy=calc_II(CMIz,MIxy)+rX[0]

                u_xy=MIxy-rZ#[0]
                u_xz=MIxz-rY#[0]
                u_yz=MIyz-rX#[0]

    #                 PID=Szxy + u_xz + u_yz + rZ[0]
                # PID

                # if j!=k:
                puc+=(u_xy/MIxy)
                # print(puc)

                    # else:
                    #     puc.append(1)
                # except:
                #     print(i,j,k)
                # print(puc)
                # pucN.append(puc)

            # if i!=j:
            cdf.append(puc)
            puc_genes.append([genes[i],genes[j]])

            
            # if (i!=j)&(i!=k)&(j!=k):
            #     cdf.append(pucN)
        with open('data/Pipeline_consolidate_220301/gcn/'+tiss+'_'+omic+'_cdf.pkl', 'wb') as f:
            pickle.dump(cdf, f)
        with open('data/Pipeline_consolidate_220301/gcn/'+tiss+'_'+omic+'_puc_genes.pkl', 'wb') as f:
            pickle.dump(puc_genes, f)
        
        # return cdf, puc_genes
        
        


def bipart(G,method,nodes):
    return pd.DataFrame.from_dict(eval('bipartite.'+method)(G,nodes),orient='index',columns=[str(method)])

def calc_bipart(Bgraphs,META,ii,i):
# for ii,i in enumerate(relgene.columns[1:]):
    # ccc=nx.convert_matrix.to_pandas_edgelist(graphs[i])
    j=(i.split('-')[1])
    i=(i.split('-')[0])
    rrr=str(META[META['id']==i].index.item())+'_'+str(META[META['id']==i]['group'].item())+'_'+str(j)
    # GG=rrr
    G=Bgraphs[ii]
    if not os.path.isfile('~/data/gcn/comp_net/'+rrr+'_post_nest_clust_BIP.txt'):

        remove = [node for node,degree in dict(G.degree()).items() if degree < 3]
        G.remove_nodes_from(remove)
        top_nodes = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
        bottom_nodes = set(G) - top_nodes
        try:
            tC=bipart(G,'clustering',top_nodes)
            tLC=bipart(G,'latapy_clustering',top_nodes)
            tNR=bipart(G,'node_redundancy',top_nodes)
            tDC=bipart(G,'degree_centrality',top_nodes)
            tBC=bipart(G,'betweenness_centrality',top_nodes)
            tCC=bipart(G,'closeness_centrality',top_nodes)
            # tMEC=bipart(G,'min_edge_cover',top_nodes)
            tSBP=bipart(G,'spectral_bipartivity',top_nodes)

            bC=bipart(G,'clustering',bottom_nodes)
            bLC=bipart(G,'latapy_clustering',bottom_nodes)
            bNR=bipart(G,'node_redundancy',bottom_nodes)
            bDC=bipart(G,'degree_centrality',bottom_nodes)
            bBC=bipart(G,'betweenness_centrality',bottom_nodes)
            bCC=bipart(G,'closeness_centrality',bottom_nodes)
            # bMEC=bipart(G,'min_edge_cover',bottom_nodes)
            bSBP=bipart(G,'spectral_bipartivity',bottom_nodes)


            data_frames = [tC,tLC,tNR,tDC,tBC,tCC,tSBP,bC,bLC,bNR,bDC,bBC,bCC,bSBP]
            df_merged = reduce(lambda  left,right: pd.merge(left,right,right_index=True,left_index=True,
                                                how='outer'), data_frames)

            df_merged.to_csv('~/data/gcn/comp_net/'+rrr+'_post_nest_clust_BIP2.txt',sep='\t',header=False)
        except:
            print([i,j,rrr])
            
def pdize_net(OTHER_00ST):
    bb=pd.DataFrame(columns=['source','target'])
    for i,ii in enumerate(OTHER_00ST):
        cc=nx.to_pandas_edgelist(OTHER_00ST[i])
        cc['weight']=1
        bb=bb.merge(cc,how='outer',right_on=['source','target'],left_on=['source','target'])
    bb.fillna(0)
    aa=pd.DataFrame(bb.set_index(['source','target']).fillna(0).mean(axis=1))
    cc=pd.DataFrame(bb.groupby('target').sum().mean(axis=1))
    dd=pd.DataFrame(bb.groupby('source').sum().mean(axis=1))
    return aa,cc,dd


def plot_time_net(partitions,i,LL):
    # if LL==0:
    #     LL='OTHER'
    # elif LL==1:
    #     LL='CLA'
    # elif LL==2:
    #     LL='LEVO'
    # i=str(i)
    # Graph.layout_bipartite()
    ig.plot(partitions[i].graph,"run/gcn/img/net/la_c"+str(i)+'_'+LL+"_net.png",layout="auto")

    G0=partitions[i].graph.get_edge_dataframe().sort_values(by=['source','target'])
    AA=G0.pivot('source',columns='target',values='weight').fillna(0)
    plt.imshow(AA, cmap='hot', interpolation='nearest',aspect='auto')
    plt.savefig("run/gcn/img/sq/la_c"+str(i)+'_'+LL+"_sq.png",dpi=300,bbox_inches = "tight")

    # ig.plot(partitions[i].graph,'LEVO_la_c0.png')
    partition = la.find_partition(eval('G_'+str(i)), la.ModularityVertexPartition)
    eval('G_'+str(i)).vs['cluster'] = partition.membership
    ig.plot(partition.graph,'run/gcn/img/out/la_c'+str(i)+'_'+LL+'_out.png',layout="auto")
    
    g = eval('G_'+str(i)).to_graph_tool()
    state = gt.minimize_blockmodel_dl(g, state=gt.PPBlockState)
    state.multiflip_mcmc_sweep(beta=np.inf, niter=100)
    state.draw(output='run/gcn/img/gt/la_c'+str(i)+'_'+LL+'_GT.png')#pos=g.vp.pos)
    state = gt.minimize_nested_blockmodel_dl(g)
    state.draw(output='run/gcn/img/nest/la_c'+str(i)+'_'+LL+'_GT.png')
    
    tree = min_spanning_tree(g)
    graph_draw(g, edge_color=tree,output="run/gcn/img/tree/min_tree"+str(i)+'_'+LL+".png")

    g.set_edge_filter(tree)
    graph_draw(g,output="run/gcn/img/tree/min_tree_filtered"+str(i)+'_'+LL+".png")
    
    bv, be = betweenness(g)
    be.a /= be.a.max() / 5
    graph_draw(g, vertex_fill_color=bv, edge_pen_width=be, output="run/gcn/img/bt_tree/filtered-bt"+str(i)+'_'+LL+".png")


   
    
    
# turn this into a function
def igraph_from_pandas(edges_table, vertices_table, source_cl='from', target_cl='to', vertex_attrs=None, vertex_id_cl='v_id', directed=False):

    import pandas as pd
    import igraph as ig
    # control parameters
    if isinstance(edges_table, pd.DataFrame):
        try:
            if source_cl and target_cl in edges_table.columns:
                id_gen = ig.UniqueIdGenerator()
                edgelist = []
                for start_edge, end_edge in edges_table[[source_cl, target_cl]].itertuples(index=False, name=None):
                    edgelist.append((id_gen[start_edge], id_gen[end_edge]))
                if directed:
                    gg = ig.Graph(edgelist, directed=True)
                else:
                    gg = ig.Graph(edgelist, directed=False)
                gg.vs["name"] = id_gen.values()
        except (KeyError, NameError):
            raise ValueError('Edges columns missing!')
    else:
        raise ValueError("edges table is required!")
    if isinstance(vertices_table, pd.DataFrame):
        if not vertex_attrs:
            raise ValueError('No attributes provided. Remove vertices table from arguments and try again.')
        else:
            try:
                # order vertices table based on edge_list
                vertices_table_ordered = pd.DataFrame(id_gen.values(), columns=['unique_id'])
                # bring previous vertices table with attributes (to be reordered)
                vertices_table_ordered = vertices_table_ordered.merge(vertices_table, left_on = 'unique_id', right_on = vertex_id_cl, how='left')
                for attr2use in vertex_attrs:
                    if attr2use in vertices_table.columns:
                        # add attributes to graph
                        gg.vs[attr2use] = vertices_table_ordered[attr2use].values.tolist()
            except (KeyError, NameError):
                raise ValueError('Vertex ID column missing!')    
    return gg

def explode_graphs(relgene,META,B2graphs):
    qqq=[]
    for ii,i in enumerate(relgene.columns[1:-2]):
        j=(i.split('-')[1])
        i=(i.split('-')[0])
        rrr=str(META[META['id']==i]['group'].item())+'_'+str(j)
        rrr=rrr.replace("03ST", "02ST")
        qqq.append(rrr)
        B2graphs[ii].name=str(rrr)
        # print(rrr,ii)

        with open('data/gcn/'+rrr+'.pkl', 'wb') as f:
                pickle.dump(B2graphs[ii], f)
        # if
    # zzz=[]
    for i in np.unique(qqq):
        zzz=[]
        for j in B2graphs:
            # print(i,j.name)
            if j.name==i:

                zzz.append(j)
            with open('data/gcn/'+i+'.pkl', 'wb') as f:
                pickle.dump(zzz, f)
                
                
def plot_comm(CC,LL,OO,palette,j): ##time, community_type
    
    c=0
    fig, ax = plt.subplots(ncols=3, nrows=3,figsize=(15,15))
    # dict_keys = [k for k in z.keys()]
    for k in [CC,LL,OO]:
        web=(k)
        # k=eval(k)
        for i in k.columns[2:]:
            
            jj0=pd.DataFrame(k[['source','target',i]])#.reset_index()
            jj0.columns=['start','end','value']
            jj0=jj0[jj0['value']>0]
            kk=pd.DataFrame(np.unique(jj0[['start','end']].melt()['value']))
            kk.columns=['name']
            kk['id']=kk['name']
            kk['types']=np.concatenate([np.ones(len(np.unique(jj0.end))),np.zeros(len(np.unique(jj0.start)))]).astype('int')
            G_00=igraph_from_pandas(edges_table=jj0, vertices_table=kk, source_cl='start', target_cl='end', vertex_attrs=list(kk.columns), vertex_id_cl='name', directed=False)

            G_00.vs['pagerank']=G_00.pagerank()
            G_00.vs['cluster'] = G_00.community_infomap().membership
                # N,Q,I,R=Parallel(n_jobs=10)(structural_analysis(ii,i,graphs,ARG_meta,rand=rand,deg_rand=deg_rand))

            # for j in [edge_betweenness,fastgreedy,infomap,label_propagation,leading_eigenvector,leiden,multilevel,_optimal_modularity,spinglass,walktrap]:
            # G_00.vs['louvain_membership']=G_00.community_multilevel().membership
            try:
                Path("run/gcn/img/"+j).mkdir(parents=True, exist_ok=True)
                G_00.vs[j+'membership']=eval('G_00.community_'+j)().membership
            except:
                G_00.vs[j+'membership']=eval('G_00.community_'+j)()
            g00e=G_00.get_edge_dataframe()#.reset_index()
            g00v=G_00.get_vertex_dataframe().reset_index()
            g00e['weight']=jj0['value']
            jj=g00e.merge(g00v,left_on='source',right_on='vertex ID').merge(g00v,left_on='target',right_on='vertex ID')
            jj.name_x, jj.name_y = np.where(jj.name_x.str.contains('UniRef'), [jj.name_y, jj.name_x], [jj.name_x, jj.name_y])
            # for j in ['louvain_','leiden_']:
            gg=jj.groupby(['name_x',j+'membership_x']).count().sort_values(['target',j+'membership_x'])
            gg=gg.reset_index()[['name_x',j+'membership_x','target']]
            # hh=hh.merge(gg,how='outer')
            
            # gg[j+'membership_x']=gg[j+'membership_x'].astype('category')
            sns.histplot(gg, x=gg[j+'membership_x'].astype('int')+1, hue='name_x',log_scale=[False,True],weights='target',multiple='stack',shrink=0.8,palette=palette,ax=ax.flat[c],legend=False).set_title(web+"_"+i)
            
            c=c+1
    # plt.figure(figsize=(10,6))
    
            # plt.set(xlabel=i+' '+j+' membership',yscale='log',ylabel='link count')
    # g.set_title(str(C)+'_'+i+'_'+j)
    # g.set_xlabel(i+' '+j+' membership')
            # g.set_yscale('log')
    plt.savefig('run/gcn/img/'+j+'.png', dpi=100)
    
    

def card_net(OTHER_00ST,uni_conv,card):
    bb=pd.DataFrame(columns=['source','target'])
    for i,ii in enumerate(OTHER_00ST):
        cc=nx.to_pandas_edgelist(OTHER_00ST[i])
        cc['weight']=1
        bb=bb.merge(cc,how='outer',right_on=['source','target'],left_on=['source','target'])
    bb.fillna(0)
    # aa=pd.DataFrame(bb.set_index(['source','target']).fillna(0).mean(axis=1))
    # ee=aa
    bb.reset_index(inplace=True)
    cc=bb.loc[:,bb.columns.str.contains('weight')]
    cc=1-(np.sum(cc,axis=0)/len(cc))
    
    bb['target']=bb['target'].str.split('_').str[1]#.to_csv('CLA_uniref.txt',sep='\t')
    bb=bb.merge(uni_conv,left_on='target',right_on='Entry')
    bb=bb.merge(card,on='Entry name')#[['source','target',0]]
    bb.fillna(0,inplace=True)
    bb=bb.loc[:,bb.columns.str.contains('weight')]
    dd=1-(np.sum(bb,axis=0)/len(bb))
    return bb,cc,dd