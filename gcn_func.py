import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
# from tqdm import tqdm
from scipy.stats import rankdata
import pandas as pd
import networkx as nx
import seaborn as sns
from joblib import delayed, wrap_non_picklable_objects
from pathlib import Path
# j=sys.argv[1]


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
    for q,ee in enumerate(np.unique(jeff.species)):
        jeff2=jeffA[jeffA['species']==ee]#.explode('target')
        dd=pd.DataFrame(jeff2['target'].to_numpy().reshape(int(len(jeff2)/tim_len),tim_len,order='F'))

        if len(dd.melt())==(tim_len*XX):
            JJ=JJ.append(dd)
            rr=np.append(rr, ee)
    jeffA=jeffA.sort_values(['species', 't'], ascending=[True, True])      
    labels,levels=pd.factorize((JJ.melt()['value']))
    cc=pd.DataFrame(np.array(labels).reshape((XX)*len(rr),tim_len,order='F'))

    for i in np.arange(0,len(cc),XX+1):
        for col in cc:
            cc.iloc[i:i+XX,col] = cc.iloc[i:i+XX,col].sort_values(ignore_index=True)
            cc.loc[i+XX]=0

    plt.figure(figsize=(10,30))
    ax=sns.heatmap(cc,cmap='rocket_r',annot=True, fmt="d",cbar=False,xticklabels=False,
                   yticklabels=False).set(ylabel='         -         '.join(rr))

    # plt.show()
    # Path("data/gcn/img/"+cc).mkdir(parents=True, exist_ok=True)
    plt.savefig("data/gcn/img/full_"+str(gg)+"_10.png",dpi=300,bbox_inches = "tight")


def plotRidge(df,typ,measur,patt):
    # import numpy as np
    # import pandas as pd
    # 
    # import matplotlib.pyplot as plt
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Create the data
    # rs = np.random.RandomState(1979)
    # x = rs.randn(XX0)
    # g = np.tile(list("ABCDEFGHIJ"), XX)
    # df = pd.DataFrame(dict(x=x, g=g))
    # m = df.g.map(ord)
    # df["x"] += m

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(df[df['type']==typ], row="gen", hue="gen", aspect=20, height=.5, palette=pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "value",
          bw_adjust=.5, clip_on=False,
          fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, "value", clip_on=False, color="w", lw=2, bw_adjust=.5)

    # passing color=None to refline() uses the hue mapping
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)


    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)


    g.map(label, "value")

    # Set the subplots to overlap
    g.figure.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)
    plt.savefig("data/gcn/"+patt+"_"+str(measur)+".png",dpi=300,bbox_inches = "tight")


    