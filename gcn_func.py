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
import plotly
from joblib import Parallel
import sklearn.utils as sku
import plotly.graph_objects as go
import plotly.express as px
# j=sys.argv[1]

sys.path.insert(1, './nestedness_modularity_in-block_nestedness_analysis-master/')
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
        ff=ff[ff['source'].str.contains('s__')]
        ff=ff[ff['target'].str.contains('UniRef')]
        ff.groupby('source').sum().transpose().to_csv('comm_'+i+'.csv')
        ff.reset_index(inplace=True)
        ff.set_index(['source', 'target'], inplace=True)
        del ff['index']
        ff.columns=(ff.columns.str.split('_').str[1]+'_'+ff.columns.str.split('_').str[2])
        gg=ff.groupby(by=ff.columns, axis=1).sum()
        traits=gg[[i]].reset_index().pivot('source','target',i).dropna(how='all',axis=1).replace(np.nan,0)
        traits.to_csv('trait_'+i+'.csv')
    

def buildNestedNess():
    C=pd.DataFrame(columns=['N','Q','I','type'f])
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
    
    a,b=pd.factorize(ccc['source']) 
    c,d=pd.factorize(ccc['target'])
    rrr=pd.DataFrame()
    rrr['from']=a
    rrr['to']=c
    rrr['value']=1
    sss=str(ARG_meta[ARG_meta['id']==i]['group'].item())+'_'+str(j)
    Path('nest/'+sss).mkdir(parents=True, exist_ok=True)
    # rrr[['from','to','value']].to_csv('~/nest/'+sss+'/'+str(ccc.columns[2])+'.csv',sep=' ',index=False,header=False)
    aa= rrr[['from','to','value']].values
    
    if rand==True: ## to randomize
        aa=pd.DataFrame(aa)
        ddd=aa.sample(frac=np.float(deg_rand), replace=False, random_state=1) ##degree randomized
        rrr=aa[~aa.isin(ddd)].dropna(how='all')
        ddd.reset_index(inplace=True)
        del ddd['index']
        sss=shuffle(ddd)
        aa=pd.concat([rrr,sss])
        aa=np.array(aa).astype(int)
    
    nodes_cols = int(max(aa[j,1] for j in range(aa.shape[0]))+1)
    nodes_rows= int(max(aa[j,0] for j in range(aa.shape[0]))+1)
    matrix=np.zeros((nodes_rows,nodes_cols),dtype='int')
    for j in range(aa.shape[0]):
        matrix[aa[j,0],aa[j,1]] = 1
    M=matrix

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
    with open("randB_ARG_nest.txt", "ab") as f:
        np.savetxt(f,np.column_stack([str(N),str(Q),str(I),sss,pww]),delimiter = '\t', fmt='%s')
    return N,Q,I,sss,pww



####################TRASH FOLLOWS#########################
   


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

    
    #!/usr/bin/env python
"""
Plot multi-graphs in 3D.
"""
# import numpy as np
# import matplotlib.pyplot as plt
# import networkx as nx

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection


class LayeredNetworkGraph(object):

    def __init__(self, graphs, node_labels=None, layout=nx.spring_layout, ax=None):
        """Given an ordered list of graphs [g1, g2, ..., gn] that represent
        different layers in a multi-layer network, plot the network in
        3D with the different layers separated along the z-axis.

        Within a layer, the corresponding graph defines the connectivity.
        Between layers, nodes in subsequent layers are connected if
        they have the same node ID.

        Arguments:
        ----------
        graphs : list of networkx.Graph objects
            List of graphs, one for each layer.

        node_labels : dict node ID : str label or None (default None)
            Dictionary mapping nodes to labels.
            If None is provided, nodes are not labelled.

        layout_func : function handle (default networkx.spring_layout)
            Function used to compute the layout.

        ax : mpl_toolkits.mplot3d.Axes3d instance or None (default None)
            The axis to plot to. If None is given, a new figure and a new axis are created.

        """

        # book-keeping
        self.graphs = graphs
        self.total_layers = len(graphs)

        self.node_labels = node_labels
        self.layout = layout

        if ax:
            self.ax = ax
        else:
            fig=plt.figure(figsize=(5,15))
            self.ax = fig.add_subplot(111, projection='3d')

        # create internal representation of nodes and edges
        self.get_nodes()
        self.get_edges_within_layers()
        self.get_edges_between_layers()

        # compute layout and plot
        self.get_node_positions()
        self.draw()


    def get_nodes(self):
        """Construct an internal representation of nodes with the format (node ID, layer)."""
        self.nodes = []
        for z, g in enumerate(self.graphs):
            self.nodes.extend([(node, z) for node in g.nodes()])


    def get_edges_within_layers(self):
        """Remap edges in the individual layers to the internal representations of the node IDs."""
        self.edges_within_layers = []
        for z, g in enumerate(self.graphs):
            self.edges_within_layers.extend([((source, z), (target, z)) for source, target in g.edges()])


    def get_edges_between_layers(self):
        """Determine edges between layers. Nodes in subsequent layers are
        thought to be connected if they have the same ID."""
        self.edges_between_layers = []
        for z1, g in enumerate(self.graphs[:-1]):
            z2 = z1 + 1
            h = self.graphs[z2]
            shared_nodes = set(g.nodes()) & set(h.nodes())
            self.edges_between_layers.extend([((node, z1), (node, z2)) for node in shared_nodes])


    def get_node_positions(self, *args, **kwargs):
        """Get the node positions in the layered layout."""
        # What we would like to do, is apply the layout function to a combined, layered network.
        # However, networkx layout functions are not implemented for the multi-dimensional case.
        # Futhermore, even if there was such a layout function, there probably would be no straightforward way to
        # specify the planarity requirement for nodes within a layer.
        # Therefor, we compute the layout for the full network in 2D, and then apply the
        # positions to the nodes in all planes.
        # For a force-directed layout, this will approximately do the right thing.
        # TODO: implement FR in 3D with layer constraints.

        composition = self.graphs[0]
        for h in self.graphs[1:]:
            composition = nx.compose(composition, h)

        pos = self.layout(composition, *args, **kwargs)

        self.node_positions = dict()
        for z, g in enumerate(self.graphs):
            self.node_positions.update({(node, z) : (*pos[node], z) for node in g.nodes()})


    def draw_nodes(self, nodes, *args, **kwargs):
        x, y, z = zip(*[self.node_positions[node] for node in nodes])
        self.ax.scatter(x, y, z, *args, **kwargs)


    def draw_edges(self, edges, *args, **kwargs):
        segments = [(self.node_positions[source], self.node_positions[target]) for source, target in edges]
        line_collection = Line3DCollection(segments, *args, **kwargs)
        self.ax.add_collection3d(line_collection)


    def get_extent(self, pad=0.1):
        xyz = np.array(list(self.node_positions.values()))
        xmin, ymin, _ = np.min(xyz, axis=0)
        xmax, ymax, _ = np.max(xyz, axis=0)
        dx = xmax - xmin
        dy = ymax - ymin
        return (xmin - pad * dx, xmax + pad * dx), \
            (ymin - pad * dy, ymax + pad * dy)


    def draw_plane(self, z, *args, **kwargs):
        (xmin, xmax), (ymin, ymax) = self.get_extent(pad=0.1)
        u = np.linspace(xmin, xmax, 10)
        v = np.linspace(ymin, ymax, 10)
        U, V = np.meshgrid(u ,v)
        W = z * np.ones_like(U)
        self.ax.plot_surface(U, V, W, *args, **kwargs)


    def draw_node_labels(self, node_labels, *args, **kwargs):
        for node, z in self.nodes:
            if node in node_labels:
                ax.text(*self.node_positions[(node, z)], node_labels[node], *args, **kwargs)


    def draw(self):

        self.draw_edges(self.edges_within_layers,  color='k', alpha=0.3, linestyle='-', zorder=2)
        self.draw_edges(self.edges_between_layers, color='k', alpha=0.3, linestyle='--', zorder=2)

        for z in range(self.total_layers):
            self.draw_plane(z, alpha=0.2, zorder=1)
            self.draw_nodes([node for node in self.nodes if node[1]==z], s=300, zorder=3)

        if self.node_labels:
            self.draw_node_labels(self.node_labels,
                                  horizontalalignment='center',
                                  verticalalignment='center',
                                  zorder=100)


# if __name__ == '__main__':




def plot_sankey(data,XX,A,B):
    # https://gist.github.com/nicolasesnis/595d34c3c7dbca2b3419332304954433
    # Working on the nodes_dict

    all_events = list(data[B].unique())

    # Create a set of colors that you'd like to use in your plot.
    palette = ['50BE97', 'E4655C', 'FBEEAC', '3E5066',
               'BFD6DE', 'FCC865', '353A3E', 'E6E6E6']#,
              # 'E4655C', 'FBEEAC', '3E5066',
               # 'BFD6DE', 'FCC865', '353A3E', 'E6E6E6']
    #  Here, I passed the colors as HEX, but we need to pass it as RGB. This loop will convert from HEX to RGB:
    for i, col in enumerate(palette):
        palette[i] = tuple(int(col[i:i+2], 16) for i in (0, 2, 4))

    # Append a Seaborn complementary palette to your palette in case you did not provide enough colors to style every event
    complementary_palette = sns.color_palette("deep", len(all_events) - len(palette))
    if len(complementary_palette) > 0:
        palette.extend(complementary_palette)

    output = dict()
    output.update({'nodes_dict': dict()})

    i = 0
    for t in data.t.unique(): # For each tri of clus...
        # Create a new key equal to the tri...
        output['nodes_dict'].update(
            {t: dict()}
        )

        # Look at all the events that were done at this step of the funnel...
        all_events_at_this_rank = data[data.t ==
                                       t][B].unique()

        # Read the colors for these events and store them in a list...
        tri_palette = []
        for clus in all_events_at_this_rank:
            tri_palette.append(palette[list(all_events).index(clus)])

        # Keep trace of the events' names, colors and indices.
        output['nodes_dict'][t].update(
            {
                'sources': list(all_events_at_this_rank),
                'color': tri_palette,
                'sources_index': list(range(i, i+len(all_events_at_this_rank)))
            }
        )
        # Finally, increment by the length of this rank's available clus to make sure next indices will not be chosen from existing ones
        i += len(output['nodes_dict'][t]['sources_index'])
        
        
    # Working on the links_dict

    output.update({'links_dict': dict()})

    # Group the DataFrame by pat_epi and tri


    grouped = data.groupby([A, 't'])
    jaws=0
    # Define a function to read the souces, targets, values species to next_clus:
    def update_source_target(user):
        try:
            # user.name[0] is the user's spec; user.name[1] is the time
            # 1st we retrieve the source and target's indices from nodes_dict
            source_index = output['nodes_dict'][user.name[1]]['sources_index'][output['nodes_dict']
                                                                               [user.name[1]]['sources'].index(user[B].values[0])]
            target_index = output['nodes_dict'][user.name[1] + 1]['sources_index'][output['nodes_dict']
                                                                                   [user.name[1] + 1]['sources'].index(user['next_clus'].values[0])]
            # print(user.name)
             # If this source is already in links_dict...
            if source_index in output['links_dict']:
                # ...and if this target is already associated to this source...
                if target_index in output['links_dict'][source_index]:
                    # ...then we increment the count of users with this source/target pair by 1
                    output['links_dict'][source_index][target_index]['unique_users'] += 1
                # ...but if the target is not already associated to this source...
                else:
                    # ...we create a new key for this target, for this source, and initiate it with 1 user and the time from source to target
                    output['links_dict'][source_index].update({target_index:
                                                               dict(
                                                                   {'unique_users': 0}
                                                                )
                                                               })
            # ...but if this source isn't already available in the links_dict, we create its key and the key of this source's target, and we initiate it with 1 user and the time from source to target
            else:
                output['links_dict'].update({source_index: dict({target_index: dict(
                    {'unique_users': 0})})})
        except Exception as e:
            pass

    # Apply the function to your grouped Pandas object:
    grouped.apply(lambda user: update_source_target(user)) 


    targets = []
    sources = []
    values = []

    for source_key, source_value in output['links_dict'].items():
        for target_key, target_value in output['links_dict'][source_key].items():
            sources.append(source_key)
            targets.append(target_key)
            values.append(target_value['unique_users'])

    labels = []
    colors = []

    for key, value in output['nodes_dict'].items():
        labels = labels + list(output['nodes_dict'][key]['sources']) 
        colors = colors + list(output['nodes_dict'][key]['color'])
    for idx, color in enumerate(colors):
        colors[idx] = "rgb" + str(color) + ""
    # [len(targets),len(sources),len(values)]
    return labels,color,sources,targets,values
    
