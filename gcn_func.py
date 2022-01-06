import matplotlib.pyplot as plt
import numpy as np
import os,glob,sys,importlib,pickle#,scipy,coolbox,pybedtools,
# from tqdm import tqdm
import pandas as pd
import networkx as nx
import seaborn as sns
from joblib import delayed, wrap_non_picklable_objects

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
    ee=nx.from_pandas_edgelist(dd,source='spec',target='gene')
    remove = [node for node,degree in dict(ee.degree()).items() if degree <10]
    ee.remove_nodes_from(remove)
    # ff.append(ee)
    
    B = nx.Graph()
    B.add_nodes_from(dd['spec'], bipartite=0)
    B.add_nodes_from(dd['gene'], bipartite=1)
    B.add_edges_from(tuple(dd[['spec','gene']].itertuples(index=False, name=None)))
    remove = [node for node,degree in dict(B.degree()).items() if degree <10]
    B.remove_nodes_from(remove)
    # C.append(B)

    # with open('data/gcn/NX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
    #     pickle.dump(ff, f)
    # with open('data/gcn/BX_'+str(patt)+'_hypert.pkl', 'ab+') as f:
    #     pickle.dump(C, f)
    return ee,B

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
    df.to_csv("data/gcn/"+patt+"_"+str(measur)+"_.txt",sep='\t')

    plt.figure(figsize=(10,30))
    sns.set_theme(style="whitegrid")

    sns.violinplot(data=df, y="gen", x="value",hue="type",
                   split=True, inner="quart", linewidth=1,
                   orient="h")
    sns.despine(left=True)
    plt.savefig("data/gcn/"+patt+"_"+str(measur)+"_violin.png",dpi=300,bbox_inches = "tight")
    return df

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

