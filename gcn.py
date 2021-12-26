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

# from joblib.externals.loky import set_loky_pickler
# from joblib import parallel_backend
# from joblib import Parallel, delayed
# from joblib import wrap_non_picklable_objects

# j=sys.argv[1]

relgene=pd.read_csv('../../groups/cgsd/gordonq/LauG_Metagenomics_CPOS-200710-CJX-3455a/50_genefamilies.tsv',sep='\t')

relgene['gene']=relgene['# Gene Family'].str.split('|').str[0]
relgene=relgene[relgene['gene']!='UniRef90_unknown']
relgene=relgene[relgene['gene']!='UNMAPPED']
relgene.index=relgene['# Gene Family']
del relgene['gene'], relgene['# Gene Family']
# relgene=relgene/relgene.sum(axis=0)
relgene=relgene/relgene.sum(axis=0)
relgene['gen']=relgene.index.str.split('|').str[1].str.split('.').str[0].tolist()
relgene['spec']=relgene.index.str.split('.').str[1]#.str.split('.').str[0].tolist()
relgene['spec'].replace('_',' ')
relgene.index=relgene.index.str.split('|').str[0]
relgene=relgene.dropna()


cc=relgene.groupby(['# Gene Family','spec']).sum()
# dd=relgene.groupby(['# Gene Family','gen']).sum()
cc=cc.reset_index()
# dd=dd.reset_index()
cc=cc.rename(columns={'# Gene Family':'gene'})#,'spec':0,'gene':1})


primary=pd.read_excel('data/Data Raw - Gut Microbiome Cohort Project Database - 300 Cohort v3.0_280921.xlsx',index_col=0,sheet_name='Primary Data')
uni_bact=primary[['Age','Hypertension Category by 24h BP w/o considering antihypertensive med']]
uni_bact=uni_bact.rename(columns={"Hypertension Category by 24h BP w/o considering antihypertensive med": "HT"})


# @delayed
# @wrap_non_picklable_objects

def bip(cc,net,ff,C):

    # print(net)
    dd=cc[['spec','gene',net]]
    dd=dd[dd[net]!=0]
    ee=nx.from_pandas_edgelist(dd,source='spec',target='gene')
    remove = [node for node,degree in dict(ee.degree()).items() if degree <10]
    ee.remove_nodes_from(remove)
    ff.append(ee)
    
    B = nx.Graph()
    B.add_nodes_from(dd['spec'], bipartite=0)
    B.add_nodes_from(dd['gene'], bipartite=1)
    B.add_edges_from(tuple(dd[['spec','gene']].itertuples(index=False, name=None)))
    remove = [node for node,degree in dict(B.degree()).items() if degree <10]
    B.remove_nodes_from(remove)
    C.append(B)

    with open('data/gcn/NX_50_hypert.pkl', 'ab') as f:
        pickle.dump(ff, f)
    with open('data/gcn/BX_50_hypert.pkl', 'ab') as f:
        pickle.dump(C, f)

def load_list_of_dicts(filename, create_using=nx.Graph):
    with open(filename, 'rb') as f:
        list_of_dicts = pickle.load(f)   
    graphs = [create_using(graph) for graph in list_of_dicts]
    return graphs


def meas(measur,uni_bact,relgene,graphs):
    HT50=uni_bact[uni_bact.index.isin(relgene.columns[:-2].str.split('-').str[0])]
    HT50['index']=np.arange(len(HT50))
    measur=eval(measur)
    S = [measur(graphs[i]) for i in HT50[HT50['HT']==0]['index'].values]
    T = [measur(graphs[i]) for i in HT50[HT50['HT']!=0]['index'].values]

    non=pd.DataFrame(S).melt()
    non['type']='NoHT'
    non.dropna(inplace=True)
    non=non[non.value!=0]
    non=non[~non['variable'].str.contains('UniRef90')]
    non.value=non.value/np.sum(non.value)
    yes=pd.DataFrame(T).melt()
    yes['type']='HT'
    yes.dropna(inplace=True)
    yes=yes[yes.value!=0]
    yes=yes[~yes['variable'].str.contains('UniRef90')]
    yes.value=yes.value/np.sum(yes.value)
    df=non.append(yes)
    # df=df.dropna()
    df['gen']=df.variable.str.split('_').str[2]
    return df

def plotRidge(df,typ,measur):
    # import numpy as np
    # import pandas as pd
    # import seaborn as sns
    # import matplotlib.pyplot as plt
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Create the data
    # rs = np.random.RandomState(1979)
    # x = rs.randn(500)
    # g = np.tile(list("ABCDEFGHIJ"), 50)
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
    plt.savefig("data/gcn/"+measur+"_"+typ+".png",dpi=300,bbox_inches = "tight")


# ff=[]
# C=[]
# for i,net in enumerate(cc.columns[2:]):
#     bip(cc,net,ff,C)

# Parallel(n_jobs=48) (bip(cc,net,ff,C) #for i,net in enumerate(cc.columns[2:])

graphs = load_list_of_dicts('data/gcn/BX_50_hypert.pkl')


measur='nx.betweenness_centrality'
df=meas(measur,uni_bact,relgene,graphs)
plotRidge(df,'NoHT',measur)
plotRidge(df,'HT',measur)

measur='nx.degree_centrality'
df=meas(measur,uni_bact,relgene,graphs)
plotRidge(df,'NoHT',measur)
plotRidge(df,'HT',measur)

measur='nx.closeness_centrality'
df=meas(measur,uni_bact,relgene,graphs)
plotRidge(df,'NoHT',measur)
plotRidge(df,'HT',measur)

measur='nx.node_redundancy'
df=meas(measur,uni_bact,relgene,graphs)
plotRidge(df,'NoHT',measur)
plotRidge(df,'HT',measur)