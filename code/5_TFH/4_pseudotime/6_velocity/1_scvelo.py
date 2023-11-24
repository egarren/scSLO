
import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import rpy2.robjects as robjects
import anndata2ri
import scvelo as scv
import re
import anndata as adata1
import glob
import pickle
from IPython.display import display
import pickle


scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view

anndata2ri.activate()
robjects.r['load']('sce.filtered.RData')
seurat = robjects.r('as(obj.sce, "SingleCellExperiment")')

#load loom data
looms=dict()
for filepath in glob.glob(os.path.join('loom', '*.loom')):
  ldata=scv.read(filepath,cache=True)
  ldata.var_names_make_unique()
  ldata.obs_names=[re.sub(":", "_", x) for x in ldata.obs_names]
  ldata.obs_names=[re.sub("x", "-1", x) for x in ldata.obs_names]
  looms[filepath]=ldata
  
#concatenate looms for all mice
ldata_all=adata1.AnnData.concatenate(list(looms.values())[0],list(looms.values())[1],list(looms.values())[2],list(looms.values())[3],
list(looms.values())[4],list(looms.values())[5],list(looms.values())[6],list(looms.values())[7],list(looms.values())[8],index_unique=None)
ldata_all.obs_names
seurat.obs_names

#add to seurat data
adata=scv.utils.merge(seurat,ldata_all)
adata
adata.obs["clusters"]=adata.obs["my.clusters3"]
adata.obs.clusters
print(sorted(adata.obs.clusters.unique()))
cols=["#F8766D" ,"#D39200" ,"#93AA00" ,"#00BA38" ,"#00C19F", "#00B9E3", "#619CFF", "#DB72FB" ,"#FF61C3"]
sc.pl.umap(adata, color='clusters', legend_loc='on data',palette=cols)
scv.pl.proportions(adata)

#velocity computation
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata,n_jobs=20)
with open('adata_1','wb') as adata_file:
  pickle.dump(adata,adata_file)

#visualization
scv.pl.velocity_embedding(adata, basis='umap',palette=cols,arrow_length=3,arrow_size=2,dpi=120)
scv.pl.velocity_embedding_grid(adata, basis='umap',palette=cols,arrow_length=2,arrow_size=1,dpi=120)
scv.pl.velocity_embedding_stream(adata, basis='umap',palette=cols)
scv.pl.velocity(adata, var_names=['Sostdc1','Pdcd1','Tcf7'],ncols=2,palette=cols)
scv.pl.scatter(adata, 'Sostdc1', color=['clusters', 'velocity'],palette=cols)#,add_outline='Ngn3 high EP, Pre-endocrine, Beta'

#identifying genes
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
kwargs = dict(frameon=False, size=10, linewidth=1.5,add_outline='Sostdc1_1')
scv.pl.scatter(adata, df['Sostdc1_1'][:5], ylabel='Sostdc1_1', **kwargs)
kwargs = dict(frameon=False, size=10, linewidth=1.5,add_outline='Sostdc1_2')
scv.pl.scatter(adata, df['Sostdc1_2'][:5], ylabel='Sostdc1_2', **kwargs)
with open('adata_2','wb') as adata_file:
  pickle.dump(adata,adata_file)

#cell cycling
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index
kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)
scv.pl.velocity(adata, ['Ung', 'Mki67'], ncols=2, add_outline=True)

#speed and coherence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
display(df)

#pseudotime
scv.pl.velocity_graph(adata,palette=cols,threshold=0.1)
#tracing single cell
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=np.flatnonzero(adata.obs['clusters']  == 'TFH_Il31ra')[10])
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
#velocity based pseudotime
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')
with open('adata_3','wb') as adata_file:
  pickle.dump(adata,adata_file)


#trajectory inference
sc.pl.umap(adata, color='clusters', legend_loc='on data')
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='clusters',minimum_spanning_tree=False)#,threshold_root_end_prior=0.5,root_key="A",end_key="E",
df = scv.get_df(adata, 'paga/transitions_confidence', precision=1).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
df
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5)
with open('adata_4','wb') as adata_file:
  pickle.dump(adata,adata_file)

#differential kinetics
var_names = ['Sostdc1']
scv.tl.differential_kinetic_test(adata, var_names=var_names, groupby='clusters')
scv.get_df(adata[:, var_names], ['fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(adata, basis=var_names, add_outline='fit_diff_kinetics', **kwargs)
diff_clusters=list(adata[:, var_names].var['fit_diff_kinetics'])
scv.pl.scatter(adata, legend_loc='right', size=40, title='diff kinetics',
               add_outline=diff_clusters, outline_width=(.8, .2))

#top genes 
scv.tl.recover_dynamics(adata,n_jobs=20)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(adata, var_names=top_genes, groupby='clusters')
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(adata, basis=top_genes[15:30], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
with open('adata_5','wb') as adata_file:
  pickle.dump(adata,adata_file)
with open('adata_5','rb') as adata_file:
  adata=pickle.load(adata_file)

#recompute velocities
scv.tl.velocity(adata, diff_kinetics=True)
scv.tl.velocity_graph(adata,n_jobs=10)
scv.pl.velocity_embedding(adata, basis='umap',palette=cols,arrow_length=3,arrow_size=2,dpi=120)
scv.pl.velocity_embedding_grid(adata, basis='umap',palette=cols,arrow_length=2,arrow_size=1,dpi=120)
scv.pl.velocity_embedding_stream(adata, basis='umap',palette=cols)

#export (for Seurat)
velodf = adata.obs[['velocity_length','velocity_confidence','velocity_pseudotime']]
velodf.to_csv("sc.velo.txt", sep="\t", float_format = '%5.10f')

with open('adata_6','wb') as adata_file:
  pickle.dump(adata,adata_file)



#dynamic modeling
adata=scv.utils.merge(seurat,ldata_all)
adata.obs["clusters"]=adata.obs["my.clusters3"]
adata.obs.clusters
print(sorted(adata.obs.clusters.unique()))
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata,n_jobs=25)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata,n_jobs=25)
scv.pl.velocity_embedding_stream(adata, basis='umap',palette=cols)

#kinetic rates
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
scv.get_df(adata, 'fit*', dropna=True).head()
#latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)
var_names = ["Sostdc1","Tcf7"]
scv.pl.scatter(adata, var_names, frameon=False)
scv.pl.scatter(adata, x='latent_time', y=var_names, frameon=False)

with open('adata_7','wb') as adata_file:
  pickle.dump(adata,adata_file)

#cluster specifics
scv.tl.rank_dynamical_genes(adata, groupby='clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)
for cluster in ["Sostdc1_1","Sostdc1_2"]:
    scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, frameon=False)




