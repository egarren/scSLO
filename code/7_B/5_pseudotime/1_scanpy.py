
import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import rpy2.robjects as robjects
import anndata2ri


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/scanpy.h5ad'
anndata2ri.activate()
robjects.r['load']('sce.filtered.RData')
adata = robjects.r('as(obj.sce, "SingleCellExperiment")')
#compute PAGA path 
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
sc.tl.paga(adata, groups='my.clusters2')
sc.pl.paga(adata,save=True)
#plot PAGA path
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    palette=["#F8766D", "#E38900" ,"#C49A00", "#99A800", "#53B400", "#00BC56" ,"#00C094","#00BFC4", "#00B6EB", "#06A4FF" ,"#A58AFF", "#DF70F8" ,"#FB61D7" ,"#FF66A8"], #"#A58AFF" ,
    legend_fontsize=0, fontsize=0, frameon=False, edges=True, save=True)
#plot pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['my.clusters2']  == 'Naive')[0]
sc.tl.dpt(adata)
dpt = adata.obs.dpt_pseudotime
dpt.to_csv("scanpy.dpt.txt", sep="\t", float_format = '%5.10f')
#Diffmap
sc.tl.diffmap(adata)
diffmap = adata.obsm['X_diffmap']
np.savetxt("scanpy.diffmap.txt", diffmap, '%5.10f',delimiter = "\t")


#gene changes along PAGA path
adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain_anno']  == '16/Stem')[0] #root cell for diffusion
sc.tl.dpt(adata)
gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
              'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
              'Irf8', 'Csf1r', 'Ctsg']                     # monocyte
adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
         ('neutrophils', [16, 0, 4, 2, 14, 19]),
         ('monocytes', [16, 0, 4, 11, 1, 9, 24])]
adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['louvain_anno']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['louvain_anno_colors']
os.mkdir("write")
_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('./write/paga_path_{}.csv'.format(descr))
pl.savefig('./figures/paga_path_paul15.pdf')
pl.show()
       

