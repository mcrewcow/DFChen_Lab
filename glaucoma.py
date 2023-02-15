#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
get_ipython().system('{sys.executable} -m pip -q install palantir fa2')


# In[28]:


import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')

sc.set_figure_params()
scf.set_figure_pubready()


# In[35]:


adata1168 = sc.read_10x_mtx(r'/mnt/c/Users/Emil/10X/1168', cache = True)


# In[37]:


adata1167 = sc.read_10x_mtx('/mnt/c/Users/Emil/10X/1167', cache = True)


# In[38]:


adata1168.obs['stage'] = '1168'
adata1167.obs['stage'] = '1167'


# In[39]:


adataall = AnnData.concatenate(adata1168, adata1167)


# In[40]:


adataall


# In[41]:


adata1168


# In[42]:


adata1167


# In[43]:


adataall.obs_names_make_unique()


# In[44]:


sc.pp.filter_genes(adataall,min_cells=3)
sc.pp.normalize_total(adataall)
sc.pp.log1p(adataall,base=10)
sc.pp.highly_variable_genes(adataall)


# In[45]:


sc.pp.pca(adataall)
pca_projections = pd.DataFrame(adataall.obsm["X_pca"],index=adataall.obs_names)


# In[46]:


dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[47]:


# generate neighbor draph in multiscale diffusion space
adataall.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adataall,n_neighbors=30,use_rep="X_palantir")


# In[48]:


adataall.obsm["X_pca2d"]=adataall.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adataall,init_pos='X_pca2d')


# In[49]:


sc.tl.leiden(adataall)


# In[53]:


sc.set_figure_params()
sc.pl.draw_graph(adataall,color=["stage",'Cx3cr1']) #PCA-based, n_neighbors = 50, n_eigs = 4


# In[58]:


sc.set_figure_params(figsize = (8,8))
sc.pl.draw_graph(adataall,color=["leiden"], legend_loc = 'on data')


# In[60]:


sc.set_figure_params()
sc.tl.rank_genes_groups(adataall, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adataall, n_genes=25, sharey=False)


# In[64]:


sc.set_figure_params(figsize = (8,8))
sc.pl.draw_graph(adataall,color=["Emcn",'Pdgfrb'], legend_loc = 'on data') #endotheliocytes, podocytes


# In[66]:


adataallsubset = adataall[adataall.obs['leiden'].isin(['18','3','10','21','12','22','26','7','0','8','2','11','1','9','5','4','6','13'])]


# In[67]:


sc.set_figure_params(figsize = (8,8))
sc.pl.draw_graph(adataallsubset,color=["leiden"], legend_loc = 'on data')


# In[69]:


sc.pp.pca(adataallsubset)
pca_projections = pd.DataFrame(adataallsubset.obsm["X_pca"],index=adataallsubset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
# generate neighbor draph in multiscale diffusion space
adataallsubset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adataallsubset,n_neighbors=30,use_rep="X_palantir")

adataallsubset.obsm["X_pca2d"]=adataallsubset.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adataallsubset,init_pos='X_pca2d')
sc.tl.leiden(adataallsubset)
sc.set_figure_params()
sc.pl.draw_graph(adataallsubset,color=["stage",'Cx3cr1']) #PCA-based


# In[75]:


sc.tl.umap(adataallsubset) #PCA-based


# In[76]:


sc.pl.umap(adataallsubset,color=["stage",'Cx3cr1']) #PCA-based


# In[77]:


sc.set_figure_params(figsize = (8,8))
sc.pl.draw_graph(adataallsubset,color=["leiden"], legend_loc = 'on data')


# In[85]:


sc.set_figure_params(figsize = (8,8))
sc.pl.draw_graph(adataallsubset,color=["Gfap",'Rlbp1'], legend_loc = 'on data')


# In[97]:


sc.set_figure_params()
sc.pl.draw_graph(adataallsubset,color=["Apoe",'stage'], legend_loc = 'on data')


# In[78]:


sc.set_figure_params(figsize = (8,8))
sc.pl.umap(adataallsubset,color=["leiden"], legend_loc = 'on data')


# In[79]:


sc.set_figure_params()
sc.tl.rank_genes_groups(adataallsubset, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adataallsubset, n_genes=25, sharey=False)


# In[98]:


adataallsubset = adataallsubset[adataallsubset.obs['leiden'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','24'])]


# In[99]:


sc.set_figure_params()
sc.pl.draw_graph(adataallsubset,color=['Cx3cr1',"Apoe",'stage'], legend_loc = 'on data')


# In[101]:


sc.pp.pca(adataallsubset)
pca_projections = pd.DataFrame(adataallsubset.obsm["X_pca"],index=adataallsubset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
# generate neighbor draph in multiscale diffusion space
adataallsubset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adataallsubset,n_neighbors=30,use_rep="X_palantir")

adataallsubset.obsm["X_pca2d"]=adataallsubset.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adataallsubset,init_pos='X_pca2d')
sc.tl.leiden(adataallsubset)
sc.set_figure_params()
sc.pl.draw_graph(adataallsubset,color=["stage",'Cx3cr1','Apoe','Mki67']) #PCA-based


# In[102]:


sc.pl.draw_graph(adataallsubset,color=["stage",'Cx3cr1','Apoe']) #PCA-based


# In[104]:


sc.pl.draw_graph(adataallsubset,color=["Mki67"], vmax = 0.5) #PCA-based


# In[105]:


sc.tl.score_genes(adataallsubset, gene_list = ['Cxcl10','Stat1','Ifit1'], score_name = 'IFN-response')


# In[106]:


sc.tl.score_genes(adataallsubset, gene_list = ['Apoe','Igf1','H2-Aa','Cd74'], score_name = 'Inflammatory')


# In[110]:


sc.tl.score_genes(adataallsubset, gene_list = ['Mki67','Cdc20'], score_name = 'Proliferating')


# In[107]:


sc.tl.score_genes(adataallsubset, gene_list = ['P2ry12', 'Ccr5', 'Cd33', 'Csf1r', 'Cx3cr1', 'Glul', 'Gpr34', 'Adgrg1', 'Tgfb1', 'Tgfbr1', 'Serinc3', 'Siglech', 'Mertk', 'Bin1', 'Tmem119', 'Spi1', 'Sall1', 'Mafb', 'Smad3', 'Mef2a', 'Egr1', 'Jun'], score_name = 'Homeostatic')


# In[108]:


sc.tl.score_genes(adataallsubset, gene_list = ['Apoe', 'Axl', 'Clec7a', 'Csf1', 'Cst7', 'Ctsb', 'Ctsd', 'Cybb', 'Fabp5', 'Fth1', 'Gnas', 'Gpnmb', 'Grn', 'Il1b', 'Lgals3', 'Lilrb4a', 'Lpl', 'Lyz2', 'Msr1', 'Spp1', 'Trem2', 'Tyrobp', 'Vegfa', 'H2-Ab1'], score_name = 'Proinflammatory')


# In[113]:


sc.pl.draw_graph(adataallsubset,color=['Homeostatic','Proinflammatory'], vmax = [0.7,1])


# In[115]:


sc.pl.draw_graph(adataallsubset,color=['Proliferating','Inflammatory','IFN-response'], vmax = [0.5,2,1])


# In[117]:


sc.tl.leiden(adataallsubset)


# In[118]:


AnnData.write_h5ad(adataallsubset, "/mnt/c/Users/Emil/10X/microglia_glaucoma.h5ad")


# In[126]:


adataSeurat = adataallsubset


# In[130]:


adataSeurat.raw = adataSeurat


# In[131]:


adataSeurat


# In[132]:


adataSeurat.raw


# In[133]:


t=adataSeurat.raw.X.toarray()
pd.DataFrame(data=t, index=adataSeurat.obs_names, columns=adataSeurat.raw.var_names).to_csv('/mnt/c/Users/Emil/10X/microglia_glaucomaraw-counts.csv')
# comma separated csv with header, no rownames
pd.DataFrame(adataSeurat.obs).to_csv("/mnt/c/Users/Emil/10X/microglia_glaucomametadata.csv")


# In[135]:


sc.pl.draw_graph(adataallsubset,color=['leiden'], legend_loc = 'on data')


# In[140]:


old_to_new = {
'5':'Inflammatory',
'14':'Inflammatory',
'12':'Inflammatory',
'21':'Inflammatory',
'7':'Inflammatory',
'23':'Proliferating',
'17':'IFN-response',
'9':'IFN-response',
'2':'IFN-response',
'16':'IFN-response',
'3':'Resident',
'22':'Resident',
'19':'Resident',
'13':'Resident',
'24':'Resident',
'11':'Resident',
'10':'Resident',
'18':'Resident',
'8':'Resident',
'15':'Resident',
'1':'Resident',
'4':'Resident',
'0':'Resident',
'6':'Resident',
'20':'Resident'}
adataallsubset.obs['annotation'] = (
adataallsubset.obs['leiden']
.map(old_to_new).astype('category')
)


# In[141]:


sc.pl.draw_graph(adataallsubset,color=['stage','annotation'], legend_loc = 'on data')


# In[153]:


adataallsubset_1167 = adataallsubset[adataallsubset.obs['stage'].isin(['1167'])]


# In[155]:


adataallsubset_1167.obs


# In[154]:


adataallsubset_1168 = adataallsubset[adataallsubset.obs['stage'].isin(['1168'])]


# In[157]:


adataallsubset_1167.obs.annotation.value_counts()


# In[158]:


adataallsubset_1168.obs.annotation.value_counts()


# In[160]:


sc.pl.draw_graph(adataallsubset,color=['Cx3cr1'], vmax = 0.5)


# In[162]:


import scvelo as scv


# In[167]:


loom67 = scv.read("/mnt/c/Users/Emil/10X/D23-1167.loom", cache = True)


# In[168]:


loom68 = scv.read("/mnt/c/Users/Emil/10X/D23-1168.loom", cache = True)


# In[169]:


adataallsubset_1167 = scv.utils.merge(adataallsubset_1167, loom67)


# In[170]:


adataallsubset_1168


# In[171]:


adataallsubset_1168 = scv.utils.merge(adataallsubset_1168, loom68)


# In[172]:


adataallsubset_1168


# In[173]:


scv.pl.proportions(adataallsubset_1167)


# In[174]:


scv.pl.proportions(adataallsubset_1168)


# In[175]:


scv.pp.filter_and_normalize(adataallsubset_1167, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adataallsubset_1167, n_pcs=30, n_neighbors=30)


# In[176]:


scv.pp.filter_and_normalize(adataallsubset_1168, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adataallsubset_1168, n_pcs=30, n_neighbors=30)


# In[177]:


scv.tl.velocity(adataallsubset_1167)
scv.tl.velocity(adataallsubset_1168)


# In[179]:


scv.pp.neighbors(adataallsubset_1167)


# In[180]:


scv.tl.velocity_graph(adataallsubset_1167)


# In[181]:


scv.pp.neighbors(adataallsubset_1168)


# In[182]:


scv.tl.velocity_graph(adataallsubset_1168)


# In[184]:


scv.pl.velocity_embedding_stream(adataallsubset_1167, basis='draw_graph_fa', color = 'annotation', dpi = 300, add_margin = 0.2, figsize = (8,8), arrow_size = 2, density = 6)


# In[185]:


scv.pl.velocity_embedding_stream(adataallsubset_1168, basis='draw_graph_fa', color = 'annotation', dpi = 300, add_margin = 0.2, figsize = (8,8), arrow_size = 2, density = 6)


# In[202]:


scv.pl.velocity_embedding_stream(adataallsubset_1168, basis='draw_graph_fa', color = 'annotation', dpi = 300, add_margin = 0.2, figsize = (8,8), arrow_size = 2, density = 6)


# In[203]:


sc.pl.dotplot(adataallsubset, var_names = ['Tlr4','Igf1','Igf1r'], groupby = ['annotation','stage'])


# In[193]:


sc.pl.dotplot(adataallsubset, var_names = ['Igf1'], groupby = ['annotation','stage'])


# In[194]:


sc.pl.dotplot(adataallsubset, var_names = ['Trem2'], groupby = ['annotation','stage'])


# In[188]:


adataallsubset.obs


# In[204]:


AnnData.write_h5ad(adataallsubset_1168, "/mnt/c/Users/Emil/10X/adataallsubset_1168.h5ad")


# In[205]:


AnnData.write_h5ad(adataallsubset_1167, "/mnt/c/Users/Emil/10X/adataallsubset_1167.h5ad")


# In[223]:


scv.pl.velocity_embedding_stream(adataallsubset_1167, basis='draw_graph_fa', color = 'annotation', palette = ['blue','green','red'], dpi = 300, density = 4, add_margin=0.2, figsize = (8,8), arrow_size = 2,save="figures/glaucoma1167.svg")


# In[224]:


scv.pl.velocity_embedding_stream(adataallsubset_1168, basis='draw_graph_fa', color = 'annotation', palette = ['blue','green','magenta','red'], dpi = 300, density = 4, add_margin=0.2, figsize = (8,8), arrow_size = 2,save="figures/glaucoma1168.svg")


# In[225]:


scv.pl.velocity_embedding_grid(adataallsubset_1167, basis='draw_graph_fa', arrow_length = 5, arrow_size = 3, dpi=120, figsize = (6,6), add_margin = 0.2,save="figures/glaucoma1167_gray.svg")


# In[226]:


scv.pl.velocity_embedding_grid(adataallsubset_1168, basis='draw_graph_fa', arrow_length = 2, arrow_size = 3, dpi=120, figsize = (6,6), add_margin = 0.2,save="figures/glaucoma1168_gray.svg")


# In[232]:


sc.pl.draw_graph(adataallsubset,color=['Homeostatic','Proinflammatory'], vmax = [1,1], save = 'glaucoma_patterns_1.pdf')


# In[233]:


sc.pl.draw_graph(adataallsubset,color=['IFN-response','Inflammatory','Proliferating'], vmax = [1,2,0.5], save = 'glaucoma_patterns_2.pdf')


# In[234]:


WT = sc.read_h5ad('/mnt/c/Users/Emil/10X/microglia_WTKO_basic_and_velocity_WT.h5ad')


# In[235]:


KO = sc.read_h5ad('/mnt/c/Users/Emil/10X/microglia_WTKO_basic_and_velocity_KO.h5ad')


# In[240]:


microgliaChen = sc.read_h5ad("/mnt/c/Users/Emil/10X/microglia.h5ad")


# In[244]:


sc.pl.draw_graph(microgliaChen,color=['leiden'], legend_loc = 'on data')


# In[245]:


microgliaChen = microgliaChen[microgliaChen.obs['leiden'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'])]


# In[246]:


sc.pl.draw_graph(microgliaChen,color=['leiden'], legend_loc = 'on data')


# In[247]:


adatasummary = AnnData.concatenate(adataallsubset, microgliaChen)


# In[248]:


adatasummary


# In[249]:


adatasummary.obs.stage.value_counts()


# In[251]:


sc.pl.dotplot(adataallsubset_1167, var_names = ['Trem2'], groupby = ['annotation','stage'])


# In[261]:


adatasummary.raw = adatasummary


# In[263]:


t=adatasummary.raw.X.toarray()
pd.DataFrame(data=t, index=adatasummary.obs_names, columns=adatasummary.raw.var_names).to_csv('/mnt/c/Users/Emil/10X/11raw-counts.csv')
# comma separated csv with header, no rownames
pd.DataFrame(adatasummary.obs).to_csv("/mnt/c/Users/Emil/10X/11metadata.csv")


# In[ ]:




