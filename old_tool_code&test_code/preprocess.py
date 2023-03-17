import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc

graph=cdt.utils.io.read_list_edges("target.csv")
od=pd.read_csv("origin_data.csv",header=0,index_col=0).T
od = pd.DataFrame(od.values, index=od.index, columns=od.columns)
od=od.drop(['MSI','AJCC Stages','Primary Site','Age'],axis=1)
adata=ad.AnnData(od)
adata.var_names_make_unique()
adata.obs_names=od.index.values.tolist()
adata.var_names=od.columns.values.tolist()

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=100)
adata = adata[:, adata.var.highly_variable]
with open("hvg.txt","w") as f:
    for i in adata.var_names:
        f.write(i+"\n")
    f.write("MSS\nMSI-L\nMSI-H")