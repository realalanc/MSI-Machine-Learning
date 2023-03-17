import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt

data=pd.read_csv("DESeq2.csv",header=0,index_col=0).T
datat=data[['CLDN11','HGF','CDH1','SEZ6','MLH1','DOCK3','SEZ6L','VTN','EPM2A','CD274','SMAD9','ACVR2A','F12','VEGFC','ALK','TCEA2','CTSF','SEZ6L2','SULT1C4','MSI','ZNF43']]
print(datat)
#glasso = cdt.independence.graph.Glasso()
#skeleton = glasso.predict(datat)
#new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
model = cdt.causality.graph.GES()
output_graph = model.predict(datat)
output_graph=model.predict(datat)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
print(nx.adjacency_matrix(output_graph).todense())