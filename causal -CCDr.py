import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt

datat=pd.read_csv("chosendata2.csv",header=0,index_col=0)

#datat=data[['CLDN11','HGF','CDH1','SEZ6','MLH1','DOCK3','SEZ6L','VTN','EPM2A','CD274','SMAD9','ACVR2A','F12','VEGFC','ALK','TCEA2','CTSF','SEZ6L2','SULT1C4','MSS-ALL','ZNF43']]
model = cdt.causality.graph.CCDr()
output_graph=model.predict(datat)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
print(nx.adjacency_matrix(output_graph).todense())