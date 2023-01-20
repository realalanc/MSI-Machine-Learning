import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt
data=pd.read_csv("chosen_data.csv")
datat=data[['CLDN11','HGF','CDH1','SEZ6','MLH1','DOCK3','SEZ6L','VTN','EPM2A','CD274','SMAD9','ACVR2A','F12','VEGFC','ALK','TCEA2','CTSF','SEZ6L2','SULT1C4','MSS','MSI-H','MSI-L','ZNF43']]
glasso = cdt.independence.graph.Glasso()
skeleton = glasso.predict(datat)
nx.draw_networkx(skeleton, font_size=5)
plt.show()