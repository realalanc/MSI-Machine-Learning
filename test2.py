import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt

od=pd.read_csv("origin_data.csv",header=0,index_col=0).T
od = pd.DataFrame(od.values, index=od.index, columns=od.columns)
od=od.drop(['MSI','AJCC Stages','Primary Site','Age'],axis=1)
model = cdt.causality.graph.GIES()
output_graph = model.predict(od)
nx.draw_networkx(output_graph, font_size=10)
plt.show()