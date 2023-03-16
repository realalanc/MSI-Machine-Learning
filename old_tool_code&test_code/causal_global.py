import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt

graph=cdt.utils.io.read_list_edges("target.csv")
od=pd.read_csv("origin_data.csv",header=0,index_col=0).T
od = pd.DataFrame(od.values, index=od.index, columns=od.columns)
#od=od.drop(['MSI','AJCC Stages','Primary Site','Age'],axis=1)
#model = cdt.causality.graph.bnlearn.GS()
t=[]
with open("hvg.txt","r") as f:
    while True:
        line=f.readline()
        if not line:
            break
        else:
            t.append(line.strip())
od=od[t]
model = cdt.causality.graph.GES()
output_graph = model.predict(od, graph)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
