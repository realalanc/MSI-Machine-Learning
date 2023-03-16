import pandas as pd
import cdt
import networkx as nx

data=pd.read_csv("chosen_data.csv")
model = cdt.causality.graph.Inter_IAMB()
output_graph = model.predict(data)
print(nx.adjacency_matrix(output_graph).todense())