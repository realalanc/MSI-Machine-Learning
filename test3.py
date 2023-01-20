import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt
graph=cdt.utils.io.read_list_edges("target.csv")
nx.draw_networkx(graph, font_size=10)
plt.show()