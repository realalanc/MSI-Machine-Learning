import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt

data=pd.read_csv("DESeq2_scored.csv",header=0,index_col=0).T

alist=['MSIscore']
with open("old_result/ttest_deseq2.txt","r") as f:
    while True:
        namet=f.readline() 
        if not namet:
            break
        name=namet.strip().split(' ')[0]
        avg=float(namet.strip().split(' ')[1][14:])
        p1=float(f.readline().strip().split('  ')[1][6:])
        p2=float(f.readline().strip().split('  ')[1][6:])   
        p3=float(f.readline().strip().split('  ')[1][6:])
        if p1=='nan' or p2=='nan' or p3=='nan': continue
        if avg>=2.0:
            alist.append(name)
datat=data[alist]
model = cdt.causality.graph.CCDr()
output_graph=model.predict(datat)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
print(nx.adjacency_matrix(output_graph).todense())
gmlgenerator=nx.generate_gml(output_graph)
gmllist=[]
for i in gmlgenerator:
    gmllist.append(i)
gmlstring="".join(gmllist)
with open("graph.txt","w") as f:
    f.write(gmlstring)