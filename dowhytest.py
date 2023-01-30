import networkx
import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt
import dowhy

data=pd.read_csv("chosen_data.csv")
datat=data[['CLDN11','HGF','CDH1','SEZ6','MLH1','DOCK3','SEZ6L','VTN','EPM2A','CD274','SMAD9','ACVR2A','F12','VEGFC','ALK','TCEA2','CTSF','SEZ6L2','SULT1C4','MSS','MSI-H','MSI-L','ZNF43']]
glasso = cdt.independence.graph.Glasso()
skeleton = glasso.predict(datat)
new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
model = cdt.causality.graph.GES()
output_graph=model.predict(datat,new_skeleton)

while not nx.is_directed_acyclic_graph(output_graph):
    print("No")
    output_graph = model.predict(datat, new_skeleton)
print("Yes")
print(output_graph)
gmlgenerator=nx.generate_gml(output_graph)
gmllist=[]
for i in gmlgenerator:
    gmllist.append(i)
gmlstring="".join(gmllist)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
#print(nx.adjacency_matrix(output_graph).todense())

model_dowhy=dowhy.causal_model.CausalModel(datat,treatment="MLH1",outcome="MSI-H",graph=gmlstring)
model_dowhy.view_model()
identified_estimand = model_dowhy.identify_effect()
print(identified_estimand)
estimate = model_dowhy.estimate_effect(identified_estimand, method_name="backdoor.linear_regression")
print(estimate)
refute1_results=model_dowhy.refute_estimate(identified_estimand, estimate,
        method_name="random_common_cause")
print(refute1_results)
refute2_results=model_dowhy.refute_estimate(identified_estimand, estimate,
        method_name="placebo_treatment_refuter")
print(refute2_results)
refute3_results=model_dowhy.refute_estimate(identified_estimand, estimate,
        method_name="data_subset_refuter")
print(refute3_results)
