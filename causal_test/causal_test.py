import cdt
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
def pseudo_data_generator(flag=True):
    if not flag:
        graph=nx.read_gml("generated_graph")
    else:
        generator = cdt.data.AcyclicGraphGenerator(causal_mechanism="gp_add", noise='gaussian',
                                                   noise_coeff=0.4,
                                                npoints=500, nodes=20, parents_max=5, expected_degree=3)
        data,graph=generator.generate()
        print(data)
        generator.to_csv('generated_graph')
        nx.write_gml(graph,"generated_graph")
    nx.draw_networkx(graph)
    plt.savefig("generated_graph.png",format="PNG")
    plt.show()
    return graph

def causal_discovery():
    data=pd.read_csv("generated_graph_data.csv")
    model = cdt.causality.graph.CGNN()  #CGNN is fxxking slow without a flamework, try other solver
    inferred_graph=model.predict(data)
    nx.draw_networkx(inferred_graph)
    plt.savefig("inferred_graph.png", format="PNG")
    plt.show()
    return inferred_graph

def distance(gt,infer): #loss used in a paper
    correct_true_true=0  #'true' indentified as 'true'
    correct_false_false=0
    for node1 in nx.nodes(gt):
        for node2 in nx.nodes(gt):
            if gt.has_edge(node1,node2) and infer.has_edge(node1,node2):
                correct_true_true+=1
            elif not(gt.has_edge(node1,node2)) and not(infer.has_edge(node1,node2)):
                correct_false_false+=1
    sensitivity=correct_true_true/nx.number_of_edges(gt)
    specificity=correct_false_false/(nx.number_of_nodes(gt)*(nx.number_of_nodes(gt)-1)-nx.number_of_edges(gt))#n*(n-1)-edges
    d=((1-sensitivity)**2+(1-specificity)**2)**0.5
    print(d)
    return d




generated_graph=pseudo_data_generator(False)# true(defalut) to re-generate data, otherwise load data
inferred_graph=causal_discovery()
distance(generated_graph,inferred_graph)
