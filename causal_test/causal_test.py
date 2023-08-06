import cdt
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split


def pseudo_data_generator(flag=True):
    if not flag:
        graph = nx.read_gml("generated_graph")
    else:
        generator = cdt.data.AcyclicGraphGenerator(causal_mechanism="gp_add", noise='gaussian',
                                                   noise_coeff=0.4,
                                                   npoints=500, nodes=20, parents_max=5, expected_degree=3)
        data, graph = generator.generate()
        print(data)
        generator.to_csv('generated_graph')
        nx.write_gml(graph, "generated_graph")
    nx.draw_networkx(graph)
    plt.savefig("generated_graph.png", format="PNG")
    plt.show()
    return graph


def causal_discovery(data_name):
    data = pd.read_csv(data_name)
    model = cdt.causality.graph.PC()
    inferred_graph = model.predict(data)
    nx.draw_networkx(inferred_graph)
    plt.savefig("inferred_graph.png", format="PNG")
    plt.show()
    return inferred_graph


def distance(gt, infer):  # loss used in a paper
    correct_true_true = 0  # 'true' identified as 'true'
    correct_false_false = 0
    for node1 in nx.nodes(gt):
        for node2 in nx.nodes(gt):
            if gt.has_edge(node1, node2) and infer.has_edge(node1, node2):
                correct_true_true += 1
            elif not (gt.has_edge(node1, node2)) and not (infer.has_edge(node1, node2)):
                correct_false_false += 1
    sensitivity = correct_true_true / nx.number_of_edges(gt)
    specificity = correct_false_false / (
            nx.number_of_nodes(gt) * (nx.number_of_nodes(gt) - 1) - nx.number_of_edges(gt))  # n*(n-1)-edges
    d = ((1 - sensitivity) ** 2 + (1 - specificity) ** 2) ** 0.5
    print(d)
    return d


def data_spliter(data_name,data1_name,data2_name,size=0.6, seed=114514):
    data = pd.read_csv(data_name)
    data1, data2 = train_test_split(data, train_size=size, random_state=seed)  # data1 is 0.6
    data1.to_csv(data1_name)
    data2.to_csv(data2_name)

def algorithm(core="core_data.csv",whole_data="generated_graph_data.csv"):
    pass   # I have generated possible framework by ChatGPT in MSInew.py, it may help


generated_graph = pseudo_data_generator(False)  # true(default) to re-generate data, otherwise load data
data_spliter("generated_graph_data.csv","visible_data.csv","hidden_data.csv")
data_spliter("visible_data.csv","core_data.csv","adding_data.csv",size=0.4/0.6)
gt_graph = causal_discovery(data_name="visible_data.csv")
inferred_graph=algorithm()
distance(gt_graph, inferred_graph)
