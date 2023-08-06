import cdt
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split


def pseudo_data_generator(num=100,rate=0.6,flag=True,is_noise=False): # rate means how many nodes are the ones to be researched
    # notice! nodes=num*rate shall be an int
    if not flag:
        graph_causal = nx.read_gml("generated_graph")
    else:
        #selected
        generator_causal = cdt.data.AcyclicGraphGenerator(causal_mechanism="gp_add", noise='gaussian',
                                                   noise_coeff=0.4,
                                                   npoints=500, nodes=int(num*rate), parents_max=5, expected_degree=3)
        data_causal, graph_causal = generator_causal.generate()
        #noise
        generator_noise = cdt.data.AcyclicGraphGenerator(causal_mechanism="gp_add", noise='gaussian',
                                                   noise_coeff=0.4,
                                                   npoints=500, nodes=int(num*(1-rate)), parents_max=5, expected_degree=3)
        data_noise,graph_noise=generator_noise.generate()

        #print(data)
        generator_causal.to_csv('generated_graph')
        generator_noise.to_csv('noise')
        data1=pd.read_csv("generated_graph_data.csv")
        data2=pd.read_csv("noise_data.csv")

        list_nodes=data2.columns.to_list()
        print(list_nodes)
        list_nodes_fixed=[]
        for i in list_nodes:
            list_nodes_fixed.append(i+"'")
        print(list_nodes_fixed)
        #data2.rename(columns=dict(zip(list_nodes, list_nodes_fixed)), inplace=True)# solve columns' name conflict
        data2.rename(columns={i:j for i,j in zip(list_nodes,list_nodes_fixed)},inplace=True)
        print(data2)
        merged_data=pd.concat([data1,data2],axis=1)
        print(merged_data)         # test
        merged_data.to_csv("merged_data.csv")  # merge the data between selected nodes and noise nodes
        nx.write_gml(graph_causal, "generated_graph") # only save the selected nodes

    nx.draw_networkx(graph_causal)
    plt.savefig("generated_graph.png", format="PNG")
    plt.show()    # print it
    return graph_causal


def causal_discovery(data_name):      # may not use
    data = pd.read_csv(data_name)
    model = cdt.causality.graph.PC()     # ok to change the solver
    inferred_graph = model.predict(data)
    nx.draw_networkx(inferred_graph)
    plt.savefig("inferred_graph.png", format="PNG")
    plt.show()
    return inferred_graph


def distance(gt, infer):  # loss used in a paper, may be improved
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


def data_spliter(data_name,data1_name,data2_name,rate=0.6, seed=114514):
    data = pd.read_csv(data_name)
    data=data.T
    print(data)
    data1, data2 = train_test_split(data, train_size=rate, random_state=seed)  # data1 is 0.6

    data1=data1.T
    data2=data2.T
    print(data1)
    print(data2)
    data1.to_csv(data1_name)
    data2.to_csv(data2_name)

def algorithm(core="core_data.csv",whole_data="generated_graph_data.csv"):
    pass   # I have generated possible framework by ChatGPT in MSInew.py, it may help, developed later

causality_rate=0.6
gt_graph = pseudo_data_generator()  # true(default) to re-generate data, otherwise load data
data_spliter("generated_graph_data.csv","core_data.csv","adding_data.csv",rate=2/3)
inferred_graph=algorithm()
#distance(gt_graph, inferred_graph)  # may need overwrite
