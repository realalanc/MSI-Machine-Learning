import pandas as pd
import torch as tr
import itertools
import numpy as np
import networkx as nx
from scipy.stats import spearmanr

def read_data():

    import pandas as pd

    # 读取 CSV 文件
    csv_file_path = 'DESeq2_scored.csv'  # 替换为实际的文件路径
    df = pd.read_csv(csv_file_path)
    # 将数据转换为二维数组形式，作为 data 输入
    data = df.to_numpy()       # return as numpy array
    return data

def method_dl(data):
    pass

def method_informed_search(data):
    # 先验因果图的定义
    prior_causal_graph = nx.DiGraph()
    prior_causal_graph.add_edges_from([(1, 3), (2, 3)])

    def conditional_independence_test(data, node1, node2, condition_nodes, alpha=0.05):
        # 提取相关的数据列
        cols = [node1, node2] + condition_nodes
        subset_data = data[:, cols]

        # 计算斯皮尔曼相关系数和相应的p-value
        correlation_coefficient, p_value = spearmanr(subset_data, axis=0)

        # 判断p-value是否小于显著性水平
        return p_value < alpha

    def compute_graph_difference(graph1, graph2):
        # 计算两个因果图之间的差距（不包括节点数量差异）
        return len(set(graph1.edges()) - set(graph2.edges()))

    def generate_causal_graph(prior_graph, data, nodes, alpha=0.05, max_additional_nodes=5, initial_temperature=1.0,
                              cooling_rate=0.95, max_iterations=1000):
        # 生成新的因果图
        new_causal_graph = prior_graph.copy()

        best_causal_graph = new_causal_graph.copy()

        temperature = initial_temperature
        iterations = 0

        while temperature > 0.1 and iterations < max_iterations:
            # 在当前温度下，随机选择一对节点，其中一个节点来自更大节点集合
            node1, node2 = np.random.choice(nodes, 2, replace=False)

            if node1 == node2:
                # 跳过相同的节点
                continue

            if node2 in new_causal_graph[node1] or node1 in new_causal_graph[node2]:
                # 跳过已存在的边和相互连接的节点对
                continue

            # 执行条件独立性检验，使用随机选择的两个节点之外的其他节点作为条件
            condition_nodes = list(set(larger_node_set) - {node1, node2})
            all_independent = True

            # 遍历所有可能的 condition_nodes 组合
            for r in range(len(condition_nodes) + 1):
                for condition_node_set in itertools.combinations(condition_nodes, r):
                    if not conditional_independence_test(data, node1, node2, condition_node_set, alpha):
                        all_independent = False
                        break
                if not all_independent:
                    break

            if all_independent:
                # 所有可能的条件组合都不满足条件独立性，则添加该边
                new_causal_graph.add_edge(node1, node2)

                # 计算因果图差距
                current_difference = compute_graph_difference(prior_causal_graph, new_causal_graph)
                best_difference = compute_graph_difference(prior_causal_graph, best_causal_graph)

                # 根据概率选择是否接受更好的解
                if current_difference <= best_difference or np.random.rand() < np.exp(
                        (best_difference - current_difference) / temperature):
                    best_causal_graph = new_causal_graph.copy()
                else:
                    # 回退，保持当前状态
                    new_causal_graph.remove_edge(node1, node2)

            # 降低温度
            temperature *= cooling_rate
            iterations += 1

        return best_causal_graph

    # 更大节点集合
    larger_node_set = [1, 2, 3, 4, 5]
    nodes = prior_causal_graph.nodes()

    new_causal_graph = generate_causal_graph(prior_causal_graph, data, larger_node_set, alpha=0.05,
                                             max_additional_nodes=5)
    print("新的因果图的边：", new_causal_graph.edges())

def causal_discovery(data,method):
    if method=="deep learning":
        graph=method_dl(data)
    elif method=="informed search":
        graph=method_informed_search(data)
    else:
        pass
    return graph




