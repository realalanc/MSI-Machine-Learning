import pandas as pd
import cdt
import networkx as nx
import matplotlib.pyplot as plt
import dowhy

#前面这里衔接的是马上会push上来的，相关性发现的部分
#大概流程是：按照相关性，用户选n=一个数，然后生成包含前n个的datat，也就是在真正的工具包中，用来跑算法的datat是从数据中生成的，而不是现在这样钦定

#以下是因果发现的部分，两个用户选项就是，用哪个solver（有CCDr和GES可以用），用不用骨架（只有GES能用骨架）
data=pd.read_csv("chosen_data.csv")
datat=data[['CLDN11','HGF','CDH1','SEZ6','MLH1','DOCK3','SEZ6L','VTN','EPM2A','CD274','SMAD9','ACVR2A','F12','VEGFC','ALK','TCEA2','CTSF','SEZ6L2','SULT1C4','MSS-ALL','ZNF43']]
glasso = cdt.independence.graph.Glasso()
skeleton = glasso.predict(datat)
new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
model = cdt.causality.graph.CCDr()
output_graph=model.predict(datat)



#以下一是画图，二是在cdt和dowhy间传递这个图，开发成工具的话应该给一个保存这个图的选项，为了防止用户看不清图，可以同时输出矩阵
gmlgenerator=nx.generate_gml(output_graph)
gmllist=[]
for i in gmlgenerator:
    gmllist.append(i)
gmlstring="".join(gmllist)
nx.draw_networkx(output_graph, font_size=5)
plt.show()
print(nx.adjacency_matrix(output_graph).todense())

#用户看完图后自己输入treatment和outcome分别是什么，也就是这两个string应该是用户输入的
model_dowhy=dowhy.causal_model.CausalModel(datat,treatment="MLH1",outcome="MSS-ALL",graph=gmlstring)

identified_estimand = model_dowhy.identify_effect()
print(identified_estimand)

#下面是后门准则的，这个method可能可以进一步改进，比如结合ecolml之类的别的包
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


#下面是工具变量法的，这里也可以搞成让用户选择的
estimate_iv=model_dowhy.estimate_effect(identified_estimand, method_name="iv.instrumental_variable")
print(estimate_iv)
refute1_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
        method_name="random_common_cause")
print(refute1_results)
'''
refute2_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
        method_name="placebo_treatment_refuter")
print(refute2_results)
'''
#如果在工具变量识别的里面用安慰剂检验就会报错,貌似是这个方法调用几层之后某个方法的参数不能是默认参数，所以我就不会了
refute3_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
        method_name="data_subset_refuter")
print(refute3_results)