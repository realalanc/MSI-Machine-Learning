import sys
import getopt
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import heapq
import pandas as pd
from os import mkdir
import scipy
import math
'''
MSItool是一个研究性状与基因的相关性与因果性相关的软件。\n\
option:\n\
-m,--name     待研究的性状名称，默认为MSI\n\
-o,--output   输出文件夹的名称，默认和输入文件相同\n\
-c,--corr     correlation:0:pass 1:Pearson 2:spearman\n \
-a,--alg      algorithm:0:pass 1:CCDr 2:PC 3:GES\n \
-e,--evl      evl:0:pass 1:backdoor 2:instrumental variable\n \
arg1:         表达矩阵文件，格式为csv\n\
arg2:         在不进行因果发现直接进行dowhy时需要传入，传入待检验的因果图，格式为gml\n\
-h,--help     帮助
-v,--version  版本号
'''

def exit_help():
    str="\
MSItool是一个研究性状与基因的相关性与因果性相关的软件。\n\
option:\n\
-m,--name     待研究的性状名称，默认为MSI\n\
-o,--output   输出文件夹的名称，默认和输入文件相同\n\
-c,--corr     correlation:0:pass 1:Pearson 2:spearman\n \
-a,--alg      algorithm:0:pass 1:CCDr 2:PC 3:GES\n \
-e,--evl      evl:0:pass 1:backdoor 2:instrumental variable\n \
arg1:         表达矩阵文件，格式为csv\n\
arg2:         在不进行因果发现直接进行dowhy时需要传入，传入待检验的因果图，格式为gml\n\
-h,--help     帮助\n\
-v,--version  版本号\
    "
    print(str)
    exit()

version="0.0.1"
correlation=0
algorithm=0
evalution=0
filename=""
data=0
msiname="MSI"
args=[]
opts=[]
try:
    opts,args=getopt.gnu_getopt(sys.argv[1:],"-m:-o:-v-h-c:-e:-a",["name=","output=","version","help","corr=","alg=","evl="])
except getopt.GetoptError as err:
    print(str(err))
    exit_help()
def readdata():
    try:
        filename=args[0]
        data=pd.read_csv("chosen_data.csv")
    except:
        print("you need to choose a csv file as input")
        exit()
#print(opts)
#print(args)
graphname=""
if len(args)==2:
    graphname=args[1]

for o,a in opts:
    if o in ("-h","--help"):
        exit_help()
    if o in ("-v","--version"):
        print(version)
        exit()
    if o in ("-a","--alg"):
        algorithm=int(a)
    if o in ("-e","--evl"):
        evalution=int(a)
    if o in ("-o","--output"):
        writedir=a
    if o in ("-m","--msi"):
        msiname=a
    if o in ("-c","--corr"):
        correlation=int(a)

readdata()
if msiname in data.index:
    data=data.T
if msiname not in data.column:
    print("invalid msi column(or index) name")
    exit()

if writedir=="":
    writedir=filename.split(".")[0]
try:
    mkdir("writedir")
except:
    pass

def getcorrelation():
    nrows=data.columns.get_loc(msiname)+1
    if correlation==1:
        dcorr=data.corr()
        dcorr.to_csv(writedir+'/corr_result.csv')
        # dcorr=pd.read_csv("corr_result.csv",header=0,index_col=0)
        a=np.triu(dcorr,1)
        a=a.flatten()
        topn=heapq.nlargest(n, range(len(a)), a.take)
        with open(writedir+'/corrdata.txt', 'w') as f:
            for i in topn:
                f.write((dcorr.index)[(i+1)//(nrows-1)]);f.write("    ")
                f.write((dcorr.index)[(i+1)%(nrows-1)-1]);f.write("    ")
                f.write(f"{dcorr.iloc[(i+1)//(nrows-1),(i+1)%(nrows-1)-1]}\n")
    elif correlation==2:
        # print(data.iloc[:,nrows-1])
        corr=[];o_corr=[];o_index=[];o_pv=[]
        for i in range(0,nrows):
            dcorr,p=scipy.stats.spearmanr(data.iloc[:, nrows-1], data.iloc[:, i])
            if not math.isnan(dcorr):
                corr.append(dcorr)
                o_corr.append(dcorr)
                o_index.append(i)
                o_pv.append(p)
        corr.sort(reverse=True)
        with open(writedir+'/corrdata.txt', 'w') as f:
            for i in range(0,len(corr)):
                f.write(data.columns[o_index[o_corr.index(corr[i])]]);f.write("    p=")
                f.write(f"{o_pv[o_corr.index(corr[i])]}");f.write("    ")
                f.write(f"{corr[i]}\n")
    with open(writedir+'/kendall_result.txt', 'w') as f:
            f.write('kendalltau  r and p\n\n')
            r, p = scipy.stats.kendalltau(data.iloc[:, nrows - 1], data.iloc[:, nrows])
            f.write('Primary Site and MSI\n');f.write(f"{r}");f.write("    ");f.write(f"{p}\n\n")
def getgraph():
    import cdt
    glasso = cdt.independence.graph.Glasso()
    skeleton = glasso.predict()
    new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
    if(algorithm==1):
        model = cdt.causality.graph.CCDr()
    elif algorithm==2:
        model=cdt.causality.graph.PC()
    elif algorithm==3:
        model=cdt.causality.graph.GES()
    output_graph=model.predict(data)
    nx.draw_networkx(output_graph, font_size=5)
    plt.savefig(writedir+"/graph.png",format="PNG")
    plt.show()
    gmlgenerator=nx.generate_gml(output_graph)
    gmllist=[]
    for i in gmlgenerator:
        gmllist.append(i)
    gmlstring="".join(gmllist)
    with open(writedir+"/graph.txt","w") as f:
        f.write(gmlstring)
    return output_graph

def dowhytest(graph):
    import dowhy
    if not nx.is_directed_acyclic_graph(graph):
        print("input graph should be dag")
        exit()
    gmlgenerator=nx.generate_gml(graph)
    gmllist=[]
    for i in gmlgenerator:
        gmllist.append(i)
    gmlstring="".join(gmllist)

    print("please input treatment and outcome to dowhy test")
    tt=input()
    oc=input()
    try:
        model_dowhy=dowhy.causal_model.CausalModel(data,treatment="tt",outcome="oc",graph=gmlstring)
        identified_estimand = model_dowhy.identify_effect()
        with open(writedir+"/dowhy_identified_estimand.txt","w") as f:
            f.write(identified_estimand)
        print("result of identify effect has been written into identified_estimand.txt")
    except:
        print("invalid input edge!")
        exit()
    if evalution==1:
        estimate = model_dowhy.estimate_effect(identified_estimand, method_name="backdoor.linear_regression")
        #print(estimate)
        refute1_results=model_dowhy.refute_estimate(identified_estimand, estimate,
            method_name="random_common_cause")
        #print(refute1_results)
        refute2_results=model_dowhy.refute_estimate(identified_estimand, estimate,
            method_name="placebo_treatment_refuter")
        #print(refute2_results)
        refute3_results=model_dowhy.refute_estimate(identified_estimand, estimate,
            method_name="data_subset_refuter")
        #print(refute3_results)
        with open(writedir+"/evalution.txt","w") as f:
            f.write("estimate:\n")
            f.write(estimate)
            f.write("\n----------------------------\n----------------------------\n")
            f.write("random common cause method result:\n")
            f.write(refute1_results)
            f.write("\n----------------------------\n----------------------------\n")
            f.write("placebo treatment refuter method result\n")
            f.write(refute2_results)
            f.write("\n----------------------------\n----------------------------\n")
            f.write("data subset refuter method result\n")
            f.write(refute3_results)
    elif evalution==2:
        estimate_iv=model_dowhy.estimate_effect(identified_estimand, method_name="iv.instrumental_variable")
        #print(estimate)
        refute1_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
            method_name="random_common_cause")
        #print(refute1_results)
        refute2_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
            method_name="data_subset_refuter")
        #print(refute3_results)
        with open(writedir+"/evalution.txt","w") as f:
            f.write("estimate:\n")
            f.write(estimate_iv)
            f.write("\n----------------------------\n----------------------------\n")
            f.write("random common cause method result:\n")
            f.write(refute1_results)
            f.write("\n----------------------------\n----------------------------\n")
            f.write("data subset refuter method result\n")
            f.write(refute2_results)

if correlation!=0:
    getcorrelation()
    exit()
if algorithm!=0:
    graph=getgraph()
else:
    if graphname=="":
        exit()
    else:
        graph=nx.read_gml(graphname,label='label')
if dowhytest!=0:
    dowhytest(graph)
