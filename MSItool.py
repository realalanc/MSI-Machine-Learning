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
writedir="" #貌似也是作用域的问题，不在这里写上的话后面就会认为writedir只在那个for里面定义
data=0
msiname="MSI"
args=[]
opts=[]
try:
    opts,args=getopt.gnu_getopt(sys.argv[1:],"-m:-o:-v-h-c:-a:-e:",["name=","output=","version","help","corr=","alg=","evl="])#这里最后那个-e没有加冒号，a和e还写反了
except getopt.GetoptError as err:
    print(str(err))
    exit_help()
def readdata():
    try:
        global filename
        filename=args[0] #还有就是函数内用全局变量要声明global，下面的data同理
        global data
        data=pd.read_csv(filename,header=0,index_col=0,low_memory=False) #这里硬编码死输入的csv名还是搞成输入参数？我这里改了是因为，原来的chosen_data有好多非数项不能跑
        #最重要的一点:数据预处理必须扔到非数项！包括sample的代号，除了index都要是数才行
        data=data.dropna()
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
if msiname not in data.columns: #这里原来写成了column
    print("invalid msi column(or index) name")
    exit()

if writedir=="":
    writedir=filename.split(".")[0]
try:
    mkdir("writedir")
except:
    pass

def getcorrelation():
    global data
    del data['symbol_mrna']
    nrows = data.columns.get_loc("MSI") + 1  # Excel中MSI所在excel中的行-1(除去了表头） 19966
    n = 100
    datat = pd.read_csv("origin_data.csv", header=0, index_col=0, nrows=nrows - 1).T


    if correlation==1:
        dcorr=datat.corr(numeric_only=True)
        print(dcorr)
        dcorr.to_csv('corr_result.csv')#路径的问题和下面一样
        # dcorr=pd.read_csv("corr_result.csv",header=0,index_col=0)
        a=np.triu(dcorr,1)
        a=a.flatten()
        topn=heapq.nlargest(n, range(len(a)), a.take)
        with open('corrdata.txt', 'w') as f:
            for i in topn:
                f.write((dcorr.index)[(i+1)//(nrows-1)]);f.write("    ")
                f.write((dcorr.index)[(i+1)%(nrows-1)-1]);f.write("    ")
                f.write(f"{dcorr.iloc[(i+1)//(nrows-1),(i+1)%(nrows-1)-1]}\n")
    elif correlation==2:
        # print(data.iloc[:,nrows-1])
        corr=[];o_corr=[];o_index=[];o_pv=[]
        for i in range(0,nrows):
            dcorr,p=scipy.stats.spearmanr(datat.iloc[:, nrows-1], datat.iloc[:, i])
            if not math.isnan(dcorr):
                corr.append(dcorr)
                o_corr.append(dcorr)
                o_index.append(i)
                o_pv.append(p)
        corr.sort(reverse=True)
        with open('corrdata.txt', 'w') as f:
            for i in range(0,len(corr)):
                f.write(datat.columns[o_index[o_corr.index(corr[i])]]);f.write("    p=")
                f.write(f"{o_pv[o_corr.index(corr[i])]}");f.write("    ")
                f.write(f"{corr[i]}\n")


    with open('kendall_result.txt', 'w') as f:
        f.write('kendalltau  r and p\n\n')
        #print(data.iloc[:,nrows-1])
        #print(data.loc[:, 'Primary Site'])
        r, p = scipy.stats.kendalltau(data.iloc[:, nrows - 1], data.iloc[:, nrows])
        f.write('Primary Site and MSI\n');f.write(f"{r}");f.write("    ");f.write(f"{p}\n\n")

def getgraph():
    import cdt
    #glasso = cdt.independence.graph.Glasso()
    #skeleton = glasso.predict(data)
    #new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
    if(algorithm==1):
        model = cdt.causality.graph.CCDr()
    elif algorithm==2:
        model=cdt.causality.graph.PC()
    elif algorithm==3:
        model=cdt.causality.graph.GES()
    output_graph=model.predict(data)
    nx.draw_networkx(output_graph, font_size=5)
    plt.savefig("graph.png",format="PNG")
    plt.show()
    gmlgenerator=nx.generate_gml(output_graph)
    gmllist=[]
    for i in gmlgenerator:
        gmllist.append(i)
    gmlstring="".join(gmllist)
    with open("graph.txt","w") as f:#这里如果用加上writedir的绝对路径的话会报no such file，反正就是在一个工作目录下，没必要
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
        model_dowhy=dowhy.causal_model.CausalModel(data,treatment=tt,outcome=oc,graph=gmlstring) #这里的tt和oc不应该打引号
        identified_estimand = model_dowhy.identify_effect()
        with open("dowhy_identified_estimand.txt","w",encoding='utf-8') as f:
            print(str(identified_estimand))

            f.write(str(identified_estimand))  #write()只能写一个str，没有重载写identified_estimand的,并且由于编码的问题需要指定为utf-8
        print("result of identify effect has been written into identified_estimand.txt")
    except:
        print("invalid input edge!")
        exit()
    if evalution==1:
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
        with open("evalution.txt","w",encoding='utf-8') as f:#还是那个问题：如果是全局路径就莫名会no such file
            f.write("estimate:\n")
            f.write(str(estimate))
            f.write("\n----------------------------\n----------------------------\n")
            f.write("random common cause method result:\n")
            f.write(str(refute1_results))
            f.write("\n----------------------------\n----------------------------\n")
            f.write("placebo treatment refuter method result\n")
            f.write(str(refute2_results))
            f.write("\n----------------------------\n----------------------------\n")
            f.write("data subset refuter method result\n")
            f.write(str(refute3_results))
    elif evalution==2:
        estimate_iv=model_dowhy.estimate_effect(identified_estimand, method_name="iv.instrumental_variable")
        print(estimate_iv)
        refute1_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
            method_name="random_common_cause")
        print(refute1_results)
        refute2_results=model_dowhy.refute_estimate(identified_estimand, estimate_iv,
            method_name="data_subset_refuter")
        print(refute2_results)
        with open("evalution.txt","w",encoding='utf-8') as f:
            f.write("estimate:\n")
            f.write(str(estimate_iv))
            f.write("\n----------------------------\n----------------------------\n")
            f.write("random common cause method result:\n")
            f.write(str(refute1_results))
            f.write("\n----------------------------\n----------------------------\n")
            f.write("data subset refuter method result\n")
            f.write(str(refute2_results))
    print("result of evalution has been written into evalution.txt")
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
