import pandas as pd
import heapq
import math
import scipy
from scipy import stats

n=1000
dic = {'MSS':0,'MSI-L':1,'MSI-H':2}
dic2 = {float('nan'):0,'Stage I':1,'Stage IA':2,'Stage IB':3,'Stage II':4,'Stage IIA':5,'Stage IIB':6,
        'Stage III':7,'Stage IIIA':8,'Stage IIIB':9,'Stage IIIC':10,'Stage IV':11} #空缺值我按无病状处理的，不知道是不是
data=pd.read_csv("origin_data.csv",header=0,index_col=0,low_memory=False).T
nrows=data.columns.get_loc("AJCC Stages")+1 #Excel中ajcc所在excel中的行-1(除去了表头）
print(nrows)
data.iloc[:,nrows-1]=data.iloc[:,nrows-1].map(dic2)
data.iloc[:,nrows-2]=data.iloc[:,nrows-2].map(dic)
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
with open('spearman_ajcc.txt', 'w') as f:
    for i in range(0,len(corr)):
        f.write(data.columns[o_index[o_corr.index(corr[i])]]);f.write("    p=")
        f.write(f"{o_pv[o_corr.index(corr[i])]}");f.write("    ")
        f.write(f"{corr[i]}\n")
f.close()