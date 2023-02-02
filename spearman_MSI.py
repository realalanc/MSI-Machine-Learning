import pandas as pd
import heapq
import math
import scipy
from scipy import stats

n=1000
dic = {'MSS':0,'MSI-L':1,'MSI-H':2}
data=pd.read_csv("origin_data.csv",header=0,index_col=0,low_memory=False).T
nrows=data.columns.get_loc("MSI")+1 #Excel中MSI所在excel中的行-1(除去了表头）
print(nrows)
data.iloc[:,nrows-1]=data.iloc[:,nrows-1].map(dic)
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
with open('spearman_MSI.txt', 'w') as f:
    for i in range(0,len(corr)):
        f.write(data.columns[o_index[o_corr.index(corr[i])]]);f.write("    p=")
        f.write(f"{o_pv[o_corr.index(corr[i])]}");f.write("    ")
        f.write(f"{corr[i]}\n")
f.close()