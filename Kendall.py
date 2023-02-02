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
with open('Kendall_result.txt', 'w') as f:
        f.write('kendalltau  r and p\n\n')
        r, p = scipy.stats.kendalltau(data.iloc[:, nrows - 2], data.iloc[:, nrows])
        f.write('Primary Site and MSI\n');f.write(f"{r}");f.write("    ");f.write(f"{p}\n\n")
        r, p = scipy.stats.kendalltau(data.iloc[:, nrows - 1], data.iloc[:, nrows])
        f.write('Primary Site and AJCC Stages\n');f.write(f"{r}");f.write("    ");f.write(f"{p}\n\n")
        r, p = scipy.stats.kendalltau(data.iloc[:, nrows - 2], data.iloc[:, nrows-1])
        f.write('MSI and AJCC Stages\n');f.write(f"{r}");f.write("    ");f.write(f"{p}\n")
f.close()