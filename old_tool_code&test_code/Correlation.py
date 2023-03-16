import pandas as pd
import numpy as np
import heapq

#皮尔逊相关系数 可能要运行个几十分钟  ，dcorr是19965*19965的相关矩阵，是7个多G的数据

data=pd.read_csv("origin_data.csv",header=0,index_col=0,low_memory=False).T
nrows=data.columns.get_loc("MSI")+1 #Excel中MSI所在excel中的行-1(除去了表头） 19966
n=5000
data=pd.read_csv("origin_data.csv",header=0,index_col=0,nrows=nrows-1).T
print(data)
dcorr=data.corr()
dcorr.to_csv(r'corr_result.csv')
# dcorr=pd.read_csv("corr_result.csv",header=0,index_col=0)
a=np.triu(dcorr,1)
a=a.flatten()
topn=heapq.nlargest(n, range(len(a)), a.take)
with open('corrdata.txt', 'w') as f:
    for i in topn:
        # print((dcorr.index)[(i+1)//(nrows-1)],(dcorr.index)[(i+1)%19965-1],dcorr.iloc[(i+1)//19965,(i+1)%19965-1])
        f.write((dcorr.index)[(i+1)//(nrows-1)]);f.write("    ")
        f.write((dcorr.index)[(i+1)%(nrows-1)-1]);f.write("    ")
        f.write(f"{dcorr.iloc[(i+1)//(nrows-1),(i+1)%(nrows-1)-1]}\n")
f.close()