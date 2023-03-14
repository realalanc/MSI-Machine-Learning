import pandas as pd
from scipy import stats
import numpy as np

#0-19968

data=pd.read_csv("DESeq2.csv",index_col=0,header=0)
for j in range(0,19963):
    list1=[]
    list2=[]
    list3=[]
    count=0
    sum0=0.0
    div=0
    for i,columns in data.items():
        columns[j]=np.log2(float(columns[j])+1)
        div+=1
        sum0+=float(columns[j])
        #columns[j]=float(columns[j])*10000000/float(columns["Sum"])
        if columns[j]<=0.05:
            count+=1
        if columns["MSI"]=="MSS":
            list1.append(columns[j])
        elif columns["MSI"]=="MSI-L":
            list2.append(columns[j])
        else: 
            list3.append(columns[j])     
    ar1=np.array(list1)
    ar2=np.array(list2)
    ar3=np.array(list3)
    #print(ar1)
    #print(ar2)
    #print(ar3)
    #t1,pval1=stats.ttest_ind(ar1,ar2)
    #t2,pval2=stats.ttest_ind(ar1,ar3)
    #t3,pval3=stats.ttest_ind(ar2,ar3)
    with open("ttest_deseq22.txt","a+") as f:
        f.write(data.index[j])
        f.write("\n")
        t1=np.mean(ar1)
        t2=np.mean(ar2)
        t3=np.mean(ar3)
        f.write("%lf %lf %lf\n"%(t1,t2,t3))
        if t1<=t2 and t2<=t3:
            f.write("y\n")
        elif t1>=t2 and t2>=t3:
            f.write("y\n")
        else:
            f.write("n\n")
    #exit()