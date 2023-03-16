import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

data=pd.read_csv("mlh1.csv",index_col=0,header=0)
gene1=[]
gene2=[]
gene3=[]
for i,columns in data.items():
    #columns[0]=float(columns[0])
    t=np.log2(float(columns[0])+1)
    if t=='nan':
        print(columns[0])
    if columns["MSI"]=="MSS":
        gene1.append(t)
    elif columns["MSI"]=="MSI-L":
        gene2.append(t)
    else: 
        gene3.append(t)
plt.grid(True)
plt.boxplot([gene1,gene2,gene3])
plt.show() 
ar1=np.array(gene1)
ar2=np.array(gene2)
ar3=np.array(gene3)
print(ar1)
print('\n\n\n\n')
print(ar2)
print('\n\n\n\n')
print(ar3)
t1,pval1=stats.ttest_ind(ar1,ar2)
t2,pval2=stats.ttest_ind(ar1,ar3)
t3,pval3=stats.ttest_ind(ar2,ar3)
print(pval1)
print(pval2)
print(pval3)