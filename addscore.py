import pandas as pd
score=pd.read_csv("sup_table3_caseinfo_csved.csv",header=0,index_col=0)
print(score.head())
olddata=pd.read_csv("DESeq2.csv",header=0,index_col=0)
scorelist=[]
for i in olddata:
    for j in score.index:
        if(i==j):
            scorelist.append(score.loc[j,"MSIsensor-pro(all)"])
scorelistdf=pd.DataFrame(scorelist).T
scorelistdf.to_csv("MSIScore.csv")