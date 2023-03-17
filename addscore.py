import pandas as pd
score=pd.read_csv("sup_table3_caseinfo_csved.csv",header=0,index_col=0)
olddata=pd.read_csv("DESeq2.csv",header=0,index_col=0)
scorelist=[]

print(len(olddata.T))
for i in olddata:
    flag=False
    for j in score.index:
        if(i==j):
            scorelist.append(score.loc[j,"MSIsensor-pro(all)"])
            flag=True
    if not flag:
        print("Y")
        scorelist.append("NA")
print(scorelist)
print(len(scorelist))
scorelistdf=pd.DataFrame(scorelist).T
scorelistdf.to_csv("MSIScore.csv")