import pandas as pd
genelist={}
star9=[]
star8=[]
star7=[]
star6=[]
star5=[]
star4=[]
with open("../old_result/ttest_deseq2.txt","r") as f:
    with open("../old_result/ttest_deseq22.txt","r") as g:
        while True:
            namet=f.readline()
            if not namet:
                break
            name=namet.strip().split(' ')[0]
            avg=float(namet.strip().split(' ')[1][14:])
            star=0
            p1=float(f.readline().strip().split('  ')[1][6:])
            p2=float(f.readline().strip().split('  ')[1][6:])   
            p3=float(f.readline().strip().split('  ')[1][6:])
            q1=g.readline().strip()
            q2=g.readline().strip().split(' ')
            q3=g.readline().strip()
            #if q3=='n': continue
            if p1=='nan' or p2=='nan' or p3=='nan': continue
            #if p1<=0.05: star+=1
            #if p1<=0.01: star+=1
            if p1<=0.1: star+=1
            #if p2<=0.05: star+=1
            #if p2<=0.01: star+=1
            if p2<=0.1: star+=1
            #if p3<=0.05: star+=1
            #if p3<=0.01: star+=1
            if p3<=0.1: star+=1
            genelist[name]=[p1,p2,p3]
            if avg<=2.0:
                continue
            if star==3:
                star9.append(name)
            elif star==2:
                star8.append(name)
            elif star==1:
                star7.append(name)
            #elif star==3:
            #    star6.append(name)
            #elif star==5:
            #    star5.append(name)

print(len(star9))
#print(star9)
print(len(star8))
#star8.append("CLDN11")
#star8.append("MSIscore")
#print(star8)
print(len(star7))
#print(star7)
#print(len(star6))
#print(star6)
#print(len(star5))

data=pd.read_csv("../DESeq2_scored.csv",index_col=0,header=0)

newdata=(data.loc[star8]).T
newdata.to_csv("../chosendata2.csv",index=True)