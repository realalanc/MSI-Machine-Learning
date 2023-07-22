import pandas as pd
data=pd.read_csv("STAD_DESeq2_scored_msisensor_selected.csv")
datat=data.T
print(datat.head())
datat_pd=pd.DataFrame(datat)
datat_pd.to_csv("data.csv")