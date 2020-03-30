import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kstest
from scipy.stats import norm
from scipy import stats
#generate normal distribution
np.random.seed(12345678)
x = np.random.normal(10, 0.25, 100)

#load dataset
path = r"C:\Users\miles\Documents\complete.csv"
df = pd.read_csv(path, header = 0, index_col= 0)
#use the normal distribution to make an artificial peak
d1=df.iloc[:, 0]
d1=d1.dropna()
index=d1.index.values.tolist()
z=pd.cut(x,index)
counts=pd.value_counts(z)
counts=counts.sort_index(axis=0)
counts=counts.to_frame()
counts = counts.append(pd.Series(0, index=counts.columns), ignore_index=True)

counts['index']=index
counts=counts.set_index('index')
statistic_store=[]
d1=d1.to_frame()
d1=d1.to_numpy()

#KS test each column in the dataset
counts=counts.to_numpy()
d1=np.ravel(d1)
counts=np.ravel(counts)
df=df.dropna(axis=1, how='all')
columns = list(df)
for i in columns:
    test=(df[i])
    name=(test.name)
    test=test.dropna()
    test=test.to_numpy()
    test=np.ravel(test)
    stat,p = stats.ks_2samp(counts,test)
    print(stat)
    statistic_store.append(stat)



#output the score as a csv
print(statistic_store)
p_df=pd.DataFrame(statistic_store,index=columns)
p_df.to_csv('ks_metric.csv',encoding='utf-8', index=True)