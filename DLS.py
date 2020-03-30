import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn as sk
from sklearn.cluster import KMeans
import seaborn as sns

nc=5  #set number of clusters
quickval: int=10   #set quick metric

#load and interpolate data
path = r"C:\Users\miles\Downloads\completeSMA.csv" ##set input data path here
df = pd.read_csv(path, header = 0, index_col= 0)
di=df.interpolate(method = 'index', limit_direction='both')

#calculate proportion of intensity at desired end point
columns=list(di)
proportions=[]
for a in columns:
    totals=di[a].sum()
    sums=di.loc[(quickval-(0.1*quickval)):quickval+(0.1*quickval),a].sum()
    proportion=(sums/totals)*100
    proportions.append(proportion)
#return end point metric for each sample
p_df=pd.DataFrame(proportions,index=columns)
p_df.to_csv('quick_metric.csv',encoding='utf-8', index=True)
#plot line graph of all data
di.plot.line()
plt.xscale('log')
plt.legend(ncol=3)
plt.show()
dt = di.transpose()
dt = dt.dropna()



#initiate clustering parameters
clustering_kmeans = KMeans(n_clusters=nc, precompute_distances="auto", n_jobs=-1)
#apply clustering
x = clustering_kmeans.fit_predict(dt)
#store clustering results in data frame
y=pd.DataFrame(x, index=dt.index,)
y.columns=['clust']
clustdf= pd.concat([dt, y], axis=1, sort=False, join='inner')
#return .csv of cluster allocation by well
clust_out=pd.DataFrame(clustdf.clust)
clust_out.to_csv('clustering.csv', encoding='utf-8', index=True)
#return line graphs for each cluster
for i in range(0,nc):
    new1=clustdf.loc[clustdf.clust == i]
    new1=new1.transpose()
    new1=new1[:-1]
    new1.plot.line()
    plt.xscale('log')
    plt.legend(ncol=2)
    plt.show()

