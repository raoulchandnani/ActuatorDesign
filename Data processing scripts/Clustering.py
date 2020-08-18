# k-means clustering
from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.cluster import KMeans
from matplotlib import pyplot
from sklearn import preprocessing
from math import sin,cos,atan2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.ndimage import  rotate

with open("visualize_deformation.plt", "r") as a:
    next(a)
    next(a)
    next(a)
    Data = [[float(num) for num in line.split(' ')] for line in a]
Data=np.array(Data)
DataF=pd.DataFrame(Data,columns=['X', 'Y', 'Z', 'dX','dY','dZ','dS'])
DataF=DataF.drop_duplicates()
Data=np.array(DataF)
DataF=DataF.set_index(['X'])


## Slice morphing data ##
count=0
Datalow=[[0.000000 for x in range(7)] for y in range(len(Data))]
Datalow=np.array(Datalow)
for i in range(len(Data)):
    if Data[i,6]!=0:
        if abs(Data[i,1])<0.15:
            Datalow[count]=Data[i]
            count=count+1
Datalow=Datalow[~np.all(Datalow == 0, axis=1)]

# define dataset
#Datalow, _ = make_classification(n_samples=1000, n_features=2, n_informative=2, n_redundant=0, n_clusters_per_class=1, random_state=4)
# define the model
model = KMeans(n_clusters=3)
# fit the model
Datalownorm=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
Datalownorm=np.array(Datalownorm)
for i in range(7):
    if i!=3:
        max=np.max(Datalow[:,i])
        min=np.min(Datalow[:,i])
        Datalownorm[:,i]=(Datalow[:,i]-min)/(max-min)
Datalownorm[:,0]=Datalownorm[:,0]
Datanew=[[0.000000 for x in range(2)] for y in range(len(Datalow))]
Datanew=np.array(Datanew)
for i in range(len(Datalownorm)):
    Datanew[i,0]=Datalownorm[i,1]*0.8
    Datanew[i,1]=Datalownorm[i,5]*0.2
wt_kmeansclus=model.fit(Datanew)
centers = wt_kmeansclus.cluster_centers_
# assign a cluster to each example
yhat = model.predict(Datanew)
# retrieve unique clusters
clusters = unique(yhat)
cluster_order=[0 for i in range(len(clusters))]
cluster_order=np.array(cluster_order)
# create scatter plot for samples from each cluster
ax = plt.axes(projection='3d')

for cluster in clusters:
	# get row indexes for samples with this cluster
	row_ix = where(yhat == cluster)
	cluster_order[cluster]=np.min(row_ix)
clusters=clusters[np.argsort(cluster_order[:])]
handles=[]
index=0
for cluster in clusters:
	# get row indexes for samples with this cluster
	row_ix = where(yhat == cluster)
	# create scatter of these samples
	#pyplot.scatter(Datalow[row_ix, 0], Datalow[row_ix, 1]
	ax.scatter(Datalow[row_ix,0],Datalow[row_ix,1],Datalow[row_ix,2],alpha=0.8,label='Max Displacement: %i mm'%(np.max(Datalow[row_ix,6])*1000))
	index+=1
ax.set_xlabel('Keelwise (m)')
ax.set_ylabel('Spanwise (m)')
ax.set_zlabel('Vertical (m)')
font = {'family' : 'normal',
        'size'   : 160}
ax.scatter(centers[:, 0], centers[:, 1], c='black', s=500, alpha=0.5);
show the plot
ax.legend(prop={'size':20})

