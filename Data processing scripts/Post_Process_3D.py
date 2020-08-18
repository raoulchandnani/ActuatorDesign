import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import sin,cos,pi

## Load aero data ##
with open("Aero_Data_I.txt", "r") as a:
    Datalow= [[float(num) for num in line.split(' ')] for line in a]
Datalow=np.array(Datalow)
for i in range(len(Datalow)):
    Datalow[i,0]=Datalow[i,0]-11.85
    Datalow[i,3]=Datalow[i,3]-11.85
## Re-arrange aero data to match with Abaqus outputs##
temp=Datalow.copy()
for i in range(len(Datalow)):
    Datalow[i,0]=temp[i,1]
    Datalow[i,1]=-temp[i,2]
    Datalow[i,2]=temp[i,0]
    Datalow[i,3]=temp[i,4]
    Datalow[i,4]=-temp[i,5]
    Datalow[i,5]=temp[i,3] 
## Slice aero data depending on start and end locations of Conformal surface##
Datalow=Datalow[1023:2945,:]
initial=np.loadtxt("initial.txt", delimiter=",")
initial=np.array(initial)
final=np.loadtxt("final.txt", delimiter=",")
final=np.array(final)
ax = plt.axes(projection='3d')

## Affine transformations to translate abaqus data to match with Aero data##
Vct=[-1.86930884303288e-16, -1.40464922179676, -0.432385577178925]
for i in range(len(initial)):
    for j in range(3):
        initial[i,j]=initial[i,j]-Vct[j]
        final[i,j]=final[i,j]-Vct[j]

for i in range(len(initial)):
    initial[i,0]=-initial[i,0]
    initial[i,1]=-initial[i,1]    
    final[i,0]=-final[i,0]
    final[i,1]=-final[i,1]  

for i in range(len(initial)):
    temp=initial[i].copy()
    initial[i,1]=temp[1]*cos(-93.9201298450921*pi/180)-temp[2]*sin(-93.9201298450921*pi/180)
    initial[i,2]=temp[1]*sin(-93.9201298450921*pi/180)+temp[2]*cos(-93.9201298450921*pi/180)
    temp=final[i].copy()
    final[i,1]=temp[1]*cos(-93.9201298450921*pi/180)-temp[2]*sin(-93.9201298450921*pi/180)
    final[i,2]=temp[1]*sin(-93.9201298450921*pi/180)+temp[2]*cos(-93.9201298450921*pi/180)

## Sort order of abaqus data points to match with aero data##
initial_final=np.hstack((initial,final))
initial_final=initial_final[np.argsort(initial_final[:, 2])]
start=0
end=31
for xloc in range(62):
    initial_final[start:end,:]=initial_final[start:end,:][np.argsort(initial_final[start:end, 0])]   
    start+=31
    end+=31
error=0
start=0
end=31
error_loc=[]
for i in range(1860):
    error+=np.sqrt((initial_final[i,3]-Datalow[i,3])**2+(initial_final[i,4]-Datalow[i,4])**2+(initial_final[i,5]-Datalow[i,5])**2)
    error_loc.append(np.sqrt((initial_final[i,3]-Datalow[i,3])**2+(initial_final[i,4]-Datalow[i,4])**2+(initial_final[i,5]-Datalow[i,5])**2))
    start+=31
    end+=31
	
## Plot error data##
ax = plt.axes(projection='3d')
fig = plt.figure(figsize = (16, 9)) 
sctt=ax.scatter3D(initial_final[:1860,3],initial_final[:1860,5]+11.85,initial_final[:1860,4],alpha=1.0,c=error_loc,cmap='YlGnBu_r')
ax.set_xlabel('Spanwise (m)', fontweight ='bold')  
ax.set_ylabel('Keelwise (m)', fontweight ='bold')  
ax.set_zlabel('Vertical (m)', fontweight ='bold')

xx, zz = np.meshgrid((-0.15,0.15), (0.5,0.7))
yy1=[[13.23,13.23],[13.23,13.23]]
ax.plot_surface(xx,yy1,zz,alpha=0.2,color='purple')
yy2=[[13.86,13.86],[13.86,13.86]]
ax.plot_surface(xx,yy2,zz,alpha=0.2,color='purple') 
yy3=[[14.46,14.46],[14.46,14.46]]
ax.plot_surface(xx,yy3,zz,alpha=0.2,color='purple')
fig.colorbar(sctt, ax = ax, shrink = 0.5, aspect = 5,label='RMSE(m)') 



