import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
with open("visualize_deformation.plt", "r") as a:
    next(a)
    next(a)
    next(a)
    Data = [[float(num) for num in line.split(' ')] for line in a]
Data=np.array(Data)
count=0
Datalow=[[0.000000 for x in range(7)] for y in range(len(Data))]
Datalow=np.array(Datalow)
for i in range(len(Data)):
    if Data[i,2]<-0.4:
        Datalow[count]=Data[i]
        count=count+1
Datalow=Datalow[~np.all(Datalow == 0, axis=1)]

Data2=[[0.000000 for x in range(3)] for y in range(len(Datalow))]
Data2=np.array(Data2)
Data2[:,0]=Datalow[:,0]+Datalow[:,3]
Data2[:,1]=Datalow[:,1]+Datalow[:,4]
Data2[:,2]=Datalow[:,2]+Datalow[:,5]

count=0
Data3=[[0.000000 for x in range(7)] for y in range(len(Data))]
Data3=np.array(Data3)
for i in range(len(Data)):
    if Data[i,0]==13.7085:
        if Data[i,2]<-0.52:
            Data3[count]=Data[i]
            count=count+1
Data3=Data3[~np.all(Data3 == 0, axis=1)]

Data4=[[0.000000 for x in range(3)] for y in range(len(Data3))]
Data4=np.array(Data4)
Data4[:,0]=Data3[:,0]+Data3[:,3]
Data4[:,1]=Data3[:,1]+Data3[:,4]
Data4[:,2]=Data3[:,2]+Data3[:,5]

#ax = plt.axes(projection='3d')
##ax.plot3D(Data3[:,0],Data3[:,1],Data3[:,2])
#count1=0
#count2=0
#for x in range(101):
#    num=11.85+0.0315*x
#    num=round(num,4)
#    Data5=[[0.000000 for x in range(7)] for y in range(len(Data))]
#    Data6=[[0.000000 for x in range(3)] for y in range(len(Data2))]    
#    Data5=np.array(Data5)
#    Data6=np.array(Data6)    
#    for i in range(len(Data)):
#        if Data[i,0]==num:
#                Data5[count1]=Data[i]
#                count1=count1+1  
#    for i in range(len(Data2)):
#        if Data2[i,0]==num:
#                Data6[count2]=Data2[i]
#                count2=count2+1  
#    Data5=Data5[~np.all(Data5 == 0, axis=1)]
#    len5=len(Data5)
#    Data6=Data6[~np.all(Data6 == 0, axis=1)]    
#    len6=len(Data6)    
#    ax.plot3D(Data5[0:len5-2,0],Data5[0:len5-2,1],Data5[0:len5-2,2],c='blue',alpha=0.2)
#    ax.plot3D(Data6[1:len6-3,0],Data6[1:len6-3,1],Data6[1:len6-3,2],c='red',alpha=0.2)

ax = plt.axes(projection='3d')
ax.plot3D(Data[:,0],Data[:,1],Data[:,2],c='orange',alpha=0.4)
#ax.plot3D(Datalow[:,0],Datalow[:,1],Datalow[:,2],c='blue',alpha=0.4)
ax.plot3D(Data2[:,0],Data2[:,1],Data2[:,2],c='red',alpha=0.4)

plt.plot(Data3[0:55,1],Data3[0:55,2],c='blue',label='Initial configuration')
plt.plot(Data4[0:55,1],Data4[0:55,2],c='red',label='Deformed configuration')
#plt.axis('equal')
plt.ylim(-0.6,-0.4)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.ylim(-0.6,-0.5)
plt.xlim(-0.3,0.3)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Spanwise location',fontsize=35)
plt.ylabel('Vertical location',fontsize=35)
font = {'family' : 'normal',
        'size'   : 20}
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('font', **font)
plt.yticks([-0.8,-0.6,-0.4,-0.2])
plt.legend()
plt.show()


plt.plot(Datalow[:,0],Datalow[:,2])
plt.plot(Data2[:,0],Data2[:,2])
plt.axis('equal')


Data6=[[0.000000 for x in range(2)] for y in range(99)]
Data6=np.array(Data6)
count=0
for x in range(99):
    num=11.913+0.0315*x
    num=round(num,4)
    Data5=[[0.000000 for x in range(7)] for y in range(len(Data))]
    Data5=np.array(Data5)
    for i in range(len(Data)):
        if Data[i,0]==num:
                Data5[count]=Data[i]
                count=count+1
    lowest=min(Data5[:,2])
    Data6[x,0]=num
    Data6[x,1]=lowest

Data7=[[0.000000 for x in range(2)] for y in range(99)]
Data7=np.array(Data7)
count=0
for x in range(99):
    num=11.913+0.0315*x
    num=round(num,4)
    Data5=[[0.000000 for x in range(3)] for y in range(len(Data2))]
    Data5=np.array(Data5)
    for i in range(len(Data2)):
        if Data2[i,0]==num:
                Data5[count]=Data2[i]
                count=count+1
    lowest=min(Data5[:,2])
    Data7[x,0]=num
    Data7[x,1]=lowest

fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2,sharey=ax1)
ax1.plot(Data3[0:55,1],Data3[0:55,2],c='blue',label='Initial configuration')
ax1.plot(Data4[0:55,1],Data4[0:55,2],c='red',label='Deformed configuration')
ax1.set_xlabel('Spanwise location',fontsize=20)
ax1.set_ylabel('Vertical location',fontsize=20) 
ax1.axis('equal')  
ax2.plot(Data6[:,0],Data6[:,1],c='blue',label='Initial configuration')
ax2.plot(Data7[:,0],Data7[:,1],c='red',label='Deformed configuration')
ax2.set_xlabel('Keelwise location',fontsize=20)
ax2.set_ylabel('Vertical location',fontsize=20)
ax2.legend()
ax1.legend()



#plt.axis('equal')


plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel('Keelwise location',fontsize=35)
plt.ylabel('Vertical location',fontsize=35)
font = {'family' : 'normal',
        'size'   : 30}
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('font', **font)
plt.legend()
plt.show()
plt.plot(Data3[:,1],Data3[:,2],label='Initial configuration')
plt.plot(Data4[:,1],Data4[:,2],label='Deformed configuration')
plt.axis('equal')
plt.ylim(-0.8,0)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Keelwise location',fontsize=35)
plt.ylabel('Vertical location',fontsize=35)
font = {'family' : 'normal',
        'size'   : 30}
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('font', **font)
plt.legend()
plt.show()
    

ax = plt.axes(projection='3d')
ax.plot3D(Data[:,0],Data[:,1],Data[:,2],c='green',alpha=0.2)
ax.plot3D(Data6[:,0],[0 for v in range(99)],Data6[:,1],c='blue',alpha=1,label='Initial configuration')
ax.plot3D(Data7[:,0],[0 for v in range(99)],Data7[:,1],c='red',alpha=1,label='Deformed configuration')

ax.plot3D(Data3[0:55,0],Data3[0:55,1],Data3[0:55,2],c='blue',alpha=1)
ax.plot3D(Data4[0:55,0],Data4[0:55,1],Data4[0:55,2],c='red',alpha=1)
ax.legend()