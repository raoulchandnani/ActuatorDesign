import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from math import atan2
from numpy import cov
from scipy import stats

Panels=5                                #W 3,4 or 5 Panel cases
DOE='Act'                               #'Full' or 'Act' 
#with open('%iPanel_%sDOE_FinalCoord.csv' % (Panels,DOE), 'r') as f:
with open('PostDataFinal.csv', 'r') as f: 
    next(f)
    Data = [[float(num) for num in line.split(',')] for line in f]
#with open('%iPanel_%sDOE_InitialCoord.csv' % (Panels,DOE), 'r') as f:
with open('PostDataInitial.csv', 'r') as f: 
    next(f)
    Initial = [[float(num) for num in line.split(',')] for line in f]
Total=1.000
Case=0
data_plot=len(Data)
if Panels==3:
    n=8
    start=7
elif Panels==4:
    n=10
    start=8
elif Panels==5:
    n=12
    start=9
Data=np.array(Data)
Initial=np.array(Initial)
sorted_Data=Data[np.argsort(Data[:, -1])]
sorted_Initial=Initial[np.argsort(Initial[:, -1])]

# Plot Curve fit data#
A=[0 for num in range(n)]
B=[0 for num in range(n)]
A[0]=-0.5*Total
A[n-1]=0.5*Total
B[0]=sorted_Initial[Case][(start+1)]
B[n-1]=sorted_Initial[Case][(start+2*n-5)]
A[1:(n-1)]=sorted_Initial[Case][start:(start+2*n-5):2]
B[1:(n-1)]=sorted_Initial[Case][(start+1):(start+2*n-4):2]

x=[0 for ind in range(n)]
z=[0 for ind in range(n)]
x[0]=-0.5
x[n-1]=0.5
z[0]=sorted_Initial[Case][(start+1)]
z[n-1]=sorted_Initial[Case][(start+2*n-5)]
x[1:(n-1)]=sorted_Data[Case][start:(start+2*n-5):2]
z[1:(n-1)]=sorted_Data[Case][(start+1):(start+2*n-4):2]

f = interpolate.interp1d(x, z)
I=interpolate.interp1d(A, B)
xnew = np.arange(-0.5,0.5, 0.01)
znew = f(xnew)
Inew=I(xnew)
curve = +0.5*(z[0]+z[(n-1)])+0.16*xnew*xnew-0.04
error=0
for ind in range(100):
    error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
    xnew[ind]=xnew[ind]
    Inew[ind]=Inew[ind]
    znew[ind]=znew[ind]
    curve[ind]=curve[ind]
for ind in range(n):
    x[ind]=x[ind]
    z[ind]=z[ind]
plt.figure(0)
plt.plot(xnew, Inew,'-',label='Initial configuration')
plt.plot(xnew, znew, '-',label='Deformed configuration')
plt.plot(xnew,curve,'-',label='Desired Geometry')
plt.title('Actuated structure and desired geometry- %i Panel' % Panels,fontsize=40)
plt.axis('equal')
plt.ylim(-0.15,0.25)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Horizontal location (%)',fontsize=35)
plt.ylabel('Vertical location (%)',fontsize=35)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
plt.rc('xtick', labelsize=30) 
plt.rc('ytick', labelsize=30) 
plt.rc('font', **font)
plt.legend()
plt.show()

# Plot scatter data#
sorted_Data=sorted_Data[0:100]
Change_P=[[0.00 for num in range(5)] for num in range(len(sorted_Data))]
Change_H=[[0.00 for num in range(6)] for num in range(len(sorted_Data))]
Change_P=np.array(Change_P)
Change_H=np.array(Change_H)
for i in range(len(sorted_Data)):
    Change_P[i,0]=atan2((sorted_Initial[i,12]-sorted_Initial[i,10]),(sorted_Initial[i,11]-sorted_Initial[i,9])) - atan2((sorted_Data[i,12]-sorted_Data[i,10]),(sorted_Data[i,11]-sorted_Data[i,9]))
    Change_P[i,1]=atan2((sorted_Initial[i,16]-sorted_Initial[i,14]),(sorted_Initial[i,15]-sorted_Initial[i,13])) - atan2((sorted_Data[i,16]-sorted_Data[i,14]),(sorted_Data[i,15]-sorted_Data[i,13]))
    Change_P[i,2]=atan2((sorted_Initial[i,20]-sorted_Initial[i,18]),(sorted_Initial[i,19]-sorted_Initial[i,17])) - atan2((sorted_Data[i,20]-sorted_Data[i,18]),(sorted_Data[i,19]-sorted_Data[i,17]))
    Change_P[i,3]=atan2((sorted_Initial[i,24]-sorted_Initial[i,22]),(sorted_Initial[i,23]-sorted_Initial[i,21])) - atan2((sorted_Data[i,24]-sorted_Data[i,22]),(sorted_Data[i,23]-sorted_Data[i,21]))
    Change_P[i,4]=atan2((sorted_Initial[i,28]-sorted_Initial[i,26]),(sorted_Initial[i,27]-sorted_Initial[i,25])) - atan2((sorted_Data[i,28]-sorted_Data[i,26]),(sorted_Data[i,27]-sorted_Data[i,25]))

    Change_H[i,0]=atan2((sorted_Initial[i,12]-sorted_Initial[i,10]),(sorted_Initial[i,11]-sorted_Initial[i,9])) - atan2((sorted_Data[i,12]-sorted_Data[i,10]),(sorted_Data[i,11]-sorted_Data[i,9]))
    Change_H[i,1]=Change_P[i,0]-Change_P[i,1]
    Change_H[i,2]=Change_P[i,1]-Change_P[i,2]
    Change_H[i,3]=Change_P[i,2]-Change_P[i,3]
    Change_H[i,4]=Change_P[i,3]-Change_P[i,4]
    Change_H[i,5]=atan2((sorted_Initial[i,28]-sorted_Initial[i,26]),(sorted_Initial[i,27]-sorted_Initial[i,25])) - atan2((sorted_Data[i,28]-sorted_Data[i,26]),(sorted_Data[i,27]-sorted_Data[i,25]))

hinge_num=2
corr1, _ = stats.spearmanr(sorted_Data[:,0], Change_P[:,hinge_num])
corr2, _ = stats.spearmanr(sorted_Data[:,1], Change_P[:,hinge_num])
corr3, _ = stats.spearmanr(sorted_Data[:,2], Change_P[:,hinge_num])
corr4, _ = stats.spearmanr(sorted_Data[:,3], Change_P[:,hinge_num])
corr5, _ = stats.spearmanr(sorted_Data[:,4], Change_P[:,hinge_num])
corr6, _ = stats.spearmanr(sorted_Data[:,5], Change_P[:,hinge_num])
corr7, _ = stats.spearmanr(sorted_Data[:,6], Change_P[:,hinge_num])
corr8, _ = stats.spearmanr(sorted_Data[:,7], Change_P[:,hinge_num])
corr9, _ = stats.spearmanr(sorted_Data[:,8], Change_P[:,hinge_num])
hinge_num=0
corr1, _ = stats.spearmanr(sorted_Data[:,0], Change_H[:,hinge_num])
corr2, _ = stats.spearmanr(sorted_Data[:,1], Change_H[:,hinge_num])
corr3, _ = stats.spearmanr(sorted_Data[:,2], Change_H[:,hinge_num])
corr4, _ = stats.spearmanr(sorted_Data[:,3], Change_H[:,hinge_num])
corr5, _ = stats.spearmanr(sorted_Data[:,4], Change_H[:,hinge_num])
corr6, _ = stats.spearmanr(sorted_Data[:,5], Change_H[:,hinge_num])
corr7, _ = stats.spearmanr(sorted_Data[:,6], Change_H[:,hinge_num])
corr8, _ = stats.spearmanr(sorted_Data[:,7], Change_H[:,hinge_num])
corr9, _ = stats.spearmanr(sorted_Data[:,8], Change_H[:,hinge_num])

#max_disp=[0.00 for num in range(len(Data))]
#for i in range(len(Data)):
#    max_disp[i]=max(sorted_Initial[i,10:29:2]-sorted_Data[i,10:29:2])
#corr1, _ = stats.spearmanr(sorted_Data[:,0], max_disp[:])

plt.figure(1)
plt.scatter(sorted_Data[:,2]*0.1,Change_H[:,1],label='Shear strain at Tube 1')
plt.scatter(sorted_Data[:,3]*0.1,Change_H[:,1],label='Shear strain at Tube 2 ')
plt.scatter(sorted_Data[:,4]*0.1,Change_H[:,1],label='Shear strain at Tube 3 ')
#plt.title('Shear Strain vs Change in angle',fontsize=40)
plt.xlabel('Shear Strain',fontsize=35)
plt.ylabel('Change in Plate angle (rad)',fontsize=35)
plt.rc('xtick', labelsize=30) 
plt.rc('ytick', labelsize=30) 
plt.rc('font', **font)
plt.legend()
plt.show()