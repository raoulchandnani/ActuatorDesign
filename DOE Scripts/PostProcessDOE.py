import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
Panels=5                                # 3,4 or 5 Panel cases
DOE='Full'                               #'Full' or 'Act' 
with open('%iPanel_%sDOE_FinalCoord.csv' % (Panels,DOE), 'r') as f:
    next(f)
    Data = [[float(num) for num in line.split(',')] for line in f]
with open('%iPanel_%sDOE_InitialCoord.csv' % (Panels,DOE), 'r') as f:
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
plt.figure(1)
plt.scatter((sorted_Data[:,1]-sorted_Data[:,0]),sorted_Data[:,-2])
plt.title('Output vs Input',fontsize=40)
plt.xlabel('Input Parameter',fontsize=35)
plt.ylabel('Output Parameter',fontsize=35)
plt.rc('xtick', labelsize=30) 
plt.rc('ytick', labelsize=30) 
plt.rc('font', **font)
plt.legend()
plt.show()