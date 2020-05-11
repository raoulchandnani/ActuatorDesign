import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
# A= [-2.429598,-1.37092,0.449002,-58.807819,0.445219,-61.051968,-2.433388,-118.5761260,-119.9489360,-119.948936,2.433388,-118.576126,-0.445219,-61.051968,-0.449002,-58.807819,2.429598,-1.37092,-119.948936,0.755812]
# y=A[1:20:2]
# x=[-1500,-1170,-830,-500,-170,170,500,830,1170,1500]
# f = interpolate.interp1d(x, y)
# xnew = np.arange(-1500,1500, 10)
# ynew = f(xnew)   # use interpolation function returned by `interp1d`
# znew = 5.33*xnew*xnew/100000 - 120
# plt.plot(x, y, 'o', xnew, ynew, '-',xnew,znew,'-')
# plt.show()
# error = 0
# for ind in range(450):
    # error=error+(znew[ind]-ynew[ind])*(znew[ind]-ynew[ind])

with open('PostData.txt', 'r') as f:
    Data = [[float(num) for num in line.split(',')] for line in f]
Data=np.array(Data)
for ind in range(len(Data)):
    Data[i][2]=phiL
    Data[i][3]=phiR
    Data[i][4]=G
    Data[i][5]=HNG
    Total=3000
    phi02=atan2((0.5*G),HNG)
    phi1=phi02-phiL
    phi2=phi02-phiR
    LHNG2=sqrt((0.5*G)*(0.5*G)+(HNG*HNG))
    H = (Total-G*(3+cos(phiL)+cos(phiR)))/(3+cos(phiL)+cos(phiR));
    x=[0 for ind in range(11)]
    z=[0 for ind in range(11)]
    x[1:10]=sorted_Data[0][6:25:2]
    z[1:10]=sorted_Data[0][7:26:2]
    x[0]=-1500
    x[11]=1500
    z[0]=(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))
    z[11]=(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))  
    f = interpolate.interp1d(x, z)
    xnew = np.arange(-1500,1500, 10)
    znew = f(xnew)   # use interpolation function returned by `interp1d`
    curve = 0.00002*(xnew*xnew-1500*1500) +z[0]
    error = 0
    for ind in range(300):
        error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
    Data[i][27]=error
sorted_Data=Data[np.argsort(Data[:, 27])]
sorted_Data[0][2]=phiL
sorted_Data[0][3]=phiR
sorted_Data[0][4]=G
sorted_Data[0][5]=HNG
Total=3000
phi02=atan2((0.5*G),HNG)
phi1=phi02-phiL
phi2=phi02-phiR
LHNG2=sqrt((0.5*G)*(0.5*G)+(HNG*HNG))
H = (Total-G*(3+cos(phiL)+cos(phiR)))/(3+cos(phiL)+cos(phiR));
A=[0 for num in range(12)]
B=[0 for num in range(12)]
A=[-1500,-(1.5*H+1.5*G+LHNG2*sin(phi1)+(H)*cos(phiL)),(1.5*H+1.5*G+LHNG2*sin(phi1)),-(1.5*H+G),-(0.5*H+G),-(0.5*H),(0.5*H),(0.5*H+G),(1.5*H+G),(1.5*H+1.5*G+LHNG2*sin(phi2)),(1.5*H+1.5*G+LHNG2*sin(phi2)+(H)*cos(phiR)),1500]
B=[(HNG-LHNG2*cos(phi1)-(H)*sin(phiL)),-(HNG-LHNG2*cos(phi1)-(H)*sin(phiL)),(HNG-LHNG2*cos(phi1)),0,0,0,0,0,0,(HNG-LHNG2*cos(phi2)),(HNG-LHNG2*cos(phi2)-(H)*sin(phiR)),(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))]
x=[0 for ind in range(11)]
z=[0 for ind in range(11)]
x[1:10]=sorted_Data[0][6:25:2]
z[1:10]=sorted_Data[0][7:26:2]
x[0]=-1500
x[11]=1500
z[0]=(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))
z[11]=(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))

f = interpolate.interp1d(x, z)
I=interpolate.interp1d(A, B)
xnew = np.arange(-1500,1500, 10)
znew = f(xnew)
Inew=I(xnew)
curve = 0.00002*(xnew*xnew-1500*1500) + z[0]
error=0
for ind in range(300):
	error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
plt.plot(x, z, 'o', xnew, Inew, '-',xnew, znew, '-',xnew,curve,'-')
plt.show()