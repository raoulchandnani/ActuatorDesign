import subprocess as sp
from math import atan, sin, cos, tan
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import random
# These are the defined design variables
DVs=[0 for x in range(5)]
DataFile = open('PostData.txt','w')
DataFile.write('X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM\n')
DataFile.close()
Total = 1.0;
Tact1=[0.1*T for T in range(14,21,2)]
Tact2=[0.1*T for T in range(1,6,2)]
Thetain = [0.15,0.18,0.21]
# Var number (nV dimension) setpoint
# You can change the dimension's size
nV = 5

# Sample number setpoint
# You can change the Sample number
nS = 350

# Initialisation:
k=1
# Creation of a list dictionnary
x = {}
x[k] = []

# Loop elements (part1)
for i in range(1,(nV+1)):
	x1 = []

	for j in range(1,(nS+1)):
		a = ((float(j)-1)/nS)
		b = ((float(j))/nS)
		listesample = random.uniform(a,b)
		x1.append(listesample)

	# Select a random number nP times between each Sample and for each Var (part2)
	for k in range(1,nS+1):
		listechoice = random.choice(x1)
		x.setdefault(k, []).append(listechoice)
		x1.remove(listechoice)

for key in x:
	T1=(x[key][0])
	T1=round(T1, 2)
	T2=(x[key][1])*0.5*T1+0.5*T1
	T2=round(T2, 2)
	T3=(x[key][2])*0.2
	T3=round(T3, 2)
	ORTT=(x[key][3])*20+10
	ORTT=round(ORTT, 0)
	t=(x[key][4])*0.65*ORTT+0.25*ORTT
	t=round(t, 0)
	DVs[0]=T1
	DVs[1]=T2
	DVs[2]=T3
	DVs[3]=ORTT/1000
	DVs[4]=(ORTT-t)/1000
	print 'inputs in optimization code: ', DVs
	# This creates the input txt file for the Abaqus script
	fileobject = open('input.txt','wb')
	for DV in DVs:
		fileobject.write('%.4f\n' % DV)
	fileobject.close()
	# Command line to run
	command_file = r'abaqus cae nogui=5Panel_Thermal.py'
	# This is the command that runs the Abaqus script
	ps = sp.Popen(command_file, shell= True)
	# Wait until the Abaqus script is done
	ps.wait()
	# read the files in the output file
	# fileobject = open('temp.txt','rb')
	# loc = []
	# for line in fileobject:
		# loc.append(float(line))
	# fileobject.close()
	# x=loc[0:12]
	# print(x)
	# z=loc[12:24]
	# print(z)
	# f = interpolate.interp1d(x, z)
	# xnew = np.arange(-0.5,0.5, 0.01)
	# znew = f(xnew)   # use interpolation function returned by `interp1d`
	# curve = +0.5*(z[0]+z[11])+0.16*xnew*xnew-0.04
	# error = 0
	# for ind in range(100):
		# error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
	# DataFile = open('PostData.txt','a')
	# DataFile.write('%10f\n' % error)
	# DataFile.close()
	# fileobject = open('output.txt','rb')
	# outputs = []
	# for line in fileobject:
		# outputs.append(float(line))
	# fileobject.close()
#M=dlmread('PostData.txt')
# fileobject = open('PostData.txt','rb')
# outputs = []
		# for line in fileobject:
			# outputs.append(float(line))
		# fileobject.close()
# with open('PostData.txt', 'r') as f:
    # Data = [[float(num) for num in line.split(',')] for line in f]
# Data=np.array(Data)
# sorted_Data=Data[np.argsort(Data[:, 27])]

# x[1:10]=loc[6:25:2]
# z[1:10]=loc[7:26:2]
# f = interpolate.interp1d(x, z)
# xnew = np.arange(-1500,1500, 10)
# znew = f(xnew)
# curve = +0.5*(z[0]+z[11])+5.33*xnew*xnew/100000 - 120
# error=0
# for ind in range(300):
	# error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
# plt.plot(x, z, 'o', xnew, znew, '-',xnew,curve,'-')
# plt.show()