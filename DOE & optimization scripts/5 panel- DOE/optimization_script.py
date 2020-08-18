import subprocess as sp
from math import atan, sin, cos, tan
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import random
# These are the defined design variables
DVs=[0 for x in range(9)]
DataFile = open('PostData.txt','w')
DataFile.write('IRTT	ORTT	T1		T2		T3		G		HNG		phiL	phiR	UX1		UY1		UX2		UY2		UX3		UY3		UX4		UY4		UX5		UY5		UX6		UY6		UX7		UY7		UX8		UY8		UX9		UY9		UX10		UY10		MaxVM\n')
DataFile.close()
DataFile = open('PostDataInitial.txt','w')
DataFile.write('IRTT	ORTT	T1		T2		T3		G		HNG		phiL	phiR	X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM\n')
DataFile.close()
DataFile = open('PostDataFinal.txt','w')
DataFile.write('IRTT	ORTT	T1		T2		T3		G		HNG		phiL	phiR	X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM		Error\n')
DataFile.close()
Total = 1.0;
Tact1=[0.1*T for T in range(14,21,2)]
Tact2=[0.1*T for T in range(1,6,2)]
Thetain = [0.15,0.18,0.21]
x1=[0 for x in range(12)]
z1=[0 for x in range(12)]
# Var number (nV dimension) setpoint
# You can change the dimension's size
nV = 9

# Sample number setpoint
# You can change the Sample number
nS = 600

# Initialisation:
k=1
# Creation of a list dictionnary
DOE_Inputs = {}
DOE_Inputs[k] = []

# Loop elements (part1)
for i in range(1,(nV+1)):
	DOE_Inputs1 = []

	for j in range(1,(nS+1)):
		a = ((float(j)-1)/nS)
		b = ((float(j))/nS)
		listesample = random.uniform(a,b)
		DOE_Inputs1.append(listesample)

	# Select a random number nP times between each Sample and for each Var (part2)
	for k in range(1,nS+1):
		listechoice = random.choice(DOE_Inputs1)
		DOE_Inputs.setdefault(k, []).append(listechoice)
		DOE_Inputs1.remove(listechoice)

for key in DOE_Inputs:
	T1=(DOE_Inputs[key][0])*1.2
	T1=round(T1, 2)
	T2=(DOE_Inputs[key][1])*1.2
	T2=round(T2, 2)
	T3=(DOE_Inputs[key][2])*0.4
	T3=round(T3, 2)
	ORTT=(DOE_Inputs[key][3])*0.2+0.1
	ORTT=round(ORTT, 1)
	t=(DOE_Inputs[key][4])*0.8*ORTT+0.1*ORTT
	t=round(t, 2)
	G=(DOE_Inputs[key][5])*0.3+0.2
	G=round(G, 1)
	HNG=(DOE_Inputs[key][6])*0.3+0.2
	HNG=round(HNG, 1)
	phiL=(DOE_Inputs[key][7])*0.2+0.1
	phiL=round(phiL, 2)	
	phiR=(DOE_Inputs[key][8])*0.2+0.1
	phiR=round(phiR, 2)	
	DVs[0]=T1
	DVs[1]=T2
	DVs[2]=T3
	DVs[3]=ORTT/10
	DVs[4]=round(((ORTT-t)/10),3)
	DVs[5]=G/10
	DVs[6]=HNG/10
	DVs[7]=phiL
	DVs[8]=phiR
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
	#read the files in the output file
	fileobject = open('temp.txt','rb')
	loc = []
	for line in fileobject:
		loc.append(float(line))
	fileobject.close()
	x=loc[0:12]
	print(x)
	z=loc[12:24]
	print(z)
	if ((x==x1) & (z==z1)):
		print('Job Failed')
	else:
		f = interpolate.interp1d(x, z)
		xnew = np.arange(-0.5,0.5, 0.01)
		znew = f(xnew)   # use interpolation function returned by `interp1d`
		curve = +0.5*(z[0]+z[11])+0.16*xnew*xnew-0.04
		error = 0
		for ind in range(100):
			error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
		DataFile = open('PostDataFinal.txt','a')
		DataFile.write('%10f\n' % error)
		DataFile.close()
		fileobject = open('output.txt','rb')
		outputs = []
		for line in fileobject:
			outputs.append(float(line))
		fileobject.close()
	x1=x
	z1=z
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