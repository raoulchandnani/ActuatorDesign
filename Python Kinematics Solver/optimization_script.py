import subprocess as sp
from math import atan, sin, cos, tan
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
# These are the defined design variables
DVs=[0 for x in range(6)]
DataFile = open('PostData.txt','w')
#DataFile.write('theta_l		theta_r		phiL		phiR		G		HNG		X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM		Error\n')
DataFile.close()
Total = 3000;
Thetact = [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40]
Thetain = [0.1,0.2,0.3,0.4]
for actuateleft in range(8):
	theta_l=Thetact[actuateleft];
	for actuateright in range(1):
		theta_r =theta_l	#Thetact[actuateright];
		for initialleft in range(4):
			phiL = Thetain[initialleft];
			for initialright in range(1):
				phiR = phiL #Thetain[initialright];
				for hingegap in range(5,65,5):
					G = hingegap;
					H = (Total-G*(3+cos(phiL)+cos(phiR)))/(3+cos(phiL)+cos(phiR));
					for hingeheight in range(15,65,5):
						HNG = hingeheight;
						DVs[0]=theta_l
						DVs[1]=theta_r
						DVs[2]=phiL
						DVs[3]=phiR
						DVs[4]=G
						DVs[5]=HNG
						print 'inputs in optimization code: ', DVs
						# This creates the input txt file for the Abaqus script
						fileobject = open('input.txt','wb')
						for DV in DVs:
							fileobject.write('%.4f\n' % DV)
						fileobject.close()
						# Command line to run
						command_file = r'abaqus cae nogui=abaqus_script.py'
						# This is the command that runs the Abaqus script
						ps = sp.Popen(command_file, shell= True)
						# Wait until the Abaqus script is done
						ps.wait()
						# read the files in the output file
						fileobject = open('temp.txt','rb')
						loc = []
						for line in fileobject:
							loc.append(float(line))
						fileobject.close()
						x=loc[0:12]
						print(x)
						z=loc[12:24]
						print(z)
						f = interpolate.interp1d(x, z)
						xnew = np.arange(-1500,1500, 10)
						znew = f(xnew)   # use interpolation function returned by `interp1d`
						curve = +0.5*(z[0]+z[11])+5.33*xnew*xnew/100000 - 50
						error = 0
						for ind in range(300):
							error=error+(curve[ind]-znew[ind])*(curve[ind]-znew[ind])
						DataFile = open('PostData.txt','a')
						DataFile.write('%10f\n' % error)
						DataFile.close()
						fileobject = open('output.txt','rb')
						outputs = []
						for line in fileobject:
							outputs.append(float(line))
						fileobject.close()
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