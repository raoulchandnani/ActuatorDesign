import subprocess as sp
from math import atan, sin, cos, tan,pi
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import random
from pyswarm import pso
# These are the defined design variables
DVs=[0 for x in range(8)]
x1=[0 for x in range(62)]
z1=[0 for x in range(62)]

# Define the objective (to be minimize)
def weight(x):
    T1, G,HNG, phiL2,phiL,H1net,H2net,min_thk = x
    DVs[0]=T1
    DVs[1]=G
    DVs[2]=HNG
    DVs[3]=phiL2
    DVs[4]=phiL
    DVs[5]=H1net
    DVs[6]=H2net
    DVs[7]=min_thk
    print 'inputs in optimization code: ', DVs
    # This creates the input txt file for the Abaqus script
    fileobject = open('input.txt','wb')
    for DV in DVs:
		fileobject.write('%.4f\n' % DV)
    fileobject.close()
    # Command line to run
    try:
		command_file = r'abaqus cae nogui=5Panel_Thermal_Curved.py'
		# This is the command that runs the Abaqus script
		ps = sp.Popen(command_file, shell= True)
		# Wait until the Abaqus script is done
		ps.wait()
		fileobject = open('temp.txt','rb')
		#read the files in the output file
		fileobject = open('temp.txt','rb')
		loc = []
		for line in fileobject:
			loc.append(float(line))
		fileobject.close()
		Volume=loc[-1]
    except:
        Volume=1E08
    print(Volume)
    return Volume

def error_check(x):	
	try:
		#read the files in the output file
		fileobject = open('temp.txt','rb')
		loc = []
		for line in fileobject:
			loc.append(float(line))
		fileobject.close()
		initial=loc[0:62]
		#	print(x)
		final=loc[62:124]
		#	print(z)
		Desired=loc[124:-1]
		f = interpolate.interp1d(final[:31], final[31:], fill_value='extrapolate')
		f2 = interpolate.interp1d(Desired[:31], Desired[31:] ,fill_value='extrapolate')
		f3= interpolate.interp1d(initial[:31], initial[31:] ,fill_value='extrapolate')
		xnew = np.arange(-0.15,0.16, 0.01)
		znew = f(xnew)   # use interpolation function returned by `interp1d`
		znew2 = f2(xnew)   # use interpolation function returned by `interp1d`
		znew2=znew2+(initial[31]-Desired[31])
		error = 0
		for ind in range(31):
			error=error+(znew2[ind]-znew[ind])*(znew2[ind]-znew[ind])
		error=np.sqrt(error/31)
	except:
		error=1E08
	print(error)
	return (1.0E-03-error)
	
def length_check(x):
    T1, G,HNG, phiL2,phiL,H1net,H2net,min_thk = x
    return (0.28-2*H1net*cos(phiL2*pi/180)-2*H2net*cos(phiL*pi/180))
	
# args = (T1, G,HNG, phiL2,phiL,H1net,H2net)
# Define the lower and upper bounds for H, d, t, respectively
lb = [0.1, 0.003, -0.01,5,0,0.02,0.02,0.000]
ub = [1.5, 0.01, -0.003,20,10,0.1,0.1,0.010]

xopt, fopt = pso(weight, lb, ub,ieqcons=[length_check,error_check],swarmsize=160,maxiter=5,minstep=1e-4,minfunc=1e-8,debug=False)
print(xopt,fopt)
Sol=np.hstack((xopt,fopt))
np.savetxt('Solution.txt', Sol) 