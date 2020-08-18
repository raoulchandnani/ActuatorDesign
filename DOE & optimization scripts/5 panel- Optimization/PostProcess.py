import subprocess as sp
from math import atan, sin, cos, tan,pi
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import random
from pyswarm import pso

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
znew3=f3(xnew)

plt.plot(xnew,-znew,label='Final configuration')
plt.plot(xnew,-znew2,label='Desired configuration')
plt.plot(xnew,-znew3,label='Initial configuration')
plt.title('Keelwise location 13.86 m')
plt.xlabel('Spanwise location (m)')
plt.ylabel('Vertical location (m)')
# plt.gca().set_aspect('equal', adjustable='box')
plt.ylim(-0.06,0.02)
plt.axis('equal') 
plt.legend()
plt.show()
