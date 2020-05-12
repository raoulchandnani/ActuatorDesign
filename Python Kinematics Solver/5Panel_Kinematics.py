import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from math import atan, sin, cos, tan,pi,asin
## Initial geometry Variables ##
Total=3000        # Total length of mechanism
phiL= 30            #Initial angle of first plate (degrees)          
phiR= 30            # Initial angle of fifth plate (degrees)
## Inputs ##
Case = 2           # Case 1 controls plates 1,2 & 3; Case 2 controls plates 1,2 & 5
theta_1= 0        # Final angle of first plate (degrees)
theta_2= -30         # Final angle of second plate (degrees)
theta_3= 0         # Final angle of third plate (Case 1) OR fifth plate (Case 2) (degrees)
## Geometry calculations ##
phiLr=phiL*pi/180   # COnverting angles to radian
phiRr=phiR*pi/180
theta=asin((sin(phiRr)-sin(phiLr))/3)
H=(Total)/(3*cos(theta)+cos(phiLr)+cos(phiRr));
## Kinematics ##
from pyslvs import (
    parse_vpoints,
    t_config,
    data_collecting,
    expr_solving,
)
if __name__ == '__main__':
    vpoints = parse_vpoints(
        "M["
        "J[R, color[green], P[%d, %d], L[ground, link_1]], " 
        "J[R, color[green], P[%d, %d], L[link_1, link_2]], "
        "J[R, color[green], P[%d, %d], L[link_2, link_3]], "
        "J[R, color[green], P[%d, %d], L[link_3, link_5]], "
        "J[R, color[green], P[%d,%d], L[link_4, link_5]], "
        "J[R, color[green], P[%d, %d], L[ground, link_4]], "
        "]"% (-0.5*Total,0,-0.5*Total+H*cos(phiLr),H*sin(phiLr),-0.5*Total+H*cos(phiLr)+H*cos(theta),H*sin(phiLr)+H*sin(theta),-0.5*Total+H*cos(phiLr)+2*H*cos(theta),H*sin(phiLr)+2*H*sin(theta),-0.5*Total+H*cos(phiLr)+3*H*cos(theta),H*sin(phiLr)+3*H*sin(theta),0.5*Total,0)) 
    if Case==1:
        exprs = t_config(vpoints, ((0, 1), (1, 2), (2, 3)))
    else:
        theta_3=180-theta_3;
        exprs = t_config(vpoints, ((0, 1), (1, 2), (5, 4)))
    mapping = {n: f'P{n}' for n in range(len(vpoints))}
    data_dict, dof = data_collecting(exprs, mapping, vpoints)
    pos = expr_solving(exprs, mapping, vpoints, [theta_1,theta_2,theta_3])
position=[[0 for i in range(2)] for j in range(6)]
for i in range(6):
    for j in range(2):
        position[i][j]=float(pos[i][j])
print(pos)
x=[position[0][0],position[1][0],position[2][0],position[3][0],position[4][0],position[5][0]];
z=[position[0][1],position[1][1],position[2][1],position[3][1],position[4][1],position[5][1]];
print(x)
print(z)
## Plot initial and final configuration ##
A=[0 for num in range(6)]
B=[0 for num in range(6)]
A=[-0.5*Total,-0.5*Total+H*cos(phiLr),-0.5*Total+H*cos(phiLr)+H*cos(theta),-0.5*Total+H*cos(phiLr)+2*H*cos(theta),-0.5*Total+H*cos(phiLr)+3*H*cos(theta),0.5*Total]
B=[0,H*sin(phiLr),H*sin(phiLr)+H*sin(theta),H*sin(phiLr)+2*H*sin(theta),H*sin(phiLr)+3*H*sin(theta),0]
f = interpolate.interp1d(x, z)
I=interpolate.interp1d(A, B)
xnew = np.arange(-0.5*Total,0.5*Total+1, (Total/100))
znew = f(xnew)
Inew=I(xnew)
for ind in range(101):
    xnew[ind]=xnew[ind]/30
    Inew[ind]=Inew[ind]/30
    znew[ind]=znew[ind]/30
plt.plot(xnew, Inew,'-',label='Initial configuration')
plt.plot(xnew, znew, '-',label='Actuated configuration')
plt.title('Actuated structure and Initial geometry',fontsize=30)
plt.axis('equal')
plt.ylim(-10,20)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Horizontal location (%)',fontsize=20)
plt.ylabel('Vertical location (%)',fontsize=20)
font = {
        'size'   : 20}
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('font', **font)
plt.legend()
plt.show()
## Printing variables and output locations ##

DataFile = open('PostData.txt','w')
DataFile.write('Initial_left    Initial_right   Theta1   Theta2   Theta3    x1   y1   x2   y2   x3   y3   x4   y4   x5   y5   x6   y6\n')
DataFile.write('%10f,%10f,%10f,%10f,%10f,' % (phiL,phiR,theta_1,theta_2,theta_3, ))
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (position[0][0],position[0][1],position[1][0],position[1][1],position[2][0],position[2][1],position[3][0],position[3][1],position[4][0],position[4][1],position[5][0],position[5][1]))
DataFile.close()

