import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from math import atan2
from numpy import cov
from scipy import stats
from mpl_toolkits import mplot3d
import pandas as pd
from scipy.interpolate import griddata
Panels=5                                #W 3,4 or 5 Panel cases
DOE='Act'                               #'Full' or 'Act' 
#with open('%iPanel_%sDOE_FinalCoord.csv' % (Panels,DOE), 'r') as f:
with open('SurfaceData.txt', 'r') as f: 
    next(f)
    Data = [[float(num) for num in line.split(',')] for line in f]
#with open('%iPanel_%sDOE_InitialCoord.csv' % (Panels,DOE), 'r') as f:

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])



Data=np.array(Data)

ax = plt.axes(projection='3d')
#ax.set_aspect('equal')
ax.scatter3D(Data[:,0],Data[:,1],-Data[:,2],c='blue',alpha=0.2,label='Initial configuration')
ax.scatter3D(Data[:,3],Data[:,4],-Data[:,5],c='red',alpha=1,label='Deformed configuration')
#ax.contour(Data[:,0],Data[:,1],-Data[:,2])
set_axes_equal(ax)


