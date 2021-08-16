# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:00:42 2020

@author: JOSE
"""


# This import registers the 3D projection, but is otherwise unused.
 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X1 = [0.3,0.6,0.9,0.3,0.6,0.9,0.3,0.6,0.9]
Y1 = [0.001,0.001,0.001,0.01,0.01,0.01,0.05,0.05,0.05]

X =np.arange(0, 1, 0.001)
Y =  np.arange(0, 0.1, 0.001)

X, Y = np.meshgrid(X, Y)

Z = (36.38*X**(2.62))*(Y**(0.53))
Z1 = [0.04,0.24,0.69,0.13,0.82,2.38,0.31,1.95,5.66]
# Plot the surface.
surf = ax.plot_surface(X, Y, Z,color = "skyblue")
ax.scatter3D(X1, Y1, Z1, color = "black");
# Customize the z axis.
# Add a color bar which maps values to colors.
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')
plt.show()