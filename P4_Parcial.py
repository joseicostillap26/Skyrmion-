# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 10:41:07 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt
from matplotlib import cm
warnings.filterwarnings("ignore")


a0_a1_a2 = [] 


def Modele(x1,x2,a0,a1,a2):

    return(a0*x1**(a1)*x2**(a2))


def LU(matriz,b,m):
   u = np.zeros([m,m])
   l = np.zeros([m,m])
   d = np.zeros([m])
   x_lu = np.zeros([m])

   for r in range (0,m):   
    for c in range(0,m):    
        u[r,c] = matriz[r,c]



   for k in range (0,m):
     for r in range (0,m):
         if (k == r ):
             l[k,r] = 1
         if ( k < r ):
             factor = (matriz[r,k]/matriz[k,k])
             l[r,k] = factor
             for c in range(0,m):
                 matriz[r,c] = matriz[r,c] - (factor*matriz[k,c])
                 u[r,c] = matriz[r,c]

 
   d[m-m] = b[m-m]/l[m-m][m-m]

   for r in range(m**2-7*m +13,m,1):# Queria que el primer valor de range 
                                     # sea siempre 1 para m = 4 o m = 3, entonces
                                     # encontre la funcion m**3-7*m+13
    suma_1 = 0
    for c in range(0,m):
        suma_1 = suma_1 + l[r,c]*d[c]
    d[r] = (b[r] - suma_1)/l[r,r]


   for p in range(m-1,-1,-1):
    s_1 = 0
    for c in range(0,m):
        s_1 = s_1 + u[p,c]*x_lu[c]
    x_lu[p] = (d[p] - s_1)/u[p,p]
   return(x_lu)


def convert_2var(y,D,S):
    n = len(y)
    A = np.zeros([n])
    B = np.zeros([n])
    C = np.zeros([n])
    for i in range(0,n):
        A[i] = np.log(y[i]) 
        B[i] = np.log(D[i])
        C[i] = np.log(S[i])
    return(A,B,C)


#flux = flujo
#E = Experimento
#D = Diamentro
#S = Inclinacion
#y = flujo
    
def R_L_2var(t,v,fluj,dia,inc):

    n = len(fluj)
    flux = np.zeros([n]) 
    D = np.zeros([n])
    S = np.zeros([n])
    a012 =  np.zeros([3])

    xd=0
    xs=0
    y=0
    yxd=0
    yxs  = 0
    xd_2 = 0
    xs_2=0
    xsxd=0

    flux,D,S = convert_2var(fluj,dia,inc)

    for i in range(0,n):

        xd = xd + D[i]
        xs = xs + S[i]
        y = y + flux[i]
        yxd = yxd + flux[i]*D[i]
        yxs = yxs + flux[i]*S[i]
        xd_2 = xd_2 + (D[i])**2
        xs_2 = xs_2 + (S[i])**2
        xsxd = xsxd + S[i]*D[i]

    A =  np.array([[n,xd,xs],[xd,xd_2,xsxd],[xs,xsxd,xs_2]])
    b = np.array([y,yxd,yxs])
    
    a012= LU(A,b,3)
    
    a012[0] = np.exp(a012[0])
    
    a0_a1_a2.append(a012[0])
    a0_a1_a2.append(a012[1])
    a0_a1_a2.append(a012[2])

    return(a012[0] + a012[1]*t + a012[2]*v) 

Dia = np.array([0.3,0.6,0.9,0.3,0.6,0.9,0.3,0.6,0.9])
Inc = np.array([0.001,0.001,0.001,0.01,0.01,0.01,0.05,0.05,0.05])
Flu = np.array([0.04,0.24,0.69,0.13,0.82,2.38,0.31,1.95,5.66])


R_L_2var(1,1,Flu,Dia,Inc)

print(a0_a1_a2)
fig = plt.figure()
ax = fig.gca(projection='3d')

X =np.arange(0, 1, 0.001)
Y =  np.arange(0, 0.1, 0.001)
X, Y = np.meshgrid(X, Y)

Z = Modele(X,Y,a0_a1_a2[0],a0_a1_a2[1],a0_a1_a2[2])

surf = ax.plot_surface(X, Y, Z,color = "skyblue")
ax.scatter3D(Dia, Inc,Flu, color = "black");
# Customize the z axis.
# Add a color bar which maps values to colors.
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')
plt.show()