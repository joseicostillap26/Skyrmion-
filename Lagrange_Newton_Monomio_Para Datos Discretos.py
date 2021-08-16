
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 23:25:04 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


lk= []
l = 0

def Monomio(t,matriz,vector):
    
    dim = len(matriz)
    x_g = np.zeros([dim])
    v = np.zeros([dim]*2)
    lk= []###
    bi = 0 
    cont = 0
    for i in range(0,dim):
        
        for j in range(0,dim):

            v[i][j] = matriz[i]**j
    

    for k in range (0,dim):
        for r in range (k+1,dim):

            factor = (v[r,k]/v[k,k])
            vector[r] = vector[r] - (factor*vector[k])
           
            for c in range(0,dim):
                v[r,c] = v[r,c]-(factor*v[k,c])
    
    x_g[dim-1] = vector[dim-1]/v[dim-1][dim-1]


    for r in range(dim-2,-1,-1):
        suma = 0
        for c in range(0,dim):
            suma = suma + v[r,c]*x_g[c]
        x_g[r] = (vector[r] - suma)/v[r,r]

    for i in range(0,dim):
        bi = t**i
        lk.append(bi)##

    for i in range(0,dim):
        z = 0
        z = z + lk[i]*x_g[i] ###
        cont = cont + z

    return(cont)


def Newton(t,matriz,vector):

    pts = len(t)
    dim = len(matriz)
    v = np.zeros([dim]*2)
    x_f = np.zeros([dim])
    z = 0
    lj = []
    lz = []

    
    for i in range(0,dim):
        v[i][0] = 1
        for j in range(0,dim-1):
          k = 0
          cont = 1
          while k < j+1:
                p = 1
                p = p*(matriz[i]-matriz[k]) 
                cont = cont*p
                k = k +1
          v[i][j+1] = cont
          
    x_f[dim-dim] = vector[dim-dim]/v[dim-dim][dim-dim]
    
    for r in range(1,dim,1):
         suma = 0
         for c in range(0,dim):
             suma = suma + v[r,c]*x_f[c]
         x_f[r] = (vector[r] - suma)/v[r,r]
    
    for m in range(0,pts):
      pi_j= np.zeros([dim])
      pi_j[0]= 1
      for j in range(0,dim-1):
          k = 0
          cont = 1
          while k < j+1:
                p = 1
                p = p*(t[m]-matriz[k]) 
                cont = cont*p
                k = k +1
          pi_j[j+1] = cont
      lj.append(pi_j)

    for i in range(0,pts):
      z = 0
      for j in range(0,dim):
          z = z + lj[i][j]*x_f[j]
      lz.append(z)

    return(lz)


def Lagrange(t, matriz,vector):

    dim = len(matriz)
    lj= []

    z = 0

    for i in range(0,dim):
        p = 1 
        for k in range(0,dim):
               
              if k != i:
                   p = p*((t - matriz[k]) / (matriz[i]-matriz[k]))        
        lj.append(p)

    for i in range(0,dim):
          z = z + lj[i]*vector[i]

    return(z)

x =  np.array([0,1,2,5.5,11,13,16,18])
y = np.array([0.5,3.134,5.3,9.9,10.2,9.35,7.2,6.2])

t = np.arange(0,20,0.1)
plt.title( "Newton-Monomio-Lagrange")
plt.plot(t,Lagrange(t,x,y),'k:',t,Newton(t,x,y),'b:',t,Monomio(t,x,y),'y')
plt.show()