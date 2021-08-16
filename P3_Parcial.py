# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 20:19:39 2020

@author: JOSE
"""

#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math

warnings.filterwarnings("ignore")

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
   
   
   print('T1, =',x_lu[0])
   print('T1, =',x_lu[0])
   
   return(x_lu)


m1=15
m2=10
m3=8
m4=5
g = 9.81


A =  np.array([[-1,0,0,-m1],[1,-1,0,-m2],[0,1,-1,-m3],[0,0,1,-m4]])
b = np.array([m1*0.8*g*np.cos(np.pi/4)-m1*g*np.sin(np.pi/4),m2*0.2*g*np.cos(np.pi/4)-m2*g*np.sin(np.pi/4),m3*g,m4*g])

LU(A,b,4)