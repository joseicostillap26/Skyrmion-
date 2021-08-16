# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 08:14:38 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


from io import open

datos1=[]
datos2=[]

veloc = [] ### = x
presion = [] #### = y

a0_a1 = []

archivo = open('datos.dat','r')

l_a = archivo.readlines()

l_a.pop(0)
l_a.pop(11)
l_a.pop(11)

for i in range (0,len(l_a),1): 
   dato1=l_a[i].find('\t')
   dato2=l_a[i].find('\t',dato1+1)
   veloc.append(float(l_a[i][0:dato1]))
   presion.append(float(l_a[i][dato1+1:dato2]))


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



def Modele(x,a0,a1):
    print('beta = ',np.exp(a0))
    print('alpha =',a1)
    return(np.exp(a0)*x**(a1))


def convert(time,y_mesure):
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]) 
        A[i] = np.log(time[i])

    return(A,B)


def R_L(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert(time,y_mesure)

    for i in range(0,n):

        xy = xy + (A[i])*B[i]
        x = x + (A[i])
        y = y + B[i]
        x_2 = x_2 + (A[i])**2
        y_2 = y_2 + (B[i])**2
         

    a1 = (n*xy-x*y)/(n*x_2-x**2)
    a0 = y/n - a1*x/n
    r = (n*xy-x*y)/((n*x_2-x**2)*(n*y_2-y**2))**(1/2)

    a0_a1.append(a0)
    a0_a1.append(a1)
    print('Coeficiente de Correlación =',r)
    return(a0+t*a1)


R_L(1,presion,veloc)

t = np.linspace(20,200,100)
print(t)
print(len(t))
plt.title( "Newton-Monomio-Lagrange")
plt.plot(presion,veloc,'ko',t,Newton(t,presion,veloc),'r:',t,Modele(t,a0_a1[0],a0_a1[1]),'k:')
plt.show()


