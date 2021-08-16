# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:37:57 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


def fun(x):
    
    return np.exp(-(x**2))*np.sin(50*x)



#### Mi polinomio interpolante es el polinomio de Lagrange #######
    
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

#### Puntos de Chebyshev ########
    
def punt_cheby(t,a,b,n):
    lj= []
    x = np.zeros([n])
    y = np.zeros([n])
    valor = 0
    z = 0

    for i in range(1,n):
        valor = a + (b-a)*(1-np.cos(i*np.pi/(n-1)))/2
        x[i] = valor
        y[i] = fun(valor)

    for i in range(0,n):
        p = 1 
        for k in range(0,n):
               
              if k != i:
                   p = p*((t - x[k]) / (x[i]-x[k]))        
        lj.append(p)

    for i in range(0,n):
          z = z + lj[i]*y[i]
          
    return(z)


#### Puntos equidistantes ########

def punt_equi(t,a,b,n):
    
    lj= []
    x = np.zeros([n])
    y = np.zeros([n])
    valor = 0
    z = 0
    for i in range(0,n):
        valor = a+(b-a)*i/(n-1)
        x[i] = valor
        y[i] = fun(valor)

    for i in range(0,n):
        p = 1 
        for k in range(0,n):
               
              if k != i:
                   p = p*((t - x[k]) / (x[i]-x[k]))        
        lj.append(p)

    for i in range(0,n):
          z = z + lj[i]*y[i]
          
    return(z)

##### Errores en la interpolacion Lagrange-Chevishev #####
def E_pe(t,a,b,m):

        y = 0
        y = abs(fun(t)- punt_equi(t,a,b,m))
        return(y)

##### Errores en la interpolacion Lagrange-Puntos equidistantes #####
def E_ch(t,a,b,m):
        y = 0
        y = abs(fun(t)- punt_cheby(t,a,b,m))
        return(y)


t = np.arange(0,1,0.01)
plt.ylim(-5,5)
plt.title( "Grafica de Polinomios interpolantes")
plt.plot(t,fun(t),'k',t,punt_equi(t,0,1,20),'r',t,punt_cheby(t,0,1,20),'b')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('f(x) = e^x^2(sin(50x))','P_Equidistantes','Chebyshev'))
plt.ylabel('Eje y')
plt.xlabel('Eje x')
plt.show()


t = np.arange(0,1+0.01,0.01)
plt.yscale("log")
plt.title( "Errores en Chebyshev - P_equidistantes ")
plt.plot(t,E_pe(t,0,1,100),'k:',t,E_ch(t,0,1,100),'b:') ## 100 puntos equidistantes
plt.rc('legend', fontsize='medium')
plt.gca().legend(('Interpolacion P_equidistantes','Interpolacion de Chebyshev',''))
plt.ylabel('Eje y')
plt.xlabel('Eje x')
plt.show()


