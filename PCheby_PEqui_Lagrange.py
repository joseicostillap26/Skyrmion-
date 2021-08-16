
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 07:29:12 2020

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
    print(lj)
    for i in range(0,dim):
          z = z + lj[i]*vector[i]

    return(z)



def fun(x):
    
    return np.tanh(20*np.sin(12*x)) + np.exp(3*x)*np.sin(300*x)/50


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

def E_pe(t,a,b,m):
        print(t)
        y = 0
        y = abs(fun(t)- punt_equi(t,a,b,m))
        return(y)

def E_ch(t,a,b,m):
        y = 0
        y = abs(fun(t)- punt_cheby(t,a,b,m))
        return(y)




x =  np.array([0,1,2,5.5,11,13,16,18])
y = np.array([0.5,3.134,5.3,9.9,10.2,9.35,7.2,6.2])


t = np.arange(0,1,0.01)
plt.ylim(-1.5,1.5)
plt.title( "punt_cheby - t,punt_equi- fun")
plt.plot(t,punt_cheby(t,0,1,100),'k:',t,punt_equi(t,0,1,100),'y:',t,fun(t),'b:')
plt.show()

t1 = np.arange(0,1,0.01)   #### 100
t2 = np.arange(0,1,0.02) #### 200
t3 = np.arange(0,1,0.03) #### 300
t4 = np.arange(0,1,0.04) #### 400
t5 = np.arange(0,1,0.05) #### 500
t6 = np.arange(0,1,0.06) #### 600
t7 = np.arange(0,1,0.07) #### 700
t8 = np.arange(0,1,0.08) #### 800
t9 = np.arange(0,1,0.09) #### 900

t = np.arange(0,1+0.01,0.01)
plt.yscale("log")
plt.title( "Errores cheby - P_equi ")
plt.plot(t,E_pe(t,0,1,100),'k',t,E_ch(t,0,1,100),'b')
plt.show()

t = np.arange(0,20,0.1)
plt.title( "Newton-Monomio-Lagrange")
plt.plot(t,Lagrange(t,x,y),'k:',t,Newton(t,x,y),'b:',t,Monomio(t,x,y),'y')
plt.show()
