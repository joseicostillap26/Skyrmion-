


# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 05:59:39 2020

@author: JOSE
"""

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

def Modele1(x):
    
    return((0.4097*x**(1/2)+1.9921)/(x**(1/2)))**2


def Modele2(x):
    
    return(np.exp(2.2681)*x*np.exp(-2.4733*x))

def convert1(time,y_mesure):
  
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        A[i] = 1/time[i]**(1/2) 
        B[i] = y_mesure[i]**(1/2)
   
    return(A,B)

def convert2(time,y_mesure):
  
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]/time[i]) 
        A[i] = time[i]
   
    return(A,B)

def R_L_1(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert1(time,y_mesure)
    for i in range(0,n):

        xy = xy + (A[i])*B[i]
        x = x + (A[i])
        y = y + B[i]
        x_2 = x_2 + (A[i])**2
        y_2 = y_2 + (B[i])**2
         
   
    a1 = (n*xy-x*y)/(n*x_2-x**2)
    a0 = y/n - a1*x/n
    r = (n*xy-x*y)/((n*x_2-x**2)*(n*y_2-y**2))**(1/2)

    return(a0+t*a1)


def R_L_2(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert2(time,y_mesure)

    for i in range(0,n):

        xy = xy + (A[i])*B[i]
        x = x + (A[i])
        y = y + B[i]
        x_2 = x_2 + (A[i])**2
        y_2 = y_2 + (B[i])**2
         
   
    a1 = (n*xy-x*y)/(n*x_2-x**2)
    a0 = y/n - a1*x/n
    r = (n*xy-x*y)/((n*x_2-x**2)*(n*y_2-y**2))**(1/2)
    print(a0,a1)
   
    return(a0+t*a1)



def R_P(t,A,B):

    dim = len(A)
    n = 3
    M = np.zeros([n,n])
    N = np.zeros([n])
    x_g = np.zeros([n])
    factor  = 0
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    x_3 = 0
    x_4 = 0
    y_3 = 0
    x_2y = 0
    a0 = 0
    a1= 0
    a2 = 0

    for i in range(0,dim):

        xy = xy + (A[i])*B[i]
        x = x + A[i]
        y = y + B[i]
        x_2 = x_2 + (A[i])**2
        y_2 = y_2 + (B[i])**2
        x_3 = x_3 + (A[i])**3
        y_3 = y_3 + (B[i])**3
        x_4 = x_4 + (A[i])**4
        x_2y = x_2y + (A[i]**2)*B[i]

   
    M[0][0] = dim
    M[0][1] = x
    M[1][0] = x
    M[0][2] = x_2
    M[1][1] = x_2
    M[2][0] = x_2
    M[2][1] = x_3
    M[1][2] = x_3
    M[2][2] = x_4

    N[0] = y
    N[1] = xy
    N[2] = x_2y
     
    for k in range (0,n):
        for r in range (k+1,n):

            factor = (M[r,k]/M[k,k])
            N[r] = N[r] - (factor*N[k])
           
            for c in range(0,n):
                M[r,c] = M[r,c]-(factor*M[k,c])
    
    x_g[n-1] = N[n-1]/M[n-1][n-1]


    for r in range(n-2,-1,-1):
        suma = 0
        for c in range(0,n):
            suma = suma + M[r,c]*x_g[c]
        x_g[r] = (N[r] - suma)/M[r,r]
    
    a0 = x_g[0]
    a1 = x_g[1]
    a2 = x_g[2]
    
    return(a0 + a1*t + a2*t**2)



Mesure2 = np.array([0.75,1.25,1.45,1.25,0.85,0.55,0.35,0.28,0.18])
temps2 = np.array([0.1,0.2,0.4,0.6,0.9,1.3,1.5,1.7,1.8])


t = np.arange(0,2,0.1)
plt.ylim(0,2)
plt.title( "Newton-Monomio-Lagrange")
plt.plot(t,R_P(t,temps2,Mesure2),'k:',temps2,Mesure2,'bo',t,Modele2(t),'r:')
plt.show()