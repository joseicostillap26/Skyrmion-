

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
from io import open
warnings.filterwarnings("ignore")

datos1=[] ## Listas de apoyo
datos2=[]  ## '''

veloc = []  #### Lista donde guardo las velocidades
presion = []  #### Lista donde guardo las presiones


a0_a1 = [] # Lista donde guardo a0 y a1 de la regresion lineal

archivo = open('datos.dat','r') ## Abro data.dat y lo guardo en archivo##

l_a = archivo.readlines() ## Leo data.dat todo junto como una lista

l_a.pop(0)  ## Elimino los valores de la lista l_a donde no hay numeros
l_a.pop(11) ## '''
l_a.pop(11) ##  '''

for i in range (0,len(l_a),1): 
   dato1=l_a[i].find('\t')            ## Separo los valores de la lista l_a
   dato2=l_a[i].find('\t',dato1+1)    ## En columnas  con respecto al separador \t
   veloc.append(float(l_a[i][0:dato1])) ## Luego agrego las columnas a veloc y presion
   presion.append(float(l_a[i][dato1+1:dato2]))


#### Mi polinomio interpolante es el polinomio de Newton #######
   
def Newton(t,matriz,vector):

    pts = len(t)
    dim = len(matriz)
    v = np.zeros([dim]*2)
    x_f = np.zeros([dim]) ### Lista donde guardo los coeficientes de newton###    
    z = 0
    lj = [] ### Lista donde guardo las funciones t(t-t1),(t-t1)(t-t2) ... para cada t ##
    lz = [] ### Lista donde guardo el valor de la funcion fn(t) para cada t ####

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
   

    ## Resuelvo el sistema Ax = b con Fordward ####
          
    x_f[dim-dim] = vector[dim-dim]/v[dim-dim][dim-dim]
    
    for r in range(1,dim,1):
         suma = 0
         for c in range(0,dim):
             suma = suma + v[r,c]*x_f[c]
         x_f[r] = (vector[r] - suma)/v[r,r]
    
    print("Coeficientes del polinomio interpolante de Newton = ",x_f)
   
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
    print('alpha = ',np.exp(a0))
    print('beta =',a1)
    return(np.exp(a0)*x**(a1))


def convert(time,y_mesure):
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]) 
        A[i] = np.log(time[i])

    return(A,B)

######### Regresion Lineal ##############
    
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

t = np.linspace(20,200,100) ### 100 puntos equidistantes entre 20 y 200 ###

plt.title( "Ajuste lineal de data.dat")
plt.rc('legend', fontsize='medium')
plt.plot(presion,veloc,'bo',t,Newton(t,presion,veloc),'r',t,Modele(t,a0_a1[0],a0_a1[1]),'k')
plt.gca().legend(('velocidad vs Caida de presion','Polinomio de Newton','Modelo v = a(P)^b'))
plt.ylabel('Velocidad (m/s)')
plt.xlabel('Caida de Presion (mm Hg)')
plt.show()


