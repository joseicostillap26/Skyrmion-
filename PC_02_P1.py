

# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 21:24:40 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


########### Sistema de ecuaciones #######################

#####  -(k1+k2)x1 + (k2)x2 = 0           #####
#####    (k2)x1 - (k2+k3)x2 + (k3)x3= 0  #####
#####   (k3)x2 - (k3+k4)x3 + (k4)x4 = 0  #####
#####   (k4)x3 - (k4)x4 = - F            #####


lista_erreur= [] ## lista donde guardo los errres
inte= []       # Lista donde guardo las interaciones
a0_a1 = []     # Lista donde guardo a0 y a1 de la regresion lineal


def Modele(x,a0,a1):        #### Modelo Err(k) = bx^a
    print('beta = ',np.exp(a0))
    print('alpha =',a1)
    return(np.exp(a0)*x**(a1))


def error(x2,x1):
    
    return (x1-x2)**2

def sumalista(lista):
    laSuma = 0
    for i in lista:
        laSuma = laSuma + i
    return laSuma  


########## Aplico Maximo Descenso con las subfunciones Grad y Alpha_k #### 
def Grad(matriz,vector,x):
    
    dim = len(matriz)
    grad = np.zeros([dim])
    
    for i in range (0,dim):
         z = 0
         for j in range(0,dim):
             z = z + matriz[i][j]*x[j]  

         grad[i] =  z - vector[i]
    
    return grad

def Alpha_k(matriz,vector,x):
   
    dim = len(matriz)
    num = 0
    den = 0
    
    for i in range (0,dim):
         z = 0
         z1 = 0
         for j in range(0,dim):
             z = z + Grad(matriz,vector,x)[j]*matriz[j][i]*Grad(matriz,vector,x)[i]
             
         den = den + z
         z1 =  z1 + Grad(matriz,vector,x)[i]*Grad(matriz,vector,x)[i]
         num = num + z1

    return num/den


##### Maximo descenso ####
    
def M_D(error,matriz,vector,tol):####

   dim = len(matriz)
   lista_v1= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0
   N = 0
   erreur = 1

   while (erreur > tol):

     for i in range (0,dim):
         x[i] = y[i] - Alpha_k(matriz,vector,y)*Grad(matriz,vector,y)[i]

     for k in range (0,dim):

         cont_error = abs(error(x[k],y[k]))####3

         lista_v1[k] = cont_error

     for k in range(0,dim):
         y[k] = x[k]
     N  = N +1
     
     erreur = (sumalista(lista_v1))**(1/2)
     lista_erreur.append(erreur)
     inte.append(N)
     

   return(x) ### Devuelve las interaciones X

############### Regrasion Lineal ############

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
    
    return(a0+t*a1) ### Retorna la ecuacion lineal ####


######## Matriz y vector correspondiente al sistema de ecuaciones hallado ###########
    #### k1 = 150, k2 = 50 , k3 = 75, k4 = 225 y F =2000 #########

A =  np.array([[-200,50,0,0],[50,-125,75,0],[0,75,-300,225],[0,0,225,-225]])
b = np.array([0,0,0,-2000])

print("Solucion del sistema de ecuaciones = ",M_D(error,A,b,1e-6))
 
R_L(1,inte,lista_erreur)


t = np.arange(1,400,0.1)
plt.ylim(-0.5,10)
plt.title( "Ajuste lIneal de Err(k) = bx^a")
plt.rc('legend', fontsize='medium')
plt.plot(inte,lista_erreur,'b.',t,Modele(t,a0_a1[0],a0_a1[1]),'k',)
plt.gca().legend(('Error |Xi- Xi+1| vs interacciones','Modelo Err(k)= bx^a',''))
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()

