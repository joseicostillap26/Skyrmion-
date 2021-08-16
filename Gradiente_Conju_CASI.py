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


lista_erreur= []


def error_relativo(x2,x1):
    
    return (x1-x2)**2



def sumalista(lista):
    laSuma = 0
    for i in lista:
        laSuma = laSuma + i
    return laSuma  



def Phi(matriz,vector,x):
    
    dim = len(matriz)
    phi = 0
    
    for i in range (0,dim):
         z = 0
         z1 = 0
         for j in range(0,dim):
             z = z + x[j]*matriz[j][i]*x[i]/2  #####Falta entre 2
              #####Falta entre 2
              #####Falta entre 2
              #####Falta entre 2
             
         z1 = x[i]*vector[i]
         phi = phi + z -z1
    
    return phi


def Grad(matriz,vector,x):
    
    dim = len(matriz)
    grad = np.zeros([dim])
    
    for i in range (0,dim):
         z = 0
         for j in range(0,dim):
             z = z + matriz[i][j]*x[j]  

         grad[i] =  z - vector[i]
    
    return grad


def M_D(error_relativo,matriz,vector,tol):

   ############################
    dim = len(matriz)
    x_0 = np.zeros([dim])
    x_1 = np.zeros([dim])
    p_0 = np.zeros([dim])
    p_1 = np.zeros([dim])
    r_0 = np.zeros([dim])
    r_1 = np.zeros([dim])

###########################
    num = 0
    den = 0
    num_1 = 0
    den_1 = 0
    a_0 = 0
    b_0 = 0
    ###########
    N = 0
    lista_v1= np.zeros([dim])

    cont_error = 0

    erreur = 1

    for i in range (0,dim):

        p_0[i] = Grad(matriz,vector,x_0)[i]
        r_0[i] = p_0[i]
        
   # y = k, x = k+1
    while (N<3):
       N = N+1
       ##########ALpha ##########
       for i in range (0,dim):
           z = 0
           z1 = 0
           for j in range(0,dim):
               
               z = z + p_0[j]*matriz[j][i]*p_0[i]
           den = den + z
           z1 =  z1 + r_0[i]*r_0[i]
           num = num + z1 

       a_0 = num/den 
       
       print('a_0 = ',a_0)

        ################## x #################
       
       for i in range (0,dim):
            x_1[i] = x_0[i] + a_0*p_0[i]
      
       print('x_1 =',x_1)
########## r1 ########### 
       for i in range (0,dim):
           z = 0
           for j in range(0,dim):
                z = z + a_0*matriz[j][i]*p_0[i]
           r_1[i] =  r_0[i] + z 
#################################
        
       print('r_1 =',r_1)
         
################ B1 #############
         
       for i in range (0,dim):
            z = 0
            z1 = 0
            z = z + r_1[i]*r_1[i]
            z1 =  z1 + r_0[i]*r_0[i]
            num_1 = num_1 + z
            den_1 = den_1 + z1
      
       b_0 = num_1/den_1

       print('b_0 =',b_0)
    
    ############ P1 ################


       for i in range (0,dim):
            z = 0

            for j in range(0,dim):
                z = z + b_0*p_0[i]

            p_1[i] =  -r_1[i] + z 
 
       print('p_1 =',p_1)


       for k in range (0,dim):

         cont_error = abs(error_relativo(x_1[k],x_0[k]))

         lista_v1[k] = cont_error
         
       
       for k in range(0,dim):
         x_0[k] = x_1[k]
         p_0[k] = p_1[k]
         r_0[k] = r_1[k]
     
       erreur = (sumalista(lista_v1))**(1/2)
       lista_erreur.append(erreur)
    
       num = 0
       den = 0
       num_1 = 0
       den_1 = 0

    
    return(x_1)



A =  np.array([[1,2],[2,1]])
b = np.array([1,1])


d = np.zeros([])

d =  M_D(error_relativo,A,b,1e-6)
