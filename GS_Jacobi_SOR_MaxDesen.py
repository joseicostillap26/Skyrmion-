

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 22:40:58 2020

@author: JOSE
"""


from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

d = np.zeros([])
l1 = np.zeros([])
l2 = np.zeros([])
l3 = np.zeros([])
l4 = np.zeros([])
u = np.zeros([])
m = np.zeros([])
d_in = np.zeros([])

lista_erreur1= []
lista_erreur2= []
lista_erreur3= []
lista_erreur4= []
###### Diagonal #######


def error_relativo(x2,x1):
    
    return abs(x1-x2)**2


def sumalista(lista):
    laSuma = 0
    for i in lista:
        laSuma = laSuma + i
    return laSuma


def Din_b(matriz,vector):
    
   dim = len(matriz)
     
   x = np.zeros([dim])
   
   for i in range (0,dim):
        
       x[i]= vector[i]/matriz[i][i]

   return(x)


def Din_U(matriz):
    
   dim = len(matriz)
     
   x = np.zeros([dim,dim])

   z = np.zeros([dim,dim])


   for i in range (0,dim):
       
       for j in range (0,dim):
           
           if i < j:
               
               x[i][j] = matriz[i][j]

   for i in range (0,dim):
      
      k = 0
      
      while k < dim:
           
           z[i][k] =  x[i][k]/matriz[i][i]
          
           k = k + 1
           
           
   return(z)


def Din_L(matriz):
    
   dim = len(matriz)
     
   y = np.zeros([dim,dim])
   z = np.zeros([dim,dim])


   for i in range (0,dim):
       
       for j in range (0,dim):
           
           if i > j:
               
               y[i][j] = matriz[i][j]

   for i in range (0,dim):
      
      k = 0
      
      while k < dim:
           
           z[i][k] =  y[i][k]/matriz[i][i]
          
           k = k + 1
           
           
   return(z)


def Din_LU(matriz):
    
   dim = len(matriz)
     
   x = np.zeros([dim,dim])
   y = np.zeros([dim,dim])
   z = np.zeros([dim,dim])
  
   for i in range (0,dim):
       
       for j in range (0,dim):
           
           if i > j:
               
               x[i][j] = matriz[i][j]


   for i in range (0,dim):
       
       for j in range (0,dim):
           
           if i < j:
               
               y[i][j] = matriz[i][j]

   for i in range (0,dim):
      
      k = 0
      
      while k < dim:
           
           z[i][k] =  (x[i][k] + y[i][k])/matriz[i][i]
          
           k = k + 1
           
           
   return(z)

def Jacobi(error_relativo,matriz,vector,tol):

   dim = len(matriz)
   lista_v1= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0
   erreur = 1
   while (erreur > tol):

     for i in range (0,dim):
         z = 0
         for j in range(0,dim):

              z = z + Din_LU(matriz)[i][j]*y[j] 

         x[i] =  Din_b(matriz,vector)[i] - z

     for k in range (0,dim):
         cont_error = abs(error_relativo(x[k],y[k]))
         lista_v1[k] = cont_error
     
     erreur = (sumalista(lista_v1))**(1/2)
     lista_erreur1.append(erreur)
     
     for k in range(0,dim):
         y[k] = x[k]

  

   return(x)

def Gauss_Seidel(error_relativo,matriz,vector,tol):

   dim = len(matriz)
   lista_v2= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0
   erreur = 1
  
   while (erreur > tol):

     for i in range (0,dim):
         z1 = 0
         z2 = 0
         for j in range(0,dim):

              z1 = z1 + Din_U(matriz)[i][j]*y[j] 
              z2 = z2 + Din_L(matriz)[i][j]*x[j]
              
              
         x[i] =  Din_b(matriz,vector)[i] - z1 -z2

     for k in range (0,dim):
         cont_error = abs(error_relativo(x[k],y[k]))
         lista_v2[k] = cont_error
         

     erreur = (sumalista(lista_v2))**(1/2)
     lista_erreur2.append(erreur)
     
     for k in range(0,dim):
         y[k] = x[k] 
         x[k] = 0 # reset


  
   return(y)

def SOR(error_relativo,matriz,vector,tol,w):

   dim = len(matriz)
   lista_v3= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0
   erreur = 1
   
   while (erreur > tol):

     for i in range (0,dim):
         z1 = 0
         z2 = 0
         for j in range(0,dim):
              # y = k , x = k+1
              z1 = z1 + ((1-w)*Din_L(matriz)[i][j]+Din_U(matriz)[i][j])*y[j] 
              z2 = z2 + (w*Din_L(matriz)[i][j])*x[j]
              
              
         x[i] =  Din_b(matriz,vector)[i] - z1 -z2

     for k in range (0,dim):
         cont_error = abs(error_relativo(x[k],y[k]))
         lista_v3[k] = cont_error
         

     erreur = (sumalista(lista_v3))**(1/2)
     lista_erreur3.append(erreur)
     
     for k in range(0,dim):
         y[k] = x[k] 
         x[k] = 0 # reset

   return(y)

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



def M_D(error_relativo,matriz,vector,tol):

   dim = len(matriz)
   lista_v1= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0

   erreur = 1

   # y = k, x = k+1
   while (erreur > tol):

     for i in range (0,dim):
         x[i] = y[i] - Alpha_k(matriz,vector,y)*Grad(matriz,vector,y)[i]

     for k in range (0,dim):

         cont_error = abs(error_relativo(x[k],y[k]))

         lista_v1[k] = cont_error

     for k in range(0,dim):
         y[k] = x[k]

     
     erreur = (sumalista(lista_v1))**(1/2)
     lista_erreur4.append(erreur)
   
   return(x)
m1=15
m2=10
m3=8
m4=5
g = 9.81


A =  np.array([[-1,0,0,-m1],[1,-1,0,-m2],[0,1,-1,-m3],[0,0,1,-m4]])
b = np.array([m1*0.8*g*np.cos(np.pi/4)-m1*g*np.sin(np.pi/4),m2*0.2*g*np.cos(np.pi/4)-m2*g*np.sin(np.pi/4),m3*g,m4*g])

    

l1 = Jacobi(error_relativo,A,b,1e-6)
l2 = Gauss_Seidel(error_relativo,A,b,1e-6)
l3 = SOR(error_relativo,A,b,1e-6,1.5)
l4 =  M_D(error_relativo,A,b,1e-6)


dim_e1 = len(lista_erreur1)
dim_e2 = len(lista_erreur2)
dim_e3 = len(lista_erreur3)
dim_e4 = len(lista_erreur4)

print("Solucion Jacobi : ",l1)
print("Solucion Gauss Seidel : ",l2)
print("Solucion SOR: ",l3)
print("Solucion Maximo descenso: ",l4)


plt.title( "Jacobi y Gauss Seidel,SOR")
plt.plot(range(0,dim_e1),lista_erreur1,'ko-',range(0,dim_e2),lista_erreur2,'ro-',range(0,dim_e3),lista_erreur3,'bo-',range(0,dim_e4),lista_erreur4,'yo-')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('Jacobi','Gauss Seidel','SOR','Maximo Descenso',''))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()