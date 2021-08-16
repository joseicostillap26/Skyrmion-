

#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


lista_b = [[1,2,1],[0,-4,1],[0,0,-2]]
vector1 = [5,2,4]
matriz_b = np.array(lista_b)
x1 = np.zeros([3])

print(" Matriz Backward :\n ",matriz_b)
print("Con vector :\n ",vector1)

for r in range(2,-1,-1):
    suma = 0
    for c in range(0,3):
        suma = suma + matriz_b[r,c]*x1[c]
    x1[r] = (vector1[r] - suma)/matriz_b[r,r]

print("Se obtiene las soluciones :\n ",x1)
print("\n ")

lista_f = [[2,0,0],[1,4,0],[4,3,3]]
vector2 = [4,2,5]
matriz_f = np.array(lista_f)
x2 = np.zeros([3])

print(" Matriz Forward :\n ",matriz_f)
print("Con vector :\n ",vector2)

x2[0] = vector2[0]/matriz_f[0][0] #un arreglo

for r in range(1,3,1):
    suma = 0
    for c in range(0,3):
        suma = suma + matriz_f[r,c]*x2[c]
    x2[r] = (vector2[r] - suma)/matriz_f[r,r]

print("Se obtiene las soluciones :\n ",x2)
print("\n ")


m = 3

matriz = np.array([[1,-2,4],[1.0,1.0,1.0],[1.0,0.0,0.0]])
vector = np.array([-27.0,0.0,-1.0])
x_g = np.zeros([m])



print(" Matriz  :\n ",matriz)
print("vector :\n ",vector)
print("\n ")

for k in range (0,m):

    for r in range (k+1,m):

        factor = (matriz[r,k]/matriz[k,k])
        vector[r] = vector[r] - (factor*vector[k])
        for c in range(0,m):
            
            matriz[r,c] = matriz[r,c]-(factor*matriz[k,c])
            

print(" Matriz Backward :\n ",matriz)
print("Con vector :\n ",vector)
print("\n ")

x_g[m-1] = vector[m-1]/matriz[m-1][m-1]


for r in range(m-2,-1,-1):
    suma = 0
    for c in range(0,m):
        suma = suma + matriz[r,c]*x_g[c]
    x_g[r] = (vector[r] - suma)/matriz[r,r]

print("Se obtiene las soluciones _Gaussiana:\n ",x_g)
print("\n ")



 ################ Solucion LU ################

m1=15
m2=10
m3=8
m4=5
g = 9.81


matriz =  np.array([[-1,0,0,-m1],[1,-1,0,-m2],[0,1,-1,-m3],[0,0,1,-m4]])
b = np.array([m1*0.8*g*np.cos(np.pi/4)-m1*g*np.sin(np.pi/4),m2*0.2*g*np.cos(np.pi/4)-m2*g*np.sin(np.pi/4),m3*g,m4*g])

u = np.zeros([m,m])
l = np.zeros([m,m])
d = np.zeros([m])
x_lu = np.zeros([m])

for r in range (0,m):   
    for c in range(0,m):    
        u[r,c] = matriz[r,c]

print ("Matriz:\n ",matriz)
print ("con vector:\n ",b)
print("\n ")

 ################ Descomposicion LU ################
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

 ################ Descomposicion LU ################
                 
print ("L:\n ",l)
print ("U:\n ",u)
print("\n ")

d[m-3] = b[m-3]/l[m-3][m-3]

for r in range(m-2,3,1):
    suma = 0
    for c in range(0,m):
        suma = suma + l[r,c]*d[c]
    d[r] = (b[r] - suma)/l[r,r]

print("Solucion d :\n ",d)
print("\n ")


for p in range(m-1,-1,-1):
    s = 0
    for c in range(0,m):
        s = s + u[p,c]*x_lu[c]
    x_lu[p] = (d[p] - s)/u[p,p]

print("Solucion LUssiana :\n ",x_lu)
