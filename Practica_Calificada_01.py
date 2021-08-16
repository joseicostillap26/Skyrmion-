# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:54:15 2020

@author: JOSE
"""


#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
from math import sin
warnings.filterwarnings("ignore")

#### FUNCIONES QUE DEFINO PARA LA PREGUNTA 1 Y 2 ####

### GAUSS ##

def Gauss(matriz,vector,m):

   x = np.zeros([m])

   for k in range (0,m):

     for r in range (k+1,m):

        factor = (matriz[r,k]/matriz[k,k])
        vector[r] = vector[r] - (factor*vector[k])
        for c in range(0,m):
            
            matriz[r,c] =matriz[r,c]-(factor*matriz[k,c])
            

   x[m-1] = vector[m-1]/matriz[m-1][m-1]


   for r in range(m-2,-1,-1):
    suma_2 = 0
    for c in range(0,m):
        suma_2 = suma_2 + matriz[r,c]*x[c]
    x[r] = (vector[r] - suma_2)/matriz[r,r]
    
   return(x)




### DESCOMPOSICION LU ##
   
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
   
   return(x_lu)




### DESCOMPOSICION DE CHOLESKY ##
   
def Cholesky(matriz,b,m):

   u = np.zeros([m,m])
   l = np.zeros([m,m])
   d = np.zeros([m])
   x = np.zeros([m])


   for j in range(m):
       for i in range(j,m):
           if i == j:
               sumk = 0
               for k in range(j):
                    sumk += l[i,k]**2
               l[i,j] = np.sqrt(matriz[i,j]-sumk)
           else:
               sumk = 0
               for k in range(j):
                    sumk += l[i,k]*l[j,k]
               l[i,j] = (matriz[i,j]-sumk)/l[j,j]

   u = np.transpose(l)
   d[m-m] = b[m-m]/l[m-m][m-m]


   for r in range(m**2-7*m +13,m,1): # Queria que el primer valor de range 
                                     # sea siempre 1 para m = 4 o m = 3, entonces
                                     # encontre sa funcion m**3-7*m+13
       suma = 0
       for c in range(0,m):
            suma = suma + l[r,c]*d[c]
       d[r] = (b[r] - suma)/l[r,r]


   for p in range(m-1,-1,-1):
     s = 0
     for c in range(0,m):
        s = s + u[p,c]*x[c]
     x[p] = (d[p] - s)/u[p,p]
   
   return(x)




#### FUNCIONES QUE DEFINO PARA LA PREGUNTA 3 Y 4 ####

#### FUNCIONES  PARA LA PREGUNTA 3 ####
   
def funcion_Fx(x):
    return v_x*x/((x**2+a**2)**(3/2)) # Funcion F_x


def funcion_Fy(x):
    return v_y/((x**2+a**2)**(3/2)) # Funcion F_y


def funcion_F(x):                   # Funcion F = (F_x**2 + F_y**2)**(1/2) 
    return ((v_x*x/((x**2+a**2)**(3/2)))**2 + (v_y/((x**2+a**2)**(3/2)))**2)**(1/2)



def funcion_P3_Fy(x): 
    return x**2 - (v_y1)**(2/3) + a**2 # Funcion para hallar raices F_y = -2.17

def f_d_P3_Fy(x): # Derivada de F_y
    return 2*x 

def funcion_P3(x): # Funcion para hallar raices cuando F_x = 1.56

    return x**2 - (v*x)**(2/3) + a**2

def f_d_P3(x): # Derivada de  F_x
    
    return 2*x - (2*v**(2/3)*x**(-1/3))/3


#### FUNCIONES  PARA LA PREGUNTA 4 ####
    
def funcion_P4(x): # Funcion que me da la altura (h = x)
     return 5*(x**3) - 15*(x**2) + 16 


def f_d_P4(x): # Derivada de la funcion F(h)
    return 15*(x**2) - 30*x



def error_relativo(x2,x1): #FUNCION PARA EL ERROR RELATIVO
    return abs(x1-x2)/x2



#### METODO NEWTON RHAPSON ####
    
def newton_rhapson(error_relativo,funcion,f_d,x1,tol,firma):# Firma es solo un
                                                            # numero que me sirve para rellenar 
                                                            # diferentes listas
    error = 10
    x2 = 0
    cont_error = 0
    n = 0
    m = 0
    p = 0
    r = 0
    t = 0
    while error > tol:
        
        x2 = x1 - funcion(x1)/f_d(x1)
        cont_error = error_relativo(x2,x1)
        x1 = x2
        error = abs(funcion(x2)) 

        if firma == 0:
            lista_v1_P3.append(cont_error)
            n = n+1
        if firma == 1:
            lista_v1_P4.append(cont_error)
            m = m+1
        if firma == 2:
            lista_v1_P3_Fy.append(cont_error)
            p = p+1  
            
        if firma == 3:
            lista_v1_P3_2.append(cont_error)
            r = r+1  
            
        if firma == 4:
            lista_v1_P3_Fy_2.append(cont_error)
            t = t+1  
            
    print("Solucion aproximada newton_rhapson: {:.4f}".format(x2))
    if firma == 0:
        print("Con {:d} pasos".format(n))
    if firma == 1:
         print("Con {:d} pasos".format(m))
    if firma == 2:
         print("Con {:d} pasos".format(p))
    if firma == 3:
         print("Con {:d} pasos".format(r))
    if firma == 4:
         print("Con {:d} pasos".format(t))


#### METODO DE LA SECANTE ####
         
def secante(error_relativo,funcion,x1,x2,tol,firma):
    
    error = 100
    x3 = 0
    cont_error = 0
    n = 0
    m = 0
    p = 0
    r = 0
    t = 0
    while error > tol:
        
        x3 = x1 - ((x2-x1))/((funcion(x2)-funcion(x1)))*funcion(x1)
        x1 = x2
        cont_error = abs(error_relativo(x3,x2))
        x2 = x3
        error = round(abs(funcion(x3)),6)
        if firma == 0:
            lista_v2_P3.append(cont_error)
            n = n+1
        if firma == 1:
            lista_v2_P4.append(cont_error)
            m = m+1
        if firma == 2:
            lista_v2_P3_Fy.append(cont_error)
            p = p+1 
            
        if firma == 3:
            lista_v2_P3_2.append(cont_error)
            r = r+1  
            
        if firma == 4:
            lista_v2_P3_Fy_2.append(cont_error)
            t = t+1
            
    print("Solucion aproximada-Secante- : {:.6f}".format(x2))
    
    if firma == 0:
        print("Con {:d} pasos".format(n))
    if firma == 1:
         print("Con {:d} pasos".format(m))
    if firma == 2:
         print("Con {:d} pasos".format(p))       
    if firma == 3:
         print("Con {:d} pasos".format(r))
    if firma == 4:
         print("Con {:d} pasos".format(t))



#### METODO DE LA FALSA POSICION ####

def falsa_posicion(error_relativo,funcion,xi,xu,tol,firma):
    
    error = 0
    
    xr = (xi*funcion(xu)-xu*funcion(xi))/(funcion(xu)-funcion(xi))                                            
    cont_error = 0
    n = 0
    m = 0
    p = 0
    r = 0
    t = 0
    while round(abs(funcion(xr)),6) > tol:
        
        xr = (xi*funcion(xu)-xu*funcion(xi))/(funcion(xu)-funcion(xi))
        if funcion(xi)*funcion(xr) < 0:
           
            xu = xr
            
        else:
           
            xi = xr

        cont_error = error_relativo(xr,error)
        
        if firma == 0:
            lista_v3_P3.append(cont_error)
            n = n+1
        if firma == 1:
            lista_v3_P4.append(cont_error)
            m = m+1       
        if firma == 2:
            lista_v3_P3_Fy.append(cont_error)
            p = p+1            
       
        if firma == 3:
            lista_v3_P3_2.append(cont_error)
            r = r+1  
            
        if firma == 4:
            lista_v3_P3_Fy_2.append(cont_error)
            t = t+1     
        
        error = xr

          
    print("Solucion aproximada - Falsa Posicion- : {:.6f}".format(xr))
   
    if firma == 0:
        print("Con {:d} pasos".format(n))
    if firma == 1:
         print("Con {:d} pasos".format(m))
    if firma == 2:
         print("Con {:d} pasos".format(p))
    if firma == 3:
         print("Con {:d} pasos".format(r))
    if firma == 4:
         print("Con {:d} pasos".format(t))



#### VALORES INTRODUCIDOS  ####

#### PREGUNTA 1 ####
         
x_g_P4= np.zeros([3])  # Listas donde se guardaras las soluciones
x_lu_P4 = np.zeros([3]) # de Gauss, Lu y cholesky
x_ch_P4 = np.zeros([3])


#Defini tres valores de  "matriz" porque uno de las def_funcion "chanca" los  
#valores de "matriz", asi que mejor defini la matriz y el vector 3 veces 

matriz_P1_1 = np.array([[4.0,-1.0,-1.0,-1.0],[-1.0,3.0,0.0,-1.0],[-1.0,0.0,3.0,-1.0],[-1.0,-1.0,-1.0,4.0]])
vector_P1_1 = np.array([5.0,0.0,5.0,0.0]) # Valores de la matriz de resistencias

matriz_P1_2 = np.array([[4.0,-1.0,-1.0,-1.0],[-1.0,3.0,0.0,-1.0],[-1.0,0.0,3.0,-1.0],[-1.0,-1.0,-1.0,4.0]])
vector_P1_2 = np.array([5.0,0.0,5.0,0.0])

matriz_P1_3 = np.array([[4.0,-1.0,-1.0,-1.0],[-1.0,3.0,0.0,-1.0],[-1.0,0.0,3.0,-1.0],[-1.0,-1.0,-1.0,4.0]])
vector_P1_3 = np.array([5.0,0.0,5.0,0.0])

x_g_P1 = Gauss(matriz_P1_1,vector_P1_1,4 )
x_lu_P1 = LU(matriz_P1_2,vector_P1_2,4)
x_ch_P1 = Cholesky(matriz_P1_3,vector_P1_3,4)


print("Pregunta 1,punto 2")
print("SOLUCION CON EL METODO DE GAUSS",x_g_P1)
print("SOLUCION CON EL METODO DE LU",x_lu_P1)
print("SOLUCION CON EL METODO DE CHOLESKY",x_ch_P1)
print("\n")




#### PREGUNTA 2 ####

x_g_P1 = np.zeros([4])
x_lu_P1 = np.zeros([4])
x_ch_P1= np.zeros([4])


#Defini tres valores de  "matriz" porque uno de las funciones "chanca" los  
#valores de "matriz", asi que mejor defini la matriz y el vector 3 veces 

matriz_P4_1 = np.array([[2.0,-1.0,0.0],[-1.0,2.0,-1.0],[0.0,-1.0,2.0]])
vector_P4_1 = np.array([-4.0,1.2,0.78])

matriz_P4_2 = np.array([[2.0,-1.0,0.0],[-1.0,2.0,-1.0],[0.0,-1.0,2.0]])
vector_P4_2 = np.array([-4.0,1.2,0.78])

matriz_P4_3 = np.array([[2.0,-1.0,0.0],[-1.0,2.0,-1.0],[0.0,-1.0,2.0]])
vector_P4_3 = np.array([-4.0,1.2,0.78])


x_g_P4 = Gauss(matriz_P4_1,vector_P4_1,3 )
x_lu_P4 = LU(matriz_P4_2,vector_P4_2,3)
x_ch_P4 = Cholesky(matriz_P4_3,vector_P4_3,3)


print("Pregunta 2,punto 2")
print("SOLUCION CON EL METODO DE GAUSS",x_g_P4)
print("SOLUCION CON EL METODO DE LU",x_lu_P4)
print("SOLUCION CON EL METODO DE CHOLESKY",x_ch_P4)
print("\n")




#### PREGUNTA 3 ####

a = 2         # Datos 
q = 1e-4
Q = 2e-5
F_x = 1.56
F_y = -2.17
e_0 = 8.854e-12
m = 0

v = (Q*q*sin(4*(np.pi)**2))/(4*np.pi*e_0*F_x) # Constantes definidas 
v_x = (Q*q*sin(4*(np.pi)**2))/(4*np.pi*e_0)
v_y = (-2*a*q*Q*(0.738))/(4*np.pi*e_0)
v_y1 = (-2*a*q*Q*(0.738))/(4*np.pi*e_0*F_y)

lista_v1_P4= []   # Listas donde se guardaran los errores
lista_v2_P4= []
lista_v3_P4= []

lista_v1_P3= []
lista_v2_P3= []
lista_v3_P3= []

lista_v1_P3_Fy= []
lista_v2_P3_Fy= []
lista_v3_P3_Fy= []

lista_v1_P3_2= []
lista_v2_P3_2= []
lista_v3_P3_2= []

lista_v1_P3_Fy_2= []
lista_v2_P3_Fy_2= []
lista_v3_P3_Fy_2= []



#### GRAFICA  PREGUNTA 3 PUNTO 2 ####

plt.title("Grafica Fx, Fy, F' VS x")
x = np.arange(-10,10,0.1)
x_1 = np.arange(0,10,0.1)
plt.ylim(-7,7)
plt.plot(x,funcion_Fx(x),'r',x,funcion_Fy(x),'b',x_1, funcion_F(x_1),'k')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('Fx','Fy','F'))
plt.xlabel('Eje X')
plt.ylabel('Fx, Fy, F')
plt.show()




#### GRAFICA  PREGUNTA 3 PUNTO 3 ####

plt.title("Fx(x) vs x")
x = np.arange(0,3,0.1)
plt.ylim(-2,7)
plt.plot(x,funcion_P3(x),'k')
plt.xlabel('x')
plt.ylabel('Funcion(x)')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('Fx(x) = xqQsen(4pi)^2/4pie_0d^3 ',))
plt.show()


print("Pregunta 3,punto 3")
print("\n")
print("Primera Solucion")
newton_rhapson(error_relativo,funcion_P3,f_d_P3,0.5,1e-6,0) #Punto inicial x = 1.5
secante(error_relativo,funcion_P3,0.5,0.4,1e-6,0)
falsa_posicion(error_relativo,funcion_P3,0.4,0.5,1e-6,0)
print("Segunda Solucion")
newton_rhapson(error_relativo,funcion_P3,f_d_P3,3,1e-6,3) #Punto inicial x = 1.5
secante(error_relativo,funcion_P3,3,2.8,1e-6,3)
falsa_posicion(error_relativo,funcion_P3,2.8,3,1e-6,3)
print("\n")




#### GRAFICA DE ERRORES PREGUNTA 3 PUNTO 3 ####

plt.title("Error relativo vs Interacciones para Fx = 1.56 :Solucion 0.9800")
plt.plot(range(0,4),lista_v1_P3,'b-o',range(0,6),lista_v2_P3,'g-o',range(0,14),lista_v3_P3,'k-o')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('newton rhapson','secante','falsa posicion'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()

plt.title("Error relativo vs Interacciones para Fx = 1.56 :Solucion 1.9855")
plt.plot(range(0,5),lista_v1_P3_2,'b-o',range(0,6),lista_v2_P3_2,'g-o',range(0,18),lista_v3_P3_2,'k-o')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('newton rhapson','secante','falsa posicion'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()


#### GRAFICA  PREGUNTA 3 PUNTO 4 ####

plt.title("Fy(x) vs x")
x = np.arange(-10,10,0.1)
plt.ylim(-10,20)
plt.plot(x,funcion_P3_Fy(x),'k-')
plt.xlabel('x')
plt.ylabel('Fy(x)')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('Fy(x) = -2Qqa(0.738)/4pie_0d^3',))
plt.show()


print("Pregunta 3,punto 4")
print("\n")
print("Primera Solucion")
newton_rhapson(error_relativo,funcion_P3_Fy,f_d_P3_Fy,0.5,1e-6,2) #Punto inicial x = 1.5
secante(error_relativo,funcion_P3_Fy,0.5,0.4,1e-6,2)
falsa_posicion(error_relativo,funcion_P3_Fy,0.4,0.5,1e-6,2)
print("Segunda Solucion")
newton_rhapson(error_relativo,funcion_P3_Fy,f_d_P3_Fy,-3,1e-6,4) #Punto inicial x = 1.5
secante(error_relativo,funcion_P3_Fy,-2.8,-3,1e-6,4)
falsa_posicion(error_relativo,funcion_P3_Fy,-2.8,-3,1e-6,4)

print("\n")




#### GRAFICA DE ERRORES PREGUNTA 3 PUNTO 4 ####

plt.title("Error relativo vs Interacciones para Fy = -2.17 : solucion 2.1035")
plt.plot(range(0,6),lista_v1_P3_Fy,'b-o',range(0,8),lista_v2_P3_Fy,'g-o',range(0,20),lista_v3_P3_Fy,'k-o')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('newton rhapson','secante','falsa posicion'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()

plt.title("Error relativo vs Interacciones para Fy = -2.17 : solucion -2.1035")
plt.plot(range(0,4),lista_v1_P3_Fy_2,'b-o',range(0,5),lista_v2_P3_Fy_2,'g-o',range(0,9),lista_v3_P3_Fy_2,'k-o')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('newton rhapson','secante','falsa posicion'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()


#### GRAFICA  PREGUNTA 4  ####

plt.title("Funcion(h) vs h")
x = np.arange(0,3,0.1)
plt.ylim(-5,20)
plt.plot(x,funcion_P4(x),'k-')
plt.xlabel('h')
plt.ylabel('Funcion(h)')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('F(h) = 5h^3-15h^2+16 ',))
plt.show()

print("Pregunta 4,punto 2")
newton_rhapson(error_relativo,funcion_P4,f_d_P4,0.5,1e-6,1) #Punto inicial x = 1.5
secante(error_relativo,funcion_P4,0.5,0.4,1e-6,1)
falsa_posicion(error_relativo,funcion_P4,0.4,0.5,1e-6,1)

print("\n")



#### GRAFICA  DE ERRORES PREGUNTA 2  ####

plt.title("Error relativo vs Interacciones -Pregunta 4")
plt.plot(range(0,5),lista_v1_P4,'b-o',range(0,6),lista_v2_P4,'g-o',range(0,8),lista_v3_P4,'k-o')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('newton rhapson','secante','falsa posicion'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()
