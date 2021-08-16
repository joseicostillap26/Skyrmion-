# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:44:42 2020

@author: JOSE
"""


#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

lista_v1= []
lista_v2= []
lista_v3= []
lista_v4= []
lista_v5= []
lista_v6= []

G = 6.674e-11
M = 5.974e24
m = 7.348e22
R = 3.844e8
w = 2.66210e-6

def funcion(x):
    return G*M/x**3-G*m/(x*(R-x)**2) - w**2

def funcion_FP(x):
    return x**2-x**(2/3) + 1


def funcion_g(x):
    return 5-5*np.exp(-x)

def f_d(x):
    return -5*np.exp(-x) + 1

def f_d_d(x):
    return 5*np.exp(-x) 

def error_relativo(x2,x1):
    return abs(x1-x2)/x2


def secante(error_relativo,funcion,x1,x2,tol):
    
    error = 100
    x3 = 0
    cont_error = 0
    n = 0
   
    while error > tol:
        
        x3 = x1 - ((x2-x1))/((funcion(x2)-funcion(x1)))*funcion(x1)
        x1 = x2
        cont_error = abs(error_relativo(x3,x2))
        x2 = x3
        error = round(abs(funcion(x3)),6)
        lista_v1.append(cont_error)
        n = n+1
    print("Solucion aproximada-Secante- : {:.6f}".format(x2))
    print("Con {:d} pasos".format(n))


def punto_fijo(error_relativo,funcion,funcion_g,x1,tol):
    
    error = 100
    x2 = funcion_g(x1)  
    x1 = x2 
    n = 0
    
    while error > tol:
        x2 = funcion_g(x1)
        cont_error = round(abs(error_relativo(x2,x1)),6)
        x1 = x2
        error = round(abs(funcion(x2)),6)
        lista_v2.append(cont_error)
        n = n+1 
    print("Solucion aproximada - Punto Fijo-: {:.6f}".format(x2))
    print("Con {:d} pasos".format(n))

def newton_rhapson(error_relativo,funcion,f_d,x1,tol):
    error = 10
    x2 = 0
    cont_error = 0
    n = 0
    
    while error > tol:
        
        x2 = x1 - funcion(x1)/f_d(x1)
        cont_error = error_relativo(x2,x1)
        x1 = x2
        error = round(abs(funcion(x2)),6)  
        lista_v3.append(cont_error)
        n = n+1
    print("Solucion aproximada-Newton Rhapson- : {:.6f}".format(x2))
    print("Con {:d} pasos".format(n))

def falsa_posicion(error_relativo,funcion,xi,xu,tol):
    
    error = 0
    
    xr = (xi*funcion(xu)-xu*funcion(xi))/(funcion(xu)-funcion(xi))                                            
    cont_error = 0
    n = 0
    
    while round(abs(funcion(xr)),6) > tol:
        
        xr = (xi*funcion(xu)-xu*funcion(xi))/(funcion(xu)-funcion(xi))
        if funcion(xi)*funcion(xr) < 0:
           
            xu = xr
            
        else:
           
            xi = xr
      
        n = n+1
       
        cont_error = error_relativo(xr,error)
        lista_v4.append(cont_error)
        error = xr
    print("Solucion aproximada - Falsa Posicion- : {:.6f}".format(xr))    
    print("Con {:d} pasos".format(n))

def bisecante(error_relativo,funcion,xi,xu):
    
    error = 0
    xr= (xi+xu)/2
    cont_error = 0
    n = 0
    
    while round(funcion(xi)*funcion(xr),6)!= 1e-6:
        
        
        xr= (xi+xu)/2
        
        if funcion(xi)*funcion(xr) < 0:
    
            xu = xr
            
        else:
           
            xi = xr
      
        n = n+1
        cont_error = error_relativo(xr,error)
        lista_v5.append(cont_error)
        error = xr
   
    print("Solucion aproximada -Bisecante- : {:.6f}".format(xr))
    print("Con {:d} pasos".format(n))

def newton_rhapson_mod(error_relativo,funcion,f_d,f_d_d,x1,tol):
    error = 10
    x2 = 0
    cont_error = 0
    n = 0
    
    while round(error,6) > tol:
        
        x2 = x1 - (funcion(x1)*f_d(x1))/((f_d(x1))**2 - funcion(x1)*f_d_d(x1))
        cont_error = error_relativo(x2,x1)
        x1 = x2
        error = round(abs(funcion(x2)),6)  
        lista_v6.append(cont_error)
        n = n+1
    print("Solucion aproximada -Newton Rhapson-Modificado : {:.6f}".format(x2))
    print("Con {:d} pasos".format(n))


secante(error_relativo,funcion,11,10,1e-6)
punto_fijo(error_relativo,funcion,funcion_g,1,1e-6)
newton_rhapson(error_relativo,funcion,f_d,6,1e-6)
falsa_posicion(error_relativo,funcion,3,5,1e-6)
bisecante(error_relativo,funcion,3,6) 
newton_rhapson_mod(error_relativo,funcion,f_d,f_d_d,6,1e-6) 

plt.title("Grafica Wien")
x = np.arange(-1,6,0.1)
plt.ylim(-4,10)
plt.plot(x,funcion(x),'k-')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('5exp(-x)+x-5 ',))
plt.xlabel('Eje X')
plt.ylabel('Eje Y')
plt.show()



plt.title("Error relativo vs Interacciones")
plt.plot(range(0,3),lista_v1,'b:',range(0,5),lista_v2,'g:',range(0,3),lista_v3,'r:',range(0,3),lista_v4,'m:',range(0,11),lista_v5,'k:',range(0,2),lista_v6,'y:')
plt.rc('legend', fontsize='medium')
plt.gca().legend(('secante','punto fijo','newton rhapson','falsa posicion','bisecante','NR modificado'))
plt.xlabel('Interacciones')
plt.ylabel('Error relativo')
plt.yscale("log")
plt.show()


