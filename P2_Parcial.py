# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 15:08:08 2020

@author: JOSE
"""


#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

lista_v1= []

G = 6.674*1e-11
M = 5.974*1e24
m = 7.348*1e22
R = 3.844*1e8
w = 2.662*1e-6



def funcion(x):
    return (G*M/x**2) - (G*m/(R-x)**2) - x*w**2

def f_d(x):
    return  -(2*G*M/x**3) - (2*G*m/(R-x)**3) - w**2

def error_relativo(x2,x1):
    return abs(x1-x2)/x2

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
        lista_v1.append(cont_error)
        n = n+1
    print("Solucion aproximada-Newton Rhapson- : {:.6f}".format(x2))
    print("Con {:d} pasos".format(n))
    

newton_rhapson(error_relativo,funcion,f_d,1,1e-6)


