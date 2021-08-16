# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:41:42 2020

@author: JOSE
"""


#! /usr/bin/python3

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore") 
 


def logistic(xn):
    
    return 1* xn * (1-xn)

lista_vi1 = [0]*1000
lista_vi2 = [0]*1000
lista_vi3 = [0]*1000
inicial_value = np.arange(0.45,0.55,0.05)

    
for p,x in enumerate(lista_vi1):

   
    if p == 0:
       x = 0.45
    else:   
       x = logistic(lista_vi1[p-1])
       
    
    lista_vi1[p] = x
    
    
for p,x in enumerate(lista_vi2):
  
    if p == 0:
       x = 0.50
    else:   
       x = logistic(lista_vi2[p-1])
       
    
    lista_vi2[p] = x    
 
    
for p,x in enumerate(lista_vi3):

   
    if p == 0:
       x = 5.55
    else:   
       x = logistic(lista_vi3[p-1])
       
    
    lista_vi3[p] = x    

       
fig = plt.figure(figsize=(15,15))
plt.title("Grafica 1")
plt.plot(range(0,1000),lista_vi1,range(0,1000),lista_vi2,'g+',range(0,1000),lista_vi3,'kx') 
plt.rc('legend', fontsize='medium')
plt.gca().legend(('r = 0.45','r = 0.50','r=0.55'))

plt.yscale("log")
plt.xlabel('Interacciones')
plt.ylabel('Evolucion de los Xn+1')
plt.show()
