# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:37:03 2020

@author: JOSE
"""

#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore") 


def logistic(xn,r):
    
    return r* xn * (1-xn)

 
imagen = plt.figure(figsize=(30,30))
plt.axis([0.4,4,0,1])

lista_equi =[]

i = 0

for r in arange(1,4,0.01):
    lista_vi= []
    x = 0.5
    lista_vi.append(x)
    
    while i < 5000:
       
       x = logistic(lista_vi[i],r)
       lista_vi.append(x) 
       i = i+1
    
    lista_equi.append(lista_vi) 
    plt.plot([r]*4000,lista_vi[1000:5000],'k.')

plt.title("Mapa Logistico")
plt.xlabel('Taza de crecimiento')
plt.ylabel('Poblacion de equilibrio')
plt.show()

