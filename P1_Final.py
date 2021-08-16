# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 18:42:33 2020

@author: JOSE
"""

#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")




def fun(x,a):

    return 8**(1/2)/(a**4 - x**4)**(1/2)

def cv_fun(x,a):

    return (a/2)*fun((a/2)*x+(a/2),a)

def C_G_L (a):
    
    wi = [128/225, (322 + 13*70**(1/2))/900,(322 + 13*70**(1/2))/900,  (322 - 13*70**(1/2))/900, (322 - 13*70**(1/2))/900]
    xi = [0 , (1/3)*(5 - 2*(10/7)**(1/2))**(1/2),-(1/3)*(5 - 2*(10/7)**(1/2))**(1/2), (1/3)*(5 + 2*(10/7)**(1/2))**(1/2),-(1/3)*(5 + 2*(10/7)**(1/2))**(1/2)]

    n = len(xi)
    
    per = 0
    
   
    for i in range(0,n):
        
        per = per  + wi[i]*cv_fun(xi[i],a) 
          
    return per

t = np.arange(0.01,2,0.01)
plt.title( "Periodo T vs a")
plt.rc('legend', fontsize='medium')
plt.yscale('log')
plt.plot(t,C_G_L(t),'k:',)
plt.xlabel('a')
plt.ylabel('Periodo T ')
plt.show()
