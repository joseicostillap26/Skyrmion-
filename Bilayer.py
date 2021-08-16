# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 17:37:00 2020

@author: JOSE
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 19:03:22 2020

@author: JOSE
"""

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")



xi = 1.0000e-7

pocpp = [1.2375e-7 - xi, 1.476e-7 - xi,1.719e-7 - xi, 1.962e-7 - xi, 2.196e-7 - xi, 2.4435e-7 - xi, 2.682e-7 - xi, 2.9295e-7 - xi, 3.1174e-7 - xi,3.4245e-7 - xi]
cucpp = [10,20,30,40,50,60,70,80,90,100]



n = 10

for i in range (0,n):
    pocpp[i] = round(pocpp[i]/0.4e-9,3)


#beta = 0.6

xi = 1.0000e-7

po06 = [1.224e-7 - xi, 1.449e-7 - xi,1.6785e-7 - xi, 1.908e-7 - xi, 2.133e-7 - xi, 2.362e-7 - xi, 2.592e-7 - xi, 2.814e-7 - xi, 3.046e-7 - xi,3.276e-7 - xi]
cu06 = [100,200,300,400,500,600,700,800,900,1000]

n = 10

for i in range (0,n):
    po06[i] = round(po06[i]/0.4e-9,3)



plt.xlim(5,1000)
plt.xscale('log')
plt.title( "Densidad de corriente vs Velocidad")
plt.rc('legend', fontsize='medium')
plt.plot(cucpp,pocpp,'ko:',cu06,po06,'ro:')
plt.gca().legend(('CPP'))
plt.xlabel('J(MA/cm2)')
plt.ylabel('Velocidad (m/s)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()