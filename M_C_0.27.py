# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 13:14:40 2021

@author: JOSE
"""

# alpha 0.27
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

corriente=[]

xi = 1.0000e-7

po1 = [1.3185e-7 - xi,1.6425e-7 - xi, 1.9665e-7 - xi,2.2905e-7 - xi,2.61e-7 - xi,
      2.9345e-7 - xi, 3.2625e-7 - xi, 3.5865e-7 - xi, 3.9105e-7 - xi, 4.239e-7 - xi]
cu1 = [10,20,30,40,50,60,70,80,90,100]


po2 =[ 1.318e-7 - xi, 1.48e-7 - xi, 1.64e-7 - xi, 1.808e-7 - xi, 1.97e-7 - xi, 2.136e-7 - xi,
      2.296e-7 - xi, 2.456e-7 - xi]
cu2 = [2,3,4,5,6,7,8,9]

po3 =[1.016e-7 - xi,1.038e-7 - xi, 1.062e-7 - xi, 1.08e-7 - xi, 1.104e-7 - xi, 1.126e-7 - xi,
      1.15e-7 - xi, 1.17e-7 - xi, 1.19e-7 - xi, 1.214e-7 - xi ]
cu3 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]


po = [1.016e-7 - xi,1.038e-7 - xi, 1.062e-7 - xi, 1.08e-7 - xi, 1.104e-7 - xi, 1.126e-7 - xi,
      1.15e-7 - xi, 1.17e-7 - xi, 1.19e-7 - xi, 1.214e-7 - xi,  1.318e-7 - xi, 1.48e-7 - xi, 1.64e-7 - xi, 1.808e-7 - xi,
      1.97e-7 - xi, 2.136e-7 - xi, 2.296e-7 - xi, 2.456e-7 - xi,1.3185e-7 - xi,1.6425e-7 - xi, 1.9665e-7 - xi,2.2905e-7 - xi,2.61e-7 - xi,
      2.9345e-7 - xi, 3.2625e-7 - xi, 3.5865e-7 - xi, 3.9105e-7 - xi, 4.239e-7 - xi]

cu = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]



for i in range (0,len(po)):
    
    if i <= len(po3):
        
       po[i] = round(po[i]/2e-9,3)
    
    if i > len(po3) and i < len(po3)+len(po2):
        
       po[i] = round(po[i]/1.5e-9,3)
  
    if i >= len(po2)+len(po1):
        
       po[i] = round(po[i]/0.3e-9,3)



#beta = 0.6


############################################################


po06 = [1.0035e-7 - xi, 1.0125e-7 - xi,1.0215e-7 - xi, 1.035e-7 - xi, 1.044e-7 - xi, 1.053e-7 - xi, 1.062e-7 - xi, 1.071e-7 - xi,1.08e-7 - xi,1.089e-7-xi,
         1.1205e-7 - xi,1.188e-7 - xi, 1.2465e-7 - xi, 1.3095e-7 - xi, 1.3725e-7 - xi, 1.44e-7 - xi, 1.503e-7 - xi, 1.566e-7 - xi,1.629e-7 - xi,
         1.503e-7 - xi,1.7595e-7 - xi, 2.007e-7 - xi, 2.259e-7 - xi, 2.5155e-7 - xi, 2.767e-7 - xi, 3.0124e-7 - xi,3.276e-7 - xi,3.528e-7 - xi]
cu06 = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]



for i in range (0,28):
    
    if i < 10:
        
       po06[i] = round(po06[i]/1.5e-9,3)
    
    if i >= 10 and i < 19:
        
       po06[i] = round(po06[i]/1e-9,3)
  
    if i >= 19:
        
       po06[i] = round(po06[i]/0.4e-9,3)


print(po06)


plt.xlim(0.1,1100)
plt.xscale('log')
plt.title( "Densidad de corriente vs Velocidad")
plt.rc('legend', fontsize='medium')
plt.plot(cu,po,'ko:',cu06,po06,'ro:')
plt.gca().legend(('Sistema CPP','b  = 0.6','b = 0.3','b = 0.15'))
plt.xlabel('J(MA/cm2)')
plt.ylabel('Velocidad (m/s)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()