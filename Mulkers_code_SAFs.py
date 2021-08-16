# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 15:57:31 2020

@author: JOSE
"""
# alpha 0.3
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


po061 = [1.449e-7 - xi,1.6785e-7 - xi, 1.908e-7 - xi, 2.133e-7 - xi, 2.362e-7 - xi, 2.592e-7 - xi, 2.814e-7 - xi, 3.046e-7 - xi,3.276e-7 - xi]
cu061 = [200,300,400,500,600,700,800,900,1000]#0.4ns



po062 = [ 1.108e-7 - xi,1.164e-7 - xi, 1.22e-7 - xi, 1.28e-7 - xi, 1.34e-7 - xi, 1.4e-7 - xi, 1.456e-7 - xi, 1.508e-7 - xi,1.568e-7 - xi]
cu062 = [20,30,40,50,60,70,80,90,100]#1ns

##########################################


po063 = [1.004e-7 - xi, 1.012e-7 - xi,1.024e-7 - xi, 1.032e-7 - xi, 1.04e-7 - xi, 1.048e-7 - xi, 1.056e-7 - xi, 1.064e-7 - xi,1.072e-7 - xi,1.08e-7-xi]
cu063 = [1,2,3,4,5,6,7,8,9,10]#1.5ns

############################################################


po06 = [1.004e-7 - xi, 1.012e-7 - xi,1.024e-7 - xi, 1.032e-7 - xi, 1.04e-7 - xi, 1.048e-7 - xi, 1.056e-7 - xi, 1.064e-7 - xi,1.072e-7 - xi,1.08e-7-xi,
         1.108e-7 - xi,1.164e-7 - xi, 1.22e-7 - xi, 1.28e-7 - xi, 1.34e-7 - xi, 1.4e-7 - xi, 1.456e-7 - xi, 1.508e-7 - xi,1.568e-7 - xi,
         1.449e-7 - xi,1.6785e-7 - xi, 1.908e-7 - xi, 2.133e-7 - xi, 2.362e-7 - xi, 2.592e-7 - xi, 2.814e-7 - xi, 3.046e-7 - xi,3.276e-7 - xi]
cu06 = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]


for i in range (0,len(po06)):
    
    if i < len(po063):
        
       po06[i] = round(po06[i]/1.5e-9,3)
    
    if i >= len(po063) and i < len(po063)+len(po062):
        
       po06[i] = round(po06[i]/1e-9,3)
  
    if i > len(po062)+len(po061):
        
       po06[i] = round(po06[i]/0.4e-9,3)



#beta = 0.3


po03 = [1.001e-7 - xi, 1.004e-7 - xi, 1.01e-7 - xi, 1.015e-7 - xi, 1.02e-7 - xi, 1.03e-7 - xi, 1.036e-7 - xi,1.04e-7 - xi,1.045e-7-xi,
         1.05e-7 - xi,1.065e-7 - xi, 1.10e-7 - xi, 1.135e-7 - xi, 1.17e-7 - xi, 1.21e-7 - xi, 1.245e-7 - xi,1.28e-7 - xi,1.31e-7 - xi,1.355e-7 - xi,
         1.71e-7 - xi, 2.064e-7 - xi, 2.416e-7 - xi, 2.776e-7 - xi, 3.132e-7 - xi, 3.488e-7 - xi,3.84e-7 - xi,4.2e-7 - xi,4.55e-7 - xi]

############################################################
cu03 = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]

for i in range (0,len(po03)):
    
    if i < len(po063):
        
       po03[i] = round(po03[i]/1.5e-9,3)
    
    if i >= len(po063) and i < len(po063)+len(po062):
        
       po03[i] = round(po03[i]/1e-9,3)
  
    if i > len(po062)+len(po061):
        
       po03[i] = round(po03[i]/1e-9,3)


#beta = 0.15

po015 = [1.001e-7 - xi, 1.003e-7 - xi, 1.005e-7 - xi, 1.007e-7 - xi, 1.0085e-7 - xi, 1.0125e-7 - xi,1.017e-7 - xi,1.0215e-7-xi,
         1.024e-7 - xi,1.026e-7 - xi, 1.030e-7 - xi, 1.053e-7 - xi, 1.071e-7 - xi, 1.089e-7 - xi, 1.115e-7 - xi,1.295e-7 - xi,1.147e-7 - xi,1.1655e-7 - xi,
         1.185e-7 - xi, 1.201e-7 - xi, 1.224e-7 - xi, 1.300e-7 - xi, 1.3725e-7 - xi, 1.4535e-7 - xi,1.53e-7 - xi,1.602e-7 - xi,1.6785e-7 - xi,1.755e-7 - xi]

############################################################
cu015 = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]


print(len(po015))

for i in range (0,len(po015)):
    
    if i < len(po063):
        
       po015[i] = round(po015[i]/1.5e-9,3)
    
    if i >= len(po063) and i < len(po063)+len(po062):
        
       po015[i] = round(po015[i]/1e-9,3)
  
    if i > len(po062)+len(po061):
        
       po015[i] = round(po015[i]/0.4e-9,3)


plt.xlim(0.1,1100)
plt.xscale('log')
plt.title( "Densidad de corriente vs Velocidad")
plt.rc('legend', fontsize='medium')
plt.plot(cu,po,'ro:',cu06,po06,'yo:',cu03,po03,'go:',cu015,po015,'bo:')
plt.gca().legend(('Sistema CPP','b  = 0.6','b = 0.3','b = 0.15'))
plt.xlabel('J(MA/cm2)')
plt.ylabel('Velocidad (m/s)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()