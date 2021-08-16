# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 18:44:23 2020

@author: JOSE
"""


from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

corriente=[]

x03 = 2.097e-7

po03 = [5.06e-8 - 4.94e-8, 5.47e-8 - 4.94e-8, 2.26e-7 - x03,2.588e-7 - x03, 2.912e-7 - x03,3.228e-7 - x03,3.528e-7 - x03,
      3.668e-7 - x03]
cu03 = [0.5,3,10,30,50,70,90,100]


x15 = 4.94e-8
po15 = [5.02e-8 - x15,5.2e-8 - x15,5.29e-8 - x15, 6.02e-8 - x15,8.02e-8 - x15,1.51e-7 - x15]
cu15 = [0.5,2,3,10,30,100]

x06 = 4.94e-8
po06 = [5.1e-8 - x15,5.7e-8 - x15, 1.212e-7 - x15, 1.72e-7 - x15]
cu06 = [0.5,3,30,50]


for i in range (0,len(po03)):
    
       po03[i] = round(po03[i]/2e-9,3)
       
       
for i in range (0,len(po15)):
    
       po15[i] = round(po15[i]/2e-9,3)   
       
       
       
for i in range (0,len(po06)):
    
       po06[i] = round(po06[i]/2e-9,3)          
       
    
plt.ylim(0.1,200)
plt.xlim(0.1,200)
plt.yscale('log')
plt.xscale('log')
plt.title( "Densidad de corriente vs Velocidad")
plt.rc('legend', fontsize='medium')
plt.plot(cu15,po15,'yo:',cu03,po03,'bo:',cu06,po06,'ro:')
plt.gca().legend(('b = 0.15','b = 0.3','b = 0.6'))
plt.xlabel('J(MA/cm2)')
plt.ylabel('Velocidad (m/s)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()