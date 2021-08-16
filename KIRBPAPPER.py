# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 18:54:07 2021

@author: JOSE
"""

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

# CPP, J = 1e11 A/m2, s = 5

y5=[3.0e-8,3.3e-8,3.8e-8,4.0e-8,4.3e-8,4.4e-8,4.5e-8,4.6e-8,4.65e-8,4.7e-8,4.7e-8]
x5 = [3.25e-8,3.325e-8,3.4e-8,3.6e-8,4.075e-8,4.9e-8,5.95e-8,7.225e-8,8.65e-8,10.2e-8,11.75e-8]
time=[0.1e-9,0.2e-9,0.3e-9,0.4e-9,0.5e-9,0.6e-9,0.7e-9,0.8e-9,0.9e-9,1e-9]





x25 = [2.75e-8,3.65e-8,4.975e-8,6.3e-8,7.625e-8,8.95e-8,5.95e-8,1.027e-7,1.116e-7,1.29e-7,1.425e-7]

x20 =[3.1e-8,4.15e-8,6e-8,8.1e-8,1.0275e-7,1.247e-7,1.47e-7,1.69e-7,1.915e-7,2.137e-7]

vel20 = (x20[9]-x20[8])/(time[9]-time[8])


print(vel20)







plt.title( "Densidad de corriente vs Velocidad")
plt.rc('legend', fontsize='medium')
plt.plot(time,x20,'ko:')
plt.gca().legend(('Sistema CPP','b  = 0.6','b = 0.3','b = 0.15'))
plt.xlabel('J(MA/cm2)')
plt.ylabel('Velocidad (m/s)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()



