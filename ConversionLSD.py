# -*- coding: utf-8 -*-
"""
Created on Wed May  5 18:03:50 2021

@author: JOSE
"""



from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


kei=[52.905,51.205,51.014,50.769,50.586,50.506,50.420,50.353, 50.321,50.222,50.203,
    50.172,50.127, 50.094, 50.078, 50.050, 50.014,49.924, 49.898, 49.884, 49.875,49.835]


x = len(kei)

cas=[47.095,48.795,48.986,49.231,49.231,49.494,49.580,49.647,49.729,49.778,49.797,
   49.873,49.906,49.922, 49.950, 49.986,  50.076, 50.102, 50.116,50.125,50.156 ,50.165]
y = len(cas)


time = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]


print(x,y)
plt.title( "Keiko vs Castle")
plt.rc('legend', fontsize='medium')
plt.plot(time,kei,'bo:',time,cas,'ro:')
plt.gca().legend(('Keiko','Castle'))
plt.xlabel('Actualizacion')
plt.ylabel('Porcentaje (%)')
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()