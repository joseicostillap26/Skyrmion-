# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 02:38:24 2020

@author: JOSE
"""
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

from io import open

datos1=[]
datos2=[]
datos3=[]
datos4=[]
datos5=[]
datos6=[]

gatos1=[]
gatos2=[]
gatos3=[]
gatos4=[]
gatos5=[]
gatos6=[]


mz = [] ### = x
r = [] #### = y

mza = []
ra = []

archivo = open('new33.txt','r')
archivoa = open('new44.txt','r')

l_a = archivo.readlines()
l_a.pop(0)

l_b = archivoa.readlines()
l_b.pop(0)

for i in range (0,len(l_a),1): 
   dato1=l_a[i].find(',')
   dato2=l_a[i].find(',',dato1+1)
   dato3=l_a[i].find(',',dato2+1)
   dato4=l_a[i].find(',',dato3+1)
   dato5=l_a[i].find(',',dato4+1)
   dato6=l_a[i].find('\n',dato5+1)
   mz.append(float(l_a[i][dato2+1:dato3]))
   r.append(float(l_a[i][dato4+1:dato5]))

dim = len(mz)
dimr = len(r)

for i in range(0,dim):
    mz[i] = - mz[i]
    
for i in range(0,dimr):
    r[i] = (r[i] - 278.4e-9)*1e9
 
archivoa = open('new44.txt','r')

l_b = archivoa.readlines()
l_b.pop(0)

for i in range (0,len(l_b),1): 
   gato1=l_b[i].find(',')
   gato2=l_b[i].find(',',gato1+1)
   gato3=l_b[i].find(',',gato2+1)
   gato4=l_b[i].find(',',gato3+1)
   gato5=l_b[i].find(',',gato4+1)
   gato6=l_b[i].find('\n',gato5+1)
   mza.append(float(l_b[i][gato2+1:gato3]))
   ra.append(float(l_b[i][gato4+1:gato5]))

dima = len(mza)
dimra = len(ra)

for i in range(0,dima):
    mza[i] = - mza[i]
    
for i in range(0,dimr):
    ra[i] = (ra[i] - 4.080e-7)*1e9
    

plt.title( "mz vs distancia")
plt.rc('legend', fontsize='medium')
plt.xlim(-250,250)
plt.plot(r,mz,'k',ra,mza,'r')
plt.xlabel('nm')
plt.ylabel('mz')
plt.gca().legend(('Sin Mask','Con Mask'))
plt.savefig("sample.jpg", dpi=200)
plt.figure(figsize=(60,80))
plt.show()    