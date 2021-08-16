# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:22:06 2020

@author: JOSE
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 20:44:19 2020

@author: JOSE
"""

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")


J = [1,2,3,4,5]
Jc = [1,2,3,4]

velo_j1 = np.zeros([10])
velo_j2 = np.zeros([10])
velo_j3 = np.zeros([10])
velo_j4 = np.zeros([10])
velo_j5 = np.zeros([9])

dmi = [2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.6]
dmic = [2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4]
ku = [0.73,0.76,0.8,0.86,0.89]

#Velocidades de 1 - 5 para DMI = 2.6 - 3.6 Para Ku = 0.8
#xnumero = [j1(xi - xf), j2(xf - xi) , ...] 
x26 = [3.67-4.6,4.11-5.75, 4.55-6.85,4.95-7.86, 5.31-8.649]
x27 = [3.852- 4.905,4.41-6.3,4.91 - 7.65, 5.45 - 8.91,5.9 - 10.02  ]
x28 = [4.05 - 5.205, 4.665 - 6.795, 5.3 - 8.37, 5.9- 9.86,6.51 - 11.145 ]
x29 = [4.25 - 5.49,4.97 - 7.27,5.7 - 9.07, 6.45 - 10.79,7.085 - 12.259   ]
x30 = [4.472 - 5.772, 5.278 - 7.748,6.15 - 9.765 , 6.96 - 11.67, 7.66- 13.35]
x31 = [4.64 - 6.015, 5.595 - 8.235, 6.51-10.395, 7.44-12.495, 8.24-14.38 ]
x32 = [4.815 - 6.27, 5.835 - 8.61, 6.88 - 10.995, 7.91 - 13.305, 8.82 - 15.36]
x33 = [4.98 - 6.465, 6.09 - 9.015, 7.25 - 11.58, 8.32-14.01,9.35 - 16.24 ]
x34 = [5.16 - 6.68, 6.34 - 9.4, 7.54 - 12.08, 8.72 - 14.74, 9.82 - 17.1]
x36 = [5.44 - 7.06, 6.8 - 10.08, 8.16 - 13.04 , 9.465 - 16.04]







#Velocidades de 1 - 5 para Ku = 0.73 - 0.93 Para DMI = 3.0


k073 = [7.92 - 13.32] #10:29
k076 = [7.53 - 12.68] #10:33
k080 = [6.96 - 11.67]
k086 = [6.11 - 10.21] #10:38
k089 = [5.74 - 9.54] #10:41


velo_ku = [k073[0],k076[0],k080[0],k086[0],k089[0]]

for i in range(0,5):
    velo_ku[i] = -round(velo_ku[i]/0.1,4) 

print(velo_ku)
n = len(x26)
m = len(x36)

for i in range (0,n):
    x26[i] = -round(x26[i]/0.1,4)
    x27[i] = -round(x27[i]/0.1,4)
    x28[i] = -round(x28[i]/0.1,4)
    x29[i] = -round(x29[i]/0.1,4)
    x30[i] = -round(x30[i]/0.1,4)
    x31[i] = -round(x31[i]/0.1,4)
    x32[i] = -round(x32[i]/0.1,4)
    x33[i] = -round(x33[i]/0.1,4)
    x34[i] = -round(x34[i]/0.1,4)



for i in range (0,m):
     x36[i] = -round(x36[i]/0.1,4)
    

velo_j1[0] = x26[0]
velo_j1[1] = x27[0]
velo_j1[2] = x28[0]
velo_j1[3] = x29[0]
velo_j1[4] = x30[0]
velo_j1[5] = x31[0]
velo_j1[6] = x32[0]
velo_j1[7] = x33[0]
velo_j1[8] = x34[0]
velo_j1[9] = x36[0]

velo_j2[0] = x26[1]
velo_j2[1] = x27[1]
velo_j2[2] = x28[1]
velo_j2[3] = x29[1]
velo_j2[4] = x30[1]
velo_j2[5] = x31[1]
velo_j2[6] = x32[1]
velo_j2[7] = x33[1]
velo_j2[8] = x34[1]
velo_j2[9] = x36[1]

velo_j3[0] = x26[2]
velo_j3[1] = x27[2]
velo_j3[2] = x28[2]
velo_j3[3] = x29[2]
velo_j3[4] = x30[2]
velo_j3[5] = x31[2]
velo_j3[6] = x32[2]
velo_j3[7] = x33[2]
velo_j3[8] = x34[2]
velo_j3[9] = x36[2]

velo_j4[0] = x26[3]
velo_j4[1] = x27[3]
velo_j4[2] = x28[3]
velo_j4[3] = x29[3]
velo_j4[4] = x30[3]
velo_j4[5] = x31[3]
velo_j4[6] = x32[3]
velo_j4[7] = x33[3]
velo_j4[8] = x34[3]
velo_j4[9] = x36[3]


velo_j5[0] = x26[4]
velo_j5[1] = x27[4]
velo_j5[2] = x28[4]
velo_j5[3] = x29[4]
velo_j5[4] = x30[4]
velo_j5[5] = x31[4]
velo_j5[6] = x32[4]
velo_j5[7] = x33[4]
velo_j5[8] = x34[4]



########## Figuras aparte #############


fig=plt.figure()
ax=fig.add_subplot(111, label="1")
ax2=fig.add_subplot(111, label="2", frame_on=False)

ax.scatter(dmi, velo_j4, color="C0")
ax.set_xlabel("DMI (mJ/m2)" )
ax.set_ylabel("Velocidad (m/s)")
ax.tick_params(axis='x')
ax.tick_params(axis='y')

ax2.scatter(ku,velo_ku, color="C1")
ax2.xaxis.tick_top()      
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.set_xlabel("Ku (MJ/m3)" )
ax2.set_ylabel("Velocidad (m/s)")
plt.show()