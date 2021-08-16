# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:35:23 2021

@author: JOSE
"""

from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")
P = 0.4
g = 2
ub = 9.274e-24
e = 1.6022e-19
Ms = 580e3
j = 100e10


u = -(P*g*ub*j)/(2*e*Ms)

print(u)
