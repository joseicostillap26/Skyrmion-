# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 21:29:12 2020

@author: JOSE
"""


from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

from scipy import integrate
x2 = lambda x: 1- x -4*x**3 + 2*x**5

i = integrate.quad(x2, -2, 4)

print(i)