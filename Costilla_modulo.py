#! /usr/bin/python3
from numpy import*
import numpy as np
import warnings
warnings.filterwarnings("ignore")


############# Pregunta 1 ###################

def C_G_L (cv_fun,a): # Cuadratura de Gauss-Legendre
    
    wi = [128/225, (322 + 13*70**(1/2))/900,(322 + 13*70**(1/2))/900,  (322 - 13*70**(1/2))/900, (322 - 13*70**(1/2))/900]
    xi = [0 , (1/3)*(5 - 2*(10/7)**(1/2))**(1/2),-(1/3)*(5 - 2*(10/7)**(1/2))**(1/2), (1/3)*(5 + 2*(10/7)**(1/2))**(1/2),-(1/3)*(5 + 2*(10/7)**(1/2))**(1/2)]

    n = len(xi)
    
    per = 0
    
   
    for i in range(0,n):
        
        per = per  + wi[i]*cv_fun(xi[i],a) 
          
    return per


############# Pregunta 2 ###################

pts_x=[] # Aqui guardo los valores de la posicion X
pts_y=[] # Aqui guardo los valores de la posicion Y



def runge_kutta4_2(pts_y,pts_x,f1,f2,f3,f4,t0,x0,y0,z0,r0,h): # Runge kutta4 para 4 variables
    
    x1 = 0
    y1 = 0
    z1 = 0
    r1 = 0

    while t0 < 10:
    
        k0 = h*f1(t0,x0,y0,z0,r0)
        l0 = h*f2(t0,x0,y0,z0,r0)
        n0 = h*f3(t0,x0,y0,z0,r0)
        m0 = h*f4(t0,x0,y0,z0,r0)
        
        k1 = f1(t0 + h/2, x0 + k0/2, y0 + l0/2, z0 + n0/2, r0 +m0/2)*h
        l1 = f2(t0 + h/2, x0 + k0/2, y0 + l0/2, z0 + n0/2, r0 +m0/2)*h
        n1 = f3(t0 + h/2, x0 + k0/2, y0 + l0/2, z0 + n0/2, r0 +m0/2)*h
        m1 = f4(t0 + h/2, x0 + k0/2, y0 + l0/2, z0 + n0/2, r0 +m0/2)*h
        
        k2 = f1(t0 + h/2, x0 + k1/2, y0 + l1/2, z0 + n1/2, r0 + m1/2)*h
        l2 = f2(t0 + h/2, x0 + k1/2, y0 + l1/2, z0 + n1/2, r0 + m1/2)*h
        n2 = f3(t0 + h/2, x0 + k1/2, y0 + l1/2, z0 + n1/2, r0 + m1/2)*h
        m2 = f4(t0 + h/2, x0 + k1/2, y0 + l1/2, z0 + n1/2, r0 + m1/2)*h
        
        k3 = f1(t0 + h, x0 + k2, y0 + l2, z0 + n2, r0 + m2)*h
        l3 = f2(t0 + h, x0 + k2, y0 + l2, z0 + n2, r0 + m2)*h
        n3 = f3(t0 + h, x0 + k2, y0 + l2, z0 + n2, r0 + m2)*h
        m3 = f4(t0 + h, x0 + k2, y0 + l2, z0 + n2, r0 + m2)*h
        
        x1 = x0 + (k0+ 2*k1+ 2*k2 + k3)/6
        y1 = y0 + (l0+ 2*l1+ 2*l2 + l3)/6
        z1 = z0 + (n0+ 2*n1+ 2*n2 + n3)/6
        r1 = r0 + (m0+ 2*m1+ 2*m2 + m3)/6
         
        pts_x.append(z0)
        pts_y.append(r0)

        t0 = t0 + h
        x0 = x1
        y0 = y1
        z0 = z1
        r0 = r1
    return(x1,y1,z1,r1)