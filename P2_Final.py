from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

###########################################################


############# PREGUNTA 2 ###########################

pts_w=[] # Listas donde guardaré w1,w2,theta1,theta2 y el tiempo
pts_g=[]
pts_x=[]
pts_y=[]

def f1(t,w,g,x,y): # Funciónes de la PC03
    G = 1
    M = 10
    L = 2
    return -G*M*(x/((x**2+y**2)*(x**2+y**2 + (L**2)/4)**(1/2)))

def f2(t,w,g,x,y):
    G = 1
    M = 10
    L = 2
    return -G*M*(y/((x**2+y**2)*(x**2+y**2 + (L**2)/4)**(1/2)))

def f3(t,w,g,x,y):

    return w

def f4(t,w,g,x,y):
       
    return g

def runge_kutta4_2(t0,x0,y0,z0,r0,h): # Runge kutta4 para el problema 2
    
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
        

        
        pts_w.append(x0)
        pts_g.append(y0)
        pts_x.append(z0)
        pts_y.append(r0)


     
        t0 = t0 + h
        x0 = x1
        y0 = y1
        z0 = z1
        r0 = r1
    return(x1,y1,z1,r1)


rk2 = runge_kutta4_2(0,0,1,1,0,0.01)

plt.title( "X VS Y")
plt.rc('legend', fontsize='medium')
plt.plot(pts_y,pts_x,'b')
plt.gca().legend(('Orbita',''))
plt.xlabel('Y')
plt.ylabel('X')
plt.show()


