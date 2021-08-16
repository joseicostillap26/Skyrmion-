#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import Costilla_modulo


####################### Pregunta 1 ####################


print('-------------PREGUNTA 1------------')
print('\n')



def fun(x,a): # Funcion Periodo

    return 8**(1/2)/(a**4 - x**4)**(1/2)

def cv_fun(x,a): # Funcion cambio de variable (cv) de la funcion periodo,(x = xi + 1)

    return (a/2)*fun((a/2)*x+(a/2),a)


print("El periodo para a = 2 es = ",Costilla_modulo.C_G_L(cv_fun,2))
t = np.arange(0.01,2,0.01)
plt.title( "Periodo T vs a")
plt.rc('legend', fontsize='medium')
plt.yscale('log')
plt.plot(t,Costilla_modulo.C_G_L(cv_fun,t),'k')
plt.xlabel('a')
plt.ylabel('Periodo T ')
plt.show()



####################### Pregunta 2 ####################



def f1(t,w,g,x,y): # Funciónes definidas del examen final
    G = 1
    M = 10
    L = 2
    return -G*M*(x/((x**2+y**2)*(x**2+y**2 + (L**2)/4)**(1/2)))

def f2(t,w,g,x,y): # '''
    G = 1
    M = 10
    L = 2
    return -G*M*(y/((x**2+y**2)*(x**2+y**2 + (L**2)/4)**(1/2)))

def f3(t,w,g,x,y): # '''

    return w

def f4(t,w,g,x,y): # '''
       
    return g


rk2 = Costilla_modulo.runge_kutta4_2( Costilla_modulo.pts_y, Costilla_modulo.pts_x,f1,f2,f3,f4,0,0,1,1,0,0.01)

plt.title( "X VS Y")
plt.rc('legend', fontsize='medium')
plt.plot( Costilla_modulo.pts_y, Costilla_modulo.pts_x,'b')
plt.gca().legend(('Orbita',''))
plt.xlabel('Y')
plt.ylabel('X')
plt.show()


