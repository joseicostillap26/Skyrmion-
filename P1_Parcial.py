
#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


lista_v1= []
lista_v2= []
lista_v3= []

inte1 = []
inte2 = []
inte3 = []
a0_a1 = [] 


def funcion_1(x):
    return  -6/x**(7) + np.exp(-x)


def f_d_1(x):
    return  42/x**(8) - np.exp(-x)

def funcion_12(x):
    return  1/x**(6) - np.exp(-x) + 0.1

def error_relativo(x2,x1):
    return abs(x1-x2)/x2


def Modele(x,a0,a1):        #### Modelo Err(k) = bx^a
    print('alpha = ',np.exp(a0))
    print('beta =',a1)
    return(np.exp(a0)*x**(a1))


def newton_rhapson(error_relativo,funcion,f_d,x1,tol):
    error = 10
    x2 = 0
    cont_error = 0
    n = 1
    lista_v1.append(1)
    inte1.append(1)
    while error > tol:
        
        x2 = x1 - funcion(x1)/f_d(x1)
        cont_error = error_relativo(x2,x1)
        x1 = x2
        error = round(abs(funcion(x2)),6)  
        lista_v1.append(cont_error)
        n = n+1
        inte1.append(n)
    return(n,x2)
    

def secante(error_relativo,funcion,x1,x2,tol,raiz): 
    
    error = 100
    x3 = 0
    cont_error = 0
    n = 1
    m = 1
    lista_v2.append(1)
    inte2.append(1)
    lista_v3.append(1)
    inte3.append(1)

    while error > tol:
        
        x3 = x1 - ((x2-x1))/((funcion(x2)-funcion(x1)))*funcion(x1)
        x1 = x2
        cont_error = abs(error_relativo(x3,x2))
        x2 = x3
        error = round(abs(funcion(x3)),6)
        
        if raiz == 1:
            lista_v2.append(cont_error)
            n = n+1
            inte2.append(n)
        if raiz == 2:
            lista_v3.append(cont_error)
            m = m+1
            inte3.append(m)

    print("Solucion aproximada-Secante- : {:.6f}".format(x2))
    if raiz == 1:
        print("Con {:d} pasos".format(n))
    if raiz == 2:
         print("Con {:d} pasos".format(m))


def convert_1var(time,y_mesure):
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]) 
        A[i] = np.log(time[i])

    return(A,B)

def R_L_1var(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert_1var(time,y_mesure)

    for i in range(0,n):

        xy = xy + (A[i])*B[i]
        x = x + (A[i])
        y = y + B[i]
        x_2 = x_2 + (A[i])**2
        y_2 = y_2 + (B[i])**2
         

    a1 = (n*xy-x*y)/(n*x_2-x**2)
    a0 = y/n - a1*x/n
    r = (n*xy-x*y)/((n*x_2-x**2)*(n*y_2-y**2))**(1/2)
    a0_a1.append(a0)
    a0_a1.append(a1)
    
    return(r) ### Retorna la ecuacion lineal ####

newton_rhapson(error_relativo,funcion_1,f_d_1,1,1e-6)
secante(error_relativo,funcion_12,1.1,1,1e-6,1)
secante(error_relativo,funcion_12,3,2.9,1e-6,2)



R_L_1var(1,inte1,lista_v1)

R_L_1var(1,inte2,lista_v2)

R_L_1var(1,inte3,lista_v3)


t = np.arange(1,10,0.1)
plt.ylim(0,1.5)
plt.title( "Ajuste lIneal de Err(k) = ax^b")
plt.rc('legend', fontsize='medium')
plt.plot(inte1,lista_v1,'bo',t,Modele(t,a0_a1[0],a0_a1[1]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 4.4668,b = -3,9170 ')
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()


t = np.arange(1,10,0.1)
plt.ylim(0,1.5)
plt.title( "Ajuste lIneal de Err(k) = ax^b")
plt.rc('legend', fontsize='medium')
plt.plot(inte2,lista_v2,'bo',t,Modele(t,a0_a1[2],a0_a1[3]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 2.4556,b = -4,0203 ')
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()


t = np.arange(1,10,0.1)
plt.ylim(0,1.5)
plt.title( "Ajuste lIneal de Err(k) = ax^b")
plt.rc('legend', fontsize='medium')
plt.plot(inte3,lista_v3,'bo',t,Modele(t,a0_a1[4],a0_a1[5]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 3.3583,b = -5,1479')
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()
