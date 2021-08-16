#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
warnings.filterwarnings("ignore")

##########                ###########
##########  PREGUNTA # 1  ###########
##########                ###########

########### Sistema de ecuaciones #######################

#####  -(k1+k2)x1 + (k2)x2 = 0           #####
#####    (k2)x1 - (k2+k3)x2 + (k3)x3= 0  #####
#####   (k3)x2 - (k3+k4)x3 + (k4)x4 = 0  #####
#####   (k4)x3 - (k4)x4 = - F            #####

####### k1 = 150, k2 = 50 , k3 = 75, k4 = 225 y F =2000 #########

#####     -200x1 + 50x2 = 0                 #####
#####    50x1 - 125x2 + 75x3= 0             #####
#####   75x2 - 300x3 + 225x4 = 0            #####
#####   225x3 - 225x4 = - 2000              #####



lista_erreur= [] ## lista donde guardo los errres
inte= []       # Lista donde guardo las interaciones
a0_a1 = []     # Lista donde guardo a0 y a1 de la regresion lineal
p_cheby = np.zeros([20])
f_cheby = np.zeros([20])

def Modele(x,a0,a1):        #### Modelo Err(k) = bx^a
    print('beta = ',np.exp(a0))
    print('alpha =',a1)
    return(np.exp(a0)*x**(a1))


def error(x2,x1):
    
    return (x1-x2)**2

def sumalista(lista):
    laSuma = 0
    for i in lista:
        laSuma = laSuma + i
    return laSuma  


########## Aplico SOR con las subfunciones Din_b y Din_U , Din_L #### 


def Din_b(matriz,vector):
   
   dim = len(matriz)  
   x = np.zeros([dim])
   for i in range (0,dim):    
       x[i]= vector[i]/matriz[i][i]
   
   return(x)

def Din_U(matriz):
    
   dim = len(matriz) 
   x = np.zeros([dim,dim])
   z = np.zeros([dim,dim])

   for i in range (0,dim):
       for j in range (0,dim):
           if i < j:
               x[i][j] = matriz[i][j]
   for i in range (0,dim):
      k = 0
      while k < dim:
           z[i][k] =  x[i][k]/matriz[i][i]
           k = k + 1 
   
   return(z)

def Din_L(matriz):
    
   dim = len(matriz)
   y = np.zeros([dim,dim])
   z = np.zeros([dim,dim])
   
   for i in range (0,dim):
       for j in range (0,dim):
           if i > j:
               y[i][j] = matriz[i][j]

   for i in range (0,dim):
      k = 0
      while k < dim:
           z[i][k] =  y[i][k]/matriz[i][i]
           k = k + 1
           
   return(z)


def SOR(error_relativo,matriz,vector,tol,w): ## W escojido w = 1.5

   dim = len(matriz)
   lista_v3= np.zeros([dim])
   x = np.zeros([dim])
   y = np.zeros([dim])
   cont_error = 0
   erreur = 1
   N = 0
   
   while (erreur > tol):
     Sx_2 = 0
     for i in range (0,dim):
         z1 = 0
         z2 = 0
         for j in range(0,dim):
              z1 = z1 + ((1-w)*Din_L(matriz)[i][j]+Din_U(matriz)[i][j])*y[j] 
              z2 = z2 + (w*Din_L(matriz)[i][j])*x[j]
         x[i] =  Din_b(matriz,vector)[i] - z1 -z2

     for k in range (0,dim):
         cont_error = abs(error_relativo(x[k],y[k]))
         Sx_2 = Sx_2 + x[k]**2
         lista_v3[k] = cont_error
     
     for k in range(0,dim):
         y[k] = x[k] 
         x[k] = 0 # reset
 
     N = N + 1
     erreur = (sumalista(lista_v3)/Sx_2)**(1/2) ## Guardo los errores relativos
     lista_erreur.append(erreur)                ## Guardo los errores relativos
     inte.append(N)                            ## Guardo las interaciones
   return(y)


############### Regrasion Lineal ############

def convert1(time,y_mesure):
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]) 
        A[i] = np.log(time[i])

    return(A,B)


def R_L(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert1(time,y_mesure)

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
    print('Coeficiente de Correlación =',r)
    print('Coeficiente de Determinacion =',r**2)
    
    return(a0+t*a1) ### Retorna la ecuacion lineal ####


######## Matriz y vector correspondiente al sistema de ecuaciones hallado ###########
    #### k1 = 150, k2 = 50 , k3 = 75, k4 = 225 y F =2000 #########

A =  np.array([[-200,50,0,0],[50,-125,75,0],[0,75,-300,225],[0,0,225,-225]])
b = np.array([0,0,0,-2000])
print('PREGUNTA 1 \n')
print("Solucion del sistema de ecuaciones = ",SOR(error,A,b,1e-6,1.5))
 
R_L(1,inte,lista_erreur)

t = np.arange(1,65,0.1)
plt.ylim(-0.1,1.1)
plt.title( "P1 Ajuste lIneal de Err(k) = bx^a")
plt.rc('legend', fontsize='medium')
plt.plot(inte,lista_erreur,'b.',t,Modele(t,a0_a1[0],a0_a1[1]),'k',)
plt.gca().legend(('Error |Xi- Xi+1|/|xi+1| vs interacciones','Modelo Err(k)= bx^a',''))
plt.xlabel('# de Interacciones')
plt.ylabel('Error relativo')
plt.show()


##########                ###########
##########  PREGUNTA # 2  ###########
##########                ###########

datos1=[] ## Listas de apoyo
datos2=[]  ## '''

veloc = []  #### Lista donde guardo las velocidades
presion = []  #### Lista donde guardo las presiones


a0_a1 = [] # Lista donde guardo a0 y a1 de la regresion lineal (a0 + a1(x))

archivo = open('datos.dat','r') ## Abro data.dat y lo guardo en archivo##
                                ## Para abrir tanto datos.dat como el script
# tiene q estar en la misma direccion.

l_a = archivo.readlines() ## Leo data.dat todo junto como una lista

l_a.pop(0)  ## Elimino los valores de la lista l_a donde no hay numeros
l_a.pop(11) ## '''
l_a.pop(11) ##  '''

for i in range (0,len(l_a),1): 
   dato1=l_a[i].find('\t')            ## Separo los valores de la lista l_a
   dato2=l_a[i].find('\t',dato1+1)    ## En columnas  con respecto al separador \t
   veloc.append(float(l_a[i][0:dato1])) ## Luego agrego las columnas a veloc y presion
   presion.append(float(l_a[i][dato1+1:dato2]))


#### Mi polinomio interpolante es el polinomio de Newton #######
   
def Newton(t,matriz,vector):

    pts = len(t)
    dim = len(matriz)
    v = np.zeros([dim]*2)
    x_f = np.zeros([dim]) ### Lista donde guardo los coeficientes de newton###    
    z = 0
    lj = [] ### Lista donde guardo las funciones t(t-t1),(t-t1)(t-t2) ... para cada t ##
    lz = [] ### Lista donde guardo el valor de la funcion fn(t) para cada t ####

    for i in range(0,dim):
        v[i][0] = 1
        for j in range(0,dim-1):
          k = 0
          cont = 1
          while k < j+1:
                p = 1
                p = p*(matriz[i]-matriz[k]) 
                cont = cont*p
                k = k +1
          v[i][j+1] = cont

    ## Resuelvo el sistema Ax = b con Fordward ####
          
    x_f[dim-dim] = vector[dim-dim]/v[dim-dim][dim-dim]
    
    for r in range(1,dim,1):
         suma = 0
         for c in range(0,dim):
             suma = suma + v[r,c]*x_f[c]
         x_f[r] = (vector[r] - suma)/v[r,r]
         
    print("Coeficientes del polinomio interpolante de Newton = ",x_f)
   
    for m in range(0,pts):
      pi_j= np.zeros([dim])
      pi_j[0]= 1
      for j in range(0,dim-1):
          k = 0
          cont = 1
          while k < j+1:
                p = 1
                p = p*(t[m]-matriz[k]) 
                cont = cont*p
                k = k +1
          pi_j[j+1] = cont
      lj.append(pi_j)

    for i in range(0,pts):
      z = 0
      for j in range(0,dim):
          z = z + lj[i][j]*x_f[j]
      lz.append(z)

    return(lz)

def Modele(x,a0,a1):
    print('alpha = ',np.exp(a0))
    print('beta =',a1)
    return(np.exp(a0)*x**(a1))

def convert2(time,y_mesure):
    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    for i in range(0,n):
        B[i] = np.log(y_mesure[i]) 
        A[i] = np.log(time[i])

    return(A,B)

######### Regresion Lineal ##############
    
def R_L(t,time,y_mesure):

    n = len(y_mesure)
    A = np.zeros([n])
    B = np.zeros([n])
    xy=0
    x=0
    y=0
    x_2=0
    y_2 = 0
    A,B = convert2(time,y_mesure)

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
    print("\nPREGUNTA 2\n")
    print('Coeficiente de Correlación =',r)
    print('Coeficiente de Determinacion =',r**2)
    return(a0+t*a1)

R_L(1,presion,veloc)
t = np.linspace(20,200,100) ### 100 puntos equidistantes entre 20 y 200 ###

plt.title( "P2 Ajuste lineal de data.dat")
plt.rc('legend', fontsize='medium')
plt.plot(presion,veloc,'bo',t,Newton(t,presion,veloc),'r:',t,Modele(t,a0_a1[0],a0_a1[1]),'k')
plt.gca().legend(('velocidad vs Caida de presion','Polinomio de Newton','Modelo v = a(P)^b'))
plt.ylabel('Velocidad (m/s)')
plt.xlabel('Caida de Presion (mm Hg)')
plt.show()


##########                ###########
##########  PREGUNTA # 3  ###########
##########                ###########

def fun(x):
    return np.exp(-(x**2))*np.sin(50*x)


#### Mi polinomio interpolante es el polinomio de Lagrange #######

def Lagrange(t, matriz,vector):

    dim = len(matriz)
    lj= []
    z = 0
    for i in range(0,dim):
        p = 1 
        for k in range(0,dim):
              if k != i:
                   p = p*((t - matriz[k]) / (matriz[i]-matriz[k]))        
        lj.append(p)

    for i in range(0,dim):
          z = z + lj[i]*vector[i]
          
   
    return(z)

#### Puntos de Chebyshev ########
    
def punt_cheby(t,a,b,n):
    lj= []
    x = np.zeros([n])
    y = np.zeros([n])
    valor = 0
    z = 0

    for i in range(1,n):
        valor = a + (b-a)*(1-np.cos(i*np.pi/(n-1)))/2
        x[i] = valor
        p_cheby[i] = valor 
        y[i] = fun(valor)
        f_cheby[i] = fun(valor)
        

    for i in range(0,n):
        p = 1 
        for k in range(0,n):
              if k != i:
                   p = p*((t - x[k]) / (x[i]-x[k]))        
        lj.append(p)

    for i in range(0,n):
          z = z + lj[i]*y[i]
          
    return(z)


#### Puntos equidistantes ########

def punt_equi(t,a,b,n):
    
    lj= []
    x = np.zeros([n])
    y = np.zeros([n])
    valor = 0
    z = 0
    
    for i in range(0,n):
        valor = a+(b-a)*i/(n-1)
        x[i] = valor
        y[i] = fun(valor)

    for i in range(0,n):
        p = 1 
        for k in range(0,n):
              if k != i:
                   p = p*((t - x[k]) / (x[i]-x[k]))        
        lj.append(p)

    for i in range(0,n):
          z = z + lj[i]*y[i]
          
    return(z)

##### Errores en la interpolacion Lagrange-Chevishev #####
def E_pe(t,a,b,m):

        y = 0
        y = abs(fun(t)- punt_equi(t,a,b,m))
        return(y)

##### Errores en la interpolacion Lagrange-Puntos equidistantes #####
def E_ch(t,a,b,m):
        y = 0
        y = abs(fun(t)- punt_cheby(t,a,b,m))
        return(y)


t = np.arange(0,1,0.01)
plt.ylim(-1.5,1.5)
plt.title( "P3 Grafica de Polinomios interpolantes")
plt.plot(t,fun(t),'k',t,punt_equi(t,0,1,20),'r:',t,punt_cheby(t,0,1,20),'b:',p_cheby,f_cheby,'go',)
plt.rc('legend', fontsize='medium')
plt.gca().legend(('f(x) = e^x^2(sin(50x))','P_Equidistantes','Chebyshev','20 puntos chebyshev'))
plt.ylabel('Eje y')
plt.xlabel('Eje x')
plt.show()


t = np.arange(0,1+0.01,0.01)
plt.yscale("log")
plt.title( "P3 Errores en Chebyshev - P_equidistantes")
plt.plot(t,E_pe(t,0,1,20),'k:',t,E_ch(t,0,1,20),'b:') ## 100 puntos equidistantes
plt.rc('legend', fontsize='small')
plt.gca().legend(('Interpolacion P_equidistantes','Interpolacion de Chebyshev',''))
plt.ylabel('Eje y')
plt.xlabel('Eje x')
plt.show()

######## FIN ###################
