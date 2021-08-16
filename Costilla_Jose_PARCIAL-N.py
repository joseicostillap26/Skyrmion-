#! /usr/bin/python3
from numpy import*
import numpy as np
import matplotlib.pyplot as plt
import warnings
import matplotlib as mpl
import Costilla_modulo                  # Importo las funciones 
from mpl_toolkits.mplot3d import Axes3D # Modulo para que pueda girar la grafica 3d de la pregunta 4
                                        # No se si pueda leerlo en su pc, en caso contrario borrar este modulo
warnings.filterwarnings("ignore")

####################### Pregunta 1 ####################


print('-------------PREGUNTA 1------------')
print('\n')

pasos = 0
sol = 0

r1 = 0  #Correlaciones
r2 = 0
r3 = 0

def funcion_1(x):       # Derivada del potencial normalizado
    return  -6/x**(7) + np.exp(-x)

def f_d_1(x):           # Derivada de la derivada del potencial normalizado
    return  42/x**(8) - np.exp(-x)

def funcion_12(x):      # Potencial noramlizado = -0.1
    return  1/x**(6) - np.exp(-x) + 0.1

def error_relativo(x2,x1):
    return abs(x1-x2)/x2


def Modele(x,a0,a1):        #### Modelo Err(k) = ax^b
    print('alpha y beta')
    print('a = ',np.exp(a0))
    print('b = ',a1)
    return(np.exp(a0)*x**(a1))



pasos,sol = Costilla_modulo.newton_rhapson(error_relativo,funcion_1,f_d_1,1,1e-6)
print("----El potencial normalizado alcanza su min valor en  : {:.6f}".format(sol))
print("Con el metodo de newton en  {:d} pasos".format(pasos))

Costilla_modulo.secante(error_relativo,funcion_12,1.1,1,1e-6,1) # Imprime la solucion 1
Costilla_modulo.secante(error_relativo,funcion_12,3,2.9,1e-6,2) # Imprime la solucion 2


r1 = Costilla_modulo.R_L_1var(1,Costilla_modulo.inte1,Costilla_modulo.lista_v1) # Imprime las correlaciones
r2 = Costilla_modulo.R_L_1var(1,Costilla_modulo.inte2,Costilla_modulo.lista_v2)
r3 = Costilla_modulo.R_L_1var(1,Costilla_modulo.inte3,Costilla_modulo.lista_v3)

print('Coeficientes de determinacion')
print(r1**2,r2**2,r3**2)
# GRAFICOS

t = np.arange(1,10,0.1)
plt.ylim(-0.1,1.5)
plt.title( "Ajuste lIneal de Err(k) = ax^b cuando V(x) es min")
plt.rc('legend', fontsize='medium')
plt.plot(Costilla_modulo.inte1,Costilla_modulo.lista_v1,'bo',t,Modele(t,Costilla_modulo.a0_a1[0],Costilla_modulo.a0_a1[1]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 4.4668,b = -3,9170, r^2 = 0.5506 ')
plt.xlabel('# de Interacciones')
plt.ylabel('Error relativo')
plt.show()


t = np.arange(1,10,0.1)
plt.ylim(-0.1,1.5)
plt.title( "Ajuste lIneal de Err(k) = bx^a, cuando V(x)=-0.1, sol 1 ")
plt.rc('legend', fontsize='medium')
plt.plot(Costilla_modulo.inte2,Costilla_modulo.lista_v2,'bo',t,Modele(t,Costilla_modulo.a0_a1[2],Costilla_modulo.a0_a1[3]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 2.4556,b = -4,0203, r^2 = 0.7359 ')
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()


t = np.arange(1,10,0.1)
plt.ylim(-0.1,1.5)
plt.title( "Ajuste lIneal de Err(k) = bx^a, cuando V(x)=-0.1, sol 2 ")
plt.rc('legend', fontsize='medium')
plt.plot(Costilla_modulo.inte3,Costilla_modulo.lista_v3,'bo',t,Modele(t,Costilla_modulo.a0_a1[4],Costilla_modulo.a0_a1[5]),'k',)
plt.gca().legend(('Error relativo vs interacciones','Modelo Err(k)= ax^b',''),title = 'a = 3.3583,b = -5,1479, r^2 = 0.7432')
plt.xlabel('# de Interacciones')
plt.ylabel('Error')
plt.show()



###################### Pregunta 2 ##################
print('\n')
print('-------------PREGUNTA 2------------')
print('\n')

G = 6.674*1e-11
M = 5.974*1e24
m = 7.348*1e22
R = 3.844*1e8
w = 2.662*1e-6

pasos = 0
sol = 0
def funcion(x):
    return (G*M/x**2) - (G*m/(R-x)**2) - x*w**2

def f_d(x):
    return  -(2*G*M/x**3) - (2*G*m/(R-x)**3) - w**2

def error_relativo(x2,x1):
    return abs(x1-x2)/x2


pasos,sol = Costilla_modulo.newton_rhapson(error_relativo,funcion,f_d,1,1e-6)

print("Punto de lagrange L1 aproximada-Newton Rhapson- : {:.6f}".format(sol))
print("Con {:d} pasos".format(pasos))


####################################################

print('\n')
print('-------------PREGUNTA 3------------')
print('\n')

###################### Pregunta 3 ##################

m1=15
m2=10
m3=8
m4=5
g = 9.81

A =  np.array([[-1,0,0,-m1],[1,-1,0,-m2],[0,1,-1,-m3],[0,0,1,-m4]])
b = np.array([m1*0.8*g*np.cos(np.pi/4)-m1*g*np.sin(np.pi/4),m2*0.2*g*np.cos(np.pi/4)-m2*g*np.sin(np.pi/4),m3*g,m4*g])

print('Tensiones y aceleracion')
print('[T1,T2,T3,a]=',Costilla_modulo.LU(A,b,4))

###################### Pregunta 4 ##################


def Modele(x1,x2,a0,a1,a2):

    return(a0*x1**(a1)*x2**(a2))

Dia = np.array([0.3,0.6,0.9,0.3,0.6,0.9,0.3,0.6,0.9])
Inc = np.array([0.001,0.001,0.001,0.01,0.01,0.01,0.05,0.05,0.05])
Flu = np.array([0.04,0.24,0.69,0.13,0.82,2.38,0.31,1.95,5.66])

Costilla_modulo.R_L_2var(1,1,Flu,Dia,Inc)

# GRAFICOS

fig = plt.figure()
ax = fig.gca(projection='3d')

X =np.arange(0, 1, 0.001)
Y =  np.arange(0, 0.1, 0.001)
X, Y = np.meshgrid(X, Y)
Z = Modele(X,Y,Costilla_modulo.a0_a1_a2[0],Costilla_modulo.a0_a1_a2[1],Costilla_modulo.a0_a1_a2[2])

surf = ax.plot_surface(X, Y, Z,color = "skyblue")
poi = ax.scatter3D(Dia, Inc,Flu, color = "black");

ax.set_xlabel('Diametro (m)')
ax.set_ylabel('Inclinacion')
ax.set_zlabel('Flujo (m3/s)')
plano = mpl.lines.Line2D([0],[0],linestyle="none", c='skyblue', marker = 's')
punto = mpl.lines.Line2D([0],[0],linestyle="none", c='k', marker = 'o')
ax.legend([plano,punto], ['Modelo Q = a0(D^a1)(S^a2)','( Diametro,Inclinacion,Flujo )'], numpoints = 1,title = 'a0 = 36.3813, a1 = 2.6279, a2 = 0.5319 ')
plt.show()