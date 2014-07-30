#!/usr/bin/python
# -*- coding: utf-8 -*-

###Módulo en el que encuentran las funciones de cálculo de ecuaciones complejas y que se utilizarán en módulos del programa

from math import exp, log, log10, sqrt, sin, pi, cos, acos

from scipy.constants import g, R, calorie, liter, atm, Btu, lb
from scipy.special import cbrt
from scipy.optimize import fsolve

from unidades import Dimensionless
factor_acentrico_octano=0.398


#######################################################################
###                                              Constantes no disponibles en scipy.constants                                                 ###
#######################################################################

R_cal=R/calorie
R_atml=R/liter/atm
R_Btu=R/Btu*lb*1000*5/9

#######################################################################
###                                     Métodos de cálculo del factor de fricción de darcy en tuberías                                 ###
#######################################################################


def f_colebrook(Re, eD):
    """
    Factor de fricción según la ecuación de Colebrook
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    Se considera la relación exacta, aunque al no ser implicita es más lenta de resolver
    """
    if eD:
        f=fsolve(lambda x: 1/x**0.5+2.0*log10(eD/3.7+2.51/Re/x**0.5), 0.005)
    else:
        f=fsolve(lambda x: 1/x**0.5-2.0*log10(Re*x**0.5)+0.8, 0.005)
    return f


def f_chen(Re, eD):
    """
    Factor de fricción según la ecuación de Chen (1979)
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    return 1/(-2*log10(eD/3.7065-5.0452/Re*log10(eD**1.1098/2.8257+5.8506/Re**0.8981)))**2


def f_romeo(Re, eD):
    """
    Factor de fricción según la ecuación de Romeo (2002)
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    return 1/(-2*log10(eD/3.7065-5.0272/Re*log10(eD/3.827-4.567/Re*log10((eD/7.7918)**0.9924+(5.3326/(208.815+Re))**0.9345))))**2


def f_goudar(Re, eD):
    """
    Factor de fricción según la ecuación de Goudar-Sonnad.
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    a=2/log(10)
    b=eD/3.7
    d=log(10)*Re/5.2
    s=b*d+log(d)
    q=s**(s/(s+1))
    g=b*d+log(d/q)
    z=log(q/g)
    Dla=g/(g+1)*z
    Dcfa=Dla*(1+z/2/((g+1)**2+z/3*(2*g-1)))

    return (a*(log(d/q)+Dcfa))**-2


def f_manadilli(Re, eD):
    """
    Factor de fricción según la ecuación de Manadilli (1997). Válida en el rango de Re de 5235-1e8
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    
    return 1/(-2*log10(eD/3.7+95./Re**0.983-96.82/Re))**2


def f_serghides(Re, eD):
    """
    Factor de fricción según la ecuación de Serghides.
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    
    A=-2*log10(eD/3.7+12/Re)
    B=-2*log10(eD/3.7+2.51*A/Re)
    C=-2*log10(eD/3.7+2.51*B/Re)
    return (A-(B-A)**2/(C-2*B+A))**-2


def f_churchill(Re, eD):
    """
    Factor de fricción según la ecuación de Churchill (1977)
    Ref. Darby Chem Eng fluid mechanics pag 180.
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """

    A=(2.457*log(1/(0.27*eD+(7./Re)**0.9)))**16
    B=(37530./Re)**16
    return 8.*((8./Re)**12+(A+B)**-1.5)**(1./12)


def f_zigrang(Re, eD):
    """
    Factor de fricción según la ecuación de Zigrang-Sylvester (1982)
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    A=log10(eD/3.7-5.02/Re*log10(eD/3.7+13./Re))
    return 1/(-2*log10(eD/3.7-5.02*A/Re))**2


def f_swamee(Re, eD):
    """
    Factor de fricción según la ecuación de Swamee y Jain (1976).
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    (for 5000<Re<107 and 0.00004<e/D<0.05)
    """
    return 1/(1.14-2*log10(eD+21.25/Re**0.9))**2


#Moody Equation (4000<Re<107 and e/D <0.01)
#f = 5.5x10-3(1+ (2x104e/D + 106/Re)1/3)


def f_wood(Re, eD):
    """
    Factor de fricción según la ecuación de Wood.
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    (Re>4000)
    """
    a=-1.62*eD**0.134
    return 0.094*eD**0.225+0.53*eD+88*eD**0.44*Re**a


def f_moody(Re, eD):
    """
    Factor de fricción según la ecuación de Wood.
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    (4000<Re<107 and e/D <0.01)
    """
    return 5.5e-3*(1+(2e4*eD+1e6/Re)**(1./3))
    

def f_haaland(Re, eD):
    """
    Factor de fricción según la ecuación de Haaland (1983).
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """
    
    return 1/(-1.8**2*log10((eD/3.7)**1.11+6.9/Re))**2
    
    
    

    
    




def f_barr(Re, eD):
    """
    Factor de fricción según la ecuación de Barr (1981)
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """

    return 1/(4*log10(eD/3.7+4.518*log10(Re/7)/Re/(1+Re**0.52/29*eD**0.7)))**2


def f_round(Re, eD):
    """
    Factor de fricción según la ecuación de Round (1980)
    
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    """

    return 1/(1.8*log10(eD/0.27+6.5/Re))**2




def f_blasius(Re):
    """
    Factor de fricción según la ecuación Blasius, en tubos lisos
    
    Parámetros necesarios:
    Re: Número de reynolds
    """
    return 0.079/Re**0.25


def f_Gnielinsky(Re):
    """
    Factor de fricción según la ecuación de Gnielinsky, para la sección anular en tubos concentricos lisos
    
    Parámetros necesarios:
    Re: Número de reynolds, basado en el radio hidráulico
    """
    return (1.8*log(Re)-1.5)**-2.

f_list=f_chen, f_colebrook, f_romeo, f_goudar, f_manadilli, f_serghides, f_churchill, f_zigrang, f_swamee


def f_friccion(Re, eD=0, metodo=0, geometria=0, adicional=0):
    """
    Método generalizado de cálculo del factor de fricción para flujo laminar o turbulento y diferentes geometrias
    Parametros necesarios:
    Re: Número de reynolds
    eD: Relacion entre porosidad y diametro interno de la tubería, e/D
    método de cálculo:
        0   -   Colebrook
        1   -   Chen (1979)
        2   -   Romeo (2002)
        3   -   Goudar-Sonnad
        4   -   Manadilli (1997)
        5   -   Serghides
        6   -   Churchill (1977)
        7   -   Zigrang-Sylvester (1982)
        8   -   Swamee-Jain (1976)

    geometria: indice que indica la geometria de la conducción Ref, Darby pag. 215
        0   -   Valor por defecto, conducción cilíndrica
        1   -   Cuadrada
        2   -   triangulo isosceles
        3   -   rectangular
        4   -   elipse
        5   -   triangulo escaleno
        6   -   Anillo
    adicional: parámetro adicional necesario para definir las geometrias no circulares
        Triangulo isosceles: angulo del vértice superior
        Rectangulo: realación de longitudes de los lados
        Elipse: Array con ambos diametros
        Triángulo rectangular: ángulo del vértice inferior
    """
    if Re<2100:
        if geometria==0:
            f_friccion=16./Re
        elif geometria==1:
            f_friccion=14.2/Re
        elif geometria==2:
            pass
        elif geometria==3:
            pass
        elif geometria==4:
            D, d=adicional[1], adicional[0]
            c=(D-d)/(D+d)
            Dh=4*d*D*(64-16*c**2)/((d+D)*(64-3*c**4))
            f_friccion=2*Dh**2*(D**2+d**2)/D**2/d**2/Re
        elif geometria==5:
            pass
        elif geometria==6:
            f_friccion=64./Re
            
    else:
        if geometria==6:
            f_friccion=f_Gnielinsky(Re)
        else:
            f_friccion=f_list[metodo](Re, eD)

    return Dimensionless(f_friccion)




#######################################################################
###                                                               Parámetros adimensionales                                                           ###
#######################################################################


def Re(D, V, rho, mu):
    """
    Número de Reynolds
    Fluidos
    """
    return Dimensionless(D*V*rho/mu)

def We(D, V, rho, sigma):
    """
    Número de Weber
    Formación de gotas
    """
    return D*V**2*rho/sigma
    
def Pr(Cp, mu, k):
    """
    Número de Prandtl
    Transferencia de calor
    """
    return Cp*mu/k
    
    
def Nu(h, D, k):
    """
    Número de Nusselt
    Transferencia de calor
    """
    return h*D/k


def Ar(D, rho_p, rho, mu):
    """Número de Arquímedes
    Circulación de sólidos en fluidos, ej, fluidización
    """
    return D**3*rho*(rho_p-rho)*g/mu**2
    

def Gz(w, cp, k, L):
    """Número de Graetz
    Flujo laminar en conducciones
    w: caudal másico"""
    return w*cp/k/L
    

def Gr(beta, T, To, L, mu):
    """Número de Grashof
    """
    return g*beta*(T-To)*L**3/mu**2
    
    
    

#########################################################################
####                                                Distribuciones de partículas                                                                         ####
#########################################################################

"""Ref Hoffmann,Gas Cyclones and Swirl Tubes, pag 63"""

def normal(media, varianza):
    from scipy.stats import norm
    return norm(media, varianza)

def lognormal(media, varianza):
    from scipy.stats import lognorm
    return lognorm(media, varianza)

def weibull(escala, start, forma):
    from scipy.stats import weibull_min
    return weibull_min(escala, start, forma)


#########################################################################
####                                                Distribuciones de partículas                                                                         ####
#########################################################################

"""Crane, Flow-of-Fluids-Through-Valve Pag 107"""

def K_contraction(tita, beta):
    """Tita: angulo de la contracción en grados
    beta: razón entre los diametros antes y despues de la contracción"""
    if tita < 45.:
        K=0.8*sin(tita*pi/360)*(1-beta**2)/beta**4
    else:
        K=0.5*sqrt(sin(tita*pi/360))*(1.-beta**2)/beta**4
    return K
        
def K_enlargement(tita, beta):
    """Tita: angulo del ensanchamiento en grados
    beta: razón entre los diametros antes y despues del ensanchamiento"""
    if tita < 45.:
        K=2.6*sin(tita*pi/360)*(1.-beta**2)**2/beta**4
    else:
        K=(1.-beta**2)**2/beta**4
    return K

def K_flush(rd):
    """Usando el ajuste exponencial de la tabla disponible"""
    if rd<=1e-4:
        K=0.5
    elif rd>=0.15:
        K=0.04
    else:
        K=0.038756579558111+0.45581466480399*exp(-rd/0.041195038092995)
    return K
    
def K_MitreBend(tita):
    """Usando el ajuste gausiano de la tabla disponible"""
    return -0.31591884532927+29830.477527796*sqrt(2/pi)/122.94894071438*exp(-2*((tita-183.88854928482)/122.94894071438)**2)
    
def K_longBend(rD):
    if rD<=1.2:
        K=20
    elif rD<=1.7:
        K=14
    elif rD<=3.5:
        K=12
    elif rD<=5:
        K=14
    elif rD<=7:
        K=17
    elif rD<=9:
        K=24
    elif rD<=11:
        K=30
    elif rD<=13:
        K=34
    elif rD<=15:
        K=38
    elif rD<=18:
        K=42
    else:
        K=50
    return K
def Ft(D):
    """Función dependiente del diametro que sirve para definir las K de valvulas y accesorios de tuberías
    D debe ser introducido en mm"""
    if D<=15:
        ft=0.027
    elif D<=20:
        ft=0.025
    elif D<=25:
        ft=0.023
    elif D<=32:
        ft=0.022
    elif D<=40:
        ft=0.021
    elif D<=50:
        ft=0.019
    elif D<=80:
        ft=0.018
    elif D<=100:
        ft=0.017
    elif D<=125:
        ft=0.016
    elif D<=150:
        ft=0.015
    elif D<=250:
        ft=0.014
    elif D<=400:
        ft=0.013
    else:
        ft=0.012
    return ft


#########################################################
###                                       Funciones matemáticas especiales                                     ###
#########################################################

def root_poly3(a1, a2, a3):
    """Raices de un polinomio de grado tres normalizada, coeficiente 1 para el término de potencia 3"""
    Q=(3*a2-a1**2)/9.
    L=(9*a1*a2-27*a3-2*a1**3)/54.
    D=Q**3+L**2
    if D<0:
        tita=acos(L/(-Q**3)**0.5)
        z1=2*(-Q)**0.5*cos(tita/3.+2.*pi/3)-a1/3.
        z2=2*(-Q)**0.5*cos(tita/3.+4.*pi/3)-a1/3.
        z3=2*(-Q)**0.5*cos(tita/3.)-a1/3.
        z=[z1, z2, z3]
    else:
        S1=cbrt((L+D**0.5))
        S2=cbrt((L-D**0.5))
        if D>0:
            z=[S1+S2-a1/3.]
        else:
            z1=cbrt(L)+cbrt(L)-a1/3.
            z2=-cbrt(L)-a1/3.
            z=[z1, z2]
#    print z[0]+z[1]+z[2], -a1
#    print z[0]*z[1]+z[1]*z[2]+z[2]*z[0], a2
#    print z[0]*z[1]*z[2], -a3
    z.sort()
    return z
    
    


if __name__ == "__main__":
    print f_colebrook(1e7, 0.0002)
    print f_chen(1e7, 0.0002)
    print f_goudar(1e7, 0.0002)
    print f_romeo(1e7, 0.0002)
    print f_manadilli(1e7, 0.0002)
    print f_serghides(1e7, 0.0002)
    print f_churchill(1e7, 0.0002)
    print f_zigrang(1e7, 0.0002)
    print f_swamee(1e7, 0.0002)
    print f_wood(1e7, 0.0002)
    print f_haaland(1e7, 0.0002)
    print f_barr(1e7, 0.0002)
    print f_round(1e7, 0.0002)
    print f_moody(1e7, 0.0002)
    
#    from pylab import arange, plot,  grid, show
#    x=arange(0, 2.5, 0.1)
#    y= weibull(0.5, 0, 1.).cdf(x)
#    z=weibull(1, 0, 1).pdf(x)
#    w=weibull(1.5, 0, 1).pdf(x)
#    v=weibull(5, 0, 1).pdf(x)
#    plot(x, y)
#    grid(True)
#    show()

#    print K_MitreBend(60)

#    import time
#    from scipy import roots
#    
#    def scipy(a, b, c):
#        z=roots([1, a, b, c])
#        return z
#        
#    inicio=time.time()
#    print scipy(-1., 0.29664, -0.026584)
#    fin=time.time()
#    print "scipy.roots:", (fin-inicio)*1000., "ms"
#
#    inicio=time.time()
#    print root_poly3(-1., 0.29664, -0.026584)
#    fin=time.time()
#    print "alg.roots:", (fin-inicio)*1000., "ms"
