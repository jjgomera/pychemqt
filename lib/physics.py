#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



###############################################################################
# Module for implement physics correlation
#   -Constants don't available in scipy.constants
#   -Friction factor for pipes
#   -Adimensional numbers
#   -Particle solid distributions
#   -Fitting K
#   -Other
#       root3poly
#       Cunninghan factor
###############################################################################

from math import exp, log, log10, sqrt, sin, pi, cos, acos

from scipy.constants import g, R, calorie, liter, atm, Btu, lb
from scipy.special import cbrt
from scipy.optimize import fsolve

from .unidades import Dimensionless


# Constants don't availables in scipy.constants
R_cal = R/calorie
R_atml = R/liter/atm
R_Btu = R/Btu*lb*1000*5/9

factor_acentrico_octano = 0.398


# Friction factor for pipes
def f_colebrook(Re, eD):
    """
    Friction factor by Colebrook equation

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    This is the best correlation, slowlest to solve
    """
    if eD:
        f = fsolve(lambda x: 1/x**0.5+2.0*log10(eD/3.7+2.51/Re/x**0.5), 0.005)
    else:
        f = fsolve(lambda x: 1/x**0.5-2.0*log10(Re*x**0.5)+0.8, 0.005)
    return f


def f_chen(Re, eD):
    """
    Friction factor by Chen equation (1979)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    f = 1/(-2*log10(eD/3.7065-5.0452/Re*log10(eD**1.1098/2.8257+5.8506/Re**0.8981)))**2
    return f


def f_romeo(Re, eD):
    """
    Friction factor by Romeo equation (2002)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(-2*log10(eD/3.7065-5.0272/Re*log10(eD/3.827-4.567/Re*log10((
        eD/7.7918)**0.9924+(5.3326/(208.815+Re))**0.9345))))**2


def f_goudar(Re, eD):
    """
    Friction factor by Goudar-Sonnad equation (2002)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    a = 2/log(10)
    b = eD/3.7
    d = log(10)*Re/5.2
    s = b*d+log(d)
    q = s**(s/(s+1))
    g = b*d+log(d/q)
    z = log(q/g)
    Dla = g/(g+1)*z
    Dcfa = Dla*(1+z/2/((g+1)**2+z/3*(2*g-1)))

    return (a*(log(d/q)+Dcfa))**-2


def f_manadilli(Re, eD):
    """
    Friction factor by Manadilli (1997). Valid for 5235<Re<1e8

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(-2*log10(eD/3.7+95./Re**0.983-96.82/Re))**2


def f_serghides(Re, eD):
    """
    Friction factor by Serguides

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    A = -2*log10(eD/3.7+12/Re)
    B = -2*log10(eD/3.7+2.51*A/Re)
    C = -2*log10(eD/3.7+2.51*B/Re)
    return (A-(B-A)**2/(C-2*B+A))**-2


def f_churchill(Re, eD):
    """
    Friction factor by Churchill (1977)
    Ref. Darby Chem Eng fluid mechanics pag 180.

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    A = (2.457*log(1/(0.27*eD+(7./Re)**0.9)))**16
    B = (37530./Re)**16
    return 8.*((8./Re)**12+(A+B)**-1.5)**(1./12)


def f_zigrang(Re, eD):
    """
    Friction factor by Zigrang-Sylverster (1982)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    A = log10(eD/3.7-5.02/Re*log10(eD/3.7+13./Re))
    return 1/(-2*log10(eD/3.7-5.02*A/Re))**2


def f_swamee(Re, eD):
    """
    Friction factor by Swamee-Jain (1976)
    for 5000<Re<107 and 0.00004<e/D<0.05

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(1.14-2*log10(eD+21.25/Re**0.9))**2


def f_wood(Re, eD):
    """
    Friction factor by Wood. Valid for Re>4000

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    a = -1.62*eD**0.134
    return 0.094*eD**0.225+0.53*eD+88*eD**0.44*Re**a


def f_moody(Re, eD):
    """
    Friction factor by Moody
    for 4000<Re<107 and e/D <0.01

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 5.5e-3*(1+(2e4*eD+1e6/Re)**(1./3))


def f_haaland(Re, eD):
    """
    Friction factor by Haaland (1983)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(-1.8**2*log10((eD/3.7)**1.11+6.9/Re))**2


def f_barr(Re, eD):
    """
    Friction factor by Barr (1981)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(4*log10(eD/3.7+4.518*log10(Re/7)/Re/(1+Re**0.52/29*eD**0.7)))**2


def f_round(Re, eD):
    """
    Friction factor by Round (1980)

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    """
    return 1/(1.8*log10(eD/0.27+6.5/Re))**2


def f_blasius(Re):
    """
    Friction factor by Blasius, for bare tubes

    Input parameters:
    Re: Reynolds number
    """
    return 0.079/Re**0.25


def f_Gnielinsky(Re):
    """
    Friction factor by Gnielinsky, for annulli section in bare tubes

    Input parameters:
    Re: Reynolds number, based in hydraulic number
    """
    return (1.8*log(Re)-1.5)**-2.


f_list = (f_chen, f_colebrook, f_romeo, f_goudar, f_manadilli, f_serghides,
          f_churchill, f_zigrang, f_swamee)


def f_friccion(Re, eD=0, metodo=0, geometria=0, adicional=0):
    """
    Generalized method for calculate friction factor for laminar or turbulent
    flux in several geometries

    Input parameters:
    Re: Reynolds number
    eD: Ratio between porosity and internal diameter of pipe, e/D
    metodo
        0   -   Colebrook (default)
        1   -   Chen (1979)
        2   -   Romeo (2002)
        3   -   Goudar-Sonnad
        4   -   Manadilli (1997)
        5   -   Serghides
        6   -   Churchill (1977)
        7   -   Zigrang-Sylvester (1982)
        8   -   Swamee-Jain (1976)
    geometria: Index for duct geometry Ref, Darby pag. 215
        0   -   Circular section (default)
        1   -   Square
        2   -   Isosceles triangle
        3   -   rectangular
        4   -   ellipse
        5   -   scalen triangle
        6   -   Anulli
    adicional: other parameter necessary for geometries not circular
        Triangulo isosceles: angulo del vértice superior
        Rectangulo: realación de longitudes de los lados
        Elipse: Array con ambos diametros
        Triángulo rectangular: ángulo del vértice inferior
    """
    if Re < 2100:
        if geometria == 0:
            f_friccion = 16./Re
        elif geometria == 1:
            f_friccion = 14.2/Re
        elif geometria == 2:
            pass
        elif geometria == 3:
            pass
        elif geometria == 4:
            D, d = adicional[1], adicional[0]
            c = (D-d)/(D+d)
            Dh = 4*d*D*(64-16*c**2)/((d+D)*(64-3*c**4))
            f_friccion = 2*Dh**2*(D**2+d**2)/D**2/d**2/Re
        elif geometria == 5:
            pass
        elif geometria == 6:
            f_friccion = 64./Re

    else:
        if geometria == 6:
            f_friccion = f_Gnielinsky(Re)
        else:
            f_friccion = f_list[metodo](Re, eD)

    return Dimensionless(f_friccion)


    """Número de Arquímedes
    Circulación de sólidos en fluidos, ej, fluidización
    """
    return D**3*rho*(rho_p-rho)*g/mu**2




# Particle solids distribution
# Ref Hoffmann,Gas Cyclones and Swirl Tubes, pag 63
def normal(media, varianza):
    from scipy.stats import norm
    return norm(media, varianza)


def lognormal(media, varianza):
    from scipy.stats import lognorm
    return lognorm(media, varianza)


def weibull(escala, start, forma):
    from scipy.stats import weibull_min
    return weibull_min(escala, start, forma)


# Fitting K
# Crane, Flow-of-Fluids-Through-Valve Pag 107
def K_contraction(tita, beta):
    """Tita: angulo de la contracción en grados
    beta: razón entre los diametros antes y despues de la contracción"""
    if tita < 45.:
        K = 0.8*sin(tita*pi/360)*(1-beta**2)/beta**4
    else:
        K = 0.5*sqrt(sin(tita*pi/360))*(1.-beta**2)/beta**4
    return K


def K_enlargement(tita, beta):
    """Tita: angulo del ensanchamiento en grados
    beta: razón entre los diametros antes y despues del ensanchamiento"""
    if tita < 45.:
        K = 2.6*sin(tita*pi/360)*(1.-beta**2)**2/beta**4
    else:
        K = (1.-beta**2)**2/beta**4
    return K


def K_flush(rd):
    """Usando el ajuste exponencial de la tabla disponible"""
    if rd <= 1e-4:
        K = 0.5
    elif rd >= 0.15:
        K = 0.04
    else:
        K = 0.038756579558111+0.45581466480399*exp(-rd/0.041195038092995)
    return K


def K_MitreBend(tita):
    """Usando el ajuste gausiano de la tabla disponible"""
    return -0.31591884532927+29830.477527796*sqrt(2/pi)/122.94894071438*exp(-2*((tita-183.88854928482)/122.94894071438)**2)


def K_longBend(rD):
    if rD <= 1.2:
        K = 20
    elif rD <= 1.7:
        K = 14
    elif rD <= 3.5:
        K = 12
    elif rD <= 5:
        K = 14
    elif rD <= 7:
        K = 17
    elif rD <= 9:
        K = 24
    elif rD <= 11:
        K = 30
    elif rD <= 13:
        K = 34
    elif rD <= 15:
        K = 38
    elif rD <= 18:
        K = 42
    else:
        K = 50
    return K


def Ft(D):
    """Función dependiente del diametro que sirve para definir las K de valvulas y accesorios de tuberías
    D debe ser introducido en mm"""
    if D <= 15:
        ft = 0.027
    elif D <= 20:
        ft = 0.025
    elif D <= 25:
        ft = 0.023
    elif D <= 32:
        ft = 0.022
    elif D <= 40:
        ft = 0.021
    elif D <= 50:
        ft = 0.019
    elif D <= 80:
        ft = 0.018
    elif D <= 100:
        ft = 0.017
    elif D <= 125:
        ft = 0.016
    elif D <= 150:
        ft = 0.015
    elif D <= 250:
        ft = 0.014
    elif D <= 400:
        ft = 0.013
    else:
        ft = 0.012
    return ft


# Mathematical special functions
def root_poly3(a1, a2, a3):
    """Roots for a cubic polinomy, x^3 + a1*x^2 + a2*x + a3"""
    Q = (3*a2-a1**2)/9.
    L = (9*a1*a2-27*a3-2*a1**3)/54.
    D = Q**3+L**2
    if D < 0:
        tita = acos(L/(-Q**3)**0.5)
        z1 = 2*(-Q)**0.5*cos(tita/3.+2.*pi/3)-a1/3.
        z2 = 2*(-Q)**0.5*cos(tita/3.+4.*pi/3)-a1/3.
        z3 = 2*(-Q)**0.5*cos(tita/3.)-a1/3.
        z = [z1, z2, z3]
    else:
        S1 = cbrt((L+D**0.5))
        S2 = cbrt((L-D**0.5))
        if D > 0:
            z = [S1+S2-a1/3.]
        else:
            z1 = cbrt(L)+cbrt(L)-a1/3.
            z2 = -cbrt(L)-a1/3.
            z = [z1, z2]
#    print z[0]+z[1]+z[2], -a1
#    print z[0]*z[1]+z[1]*z[2]+z[2]*z[0], a2
#    print z[0]*z[1]*z[2], -a3
    z.sort()
    return z


# Other
def Cunningham(l, Kn, method=0):
    """Cunningham slip correction factor for air
        l: Mean free path
        kn: Knudsen dimensionless number
        method: reference procedure
            0 - Jennings (1987)
            1 - Allen & Raabe (1982)
            2 - Fuchs (1964)
            3 - Davies (1945)
    """
    if method == 3:
        C = 1+Kn*(1.257+0.4*exp(-1.1/Kn))
    elif method == 2:
        C = 1+Kn*(1.246+0.418*exp(-0.867/Kn))
    elif method == 1:
        C = 1+Kn*(1.155+0.471*exp(-0.596/Kn))
    else:
        C = 1+Kn*(1.252+0.399*exp(-1.1/Kn))
    return C


if __name__ == "__main__":
    print(f_colebrook(1e7, 0.0002))
    print(f_chen(1e7, 0.0002))
    print(f_goudar(1e7, 0.0002))
    print(f_romeo(1e7, 0.0002))
    print(f_manadilli(1e7, 0.0002))
    print(f_serghides(1e7, 0.0002))
    print(f_churchill(1e7, 0.0002))
    print(f_zigrang(1e7, 0.0002))
    print(f_swamee(1e7, 0.0002))
    print(f_wood(1e7, 0.0002))
    print(f_haaland(1e7, 0.0002))
    print(f_barr(1e7, 0.0002))
    print(f_round(1e7, 0.0002))
    print(f_moody(1e7, 0.0002))

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
