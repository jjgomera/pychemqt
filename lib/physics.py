#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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
#   -Particle solid distributions
#   -Other
#       root3poly
#       Cunninghan factor
###############################################################################

from math import exp, pi, cos, acos

from scipy.constants import R, calorie, liter, atm, Btu, lb
from scipy.special import cbrt


# Constants don't availables in scipy.constants
R_cal = R/calorie
R_atml = R/liter/atm
R_Btu = R/Btu*lb*1000*5/9

factor_acentrico_octano = 0.398


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

    from pylab import arange, plot, grid, show
    x = arange(0, 2.5, 0.01)
    y = weibull(0.5, 0, 1.).cdf(x)
    z = weibull(1, 0, 1).pdf(x)
    w = weibull(1.5, 0, 1).pdf(x)
    v = weibull(5, 0, 1).pdf(x)
    plot(x, y)
    plot(x, z)
    plot(x, w)
    plot(x, v)
    grid(True)
    show()

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
