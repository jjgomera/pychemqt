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
along with this program.  If not, see <http://www.gnu.org/licenses/>.



Module for implement physics correlation:

  * Constants don't available in scipy.constants
  * Particle solid distributions
  * Other
      * root3poly
      * Cunninghan factor
      * :func:`Collision_Neufeld`: Neufeld Collision integral

'''


from math import exp, sin, cos, acos

from scipy.constants import R, calorie, liter, atm, Btu, lb

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Neufeld, P.D., Janzen, A.R., Aziz, R.A.",
         "title": "Empirical Equations to Calculate 16 of the Transport "
                  "Collision Integrals Ω for the Lennard-Jones Potential",
         "ref": "J. Chem. Phys. 57(3) (1972) 1100-1102",
         "doi": "10.1063/1.1678363"},


    2:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
        }

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
def cubicCardano(a, b, c, d):
    """
    Implementation of Cardano formula to find all roots in a cubic equation
    faster than with scipy.roots

    .. math::
       ax^3 + bx^2 + cx + d = 0

    Algorithms referenced at http://www.1728.org/cubic2.htm

    Parameters
    ----------
    a : float
        Third grade coefficient, [-]
    b : float
        Second grade coefficient, [-]
    c : float
        First grade coefficient, [-]
    d : float
        Zero grade coefficient, [-]

    Returns
    -------
    root : list
        Solution list, [-]

    Examples
    --------
    # (x+2)^3

    >>> "%i %i %i" % (cubicCardano(1, 6, 12, 8))
    '2 2 2'

    # (x-2)^3

    >>> "%i %i %i" % (cubicCardano(1, -6, 12, -8))
    '-2 -2 -2'

    >>> "%.1f %.1f %.1f" % (cubicCardano(2, -4, -22, 24))
    '-3.0 1.0 4.0'

    >>> "{:.4f} {:.4f} {:.4f}".format(*cubicCardano(3, -10, 14, 27))
    '-1.0000 2.1667+2.0750j 2.1667-2.0750j'
    """

    f = ((3*c/a) - (b**2/a**2))/3
    g = ((2*b**3/a**3) - (9*b*c/a**2) + (27*d/a))/27
    h = g**2/4 + f**3/27

    if f == 0 and g == 0 and h == 0:
        # All 3 roots are real and equal
        x1 = abs((d/a))**(1/3) * (1, -1)[d/a < 0]
        return x1, x1, x1

    elif h <= 0:
        # All 3 roots are real

        i = (g**2/4 - h)**0.5
        j = i**(1/3)
        k = acos(-g/(2*i))
        L = -j
        M = cos(k/3)
        N = 3**0.5*sin(k/3)
        P = -b/(3*a)

        x1 = 2*j*cos(k/3) - b/(3*a)
        x2 = L*(M+N) + P
        x3 = L*(M-N) + P

        return tuple(sorted((x1, x2, x3)))

    elif h > 0:
        # One real root and two conjugate complex roots

        R = -g/2 + h**0.5
        S = abs(R)**(1/3) * (1, -1)[R < 0]
        T = -g/2 - h**0.5
        U = abs(T)**(1/3) * (1, -1)[T < 0]

        x1 = S + U - b/(3*a)
        x2 = -(S+U)/2 - b/(3*a) + 1j*(S-U)*3**0.5/2
        x3 = -(S+U)/2 - b/(3*a) - 1j*(S-U)*3**0.5/2

        return x1, x2, x3


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


@refDoc(__doi__, [1])
def Collision_Neufeld(T, l=2, s=2):
    r"""Calculate the collision integral using the Neufeld correlation

    .. math::
        \varOmega^{(l,s)}=A/T^{B}+C/\exp\left(DT\right)+E/\exp\left(FT\right)+
            G/\exp\left(HT\right)+RT^{B}\sin\left(ST^{w}-P\right)

    A,B,C,D,E,F,G,H,R,S,W,P are constants for each collison order

    Parameters
    ----------
    T : float
        Reduced temperature, [-]
    l: int, optional
        Collision integral first term order, default 2
    s : int, optional
        Collision integral second term order, default 2

    Returns
    -------
    omega : float
        Transport collision integral, [-]
    """
    # Table I
    dat = {
        (1, 1): (1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996, 1.76474,
                 3.89411, 0, 0, 0, 0),
        (1, 2): (1.00220, 0.15530, 0.16105, 0.72751, 0.86125, 2.06848, 1.95162,
                 4.84492, 0, 0, 0, 0),
        (1, 3): (0.96573, 0.15611, 0.44067, 1.52420, 2.38981, 5.08063, 0, 0,
                 -5.373, 19.2866, -1.30775, 6.58711),
        (1, 4): (0.93447, 0.15578, 0.39478, 1.85761, 2.45988, 6.15727, 0, 0,
                 4.246, 12.9880, -1.36399, 3.33290),
        (1, 5): (0.90972, 0.15565, 0.35967, 2.18528, 2.45169, 7.17936, 0, 0,
                 -3.814, 9.38191, 0.14025, 9.93802),
        (1, 6): (0.88928, 0.15562, 0.33305, 2.51303, 2.36298, 8.11690, 0, 0,
                 -4.649, 9.86928, 0.12851, 9.82414),
        (1, 7): (0.87208, 0.15568, 0.36583, 3.01399, 2.70659, 9.92310, 0, 0,
                 -4.902, 10.2274, 0.12306, 9.97712),
        (2, 2): (1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787, 0, 0,
                 -6.435, 18.0323, -0.76830, 7.27371),
        (2, 3): (1.11521, 0.14796, 0.44844, 0.99548, 2.30009, 3.06031, 0, 0,
                 4.565, 38.5868, -0.69403, 2.56375),
        (2, 4): (1.08228, 0.14807, 0.47128, 1.31596, 2.42738, 3.90018, 0, 0,
                 -5.623, 3.08449, 0.28271, 3.22871),
        (2, 5): (1.05581, 0.14822, 0.51203, 1.67007, 2.57317, 4.85939, 0, 0,
                 -7.120, 4.71210, 0.21730, 4.73530),
        (2, 6): (1.03358, 0.14834, 0.53928, 2.01942, 2.72350, 5.84817, 0, 0,
                 -8.576, 7.66012, 0.15493, 7.60110),
        (3, 3): (1.05567, 0.14980, 0.30887, 0.86437, 1.35766, 2.44123, 1.29030,
                 5.55734, 2.339, 57.7757, -1.08980, 6.94750),
        (3, 4): (1.02621, 0.15050, 0.55381, 1.40070, 2.06176, 4.26234, 0, 0,
                 5.227, 11.3331, -0.82090, 3.87185),
        (3, 5): (0.99958, 0.15029, 0.50441, 1.64304, 2.06947, 4.87712, 0, 0,
                 -5.184, 3.45031, 0.26821, 3.73348),
        (4, 4): (1.12007, 0.14578, 0.53347, 1.11986, 2.28803, 3.27567, 0, 0,
                 7.427, 21.0480, -0.28759, 6.69149)}

    A, B, C, D, E, F, G, H, R, S, W, P = dat[(l, s)]

    # Eq 2
    omega = A/T**B + C/exp(D*T) + E/exp(F*T) + G/exp(G*T) + \
        R*1e-4*T**B*sin(S*T**W-P)
    return omega


if __name__ == "__main__":

    # from pylab import arange, plot, grid, show
    # x = arange(0, 2.5, 0.01)
    # y = weibull(0.5, 0, 1.).cdf(x)
    # z = weibull(1, 0, 1).pdf(x)
    # w = weibull(1.5, 0, 1).pdf(x)
    # v = weibull(5, 0, 1).pdf(x)
    # plot(x, y)
    # plot(x, z)
    # plot(x, w)
    # plot(x, v)
    # grid(True)
    # show()

    import time
    from pprint import pprint
    from scipy import roots

    args = [1, -1.0, 0.14358582833553338, -0.004243577526503186]

    # # args = [1, -1, 0.6269002858464162, -0.07953977078483085]

    inicio=time.time()
    pprint(roots(args))
    fin=time.time()
    print("scipy.roots:", (fin-inicio)*1000., "ms")

    inicio=time.time()
    pprint(cubicCardano(*args))
    fin=time.time()
    print("alg.roots:", (fin-inicio)*1000., "ms")





