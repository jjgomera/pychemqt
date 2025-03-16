#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


**Fitting accesories K**

TODO
'''


from math import atan, exp, sqrt, sin, pi

from lib.unidades import Dimensionless
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Crane",
         "title": "Flow of Fluids Through Valves, Fittings, and Pipe",
         "ref": "Crane CO, 1982",
         "doi": ""},

    # 2:
    #     {"autor": "",
    #      "title": "",
    #      "ref": "",
    #      "doi": ""},
}


# Crane, Flow-of-Fluids-Through-Valve Pag 107
@refDoc(__doi__, [1])
def K_contraction_crane(D1, D2, tita=None, L=None):
    r"""Returns loss coefficient for a pipe contraction as shown in [1]_, pag
    A-26

    If Θ < 45º:

    .. math::
        K = 0.8 \sin \frac{\theta}{2}\left(1 - \beta^2\right)

    and for Θ > 45º:

    .. math::
        K = 0.5 \right(1 - \beta^2\left) \sqrt{\sin \frac{\theta}{2}}

    Parameters
    ----------
    D1 : float
        Pipe diameter of entrance of contraction, [m]
    D2 : float
        Pipe diameter of the exit of contraction, [m]
    tita : float, optional
        Angle of contraction, [degrees]
    L : float, optional
        Length of the contraction, [m]

    Notes
    -----
    For gradual contraction we can define with the angle of the contraction Θ
    or with the length of contraction, L.

    Returns
    -------
    K : float
        Loss coefficient, [-]
    """

    # Ratio between both diameters after and before conraction
    beta = D2/D1

    if L is not None:
        if L == 0.0:
            tita = 180
        else:
            tita = 360/pi*2.0*atan((D1-D2)/(2*L))

    if tita < 45:
        # Formula 1
        K = 0.8*sin(tita*pi/360)*(1-beta**2)
    else:
        # Formula 2
        K = 0.5*sqrt(sin(tita*pi/360))*(1.-beta**2)

    return Dimensionless(K)


@refDoc(__doi__, [1])
def K_enlargement_crane(D1, D2, tita=None, L=None):
    r"""Returns loss coefficient for a pipe expansion as shown in [1]_, pag
    A-26

    If Θ < 45º:

    .. math::
        K = 2.6 \sin \frac{\theta}{2}\left(1 - \beta^2\right)

    and for Θ > 45º:

    .. math::
        K = \right(1 - \beta^2\left)^2

    Parameters
    ----------
    D1 : float
        Pipe diameter of entrance of contraction, [m]
    D2 : float
        Pipe diameter of the exit of contraction, [m]
    tita : float, optional
        Angle of contraction, [degrees]
    L : float, optional
        Length of the contraction, [m]

    Notes
    -----
    For gradual contraction we can define with the angle of the contraction Θ
    or with the length of contraction, L.

    Returns
    -------
    K : float
        Loss coefficient, [-]
    """

    # Ratio between both diameters after and before conraction
    beta = D1/D2

    if L is not None:
        if L == 0.0:
            tita = 180
        else:
            tita = 360/pi*2.0*atan((D1-D2)/(2*L))

    if tita < 45:
        # Formula 3
        K = 2.6*sin(tita*pi/360)*(1.-beta**2)**2
    else:
        # Formula 4
        K = (1-beta**2)**2

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
    return -0.31591884532927+29830.477527796*sqrt(2/pi)/122.94894071438 * \
        exp(-2*((tita-183.88854928482)/122.94894071438)**2)


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
    """Función dependiente del diametro que sirve para definir las K de
    valvulas y accesorios de tuberías D debe ser introducido en mm"""
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


if __name__ == "__main__":
    print(K_contraction_crane(3, 2, tita=None, L=0))
    print(K_contraction_crane(3, 2, tita=180))
    print(K_flush(0.10))
