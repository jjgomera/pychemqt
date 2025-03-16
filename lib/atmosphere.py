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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Library for air dispersion modeling
###############################################################################


from math import pi, exp

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "de Visscher, A.",
         "title": "Air Dispersion Modeling: Foundations and Applications",
         "ref": "John Wiley & Sons, 2014",
         "doi": ""},

    2:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
        }


@refDoc(__doi__, [1])
def Gaussian(y, z, h, Q, u, Sy, Sz):
    r"""This function evaluate the Gaussian dispersion model for air dispersion
    of contaminants

    Parameters
    ----------
    y : float
        Crosswind distance from the plume axis [m]
    z : float
        Vertical distance from the surface [m]
    h : float
        effective source height [m]
    Q : float
        Emmision rate [g/s]
    u : float
        Wind speed in the x direction [m/s]
    Sy : float
        Dispersion parameter in the horizontal direction [m]
    Sz : float
        Dispersion parameter in the vertical direction [m]

    Returns
    -------
    c : float
        Concentration of contaminants at a given point [g/m³]

    Examples
    --------
    Example 2.1.

    >>> Sy, Sz = BriggsDispersion(1500, 'D')
    >>> "%0.3e" % Gaussian(0, 0, 90, 100, 7, Sy, Sz)
    '1.603e-04'
    >>> "%0.3e" % Gaussian(100, 0, 90, 100, 7, Sy, Sz)
    '1.075e-04'
    """

    c = Q/2/pi/u/Sy/Sz*exp(-y**2/2/Sy**2)*(
        exp(-(z-h)**2/2/Sz**2)+exp(-(z+h)**2/2/Sz**2))
    return c


@refDoc(__doi__, [1])
def Pasquill_Gifford_Stability(u, day=True, solar=60, cloudiness=0.5):
    r"""This function evaluate the Pasquill-Gifford stability class

    Parameters
    ----------
    u : float
        Wind speed in the x direction measured at 10m height [m/s]
    day : boolean
        Boolean to specify day or night
    solar : float, optional
        Sun inclination above the horizon [º]
        Range between 0 (sunrise) and 90º (mediodía)
        Neccesary only is day is True
    cloudiness : float, optional
        Cloud cover in sky in range 0...1 [-]
        Neccesary only is day is False

    Returns
    -------
    class_ : string
        Pasquill-Gifford stability class name
            A : Very unstable
            B : Moderately unstable
            C : Slightly unstable
            D : Neutral
            E : Stable
    """
    if day:
        if u < 2:
            if solar > 60:
                class_ = 'A'
            elif solar > 35:
                class_ = 'AB'
            else:
                class_ = 'B'
        elif u < 3:
            if solar > 60:
                class_ = 'AB'
            elif solar > 35:
                class_ = 'B'
            else:
                class_ = 'C'
        elif u < 5:
            if solar > 60:
                class_ = 'B'
            elif solar > 35:
                class_ = 'BC'
            else:
                class_ = 'C'
        elif u < 6:
            if solar > 60:
                class_ = 'C'
            elif solar > 35:
                class_ = 'CD'
            else:
                class_ = 'D'
        else:
            if solar > 60:
                class_ = 'C'
            elif solar > 35:
                class_ = 'D'
            else:
                class_ = 'D'
    else:
        if u < 3:
            if cloudiness >= 50:
                class_ = 'E'
            else:
                class_ = 'F'
        elif u < 5:
            if cloudiness >= 50:
                class_ = 'D'
            else:
                class_ = 'E'
        else:
            class_ = 'D'

    return class_


@refDoc(__doi__, [1])
def BriggsDispersion(x, class_, urban=False):
    r"""This function evaluate the dispersion parameters

    Parameters
    ----------
    x : float
        Distance to the source [m]
    class_ : string
        Pasquill-Gifford stability class name
        Accept to mixed classes as AB for intermediate class
    urban : Boolean
        Boolean to specify urban terrain

    Returns
    -------
    Sy : float
        Dispersion parameter in the horizontal direction [m]
    Sz : float
        Dispersion parameter in the vertical direction [m]

    Examples
    --------
    Example 2.1.

    >>> "%0.1f, %0.1f" %BriggsDispersion(1500, 'D')
    '111.9, 49.9'
    """
    Sy = 0
    Sz = 0
    if urban:
        for cls in class_:
            if cls == "A" or cls == "B":
                Sy += 0.32*x/(1+0.0004*x)**0.5
                Sz += 0.24*x/(1+0.0001*x)**0.5
            elif cls == "C":
                Sy += 0.22*x/(1+0.0004*x)**0.5
                Sz += 0.2*x
            elif cls == "D":
                Sy += 0.16*x/(1+0.0004*x)**0.5
                Sz += 0.14*x/(1+0.0003*x)**0.5
            elif cls == "E" or cls == "F":
                Sy += 0.11*x/(1+0.0004*x)**0.5
                Sz += 0.08*x/(1+0.0015*x)**0.5
    else:
        for cls in class_:
            if cls == "A":
                Sy += 0.22*x/(1+0.0001*x)**0.5
                Sz += 0.2*x
            elif cls == "B":
                Sy += 0.16*x/(1+0.0001*x)**0.5
                Sz += 0.12*x
            elif cls == "C":
                Sy += 0.11*x/(1+0.0001*x)**0.5
                Sz += 0.08*x/(1+0.0002*x)**0.5
            elif cls == "D":
                Sy += 0.08*x/(1+0.0001*x)**0.5
                Sz += 0.06*x/(1+0.0015*x)**0.5
            elif cls == "E":
                Sy += 0.06*x/(1+0.0001*x)**0.5
                Sz += 0.03*x/(1+0.0003*x)
            elif cls == "F":
                Sy += 0.04*x/(1+0.0001*x)**0.5
                Sz += 0.016*x/(1+0.0003*x)

    # Sy /= len(class_)
    # Sz /= len(class_)
    return Sy, Sz


if __name__ == "__main__":
    import doctest
    doctest.testmod()
