#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Module with neumatic conveying equipment definition
###############################################################################


from scipy.constants import g, pi

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Rhodes, M.",
         "title": "Introduction to Particle Technology 2Ed",
         "ref": "(John Wiley & Sons) 2008",
         "doi": "10.1002/9780470727102"},
    2:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
}


@refDoc(__doi__, [1])
def Rizk(M, dp, rhog, D):
    r"""Calculates saltation velocity of gas for pneumatic conveying, using
    the correlation of Rizk (1973) as described in [1]_

    .. math::
        \frac{M_p}{\rho_g V_{salt} A} = \frac{1}{10^{\delta}}
        \left(\frac{V_{salt}}{\sqrt{gD}}\right)^{\Xi}

        \delta = 1440d_p + 1.96

        \Xi = 1100d_p + 2.5

        A = \frac{pi}{4} D^2

    Rearanged saltation velocity

    .. math::
        V{salt} = \left(M*10**\delta*(g*D)**(0.5*\Xi)/\rho_g/A\right)^
        {\frac{1}{Xi+1}}

    Parameters
    ----------
    M : float
        Solid mass flow rate, [kg/s]
    dp : float
        Particle diameter, [m]
    rhog : float
        Gas density, [kg/m³]
    D : float
        Diameter of pipe, [m]

    Returns
    -------
    Vs : float
        Saltation velocity of gas, [m/s]

    Examples
    --------
    Example 8.1 from [2]_, pag. 238

    >>> "%0.2f" % (Rizk(M=0.25, dp=1e-4, rhog=1.2, D=0.078))
    '9.88'
    """

    # Using dp in mm
    dp *= 1e3

    delta = 1.44*dp + 1.96
    Xi = 1.1*dp + 2.5
    A = pi/4*D**2  # Cross section of pipe

    # Solving saltation velocity
    V = (M * 10**delta * (g*D)**(0.5*Xi) / rhog / A)**(1/(Xi+1))
    return V
