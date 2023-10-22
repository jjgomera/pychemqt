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

    3:
        {"autor": "Matsumoto, S., Hara, M., Saito, S., Maeda, S.",
         "title": "Minimum Transport Velocity for Horizontal Pneumatic "
                  "Conveying",
         "ref": "J. Chem. Eng. Japan 7(6) (1974) 425-430",
         "doi": "10.1252/jcej.7.425"},
    4:
        {"autor": "Matsumoto, S., Harada, S., Saito, S., Maeda, S.",
         "title": "Saltation Velocity for Horizontal Pneumatic Conveying",
         "ref": "J. Chem. Eng. Japan 8(4) (1975) 331-333",
         "doi": "10.1252/jcej.8.331"},
    5:
        {"autor": "Matsumoto, S., Kikuta, M., Maeda, S.",
         "title": "Effect of Particle Size on the Minimum Transport Velocity "
                  "for Horizontal Pneumatic Conveying of Solids",
         "ref": "J. Chem. Eng. Japan 10(4) (1977) 273-279",
         "doi": "10.1252/jcej.10.273"},
    6:
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
        \begin{array}[t]{c}
        \frac{M}{\rho_g V_{salt} A} = \frac{1}{10^{\delta}}
        \left(\frac{V_{salt}}{\sqrt{gD}}\right)^{\chi}\\
        \delta = 1440d_p + 1.96\\
        \chi = 1100d_p + 2.5\\
        A = \frac{pi}{4} D^2\\
        \end{array}

    Rearanged saltation velocity

    .. math::
        V{salt} = \left(\frac{M 10^\delta (g D)^{0.5 \chi}}{\rho_g A}\right)^
        {\frac{1}{\chi+1}}

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
    Example 8.1 from [1]_, pag. 238

    >>> "%0.2f" % (Rizk(M=0.25, dp=1e-4, rhog=1.2, D=0.078))
    '9.88'
    """
    # Using dp in mm
    dp *= 1e3

    delta = 1.44*dp + 1.96
    Xi = 1.1*dp + 2.5
    A = pi/4*D**2  # Cross section of pipe

    # Solving saltation velocity
    Vs = (M * 10**delta * (g*D)**(0.5*Xi) / rhog / A)**(1/(Xi+1))
    return Vs


@refDoc(__doi__, [3, 4, 5])
def Matsumoto(M, rhop, dp, rhog, D, Vt, method="1977"):
    r"""Calculates saltation velocity of the gas for pneumatic conveying,
    according to any of method defined in Matsumoto papers.

    * 1974 method, [3]_.

    .. math::
        \begin{array}[t]{c}
        \frac{M}{\rho_g V_{salt} A} = 0.488 \left(\frac{\rho_p}{\rho_f}\right)
        ^{0.5} \left(\frac{Fr_t}{10}\right)^{-1.75}
        \left(\frac{Fr_s}{10}\right)^{3}\\
        Fr_s = \frac{V_{salt}}{\sqrt{g D}}\\
        Fr_t = \frac{V_{t}}{\sqrt{g d_p}}\\
        A = \frac{pi}{4} D^2\\
        \end{array}

    * 1975 method changing only the parameters of equation, [4]_.

    .. math::
        \frac{M}{\rho_g V_{salt} A} = 1.11 \left(\frac{\rho_p}{\rho_f}\right)
        ^{0.55} \left(\frac{Fr_t}{10}\right)^{-2.3}
        \left(\frac{Fr_s}{10}\right)^{3}

    * 1977, [5]_. In this case the correlations have two step formulation
    definning a critical particle diameter

    .. math::
        \frac{d_c}{D} = 1.39\left(\frac{\rho_p}{\rho_g}\right)^{0.74}

    For :math:`d_p < d_c`:

    .. math::
        \frac{M}{\rho_g V_{salt} A} = 5560 \left(\frac{d_p}{D}\right)^{1.43}
        \left(\frac{Fr_s}{10}\right)^{4}

    For :math:`d_p > d_c`:

    .. math::
        \frac{M}{\rho_g V_{salt} A} = 0.373 \left(\frac{\rho_p}{\rho_f}\right)
        ^{1.06} \left(\frac{Fr_t}{10}\right)^{-3.7}
        \left(\frac{Fr_s}{10}\right)^{3.61}

    Parameters
    ----------
    M : float
        Solid mass flow rate, [kg/s]
    rhop : float
        Particle density, [kg/m^3]
    dp : float
        Particle diameter, [m]
    rhog : float
        Gas density, [kg/m^3]
    D : float
        Diameter of pipe, [m]
    Vt : float
        Terminal velocity of particle settling in gas, [m/s]
    method : str
        Specified the method from the year of paper, 1974|1975|1977

    Returns
    -------
    Vs : float
        Saltation velocity of gas, [m/s]

    Examples
    --------
    >>> "%0.2f" % (Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, rhog=1.2, \
    ... D=0.0026, Vt=3.16, method="1974"))
    '18.03'
    >>> "%0.2f" % (Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, rhog=1.2, \
    ... D=0.0026, Vt=3.16, method="1975"))
    '16.49'
    >>> "%0.2f" % (Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, rhog=1.2, \
    ... D=0.0026, Vt=3.16))
    '10.51'
    """
    if method == "1977":
        # Defining critical diameter, Eq 1.
        dc = 1.39*D*(rhop/rhog)**-0.74

        if dp <= dc:
            # Fine particles, Eq 2.
            rhs = 5.56e3*(dp/D)**1.43 * (1/(g*D)**0.5/10)**4
            A = pi/4*D**2
            Vs = (M/rhog/A/rhs)**0.2
            return Vs

        # Coarse particles, Eq 3
        alpha, a, b, c = 0.373, 1.06, -3.7, 3.61

    elif method == "1974":
        alpha, a, b, c = 0.488, 0.5, -1.75, 3
    elif method == "1975":
        alpha, a, b, c = 1.11, 0.55, -2.3, 3

    A = pi/4*D**2
    Frt = Vt/(g*dp)**0.5
    rhs = alpha*(rhop/rhog)**a * (Frt/10.)**b * (1/(g*D)**0.5/10.)**c
    Vs = (M/rhog/A / rhs)**(1/(c+1))
    return Vs


def V_saltation(*args, **kw):
    """List with all saltation velocity correlations availables

    * 0 - Rizk
    * 1 - Matsumoto_1977
    * 2 - Matsumoto_1975
    * 3 - Matsumoto_1974
    """
    method = kw.get("method", 0)
    if method == 0:
        return Rizk(*args, **kw)
    elif method == 1:
        return Matsumoto(*args, method="1977", **kw)
    elif method == 2:
        return Matsumoto(*args, method="1975", **kw)
    elif method == 3:
        return Matsumoto(*args, method="1974", **kw)
