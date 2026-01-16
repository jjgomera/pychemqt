#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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


Module with neumatic conveying equipment library
"""


from scipy.constants import g, pi

from lib.unidades import Dimensionless, Speed
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Rhodes, M.",
         "title": "Introduction to Particle Technology 2Ed",
         "ref": "(John Wiley & Sons) 2008",
         "doi": "10.1002/9780470727102"},
    2:
        {"autor": "Rizk, F.",
         "title": "Pneumatic conveying at optimal operation conditions and a "
                  "solution of Bath's equation",
         "ref": "Proceedings of Pneumotransport 3, paper D4. BHRA Fluid "
                "Engineering, Cranfield, England (1973)",
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
        {"autor": "Ochi, M.",
         "title": "Saltation Velocity of the Gas-Solid Two-Phase Flow in a "
                  "Horizontal Pipe",
         "ref": "Trans. JSME B. 59(564) (1993) 2416-2421",
         "doi": "10.1299/kikaib.59.2416"},
    7:
        {"autor": "Stemerding, S.",
         "title": "The pneumatic transport of cracking catalyst in vertical "
                  "risers",
         "ref": "Chem. Eng. Sci. 17(8) (1962) 599-608",
         "doi": "10.1016/0009-2509(62)80053-0"},
    8:
        {"autor": "Reddy, K.V.S., Pei, D.C.T.",
         "title": "Particle Dynamics in Solids-Gas Flow in a Vertical Pipe",
         "ref": "Ind. Eng. Chem. Fundamen. 8(3) (1969) 490-497",
         "doi": "10.1021/i160031a020"},
    9:
        {"autor": "van Swaaij, W.P.M., Buurman, C., van Breugel, J.W.",
         "title": "Shear Stresses on the Wall of a Dense Gas-Solids Riser",
         "ref": "Chem. Eng. Sci. 25(11) (1970) 1818-1820",
         "doi": "10.1016/0009-2509(70)80072-0"},
    10:
        {"autor": "Capes, C.E., Nakamura, K.",
         "title": "Vertical Pneumatic Conveying: An Experimental Study with "
                  "Particles in the Intermediate and Turbulent Flow Regimes",
         "ref": "Can. J. Chem. Eng. 51(1) (1973) 31-38",
         "doi": "10.1002/cjce.5450510106"},
    11:
        {"autor": "Konno, H., Saito, S.",
         "title": "Pneumatic Conveying of Solids Through Straight Pipes",
         "ref": "J. Chem. Eng. Japan 2(2) (1969) 211-217",
         "doi": "10.1252/jcej.2.211"},
    12:
        {"autor": "Yang, W.-C.",
         "title": "A Correlation for Solid Friction Factor in Vertical "
                  "Pneumatic Conveying Lines",
         "ref": "AIChE J. 24(3) (1978) 548-552",
         "doi": "10.1002/aic.690240326"},
    13:
        {"autor": "Yang, W.-C.",
         "title": "Correlations for Solid Friction Factors in Vertical and "
                  "Horizontal Pneumatic Conveying",
         "ref": "AIChE J. 20(3) (1974) 605-607",
         "doi": "10.1002/aic.690200327"},
    14:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
}


# Saltation speed correlations
@refDoc(__doi__, [1, 2])
def saltation_Rizk(M, dp, rhog, D):
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

    >>> "%0.2f" % (saltation_Rizk(M=0.25, dp=1e-4, rhog=1.2, D=0.078))
    '9.88'
    """
    # Using dp in mm
    dp *= 1e3

    delta = 1.44*dp + 1.96
    Xi = 1.1*dp + 2.5
    A = pi/4*D**2  # Cross section of pipe

    # Solving saltation velocity
    Vs = (M * 10**delta * (g*D)**(0.5*Xi) / rhog / A)**(1/(Xi+1))
    return Speed(Vs)


@refDoc(__doi__, [3, 4, 5])
def saltation_Matsumoto(M, rhop, dp, rhog, D, Vt, method="1977"):
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
    >>> "%0.2f" % (saltation_Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, \
    ... rhog=1.2, D=0.0026, Vt=3.16, method="1974"))
    '18.03'
    >>> "%0.2f" % (saltation_Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, \
    ... rhog=1.2, D=0.0026, Vt=3.16, method="1975"))
    '16.49'
    >>> "%0.2f" % (saltation_Matsumoto(M=0.25, rhop=2500, dp=4.7e-4, \
    ... rhog=1.2, D=0.0026, Vt=3.16))
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
    return Speed(Vs)


@refDoc(__doi__, [6])
def saltation_Ochi(M, dp, rhog, D, Vt, fp=0.4):
    r"""Calculates saltation velocity of the gas for pneumatic conveying,
    according to Ochi correlation, [6]_.

    .. math::
        \begin{array}[t]{c}
        Fr_s = 1.05 f_d^{0.47}Fr_t^{0.82} \mu_s^{0.25}\\
        Fr_t = \frac{V_{t}}{\sqrt{g d_p}}\\
        \mu_s = \frac{M}{\rho_g V_{salt} A}\\
        A = \frac{pi}{4} D^2\\
        \end{array}

    Rearanged saltation velocity

    .. math::
        V_{salt} = \left(\frac{1.05 M^{0.25} f^{0.47} Fr_t^{0.82} \sqrt{g d_p}}
        {\left(\rho_g A\right)^{0.25}}\right)^{0.8}

    Parameters
    ----------
    M : float
        Solid mass flow rate, [kg/s]
    dp : float
        Particle diameter, [m]
    rhog : float
        Gas density, [kg/m^3]
    D : float
        Diameter of pipe, [m]
    Vt : float
        Terminal velocity of particle settling in gas, [m/s]
    fp : float
        Coefficient of friction between particles and surface, [-]

    Returns
    -------
    Vs : float
        Saltation velocity of gas, [m/s]

    Examples
    --------
    >>> "%0.2f" % (saltation_Ochi(M=0.25, dp=4.7e-4, rhog=1.2, D=0.0026, Vt=3.16))
    '8.82'
    """

    A = pi/4*D**2
    Frt = Vt/(g*dp)**0.5
    Vs = (1.05*fp**0.47*Frt**0.82*M**0.25*(g*dp)**0.5/(rhog*A)**0.25)**0.8
    return Speed(Vs)


def V_saltation(method=0, *args):
    """List with all saltation velocity correlations availables

    * 0 - Rizk (1973)
    * 1 - Matsumoto (1977)
    * 2 - Matsumoto (1975)
    * 3 - Matsumoto (1974)
    * 4 - Ochi (1991)
    """

    if method == 0:
        return saltation_Rizk(*args)
    elif method == 1:
        return saltation_Matsumoto(*args, method="1977")
    elif method == 2:
        return saltation_Matsumoto(*args, method="1975")
    elif method == 3:
        return saltation_Matsumoto(*args, method="1974")
    elif method == 4:
        return saltation_Ochi(*args)


# Solid friction factor correlations
@refDoc(__doi__, [7])
def fs_Stemerding():
    """Calculate solid friction factor as define in [7]_.

    This method give a constant value withoud dependent variables based in
    experiment in vertical neumatic conveying

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    # 0.012/4
    return Dimensionless(0.003)


@refDoc(__doi__, [8])
def fs_Reddy(Vs):
    """Calculate solid friction factor as define in [8]_.

    Parameters
    ----------
    Vs : float
        Solid velocity of particle, [m/s]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    return Dimensionless(0.045/Vs)


@refDoc(__doi__, [9])
def fs_Swaaij(Vs):
    """Calculate solid friction factor as define in [9]_.

    Parameters
    ----------
    Vs : float
        Solid velocity of particle, [m/s]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    # Using the graphical slope in Fig. 1
    return Dimensionless(0.08/Vs)


@refDoc(__doi__, [10])
def fs_Capes(Vs):
    """Calculate solid friction factor as define in [10]_.

    Parameters
    ----------
    Vs : float
        Solid velocity of particle, [m/s]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    # Converting parameter to Si units
    # >>> 0.206*k.foot**1.22
    # 0.04834654849325329
    return Dimensionless(0.0483/Vs)


@refDoc(__doi__, [11])
def fs_Konno(Vs, D):
    """Calculate solid friction factor as define in [11]_.

    Parameters
    ----------
    Vs : float
        Solid velocity of particle, [m/s]
    D : float
        Diameter of pipe, [m]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    return Dimensionless(0.0285/Vs*(g*D)**0.5)


@refDoc(__doi__, [12])
def fs_Yang(Vs, eps, Vt, V):
    """Calculate solid friction factor as define in [12]_.

    Parameters
    ----------
    Vs : float
        Solid velocity of particle, [m/s]
    eps : float
        Voidage in transporting line, [-]
    Vt : float
        Terminal velocity of faling particle, [m/s]
    V: float
        Fluid velocity, [m/s]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    # Eq 15
    # Converting factor /4 to fanning friction factor
    return Dimensionless(0.00315*(1-eps)/eps**3*((1-eps)*Vt/(V-Vs))**-0.979)


def f_solid(method=0, *args):
    """List with all solid friction factor correlations availables for vertical
    conveying

    * 0 - Stemerding (1962)
    * 1 - Reddy-Pei (1969)
    * 2 - Swaaij-Buurman-Breugel (1970)
    * 3 - Capes-Nakamura (1973)
    * 4 - Konno-Saito (1969)
    * 5 - Yang (1978)
    """

    if method == 0:
        return fs_Stemerding()
    elif method == 1:
        return fs_Reddy(*args)
    elif method == 2:
        return fs_Swaaij(*args)
    elif method == 3:
        return fs_Capes(*args)
    elif method == 4:
        return fs_Konno(*args)
    else:
        return fs_Yang(*args)


@refDoc(__doi__, [13])
def fs_Yang_Horizontal(eps, V, D):
    """Calculate solid friction factor as define in [13]_.

    Parameters
    ----------
    eps : float
        Voidage in transporting line, [-]
    V: float
        Fluid velocity, [m/s]
    D : float
        Diameter of pipe, [m]

    Returns
    -------
    fs : float
        Solid friction factor, [-]
    """
    # Eq 5
    # Converting factor /4 to fanning friction factor
    return Dimensionless(0.02925*(1-eps)/eps**3*((1-eps)*V/(g*D)**0.5)**-1.15)


