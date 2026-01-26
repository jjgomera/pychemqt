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


from math import pi

from lib.unidades import Dimensionless, Area, Length
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "du Plessis, J.P., Kröger, D.G.",
         "title": "Friction factor prediction for fully developed laminar "
                  "twisted-tape flow",
         "ref": "Int. J. Heat Mass Transfer 27(11) (1984) 2095-2100",
         "doi": "10.1016/0017-9310(84)90196-0"},
    2:
        {"autor": "du Plessis, J.P., Kröger, D.G.",
         "title": "Heat transfer correlation for thermally developing laminar "
                  "flow in a smooth tube with a twisted-tape insert",
         "ref": "Int. J. Heat Mass Transfer 30(3) (1987) 509-515",
         "doi": "10.1016/0017-9310(87)90265-1"},
    3:
        {"autor": "Shah, R.K., London, A.L.",
         "title": "Laminar Flow Forced Convection in Ducts: A Source Book for "
                  "Compact Heat Exchanger Analytical Data",
         "ref": "Academic Press 1978",
         "doi": ""},
    4:
        {"autor": "Manglik, R.M., Bergles, A.E.",
         "title": "Heat Transfer and Pressure Drop Correlations for "
                  "Twisted-Tape Inserts in Isothermal Tubes: Part I - Laminar "
                  "Flows",
         "ref": "J. Heat Transfer 115(4) (1993) 881-889",
         "doi": "10.1115/1.2911383"},
    5:
        {"autor": "Manglik, R.M., Bergles, A.E.",
         "title": "Heat Transfer and Pressure Drop Correlations for "
                  "Twisted-Tape Inserts in Isothermal Tubes: Part II - "
                  "Transition and Turbulent Flows",
         "ref": "J. Heat Transfer 115(4) (1993) 890-896",
         "doi": "10.1115/1.2911384"},

        }


class TwistedTape():
    """Twisted-tape insert used in heat exchanger to improve efficiency.
    This tape, generally a thin metal strip, is twisted about its longitudinal
    axis"""

    def __init__(self, H, D, delta):
        """
        Definition of twisted tape accesory

        Parameters
        ----------
        H : float
            Tape pitch for twist of π radians (180º), [m]
        D : float
            Internal diameter of tube, [m]
        delta : float
            Tape thickness, [m]
        """
        # Geometrical definition of parameters in [1]_

        # Helical factor, Eq 2
        self.G = Dimensionless((1+pi**2/D**2/4/H**2)**0.5)

        # Effective cross-sectional flow area, Eq 3
        self.Ae = Area(2*H**2/pi*(self.G-1) - D*delta)

        # Effective wetted perimeter, Eq 4
        self.Pe = Length(2 * (D - delta + pi*D/2/self.G))

        # Effective hydraulic diameter, Eq 7
        self.De = Length(4*self.Ae/self.Pe)

        # Area tube without tape
        self.A = pi*D**2/4

        # Tape twist parameter
        self.y = Dimensionless(H/D)


# Friction factor correlations
@refDoc(__doi__, [1])
def f_twisted_Plessis(Re, D, H, delta, Ae, De):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Plessis and Kröger correlation

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    delta : float
        Tape thickness, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    Ae : float
        Effective flow area, [m²]
    De : float
        Effective hydraulic diameter, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    Examples
    --------
    Selected point from Table 2 in [1]_

    >>> st = TwistedTape(10, 1, 0)
    >>> print("%0.3f" % f_twisted_Plessis(50, 1, 10, 0, st.Ae, st.De))
    0.849

    >>> print("%0.4f" % f_twisted_Plessis(2000, 1, 10, 0, st.Ae, st.De))
    0.0296
    """
    A = pi*D**2/4
    y = H/D

    # Eq 15
    Deltae = A*D**2/Ae/De**2

    # Eq 16
    fe = Deltae/Re*(15.767-0.14706*delta/D)

    # Eq 17
    f = fe * (1 + (Re/(70*y**1.3))**1.5)**(1/3)
    return f


@refDoc(__doi__, [3])
def f_twisted_Shah(Re, D, H, delta):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Plessis and Kröger correlation

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]

    Returns
    -------
    f : float
        Friction factor, [-]

    """
    # Chapter XVI: Longitudinal Fins and Twisted Tapes within Ducts
    # F: Circular Duct with a Twisted Tape, Pag 379 and on next

    Xl = H/D

    # Eq 559
    C = 8.8201*Xl - 2.1193*Xl**2 + 0.2108*Xl**3 - 0.0069*Xl**4

    # Eq 563
    Xi = (pi/(pi+2))**2 * ((pi+2-2*delta/D)/(pi-4*delta/D))**2 \
        * (pi/(pi-4*delta/D))

    if Re/Xl < 6.7:
        # Eq 560
        fRe = 42.23*Xi
    elif Re/Xl > 100:
        # Eq 562
        fRe = C*(Re/Xl)**0.3*Xi
    else:
        # Eq 561
        fRe = 38.4*(Re/Xl)**0.05*Xi

    return fRe/Re


@refDoc(__doi__, [4, 5])
def f_twisted_Manglik(Re, D, H, delta, Dh):
    """Calculate friction factor for a pipe with a twisted-tape insert using
    the Manglik and Bergles correlation

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal pipe diameter, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]
    Dh : float
        Hydraulic diameter, [m]
    """

    if Re <= 2000:
        # Laminar flow

        Resw = Re*(pi/(pi-4*delta/D)) * ((pi*D/H)**2)**0.5
        Sw = Resw/(H/2/D)**0.5

        fsw = 15.767/Resw * ((pi+2-2*delta/D)/(pi-4*delta/D))**2 \
            * (1+1e-6*Sw**2.55)**(1/6)

        f = fsw * Dh/D*(1+(pi*D/H)**2)**1.5

    elif Re > 1e4:
        # Turbulent flow
        f = 0.0791/Re**0.25 * (pi/(pi-4*delta/D))**1.75 \
            * ((pi+2-2*delta/D)/(pi-4*delta/D))**1.25 * (1+2.752/(H/D)**1.29)

    else:
        # Transition flow
        fl = f_twisted_Manglik(2000, D, H, delta, Dh)
        ft = f_twisted_Manglik(1e4, D, H, delta, Dh)

        # Eq 11 in 5
        f = (fl**10 + ft**10)**0.1

    return f


# Heat Transfer coefficient correlations
@refDoc(__doi__, [2])
def Nu_twisted_Plessis(Re, Pr, D, H, delta, Ae, De, x=None):
    """Calculate Nusselt number for a pipe with a twisted-tape insert using
    the Plessis and Kröger correlation

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    D : float
        Internal diameter of tube, [m]
    H : float
        Tape pitch for twist of π radians (180º), [m]
    delta : float
        Tape thickness, [m]
    Ae : float
        Effective flow area, [m²]
    De : float
        Effective hydraulic diameter, [m]
    x : float, optional
        Length in axial flow, [m]

    Returns
    -------
    Nu : float
        Nusselt number, [-]

    """
    y = H/D
    A = pi*D**2/4
    Ac = A-D*delta

    # Eq 2
    Psye = (D/De)**2*Ae/A

    Ree = Re*A/Ae*De/D
    Omge = Ree/y

    if x is not None:
        # Correction for flow don't fully developed
        xe = x*Ac/Ae
        xe_ = xe/(Ree*Pr*De)
        xinf = xe_ * (1+0.04*(Omge*Pr)**3)**(1/3)
        Psye *= (1+0.153*xinf**-1.05)**(1/3)

    # Eq 14
    Nu = 1.58*Psye*(1+6.4e-5*(Omge*Pr)**3)**0.117*(1+0.002*Omge**1.4)**(1/7)

    return Nu


if __name__ == "__main__":
    st = TwistedTape(10, 1, 0)
    print(f_twisted_Plessis(50, 1, 10, 0, st.Ae, st.De))
