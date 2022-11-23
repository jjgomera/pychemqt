#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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

Psychrometry is the field of science concerned with the thermodynamic
properties of any gas-vapor mixture, nevertheless it's mainly used for moist
air, the mixture of dry air and water vapor. This module try to implement
the procedures to calculate the properties for humid air because of its great
practical importance in the simulation of air-conditioning systems, solid
dryers or scrubbers.


The module include several functions with the definition of stantdard
atmosphere:
    * :func:`_Pbar`: Calculate barometric pressure for a specified altitude
    * :func:`_height`: Calculate the altitude for a specified pressure
    * :func:`_Tbar`: Calculate standard temperature for a specified altitude

Saturation state properties:
    * :func:`_Psat`: Calculate saturation pressure for a specified Tdb
    * :func:`_Tsat`: Calculate saturation temperature for a specified pressure
    * :func:`_Ws`: Saturation humidity calculation procedure

Calculation procedures using a perfect gas definition for humid air:
    * :func:`_h`: Specific enthalpy calculation procedure
    * :func:`_v`: Specific volume calculation procedure
    * :func:`_W_twb`: Humidity ratio calculation procedure
    * :func:`_tdp`: Dew point temperature calculation procedure
    * :func:`_twb`: Wet bulb temperature calculation procedure

Calculation procedure used in the plot procedure:
    * :func:`_Tdb`: Dry bulb temperature calculation procedure
    * :func:`_Tdb_V`: Tdb calculation procedure from specified volume
    * :func:`_W_V`: Humidity ratio calculation procedure from specified volume

Finally for grouping all functionality and integrate in main program with a
OOP scheme it's define the class:
    * :class:`PsyState`: Psychrometric state general class with common
    functionality
    * :class:`PsyIdeal`: Psychrometric state model using idial gas equation
    * :class:`PsyVirial`: Unimplemented
    * :class:`PsyCoolprop`: Psychrometric state model using CoolProp library
    * :class:`PsyRefprop`: Unimplemented
    * :class:`PsychroState`: PsyState subclass used as define in preferences
'''


from configparser import ConfigParser
import logging
import os

from iapws import _Sublimation_Pressure, _Ice, IAPWS95
from iapws._iapws import _Henry
from iapws.humidAir import _virial
from iapws.iapws97 import _PSat_T, _Region1, prop0

from numpy import exp, roots, linspace, arange, concatenate
from numpy.lib.scimath import log
from tools.qt import QtWidgets
from scipy.optimize import fsolve

# Avoid raise error at import this module if the the optional dependence isn't
# meet
try:
    from CoolProp.HumidAirProp import HAProps, HAProps_Aux
except ImportError:
    pass

from lib import unidades
from lib.config import conf_dir
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "",
         "title": "2013 ASHRAE Handook. Fundamentals (SI Edition)",
         "ref": "",
         "doi": ""},
    2:
        {"autor": "Herrmann, S., Kretzschmar, H.-J., Gatley, D.P.",
         "title": "Thermodynamic Properties of Real Moist Air, Dry Air, Steam"
                  "Water, and Ice",
         "ref": "ASHRAE RP-1485",
         "doi": ""},
    3:
        {"autor": "Nelson, H.F., Sauer, H.J.",
         "title": "Formulation of High-Temperature Properties for Moist Air",
         "ref": "HVAC&R Research 8(2) (2002) 311-334",
         "doi": "10.1080/10789669.2002.10391444"},
    4:
        {"autor": "Harvey, A.H., Huang, P.H.",
         "title": "First-Principles Calculation of the Air-Water Second Virial"
                  " Coefficient",
         "ref": "Int. J. Thermophisics 28(2) (2007) 556-565",
         "doi": "10.1007/s10765-007-0197-8"}
        }


R = 8.314472
Ma = 28.966
Mw = 18.015268
e = Mw/Ma


@refDoc(__doi__, [1])
def _Pbar(Z):
    """
    Standard atmosphere pressure as a function of altitude as explain in [1]_
    pag 1.1 Eq 3

    Parameters
    ----------
    Z : float
        Altitude, [m]

    Returns
    -------
    P : float
        Standard barometric pressure, [Pa]

    Examples
    --------
    Selected point from Table 1 in [1]_

    >>> "%0.3f" % _Pbar(-500).kPa
    '107.478'
    >>> "%0.3f" % _Pbar(8000).kPa
    '35.600'
    """
    P = (1-2.25577e-5*Z) ** 5.2559
    return unidades.Pressure(P, "atm")


@refDoc(__doi__, [1])
def _height(P):
    """
    Inverted _Pbar function

    Parameters
    ----------
    P : float
        Standard barometric pressure, [Pa]

    Returns
    -------
    Z : float
        Altitude, [m]

    Examples
    --------
    Selected point from Table 1 in [1]_

    >>> "%0.0f" % _height(107478)
    '-500'
    """
    P_atm = P/101325.
    Z = 1/2.25577e-5*(1-exp(log(P_atm)/5.2559))
    return unidades.Length(Z)


@refDoc(__doi__, [1])
def _Tbar(Z):
    """
    Standard temperature as a function of altitude as explain in [1]_
    pag 1.1 Eq 4

    Parameters
    ----------
    Z : float
        Altitude, [m]

    Returns
    -------
    T : float
        Temperature, [K]

    Examples
    --------
    Selected point from Table 1 in [1]_

    >>> "%0.1f" % _Tbar(-500).C
    '18.2'
    >>> "%0.1f" % _Tbar(8000).C
    '-37.0'
    """
    t = 15-0.0065*Z
    return unidades.Temperature(t, "C")


@refDoc(__doi__, [1])
def _Psat(T):
    """
    Water vapor saturation pressure calculation as explain in [1]_
    pag 1.2, Eq 5-6

    Parameters
    ----------
    T : float
        Temperature, [K]

    Returns
    -------
    P : float
        Saturation pressure, [Pa]
    """
    if 173.15 <= T < 273.15:
        # Saturation pressure over ice, Eq 5
        C = [-5674.5359, 6.3925247, -0.009677843, 0.00000062215701,
             2.0747825E-09, -9.484024E-13, 4.1635019]
        pws = exp(C[0]/T + C[1] + C[2]*T + C[3]*T**2 + C[4]*T**3 + C[5]*T**4
                  + C[6]*log(T))
    elif 273.15 <= T <= 473.15:
        # Saturation pressure over liquid water, Eq 6
        C = [-5800.2206, 1.3914993, -0.048640239, 0.000041764768,
             -0.000000014452093, 6.5459673]
        pws = exp(C[0]/T + C[1] + C[2]*T + C[3]*T**2 + C[4]*T**3 + C[5]*log(T))
    else:
        raise NotImplementedError("Incoming out of bound")

    return unidades.Pressure(pws)


@refDoc(__doi__, [1])
def _Tsat(Pv):
    """
    Water vapor saturation temperature calculation as explain in [1]_, inverted
    from pag 1.2, Eq 5-6 to calculate Tdb

    Parameters
    ----------
    P : float
        Saturation pressure, [Pa]

    Returns
    -------
    T : float
        Temperature, [K]
    """
    Pv_min = _Psat(173.15)
    Pv_lim = _Psat(273.15)
    Pv_max = _Psat(473.15)
    if Pv_min <= Pv < Pv_lim:

        C = [-5674.5359, 6.3925247, -0.009677843, 0.00000062215701,
             2.0747825E-09, -9.484024E-13, 4.1635019]

        def f(T):
            pws = exp(C[0]/T + C[1] + C[2]*T + C[3]*T**2 + C[4]*T**3 +
                      C[5]*T**4 + C[6]*log(T))
            return Pv - pws
        T = fsolve(f, 250)[0]

    elif Pv_lim <= Pv <= Pv_max:

        C = [-5800.2206, 1.3914993, -0.048640239, 0.000041764768,
             -0.000000014452093, 6.5459673]

        def f(T):
            pws = exp(C[0]/T + C[1] + C[2]*T + C[3]*T**2 + C[4]*T**3 +
                      C[5]*log(T))
            return Pv - pws
        T = fsolve(f, 300)[0]

    else:
        raise NotImplementedError("Incoming out of bound")

    return unidades.Temperature(T)


@refDoc(__doi__, [1])
def _Ws(P, Tdb):
    """
    Calculate the saturated humidity ratio of a humid air as a perfet gas, Ws,
    as explain in [1]_, pag 1.8, Eq 22

    Parameters
    ----------
    P : float
        Pressure, [Pa]
    Tdb: float
        Dry bulb temperature, [K]

    Returns
    -------
    Ws : float
       Saturated Humidity ratio, [kgw/kgda]
    """
    pv = _Psat(Tdb)
    return 0.621945*pv/(P-pv)


# Humid air as a perfect gas correlations
@refDoc(__doi__, [1])
def _h(Tdb, W):
    """
    Calculate the enthalpy of a humid air as a perfect gas, h, as explain
    in [1]_, pag 1.8, Eq 22

    Parameters
    ----------
    Tdb: float
       Dry bulb temperature, [K]
    W : float
       Humidity ratio, [kgw/kgda]

    Returns
    -------
    h : float
       Enthalpy of humid air, [kJ/kgda]
    """
    # Temperature en celsius
    tc = Tdb-273.15
    h = 1.006*tc + W*(2501 + 1.86*tc)
    return unidades.Enthalpy(h)


@refDoc(__doi__, [1])
def _v(P, Tdb, W):
    """
    Calculate the specific volume of a humid air as a perfect gas, v, as
    explain in [1]_, pag 1.8, Eq 28

    Parameters
    ----------
    P : float
        Pressure, [Pa]
    Tdb: float
       Dry bulb temperature, [K]
    W : float
       Humidity ratio, [kgw/kgda]

    Returns
    -------
    v : float
       Specific volume, [m³/kgda]
    """
    # Pressure in kPa
    P_kpa = P/1000
    v = 0.287042*Tdb*(1+1.607858*W)/P_kpa
    return unidades.SpecificVolume(v)


@refDoc(__doi__, [1])
def _W_twb(tdb, twb, P):
    """
    Calculate the humidity ratio of a humid air as a perfet gas, W, as explain
    in [1]_, pag 1.9, Eq 35-37

    Parameters
    ----------
    Tdb: float
        Dry bulb temperature, [K]
    Twb: float
        Wet bulb temperature, [K]
    P : float
        Pressure, [Pa]

    Returns
    -------
    W : float
       Humidity ratio, [kgw/kgda]
    """
    tdb_C = tdb - 273.15
    twb_C = twb - 273.15
    Ws = _Ws(P, tdb)
    if tdb_C >= 0:
        # Eq 35 for liquid water
        w = ((2501-2.326*twb_C)*Ws-1.006*(tdb-twb)) / \
            (2501+1.86*tdb_C-4.186*twb_C)
    else:
        # Eq 37 for ice water
        w = ((2830-0.24*twb_C)*Ws-1.006*(tdb-twb))/(2830+1.86*tdb_C-2.1*twb_C)
    return w


@refDoc(__doi__, [1])
def _tdp(Pw):
    """
    Calculate the dew-point temperature of a humid air as a perfet gas, as
    explain in [1]_, pag 1.9, Eq 39-40

    Parameters
    ----------
    P : float
        Water vapor partial pressure, [Pa]

    Returns
    -------
    Tdp: float
        Dew-point temperature, [K]
    """
    C = [6.54, 14.526, 0.7389, 0.09486, 0.4569, 0.1984]
    D = [6.09, 12.608, 0.4959]

    # Pw used in kPa
    a = log(Pw/1000.)

    Tdp1 = C[0] + C[1]*a + C[2]*a**2 + C[3]*a**3 + C[4]*(Pw/1000.)**C[5]
    Tdp2 = D[0] + D[1]*a + D[2]*a**2

    if 0 <= Tdp1 <= 93:
        t = Tdp1
    elif Tdp2 < 0:
        t = Tdp2
    else:
        raise NotImplementedError("Incoming out of bound")

    return unidades.Temperature(t, "C")


@refDoc(__doi__, [1])
def _twb(tdb, W, P):
    """
    Calculate the wet-bulb temperature of a humid air as a perfet gas, as
    explain in [1]_, pag 1.9, inveted Eq 35-37

    Parameters
    ----------
    Tdb: float
        Dry bulb temperature, [K]
    W : float
       Humidity ratio, [kgw/kgda]
    P : float
        Pressure, [Pa]

    Returns
    -------
    Twb: float
        Wet bulb temperature, [K]
    """
    tdb_C = tdb - 273.15

    def f(twb):
        Pvs = _Psat(twb)
        Ws = 0.62198*Pvs/(P-Pvs)
        twb_C = twb - 273.15
        if tdb >= 0:
            w = ((2501.-2.326*twb_C)*Ws-1.006*(tdb_C-twb_C)) / \
                (2501.+1.86*tdb_C-4.186*twb_C)-W
        else:
            w = ((2830-0.24*twb_C)*Ws-1.006*(tdb-twb)) / \
                (2830+1.86*tdb_C-2.1*twb_C)
        return w

    twb = fsolve(f, tdb)[0]
    return twb


# Procedures only used in plot
@refDoc(__doi__, [1])
def _Tdb(twb, w, P):
    """
    Calculate the dry-bulb temperature of a humid air as a perfet gas, as
    explain in [1]_, pag 1.9, inveted Eq 35-37

    Parameters
    ----------
    Twb: float
        Wet bulb temperature, [K]
    W : float
       Humidity ratio, [kgw/kgda]
    P : float
        Pressure, [Pa]

    Returns
    -------
    Tdb: float
        Dry bulb temperature, [K]
    """
    tw = twb-273.15
    ws = _Ws(P, twb)
    td = ((2501-2.381*tw)*ws+1.006*tw-w*(2501-4.186*tw))/(w*1.805+1.006)
    return td+273.15


@refDoc(__doi__, [1])
def _Tdb_V(v, P):
    """
    Calculate the specific volume of a humid air as a perfect gas, v, as
    explain in [1]_, pag 1.8, inverted Eq 28, used to calculate Tdb in a
    isochor line in plots

    Parameters
    ----------
    v : float
       Specific volume, [m³/kgda]
    P : float
        Pressure, [Pa]

    Returns
    -------
    Tdb: float
       Dry bulb temperature, [K]
    """
    P_kpa = P/1000

    def f(Tdb):
        w = _Ws(P, Tdb)
        return v-0.287042*Tdb*(1+1.607858*w)/P_kpa
    ts = fsolve(f, 300)[0]
    return ts


def _W_V(Tdb, P, v):
    """
    Calculate the humidity ratio of a humid air as a perfect gas, as
    explain in [1]_, pag 1.8, inverted Eq 28, used to calculate W in a
    isochor line in plots

    Parameters
    ----------
    Tdb: float
       Dry bulb temperature, [K]
    P : float
        Pressure, [Pa]
    v : float
       Specific volume, [m³/kgda]

    Returns
    -------

    Function to calculate isochor line
    input:
        dry bulb temperature, K
        barometric pressure, Pa
        specified volume, m3/kg air
    return
        humidity ratio, kg H2O/kg air
    """
    P_kpa = P/1000
    return (v*P_kpa-0.287042*Tdb)/(0.287042*1.607858*Tdb)


class PsyState(object):
    """
    Class to model a psychrometric state with properties
    kwargs definition parameters:
        P: Pressure, Pa
        z: altitude, m

        tdp: dew-point temperature
        tdb: dry-bulb temperature
        twb: web-bulb temperature
        w: Humidity Ratio [kg water/kg dry air]
        HR: Relative humidity
        h: Mixture enthalpy
        v: Mixture specified volume

    P: mandatory input for barometric pressure, z as an alternative P input
    it needs other two input parameters:
        0 - tdb, w
        1 - tdb, HR
        2 - tdb, twb
        3 - tdb, tdp
        4 - tdp, HR
        5 - tdp, twb
        6 - twb, w
    """
    kwargs = {"z": 0.0,
              "P": 0.0,

              "tdb": 0.0,
              "tdb": 0.0,
              "twb": 0.0,
              "w": None,
              "HR": None,
              "h": None,
              "v": 0.0}
    status = 0
    msg = "Unknown variables"

    TEXT_MODE = [
        QtWidgets.QApplication.translate("pychemqt", "T dry bulb, Humidity Ratio"),
        QtWidgets.QApplication.translate("pychemqt", "T dry bulb, Relative humidity"),
        QtWidgets.QApplication.translate("pychemqt", "T dry bulb, T wet bulb"),
        QtWidgets.QApplication.translate("pychemqt", "T dry bulb, T dew point"),
        QtWidgets.QApplication.translate("pychemqt", "T dew point, Relative humidity"),
        QtWidgets.QApplication.translate("pychemqt", "T wet bulb, Relative humidity")
        ]
    VAR_NAME = [
        ("tdb", "w"),
        ("tdb", "HR"),
        ("tdb", "twb"),
        ("tdb", "tdp"),
        ("tdp", "HR"),
        ("twb", "HR")
        ]

#        QtWidgets.QApplication.translate("pychemqt", "T dry bulb, Enthalpy"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª bulbo seco, Densidad"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª bulbo húmedo, H absoluta"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª bulbo húmedo, Entalpia"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª bulbo húmedo, Densidad"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª bulbo húmedo, Tª rocio"))
#        QtWidgets.QApplication.translate("pychemqt", "H absoluta, entalpía"))
#        QtWidgets.QApplication.translate("pychemqt", "H relativa, entalpía"))
#        QtWidgets.QApplication.translate("pychemqt", "H absoluta, densidad"))
#        QtWidgets.QApplication.translate("pychemqt", "H relativa, densidad"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª rocio, entalpía"))
#        QtWidgets.QApplication.translate("pychemqt", "Tª rocio, densidad"))

    def __init__(self, **kwargs):
        self.kwargs = self.__class__.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()
            logging.debug(QtWidgets.QApplication.translate(
                "pychemqt", "Calculate psychrometric point"))
            logging.debug(self.kwargs)
            self.msg = "Solved"

    @property
    def calculable(self):
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        # h = self.kwargs.get("h", None)
        # v = self.kwargs.get("v", 0)

        self.mode = -1
        if tdb and w is not None:
            self.mode = 0
        elif tdb and HR is not None:
            self.mode = 1
        elif tdb and twb:
            self.mode = 2
        elif tdb and tdp:
            self.mode = 3
        elif tdp and HR is not None:
            self.mode = 4

        return bool(self.mode+1)

    def _P(self):
        """Barometric pressure calculation, Pa"""
        if self.kwargs["P"]:
            P = self.kwargs["P"]
        elif self.kwargs["z"]:
            P = _Pbar(self.kwargs["z"])
        else:
            P = 101325.
        return P

    def _lib(self):
        """Properties calculate library, customize in each subclass"""
        pass

    def calculo(self):
        tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h = self._lib()
        self.tdp = unidades.Temperature(tdp)
        self.tdb = unidades.Temperature(tdb)
        self.twb = unidades.Temperature(twb)
        self.P = unidades.Pressure(P)
        self.Pvs = unidades.Pressure(Pvs)
        self.Pv = unidades.Pressure(Pv)
        self.ws = unidades.Dimensionless(ws, txt="kgw/kgda")
        self.w = unidades.Dimensionless(w, txt="kgw/kgda")
        self.HR = unidades.Dimensionless(HR, txt="%")
        self.mu = unidades.Dimensionless(w/ws*100)
        self.v = unidades.SpecificVolume(v)
        self.rho = unidades.Density(1/v)
        self.h = unidades.Enthalpy(h, "kJkg")
        self.Xa = 1/(1+self.w/0.62198)
        self.Xw = 1-self.Xa

    @classmethod
    def calculatePlot(cls):
        """Funtion to calculate point in chart, each child class must define
        it, as default use ideal gas equation of state"""
        return PsyIdeal.calculatePlot()

    @staticmethod
    def LineList(name, Preferences):
        """Return a list with the values of isoline name to plot"""
        if Preferences.getboolean("Psychr", name+"Custom"):
            t = []
            for i in Preferences.get("Psychr", name+'List').split(','):
                if i:
                    t.append(float(i))
        else:
            start = Preferences.getfloat("Psychr", name+"Start")
            end = Preferences.getfloat("Psychr", name+"End")
            step = Preferences.getfloat("Psychr", name+"Step")
            t = list(arange(start, end, step))
        return t


class PsyIdeal(PsyState):
    """Psychrometric state using ideal gas equation"""
    def _lib(self):
        """Properties calculate library"""
        P = self._P()
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

        if self.mode == 0:
            # Tdb and w
            Pvs = _Psat(tdb)
            ws = 0.621945*Pvs/(P-Pvs)
            Pv = w*P/(0.621945+w)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)

        elif self.mode == 1:
            # Tdb and HR
            Pvs = _Psat(tdb)
            ws = 0.621945*Pvs/(P-Pvs)
            Pv = Pvs*HR/100
            w = 0.621945*Pv/(P-Pv)
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)

        elif self.mode == 2:
            # Tdb and Twb
            Pvs = _Psat(tdb)
            ws = 0.621945*Pvs/(P-Pvs)
            w = _W_twb(tdb, twb, P)
            Pv = w*P/(0.621945+w)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)

        elif self.mode == 3:
            # Tdb and Tdp
            Pv = _Psat(tdp)
            w = 0.621945*Pv/(P-Pv)
            Pvs = _Psat(tdb)
            ws = 0.621945*Pvs/(P-Pvs)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            twb = _twb(tdb, w, P)

        elif self.mode == 4:
            # Tdp and HR
            Pv = _Psat(tdp)
            if HR:
                w = 0.621945*Pv/(P-Pv)
                Pvs = Pv/HR*100
            else:
                w = 0
                Pvs = Pv
            ws = 0.621945*Pvs/(P-Pvs)
            tdb = _Tsat(Pvs)
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            twb = _twb(tdb, w, P)

        elif self.mode == 5:
            # Tdp and Twb
            # Pv = _Psat(tdp)
            pass

        elif self.mode == 6:
            # Tdb and h
            pass
#            self.Tdb = tdb
#            self.h = unidades.Enthalpy(h, "kJkg")
#            f = lambda w: self.Entalpia(self.Tdb, w).kJkg-h
#            self.w = fsolve(f, 0.001)
#            self.Hs = self.Humedad_Absoluta(tdb)
#            self.twb = self.Tw(tdb, self.w)
#            self.HR = self.w/self.Hs*100
#            self.Xa = 1/(1+self.w*self.aire.M/self.agua.M)
#            self.Xw = 1-self.Xa
#            self.V = self.Volumen(tdb, self.Xa)
#            self.rho = unidades.Density(1/self.V)

        return tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h

    @classmethod
    def calculatePlot(cls, parent):
        """Funtion to calculate point in chart"""
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        parent.setProgressValue(0)

        data = {}
        P = parent.inputs.P.value
        t = cls.LineList("isotdb", Preferences)

        # Saturation line
        Hs = []
        Pvs = []
        for ti in t:
            Pv = _Psat(ti)
            Pvs.append(Pv)
            Hs.append(0.62198*Pv/(P-Pv))
            parent.setProgressValue(5*len(Hs)/len(t))
        data["t"] = t
        data["Hs"] = Hs

        # left limit of isow lines
        H = cls.LineList("isow", Preferences)
        th = []
        for w in H:
            if w:
                Pv = w*P/(0.62198+w)
                th.append(unidades.Temperature(_tdp(Pv)))
            else:
                tmin = Preferences.getfloat("Psychr", "isotdbStart")
                th.append(unidades.Temperature(tmin))
        data["H"] = H
        data["th"] = th

        # Humidity ratio lines
        hr = cls.LineList("isohr", Preferences)
        Hr = {}
        cont = 0
        for i in hr:
            Hr[i] = []
            for pvs in Pvs:
                pv = pvs*i/100
                Hr[i].append(0.62198*pv/(P-pv))
                cont += 1
                parent.setProgressValue(5+10*cont/len(hr)/len(Hs))
        data["Hr"] = Hr

        # Twb
        lines = cls.LineList("isotwb", Preferences)
        Twb = {}
        cont = 0
        for T in lines:
            H = concatenate((arange(_Ws(P, T), 0, -0.001), [0.]))
            Tw = []
            for h in H:
                Tw.append(unidades.Temperature(_Tdb(T, h, P)))
            cont += 1
            parent.setProgressValue(15+75*cont/len(lines))
            Twb[T] = (list(H), Tw)
        data["Twb"] = Twb

        # v
        lines = cls.LineList("isochor", Preferences)
        V = {}
        for cont, v in enumerate(lines):
            ts = _Tdb_V(v, P)
            T = linspace(ts, v*P/287.055, 50)
            Td = [unidades.Temperature(ti) for ti in T]
            H = [_W_V(ti, P, v) for ti in T]
            parent.setProgressValue(90+10*cont/len(lines))
            V[v] = (Td, H)
        data["v"] = V

        return data


class PsyVirial(PsyState):
    """Psychrometric state using virial equation of state"""
    def _lib(self):
        """Properties calculate library"""
        P = self._P()/1e6  # Convert to MPa
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

        # f = self._f(tdb, P)
        # print(f)
        # self.f = f

        # from numpy import roots

        # vir = self._Virial(tdb, w)
        # Bm = vir["Bm"]*1e6
        # Cm = vir["Bm"]*1e6
        # vm = roots([1, -R*tdb/P, -R*tdb*Bm/P, -R*tdb*Cm/P])
        # print(vm)
        # # if vm[0].imag==0.0:
            # v=vm[0].real
        # else:
            # v=vm[2].real
         # return unidades.SpecificVolume(v/self.aire.M/Xa)

        if self.mode == 0:
            # Tdb and w
            Pvs = self._Pwsat(tdb)
            f = self._f(tdb, P)
            phiws = f*Pvs/P
            ws = e*phiws/(1-phiws)
            phiw = w/(e+w)
            M = (1-phiw)*Ma+phiw*Mw
            Pv = phiw*P/f
            HR = Pv/Pvs*100
            v_ = self._v(P, tdb, w)
            v = (1+w)*v_/M*1e-3
            tdp = self._tdp(P, phiw, tdb)
            h_ = self._h(v_, tdb, phiw)
            h = h_/M*(1+w)
            twb = self._twb(tdb, w, P)
            print(Pvs, Pv, phiw/phiws*100)
            print(w, phiw, HR, v, h)

        elif self.mode == 1:
            # Tdb and HR
            Pvs = self._Pwsat(tdb)
            f = self._f(tdb, P)
            phiws = f*Pvs/P
            ws = e*phiws/(1-phiws)
            Pv = HR*Pvs*100
            w = HR/100*ws
            phiw = w/(e+w)
            M = (1-phiw)*Ma+phiw*Mw
            Pv = phiw*P/f
            v_ = self._v(P, tdb, w)
            v = (1+w)*v_/M*1e-3
            tdp = self._tdp(P, phiw, tdb)
            h_ = self._h(v_, tdb, phiw)
            h = h_/M*(1+w)
            twb = self._twb(tdb, w, P)
            print(Pvs, Pv, phiw/phiws*100)
            print(w, v, h, phiw, HR)
            # Pvs = _Psat(tdb)
            # ws = 0.621945*Pvs/(P-Pvs)
            # Pv = Pvs*HR/100
            # w = 0.621945*Pv/(P-Pv)
            # v = _v(P, tdb, w)
            # h = _h(tdb, w)
            # tdp = _tdp(Pv)
            # twb = _twb(tdb, w, P)

        # elif self.mode == 2:
            # # Tdb and Twb
            # Pvs = _Psat(tdb)
            # ws = 0.621945*Pvs/(P-Pvs)
            # w = _W_twb(tdb, twb, P)
            # Pv = w*P/(0.621945+w)
            # HR = Pv/Pvs*100
            # v = _v(P, tdb, w)
            # h = _h(tdb, w)
            # tdp = _tdp(Pv)

        # elif self.mode == 3:
            # # Tdb and Tdp
            # Pv = _Psat(tdp)
            # w = 0.621945*Pv/(P-Pv)
            # Pvs = _Psat(tdb)
            # ws = 0.621945*Pvs/(P-Pvs)
            # HR = Pv/Pvs*100
            # v = _v(P, tdb, w)
            # h = _h(tdb, w)
            # twb = _twb(tdb, w, P)

        # elif self.mode == 4:
            # # Tdp and HR
            # Pv = _Psat(tdp)
            # if HR:
                # w = 0.621945*Pv/(P-Pv)
                # Pvs = Pv/HR*100
            # else:
                # w = 0
                # Pvs = Pv
            # ws = 0.621945*Pvs/(P-Pvs)
            # tdb = _Tsat(Pvs)
            # v = _v(P, tdb, w)
            # h = _h(tdb, w)
            # twb = _twb(tdb, w, P)

        # elif self.mode == 5:
            # # Tdp and Twb
            # # Pv = _Psat(tdp)
            # pass

        # elif self.mode == 6:
            # # Tdb and h
            # pass
# #            self.Tdb = tdb
# #            self.h = unidades.Enthalpy(h, "kJkg")
# #            f = lambda w: self.Entalpia(self.Tdb, w).kJkg-h
# #            self.w = fsolve(f, 0.001)
# #            self.Hs = self.Humedad_Absoluta(tdb)
# #            self.twb = self.Tw(tdb, self.w)
# #            self.HR = self.w/self.Hs*100
# #            self.Xa = 1/(1+self.w*self.aire.M/self.agua.M)
# #            self.Xw = 1-self.Xa
# #            self.V = self.Volumen(tdb, self.Xa)
# #            self.rho = unidades.Density(1/self.V)

        return tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h

    def _virial(self, T):
        """Calculate the humid-air virial coefficient.

        Parameters
        ----------
        T : float
            Temperature [K]

        Returns
        -------
        prop : dict
            Dictionary with critical coefficient:

            * Baa: Second virial coefficient of dry air, [cm³/mol]
            * Baw: Second air-water cross virial coefficient, [cm³/mol]
            * Bww: Second virial coefficient of water, [cm³/mol]
            * Caaa: Third virial coefficient of dry air, [cm⁶/mol]
            * Caaw: Third air-water cross virial coefficient, [cm⁶/mol]
            * Caww: Third air-water cross virial coefficient, [cm⁶/mol]
            * Cwww: Third virial coefficient of dry air, [cm⁶/mol]
        """
        vir = _virial(T)
        vir["Baa"] *= 1e6
        vir["Baw"] *= 1e6
        vir["Bww"] *= 1e6
        vir["dBaaT"] *= 1e6
        vir["dBawT"] *= 1e6
        vir["dBwwT"] *= 1e6
        vir["Caaa"] *= 1e12
        vir["Caaw"] *= 1e12
        vir["Caww"] *= 1e12
        vir["Cwww"] *= 1e12
        vir["dCaaaT"] *= 1e12
        vir["dCaawT"] *= 1e12
        vir["dCawwT"] *= 1e12
        vir["dCwwwT"] *= 1e12

        # Caww slight differences in Nelson paper to iapws TEOS-10 model
        # Using Eq 11
        T_ = T/100
        bi = [-10.72887, 34.7804, -38.3383, 33.406]
        vir["Caww"] = -1e6*exp(sum([b/T_**i for i, b in enumerate(bi)]))
        vir["dCawwT"] = 1e6*T_/T * \
            sum([i*b*T_**(-i-1) for i, b in enumerate(bi)]) * \
            exp(sum([b/T_**i for i, b in enumerate(bi)]))

        return vir

    def _virialMixture(self, T, phi_w):
        """Calculate the mixture humid air virial coefficient

        Parameters
        ----------
        T : float
            Temperature, [K]
        phi_w : float
            Molar fraction of water, [-]

        Returns
        -------
        """

        vir = self._virial(T)
        Baa = vir["Baa"]
        Caaa = vir["Caaa"]
        Bww = vir["Bww"]
        Cwww = vir["Cwww"]
        Baw = vir["Baw"]
        Caaw = vir["Caaw"]
        Caww = vir["Caww"]

        dBaaT = vir["dBaaT"]
        dBawT = vir["dBawT"]
        dBwwT = vir["dBwwT"]
        dCaaaT = vir["dCaaaT"]
        dCaawT = vir["dCaawT"]
        dCawwT = vir["dCawwT"]
        dCwwwT = vir["dCwwwT"]

        Bm = (1-phi_w)**2*Baa + 2*(1-phi_w)*phi_w*Baw + phi_w**2*Bww
        Cm = (1-phi_w)**3*Caaa + 3*(1-phi_w)**2*phi_w*Caaw + \
            3*(1-phi_w)*phi_w**2*Caww + phi_w**3*Cwww
        dBmT = (1-phi_w)**2*dBaaT + 2*(1-phi_w)*phi_w*dBawT + phi_w**2*dBwwT
        dCmT = (1-phi_w)**3*dCaaaT + 3*(1-phi_w)**2*phi_w*dCaawT + \
            3*(1-phi_w)*phi_w**2*dCawwT + phi_w**3*dCwwwT

        prop = {}
        prop["Bm"] = Bm
        prop["Cm"] = Cm
        prop["dBmT"] = dBmT
        prop["dCmT"] = dCmT
        return prop

    def _Pwsat(self, T):
        """Calculate the partial pressure of water in saturated moint air"""
        if T < 273.15:
            Pws = _Sublimation_Pressure(T)
        else:
            Pws = _PSat_T(T)
        return Pws

    def _f(self, T, P):
        P /= 1e6

        vir = self._virial(T)
        Baa = vir["Baa"]
        Bww = vir["Bww"]
        Baw = vir["Baw"]
        Caaa = vir["Caaa"]
        Caaw = vir["Caaw"]
        Caww = vir["Caww"]
        Cwww = vir["Cwww"]

        Pws = self._Pwsat(T)
        if T < 273.15:
            kt = _Ice(T, P)["xkappa"]
            vws = _Ice(T, Pws)["v"]
        else:
            kt = _Region1(T, P)["kt"]
            vws = _Region1(T, Pws)["v"]

        vws *= Mw*1000  # Convert from m³/kg to cm³/mol

        # Henry constant
        H_N2 = _Henry(T, "N2")
        H_O2 = _Henry(T, "O2")
        H_Ar = _Henry(T, "Ar")
        H = 1/1.01325*(0.7812/H_N2+0.2095/H_O2+0.0093/H_Ar)

        if Pws > P:
            kt = 0
            H = 0

        def _f(f):
            ws = f*Pws/P
            lnf = ((1+kt*Pws)*(P-Pws)-kt*(P**2-Pws**2)/2)/R/T*vws \
                + log(1-H*(1-ws)*P) \
                + (1-ws)**2*P/R/T*Baa \
                - 2*(1-ws)**2*P/R/T*Baw \
                - (P-Pws-(1-ws)**2*P)/R/T*Bww \
                + (1-ws)**3*P**2/(R*T)**2*Caaa \
                + 3*(1-ws)**2*(1-2*(1-ws))*P**2/2/(R*T)**2*Caaw \
                - 3*(1-ws)**2*ws*P**2/(R*T)**2*Caww \
                - ((3-2*ws)*ws**2*P**2-Pws**2)/2/(R*T)**2*Cwww \
                - (1-ws)**2*(3*ws-2)*ws*P**2/(R*T)**2*Baa*Bww \
                - 2*(1-ws)**3*(3*ws-1)*P**2/(R*T)**2*Baa*Baw \
                + 6*(1-ws)**2*ws**2*P**2/(R*T)**2*Bww*Baw \
                - 3*(1-ws)**4*P**2/2/(R*T)**2*Baa**2 \
                - 2*(1-ws)**2*ws*(3*ws-2)*P**2/(R*T)**2*Baw**2 \
                - (Pws**2-(4-3*ws)*ws**3*P**2)/2/(R*T)**2*Bww**2
            return log(f) - lnf

        f = fsolve(_f, 1)[0]

        if f < 1:
            f = 1
        return f

    def _v(self, P, T, phi_w):
        """Calculate the real volume of a humid air using the virial equation
        of state

        Parameters
        ----------
        T : float
            Temperature, [K]
        phi_w : float
            Molar fraction of water, [-]

        Returns
        -------
        """
        vir = self._virialMixture(T, phi_w)
        Bm = vir["Bm"]
        Cm = vir["Cm"]

        vm = roots([1, -R*T/P, -R*T*Bm/P, -R*T*Cm/P])
        if vm[0].imag == 0.0:
            v = vm[0].real
        else:
            v = vm[2].real

        return v

    def _h(self, v, T, phi_w):


        ho = 2.924425468

        # Air ideal-gas enthalpy
        from iapws.humidAir import Air
        ho_lem = -7914.149298 + 7906.10273632
        air = Air()
        st0 = air._prop0(1e5, T)
        ha = ho_lem + st0.h*Ma                                        # Eq 3.46

        # Water ideal-gas enthalpy
        if T > 273.15:
            ho_97 = -0.01102142797 - 45064.07018385615 + 0.000611782*Mw
            st0 = prop0(T, 1e5)
            hw = ho_97 + st0["h"]*Mw                                  # Eq 3.47
        else:
            ho_95 = -0.01102303806 - 45064.4896716 + 0.000611782*Mw
            wt = IAPWS95()
            st0 = wt._prop0(1e5, T)
            hw = ho_95 + st0.h*Mw                                     # Eq 3.48

        print("ha", ha)
        print("hw", hw)

        vir = self._virialMixture(T, phi_w)
        Bm = vir["Bm"]
        Cm = vir["Cm"]
        dBmT = vir["dBmT"]
        dCmT = vir["dCmT"]
        print(T, v/Ma, phi_w, Bm, Cm, dBmT, dCmT)

        # Eq 3.45
        v *= 1e-0
        h = ho + (1-phi_w)*ha + phi_w*hw + \
            R*T*((Bm-T*dBmT)/v+(Cm-T/2*dCmT)/v**2)

        return h

    def _tdp(self, P, phi_w, tdb):
        """Calculation of dew-point temperature"""
        if not phi_w:
            return None

        def f(t):
            f = self._f(t, P)
            pws = self._Pwsat(t)
            pw = phi_w*P
            return pw-f*pws

        tdp = fsolve(f, tdb)
        return tdp

    def _twb(self, P, phi_w, tdb):
        return None

    # def _V_Virial(self, T, P, phi_w):
#        """volumen por unidad de masa de aire seco"""
#        # FIXME: Don't work, for now use ideal gas equation
#        Bm, Cm=self.Virial(T, Xa)
#        vm=roots([1, -R_atml*T/self.P.atm, -R_atml*T*Bm/self.P.atm,-R_atml*T*Cm/self.P.atm])
#        if vm[0].imag==0.0:
#            v=vm[0].real
#        else:
#            v=vm[2].real
#        return unidades.SpecificVolume(v/self.aire.M/Xa)
#
#    def _h(self, Td, w):
#        """Enthalpy calculation procedure"""
#        cp_air = self.Air.Cp_Gas_DIPPR(Td)
#        cp_water = self.Water.Cp_Gas_DIPPR(Td)
#        h = (Td-273.15)*(cp_air+cp_water*w)
#        return unidades.Enthalpy(h, "kJkg")
#
#    def _HS(self, Td):
#        """Saturation humidity calculation procedure"""
#        pv = self._Ps(Td)
#        return self.Water.M*pv/(self.Air.M*(self.P-pv))
#
#    def _Ps(self, Td):
#        """Saturation pressure calculation procedure"""
#        if Td < 273.15:
#            return _Sublimation_Pressure(Td)
#        else:
#            return _PSat_T(Td)
#
#    def _Tw(self, Td, w):
#        Hv = self.agua.Hv_DIPPR(Td)
#        Cs = self.Calor_Especifico_Humedo(Td, w)
#
#        def f(Tw):
#            return self.Humedad_Absoluta(Tw)-w-Cs/Hv*(Td-Tw)
#        Tw = fsolve(f, Td)
#        return unidades.Temperature(Tw)



class PsyCoolprop(PsyState):
    """Psychrometric state using coolprop external library"""

    @property
    def calculable(self):
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        # h = self.kwargs.get("h", None)
        # v = self.kwargs.get("v", 0)

        self._mode = 0
        if tdb and w is not None:
            self._mode = ("Tdb", "W")
        elif tdb and HR is not None:
            self._mode = ("Tdb", "RH")
        elif tdb and twb:
            self._mode = ("Tdb", "Twb")
        elif tdb and tdp:
            self._mode = ("Tdb", "Tdp")
        elif tdp and HR is not None:
            self._mode = ("Tdp", "RH")

        return bool(self._mode)

    def args(self):
        # Correct coolprop custom namespace versus pychemqt namespace
        if "Tdb" in self._mode:
            self.kwargs["Tdb"] = self.kwargs["tdb"]
        if "Twb" in self._mode:
            self.kwargs["Twb"] = self.kwargs["twb"]
        if "Tdp" in self._mode:
            self.kwargs["Tdp"] = self.kwargs["tdp"]
        if "RH" in self._mode:
            self.kwargs["RH"] = self.kwargs["HR"]
        if "W" in self._mode:
            self.kwargs["W"] = self.kwargs["w"]

        var1 = self.kwargs[self._mode[0]]
        var2 = self.kwargs[self._mode[1]]

        # units conversion to coolprop expected unit:
        # HR in 0-1, H in kJ/kg, S in kJ/kgK
#        if self._mode[0] == "P":
#            var1 /= 1000.
        if "RH" in self._mode[0]:
            var1 /= 100.
        if "RH" in self._mode[1]:
            var2 /= 100.
#        if self._mode == "HS":
#            var1 /= 1000.
#        if "S" in self._mode:
#            var2 /= 1000.

        args = ("P", self._P_kPa, self._mode[0], var1, self._mode[1], var2)
        return args

    @property
    def _P_kPa(self):
        """Property for ease access to pressure in kPa"""
        P = self._P()
        return P/1000.

    def _lib(self):
        args = self.args()
        P = self._P()

        if "Tdb" in self._mode:
            tdb = self.kwargs["Tdb"]
        else:
            tdb = HAProps("Tdb", *args)
        tdp = HAProps("Tdp", *args)
        twb = HAProps("Twb", *args)
        w = HAProps("W", *args)
        HR = HAProps("RH", *args)*100
        Pvs = HAProps_Aux("p_ws", tdb, self._P_kPa, w)[0]*1000
        Pv = Pvs*HR/100
        ws = HAProps("W", "P", self._P_kPa, "Tdb", tdb, "RH", 1)
        v = HAProps("V", *args)
        h = HAProps("H", *args)

        return tdp, tdb, twb, P, Pvs, Pv, ws, w, HR, v, h

    @classmethod
    def calculatePlot(cls, parent):
        """Funtion to calculate point in chart"""
        Preferences = ConfigParser()
        Preferences.read(conf_dir+"pychemqtrc")
        parent.setProgressValue(0)

        data = {}
        P = parent.inputs.P.value
        P_kPa = P/1000
        t = cls.LineList("isotdb", Preferences)

        # Saturation line
        Hs = []
        for tdb in t:
            Hs.append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", 1))
            parent.setProgressValue(5*len(Hs)/len(t))
        data["t"] = t
        data["Hs"] = Hs

        # left limit of isow lines
        H = cls.LineList("isow", Preferences)
        th = []
        for w in H:
            if w:
                tdp = HAProps("Tdp", "P", 101.325, "W", w, "RH", 1)
                th.append(unidades.Temperature(tdp))
            else:
                tmin = Preferences.getfloat("Psychr", "isotdbStart")
                th.append(unidades.Temperature(tmin))
        data["H"] = H
        data["th"] = th

        # Humidity ratio lines
        hr = cls.LineList("isohr", Preferences)
        Hr = {}
        cont = 0
        for i in hr:
            Hr[i] = []
            for tdb in t:
                Hr[i].append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", i/100))
                cont += 1
                parent.progressBar.setValue(5+10*cont/len(hr)/len(Hs))
        data["Hr"] = Hr

        # Twb
        lines = cls.LineList("isotwb", Preferences)
        Twb = {}
        cont = 0
        for T in lines:
            ws = HAProps("W", "P", P_kPa, "RH", 1, "Tdb", T)
            H = [ws, 0]
            _td = HAProps("Tdb", "P", P_kPa, "Twb", T, "RH", 0)
            Tw = [unidades.Temperature(T), unidades.Temperature(_td)]
            cont += 1
            parent.progressBar.setValue(15+75*cont/len(lines))
            Twb[T] = (H, Tw)
        data["Twb"] = Twb

        # v
        lines = cls.LineList("isochor", Preferences)
        V = {}
        # rh = arange(1, -0.05, -0.05)
        rh = arange(0.00, 1.05, 0.05)
        for cont, v in enumerate(lines):
            w = []
            Td = []
            for r in rh:
                try:
                    w.append(HAProps("W", "P", P_kPa, "RH", r, "V", v))
                    _td = HAProps("Tdb", "P", P_kPa, "RH", r, "V", v)
                    Td.append(unidades.Temperature(_td))
                except ValueError:
                    pass
            parent.progressBar.setValue(90+10*cont/len(lines))
            V[v] = (Td, w)
        data["v"] = V

        return data


class PsyRefprop(PsyState):
    """Psychrometric state using refprop external library"""
    pass


Preferences = ConfigParser()
Preferences.read(conf_dir+"pychemqtrc")

if Preferences.getboolean("Psychr", "virial"):
    if Preferences.getboolean("Psychr", "coolprop") and \
       os.environ["CoolProp"] == "True":
        PsychroState = PsyCoolprop
    elif Preferences.getboolean("Psychr", "refprop") and \
            os.environ["refprop"] == "True":
        PsychroState = PsyRefprop
    else:
        PsychroState = PsyVirial

    # TODO: Enable other option when availables
    if PsychroState == PsyRefprop or PsychroState == PsyVirial:
        if os.environ["CoolProp"] == "True":
            PsychroState = PsyCoolprop
        else:
            PsychroState = PsyIdeal

else:
    PsychroState = PsyIdeal


if __name__ == '__main__':
    aire = PsyIdeal(tdb=40+273.15, w=0.001)
    print(aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws)

    # aire = PsyCoolprop(tdb=40+273.15, w=0.001)
    # print(aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws)

#    aire = PsyIdeal(tdb=40+273.15, HR=10)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa
#
#    aire = PsyCoolprop(tdb=40+273.15, HR=10)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa

#    aire = PsyIdeal(tdb=40+273.15, twb=20+273.15)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws
#
#    aire = PsyCoolprop(tdb=40+273.15, tdp=20+273.15)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws

#    aire = PsyCoolprop(HR=28.92, tdp=10+273.15)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws
#
#    aire = PsyIdeal(HR=28.92, tdp=10+273.15)
#    print aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws

    # from matplotlib import pyplot
    # ti = range(0, 300, 10)
    # for p in [1e5, 5e5, 1e6, 2e6, 4e6, 6e6, 8e6, 1e7]:
        # T = []
        # f = []
        # for t in ti:
            # try:
                # aire = PsyVirial(tdb=t+273.15, w=0.001, P=p)
                # T.append(t)
                # f.append(aire.f)
            # except:
                # pass

        # pyplot.plot(T, f)
    # pyplot.ylim(1, 1.5)
    # pyplot.show()

    # print(aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws)
    # aire = PsyVirial(tdb=-15+273.15, HR=100, P=101325)
    # air = IAPWS95(T=200+273.15, P=1.01325)
    # print(air.v, 1/air.rhoM)
