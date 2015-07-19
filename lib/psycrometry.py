#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module for psychrometry calculation
###############################################################################

import os
from configparser import ConfigParser

from PyQt5.QtWidgets import QApplication
from scipy.optimize import fsolve
from scipy import log, exp, arange, concatenate, linspace

try:
    from CoolProp.HumidAirProp import HAProps, HAProps_Aux
except:
    pass

from lib.config import conf_dir
from lib import unidades


def _Pbar(Z):
    """
    ASHRAE Fundamentals Handbook pag 1.1 eq. 3
    input:
        Z: altitude, m
    return
        standard atmosphere barometric pressure, Pa
    """
    return 101325.*(1-2.25577e-5*Z)**5.256


def _height(P):
    """
    Inverted _Pbar function
    input:
        standard atmosphere barometric pressure, Pa
    return
        Z: altitude, m
    """
    P_atm = P/101325.
    return 1/2.25577e-5*(1-exp(log(P_atm)/5.2559))


def _Tbar(self, Z):
    """
    ASHRAE Fundamentals Handbook pag 1.2 eq. 4
    input:
        Z: altitude, m
    return
        standard atmosphere dry bulb temperature, K
    """
    return 288.15-0.0065*Z


def _Psat(Tdb):
    """
    ASHRAE Fundamentals Handbook pag 1.2 eq. 4
    input:
        Dry bulb temperature, K
    return:
        Saturation pressure, Pa
    """
    if 173.15 <= Tdb < 273.15:
        C1 = -5674.5359
        C2 = 6.3925247
        C3 = -0.009677843
        C4 = 0.00000062215701
        C5 = 2.0747825E-09
        C6 = -9.484024E-13
        C7 = 4.1635019
        pws = exp(C1/Tdb + C2 + C3*Tdb + C4*Tdb**2 + C5*Tdb**3 + C6*Tdb**4 +
                  C7*log(Tdb))
    elif 273.15 <= Tdb <= 473.15:
        C8 = -5800.2206
        C9 = 1.3914993
        C10 = -0.048640239
        C11 = 0.000041764768
        C12 = -0.000000014452093
        C13 = 6.5459673
        pws = exp(C8/Tdb + C9 + C10*Tdb + C11*Tdb**2 + C12*Tdb**3 + C13*log(Tdb))
    else:
        raise NotImplementedError("Incoming out of bound")

    return pws


def _Tsat(Pv):
    """
    ASHRAE Fundamentals Handbook pag 1.2 eq. 4, inverted for calculate Tdb
    input:
        Saturation pressure, Pa
    return:
        Dry bulb temperature, K
    """
    Pv_min = _Psat(173.15)
    Pv_lim = _Psat(273.15)
    Pv_max = _Psat(473.15)
    if Pv_min <= Pv < Pv_lim:
        C1 = -5674.5359
        C2 = 6.3925247
        C3 = -0.009677843
        C4 = 0.00000062215701
        C5 = 2.0747825E-09
        C6 = -9.484024E-13
        C7 = 4.1635019

        def f(T):
            return Pv - exp(C1/T+C2+C3*T+C4*T**2+C5*T**3+C6*T**4+C7*log(T))
        t = fsolve(f, 250)[0]

    elif Pv_lim <= Pv <= Pv_max:
        C8 = -5800.2206
        C9 = 1.3914993
        C10 = -0.048640239
        C11 = 0.000041764768
        C12 = -0.000000014452093
        C13 = 6.5459673

        def f(T):
            return Pv - exp(C8/T+C9+C10*T+C11*T**2+C12*T**3+C13*log(T))
        t = fsolve(f, 300)[0]

    else:
        raise NotImplementedError("Incoming out of bound")

    return t


def _W(P, Tdb):
    """
    ASHRAE Fundamentals Handbook pag 1.12 eq. 22
    input:
        Dry bulb temperature, K
    return:
        Saturation pressure, Pa
    Saturation humidity calculation procedure"""
    pv = _Psat(Tdb)
    return 0.62198*pv/(P-pv)


def _h(Tdb, W):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 32
    input:
        Dry bulb temperature, K
        Humidity ratio, kg water/kg dry air
    return:
        Specific enthalpy, kJ/kg (dry air)
    """
    T_c = Tdb-273.15
    return 1.006*T_c + W*(2501 + 1.805*T_c)


def _v(P, Tdb, W):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 28
    input:
        Dry bulb temperature, K
        Humidity ratio, kg water/kg dry air
        Barometric pressure, Pa
    return:
        Specific volume, m3/kg (dry air)
    """
    P_kpa = P/1000
    return 0.2871*Tdb*(1+1.6078*W)/P_kpa


def _W_twb(tdb, twb, P):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 35-37
    input:
        dry bulb temperature, K
        wet bulb temperature, K
        barometric Pressure, Pa
    return:
        humidity ratio, kg H2O/kg air
    """
    tdb_C = tdb - 273.15
    twb_C = twb - 273.15
    Pvs = _Psat(twb)
    Ws = 0.62198*Pvs/(P-Pvs)
    if tdb >= 0:
        w = ((2501-2.381*twb_C)*Ws-1.006*(tdb-twb))/(2501+1.805*tdb_C-4.186*twb_C)
    else:
        w = ((2830-0.24*twb_C)*Ws-1.006*(tdb-twb))/(2830+1.86*tdb_C-2.1*twb_C)
    return w


def _Tdb(twb, w, P):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 35-37 inverted to calculate Tdb
    input:
        wet bulb temperature, K
        humidity ratio, kg H2O/kg air
        saturated humidity ratio, kg H2O/kg air
    return:
        dry bulb temperature, K
    """
    tw = twb-273.15
    ws = _W(P, twb)
    td = ((2501-2.381*tw)*ws+1.006*tw-w*(2501-4.186*tw))/(w*1.805+1.006)
    return td+273.15


def _Tdb_V(v, P):
    """
    Function to calculate isochor line
    input:
        specified volume, m3/kg air
        barometric pressure, Pa
    return
        dry bulb temperature, K
    """
    P_kpa = P/1000

    def f(Tdb):
        w = _W(P, Tdb)
        return v-0.2871*Tdb*(1+1.6078*w)/P_kpa
    ts = fsolve(f, 300)
    return ts


def _W_V(Td, P, v):
    """
    Function to calculate isochor line
    input:
        dry bulb temperature, K
        barometric pressure, Pa
        specified volume, m3/kg air
    return
        humidity ratio, kg H2O/kg air
    """
    P_kpa = P/1000
    return (v*P_kpa-0.2871*Td)/(0.2871*1.6078*Td)


def _tdp(Pw):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 39-40
    input:
        water vapor partial pressure, Pa
    return:
        dew point temperature, K
    """
    C14 = 6.54
    C15 = 14.526
    C16 = 0.7389
    C17 = 0.09486
    C18 = 0.4569

    alpha = log(Pw/1000.)
    Tdp1 = C14 + C15*alpha + C16*alpha**2 + C17*alpha**3 + C18*(Pw/1000.)**0.1984
    Tdp2 = 6.09 + 12.608*alpha + 0.4959*alpha**2
    if 0 <= Tdp1 <= 93:
        t = Tdp1
    elif Tdp2 < 0:
        t = Tdp2
    else:
        print(Pw, Tdp1, Tdp2)
        raise NotImplementedError("Incoming out of bound")

    return t+273.15


def _twb(tdb, W, P):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 35-37
    input:
        dry bulb temperature, K
        humidity ratio, kg H2O/kg air
        saturated humidity ratio, kg H2O/kg air
    return:
        wet bulb temperature, K
    """
    tdb_C = tdb - 273.15

    def f(twb):
        Pvs = _Psat(twb)
        Ws = 0.62198*Pvs/(P-Pvs)
        twb_C = twb - 273.15
        if tdb >= 0:
            w = ((2501.-2.326*twb_C)*Ws-1.006*(tdb_C-twb_C))/(2501.+1.86*tdb_C-4.186*twb_C)-W
        else:
            w = ((2830-0.24*twb_C)*Ws-1.006*(tdb-twb))/(2830+1.86*tdb_C-2.1*twb_C)
        return w

    twb = fsolve(f, tdb)[0]
    return twb


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

    P: mandatory input for barometric pressure, z is an alternate pressure input
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
        QApplication.translate("pychemqt", "T dry bulb, Humidity Ratio"),
        QApplication.translate("pychemqt", "T dry bulb, Relative humidity"),
        QApplication.translate("pychemqt", "T dry bulb, T wet bulb"),
        QApplication.translate("pychemqt", "T dry bulb, T dew point"),
        QApplication.translate("pychemqt", "T dew point, Relative humidity"), 
        QApplication.translate("pychemqt", "T wet bulb, Relative humidity")
        ]
    VAR_NAME = [
        ("tdb", "w"),
        ("tdb", "HR"),
        ("tdb", "twb"),
        ("tdb", "tdp"),
        ("tdp", "HR"), 
        ("twb", "HR")
        ]

#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, Enthalpy"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo seco, Densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, H absoluta"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Entalpia"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Tª rocio"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, densidad"))

    def __init__(self, **kwargs):
        self.kwargs = self.__class__.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()
            self.msg = "Solved"

    @property
    def calculable(self):
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

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
            t = arange(start, end, step)
        return t




#    def _Volume(self, T, Xa, eos=0):
#        if eos:
#            return _V_Virial(T, Xa)
#        else:
#            return _V_Ideal(T, Xa)
#
#    def _V_Ideal(self, T, Xa):
#        """volumen por unidad de masa de aire seco"""
#        return unidades.SpecificVolume(R_atml*T/self.P.atm/self.aire.M/Xa)
#
#    def Virial(self, T, Xa):
#        """Método que devuelve los coeficientes de la ecuación del virial
#        Temperatura en kelvin
#        X fracción molar de aire en la mezcla"""
#        Baa = 0.349568e2-0.668772e4/T-0.210141e7/T**2+0.924746e8/T**3
#        Caaa = 0.125975e4-0.190905e6/T+0.632467e8/T**2
#        Bww = R_atml*T*(0.7e-8-0.147184e-8*exp(1734.29/T))
#        Cwww = R_atml**2*T**2*(0.104e-14-0.335297e-17*exp(3645.09/T))**2*Bww**2
#
#        Xw = 1-Xa
#        Baw = 0.32366097e2-0.141138e5/T-0.1244535e7/T**2-0.2348789e10/T**4
#        Caaw = 0.482737e3+0.105678e6/T-0.656394e8/T**2+0.299444e11/T**3 - \
#            0.319317e13/T**4
#        Caww = -1e-6*exp(
#            -0.10728876e2+0.347802e4/T-0.383383e6/T**2+0.33406e8/T**3)
#        Bm = Xa**2*Baa+2*Xa*Xw*Baw+Xw**2*Bww
#        Cm = Xa**3*Caaa+3*Xa**2*Xw*Caaw+3*Xa*Xw**2*Caww+Xw**3*Cwww
#        return Bm, Cm
#
#    def _V_Virial(self, T, Xa):
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
            ws = 0.62198*Pvs/(P-Pvs)
            Pv = w*P/(0.62198+w)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)

        elif self.mode == 1:
            # Tdb and HR
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            Pv = Pvs*HR/100
            w = 0.62198*Pv/(P-Pv)
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)

        elif self.mode == 2:
            # Tdb and Twb
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            w = _W_twb(tdb, twb, P)
            Pv = w*P/(0.62198+w)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)

        elif self.mode == 3:
            # Tdb and Tdp
            Pv = _Psat(tdp)
            w = 0.62198*Pv/(P-Pv)
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            HR = Pv/Pvs*100
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            twb = _twb(tdb, w, P)

        elif self.mode == 4:
            # Tdp and HR
            Pv = _Psat(tdp)
            if HR:
                w = 0.62198*Pv/(P-Pv)
                Pvs = Pv/HR*100
            else:
                w = 0
                Pvs = Pv
            ws = 0.62198*Pvs/(P-Pvs)
            tdb = _Tsat(Pvs)
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            twb = _twb(tdb, w, P)

        elif self.mode == 5:
            # Tdp and Twb
#            Pv = _Psat(tdp)
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
            parent.progressBar.setValue(5*len(Hs)/len(t))
        data["t"] = t
        data["Hs"] = Hs

        # left limit of isow lines
        H = cls.LineList("isow", Preferences)
        th = []
        for w in H:
            if w:
                Pv = w*P/(0.62198+w)
                th.append(unidades.Temperature(_tdp(Pv)).config())
            else:
                tmin = Preferences.getfloat("Psychr", "isotdbStart")
                th.append(unidades.Temperature(tmin).config())
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
                parent.progressBar.setValue(5+10*cont/len(hr)/len(Hs))
        data["Hr"] = Hr

        # Twb
        lines = cls.LineList("isotwb", Preferences)
        Twb = {}
        cont = 0
        for T in lines:
            H = concatenate((arange(_W(P, T), 0, -0.001), [0.]))
            Tw = []
            for h in H:
                Tw.append(unidades.Temperature(_Tdb(T, h, P)).config())
            cont += 1
            parent.progressBar.setValue(15+75*cont/len(lines))
            Twb[T] = (H, Tw)
        data["Twb"] = Twb

        # v
        lines = cls.LineList("isochor", Preferences)
        V = {}
        for cont, v in enumerate(lines):
            ts = _Tdb_V(v, P)
            T = linspace(ts, v*P/287.055, 50)
            Td = [unidades.Temperature(ti).config() for ti in T]
            _W_V
            H = [_W_V(ti, P, v) for ti in T]
            parent.progressBar.setValue(90+10*cont/len(lines))
            V[v] = (Td, H)
        data["v"] = V

        return data


class PsyVirial(PsyState):
    """Psychrometric state using virial equation of state"""
    pass


class PsyCoolprop(PsyState):
    """Psychrometric state using coolprop external library"""

    @property
    def calculable(self):
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

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

        data = {}
        P = parent.inputs.P.value
        P_kPa = P/1000
        t = cls.LineList("isotdb", Preferences)

        # Saturation line
        Hs = []
        for tdb in t:
            Hs.append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", 1))
            parent.progressBar.setValue(5*len(Hs)/len(t))
        data["t"] = t
        data["Hs"] = Hs

        # left limit of isow lines
        H = cls.LineList("isow", Preferences)
        th = []
        for w in H:
            if w:
                tdp = HAProps("Tdp", "P", 101.325, "W", w, "RH", 1)
                th.append(unidades.Temperature(tdp).config())
            else:
                tmin = Preferences.getfloat("Psychr", "isotdbStart")
                th.append(unidades.Temperature(tmin).config())
        data["H"] = H
        data["th"] = th

        # Humidity ratio lines
        hr = cls.LineList("isohr", Preferences)
        Hr = {}
        cont = 0
        for i in hr:
            Hr[i] = []
            for tdb in t:
                Hr[i].append(HAProps("W", "P", P_kPa, "Tdb", tdb, "RH", i/100.))
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
            Tw = [unidades.Temperature(T).config(),
                  unidades.Temperature(HAProps("Tdb", "P", P_kPa, "Twb", T, "RH", 0)).config()]
            cont += 1
            parent.progressBar.setValue(15+75*cont/len(lines))
            Twb[T] = (H, Tw)
        data["Twb"] = Twb

        # v
        lines = cls.LineList("isochor", Preferences)
        V = {}
        rh = arange(1,-0.05,-0.05)
        for cont, v in enumerate(lines):
            w = []
            Td = []
            for r in rh:
                w.append(HAProps("W", "P", P_kPa, "RH", r, "V", v))
                Td.append(unidades.Temperature(HAProps("Tdb", "P", P_kPa, "RH", r, "V", v)).config())
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
else:
    PsychroState = PsyIdeal


if __name__ == '__main__':
    aire = PsyIdeal(tdb=40+273.15, w=0.001)
    print(aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws)

    aire = PsyCoolprop(tdb=40+273.15, w=0.001)
    print(aire.tdb.C, aire.twb.C, aire.tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa, aire.ws)

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
