#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module for psychrometry calculation
###############################################################################

from scipy.optimize import fsolve
from scipy import log, exp

from compuestos import Componente
from physics import R_atml
from lib import unidades
from lib.iapws import _PSat_T, _TSat_P, _Sublimation_Pressure, _Melting_Pressure
#from lib.corriente import Corriente



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


def _h(Tdb, W):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 32
    input:
        Dry bulb temperature, K
        Humidity ratio, kg water/kg dry air
    return:
        Specific enthalpy, kJ/kg (dry air)
    """
    T_c= Tdb-273.15
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


def _W_twb(tdb, twb, Ws):
    """
    ASHRAE Fundamentals Handbook pag 1.13 eq. 35-37
    input:
        dry bulb temperature, K
        wet bulb temperature, K
        barometric pressure, Pa
    return:
        humidity ratio, kg H2O/kg air
    """
    tdb_C = tdb - 273.15
    twb_C = twb - 273.15
    if tdb >= 0:
        w = ((2501-2.381*twb_C)*Ws-1.006*(tdb-twb))/(2501+1.805*tdb_C-4.186*twb_C)
    else:                                   # Equation 37, p6.9
        w = ((2830-0.24*twb_C)*Ws-1.006*(tdb-twb))/(2830+1.86*tdb_C-2.1*twb_C)
    return w


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



class Psychrometry(object):
    def __init__(self, P=1):
        self.agua = Componente(62)
        self.aire = Componente(475)
        self.P = unidades.Pressure(P, "atm")

        self.Volumen(300, 0.999)

    def Virial(self, T, Xa):
        """Método que devuelve los coeficientes de la ecuación del virial
        Temperatura en kelvin
        X fracción molar de aire en la mezcla"""
        Baa = 0.349568e2-0.668772e4/T-0.210141e7/T**2+0.924746e8/T**3
        Caaa = 0.125975e4-0.190905e6/T+0.632467e8/T**2
        Bww = R_atml*T*(0.7e-8-0.147184e-8*exp(1734.29/T))
        Cwww = R_atml**2*T**2*(0.104e-14-0.335297e-17*exp(3645.09/T))**2*Bww**2

        Xw = 1-Xa
        Baw = 0.32366097e2-0.141138e5/T-0.1244535e7/T**2-0.2348789e10/T**4
        Caaw = 0.482737e3+0.105678e6/T-0.656394e8/T**2+0.299444e11/T**3 - \
            0.319317e13/T**4
        Caww = -1e-6*exp(
            -0.10728876e2+0.347802e4/T-0.383383e6/T**2+0.33406e8/T**3)
        Bm = Xa**2*Baa+2*Xa*Xw*Baw+Xw**2*Bww
        Cm = Xa**3*Caaa+3*Xa**2*Xw*Caaw+3*Xa*Xw**2*Caww+Xw**3*Cwww
        return Bm, Cm

    def Volumen(self, T, Xa):
        """volumen por unidad de masa de aire seco"""
        # FIXME: Don't work, for now use ideal gas equation
#        Bm, Cm=self.Virial(T, Xa)
#        vm=roots([1, -R_atml*T/self.P.atm, -R_atml*T*Bm/self.P.atm,-R_atml*T*Cm/self.P.atm])
#        if vm[0].imag==0.0:
#            v=vm[0].real
#        else:
#            v=vm[2].real
#        return unidades.SpecificVolume(v/self.aire.M/Xa)
        return unidades.SpecificVolume(R_atml*T/self.P.atm/self.aire.M/Xa)


    def Humedad_Absoluta(self, T):
        pv = self.agua.Pv(T)
        return self.agua.M*pv/(self.aire.M*(self.P-pv))

    def Tdp(self, w):
        Pv = w*self.aire.M*self.P/(w*self.aire.M+self.agua.M)
        if Pv:
            T = 19.5322+13.6626*log(0.00145*Pv)+1.17678*log(0.00145*Pv)**2 - \
                0.189693*log(0.00145*Pv)**3+0.087453*log(0.00145*Pv)**4 - \
                0.0174053*log(0.00145*Pv)**5+0.00214768*log(0.00145*Pv)**6 - \
                0.138343e-3*log(0.00145*Pv)**7+0.38e-5*log(0.00145*Pv)**8+255.38
        else:
            T = 19.5322+255.38
        return unidades.Temperature(T)

    def Calor_Especifico_Humedo(self, T, w):
        return unidades.SpecificHeat(self.aire.Cp_Gas_DIPPR(T)+self.agua.Cp_Gas_DIPPR(T)*w)

    def Isoentalpica(self, T, w):
        Hv = self.agua.Hv_DIPPR(T)
        Hs = self.Humedad_Absoluta(T)
        Cs = self.Calor_Especifico_Humedo(T, w)
        return unidades.Temperature(T+(Hs-w)*Hv/Cs)

    def Tw(self, Td, w):
        Hv = self.agua.Hv_DIPPR(Td)
        Cs = self.Calor_Especifico_Humedo(Td, w)

        def f(Tw):
            return self.Humedad_Absoluta(Tw)-w-Cs/Hv*(Td-Tw)
        Tw = fsolve(f, Td)[0]
        return unidades.Temperature(Tw)

    def Entalpia(self, Td, w):
        h = 1.006*Td.C+w*(2501+1.805*Td.C)
        return unidades.Enthalpy(h, "kJkg")

    def Isocora_H(self, T, v):
        V = unidades.SpecificVolume(v)
        return (self.P.atm*V*self.aire.M/R_atml/T-1)*self.agua.M/self.aire.M

    def Isocora_T(self, w, v):
        V = unidades.SpecificVolume(v)
        return unidades.Temperature(self.P.atm*V.lg*self.aire.M/R_atml /
                                    (1+w*self.aire.M/self.agua.M))

    def definirPunto(self, modo, tdb=None, twb=None, tdp=None, w=None,
                     HR=None, V=None, h=None):
        """
        modo: datos conocidos
            0   -   tdb, w
            1   -   tdb, HR
            2   -   tdb, tdp
            3   -   tdb, twb
            4   -   tdp, HR"""
        punto = Punto_Psicrometrico()
        if modo == 0:
            punto.P = self.P
            punto.Tdb = unidades.Temperature(tdb)
            punto.w = unidades.Dimensionless(w, txt="kgw/kgda")
            punto.Hs = self.Humedad_Absoluta(tdb)
            punto.Twb = self.Tw(tdb, w)
            punto.HR = unidades.Dimensionless(w/punto.Hs*100, txt="%")
            punto.Tdp = self.Tdp(w)
            punto.h = self.Entalpia(tdb, w)
            punto.Xa = 1/(1+w*self.aire.M/self.agua.M)
            punto.Xw = 1-punto.Xa
            punto.V = self.Volumen(tdb, punto.Xa)
            punto.rho = unidades.Density(1/punto.V)
        elif modo == 1:
            punto.P = self.P
            punto.Tdb = tdb
            punto.Hs = self.Humedad_Absoluta(tdb)
            punto.HR = HR
            punto.w = punto.Hs*HR/100
            punto.Twb = self.Tw(tdb, punto.w)
            punto.Tdp = self.Tdp(punto.w)
            punto.h = self.Entalpia(tdb, punto.w)
            punto.Xa = 1/(1+punto.w*self.aire.M/self.agua.M)
            punto.Xw = 1-punto.Xa
            punto.V = self.Volumen(tdb, punto.Xa)
            punto.rho = unidades.Density(1/punto.V)
        elif modo == 2:
            punto.P = self.P
            punto.Tdb = tdb
            punto.Tdp = tdp
            punto.w = self.Humedad_Absoluta(tdp)
            punto.Hs = self.Humedad_Absoluta(tdb)
            punto.Twb = self.Tw(tdb, punto.w)
            punto.HR = punto.w/punto.Hs*100
            punto.h = self.Entalpia(tdb, punto.w)
            punto.Xa = 1/(1+punto.w*self.aire.M/self.agua.M)
            punto.Xw = 1-punto.Xa
            punto.V = self.Volumen(tdb, punto.Xa)
            punto.rho = unidades.Density(1/punto.V)
        elif modo == 3:
            punto.P = self.P
            punto.Tdb = tdb
            punto.Twb = twb
            punto.Hs = self.Humedad_Absoluta(tdb)
            Hv = self.agua.Hv_DIPPR(tdb)

            def f(w):
                return self.Humedad_Absoluta(twb)-w - \
                    self.Calor_Especifico_Humedo(tdb, w)/Hv*(tdb-twb)
            punto.w = fsolve(f, 0.001)
            punto.Tdp = self.Tdp(punto.w)
            punto.HR = punto.w/punto.Hs*100
            punto.h = self.Entalpia(tdb, punto.w)
            punto.Xa = 1/(1+punto.w*self.aire.M/self.agua.M)
            punto.Xw = 1-punto.Xa
            punto.V = self.Volumen(tdb, punto.Xa)
            punto.rho = unidades.Density(1/punto.V)
        elif modo == 4:
            punto.P = self.P
            punto.Tdp = tdp
            punto.HR = HR
            punto.w = self.Humedad_Absoluta(tdp)
            punto.Hs = punto.w/punto.HR*100
            punto.Tdb = self.Tdp(punto.Hs)
            punto.Twb = self.Tw(punto.Tdb, punto.w)
            punto.h = self.Entalpia(punto.Tdb, punto.w)
            punto.Xa = 1/(1+punto.w*self.aire.M/self.agua.M)
            punto.Xw = 1-punto.Xa
            punto.V = self.Volumen(punto.Tdb, punto.Xa)
            punto.rho = unidades.Density(1/punto.V)
        elif modo == 5:
            self.Tdb = tdb
            self.h = unidades.Enthalpy(h, "kJkg")
            f = lambda w: self.Entalpia(self.Tdb, w).kJkg-h
            self.w = fsolve(f, 0.001)
            self.Hs = self.Humedad_Absoluta(tdb)
            self.Twb = self.Tw(tdb, self.w)
            self.HR = self.w/self.Hs*100
            self.Xa = 1/(1+self.w*self.aire.M/self.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        return punto


class Psy_state(object):
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
        1 - tdb, w
        2 - tdb, HR
        3 - tdb, twb
        4 - tdb, tdp
        5 - tdp, HR
        6 - tdp, twb
        7 - twb, w
        
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

    def __init__(self, **kwargs):
        self.kwargs = Psy_state.kwargs.copy()
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
        
        self.mode = 0
        if tdb and w is not None:
            self.mode = 1
        elif tdb and HR is not None:
            self.mode = 2
        elif tdb and twb:
            self.mode = 3
        elif tdb and tdp:
            self.mode = 4
        elif tdp and HR is not None:
            self.mode = 5
            
        return bool(self.mode)
        
    def calculo(self):
        # Barometric pressure calculation
        if self.kwargs["P"]:
            P = self.kwargs["P"]
        elif self.kwargs["z"]:
            P = _Pbar(self.kwargs["z"])
        else:
            P = 101325.
        self.P = unidades.Pressure(P)

        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

        if self.mode == 1:
            # Tdb and w
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            Pv = w*P/(0.62198+w)
            HR = Pv/Pvs*100
            mu = w/ws
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)

        elif self.mode == 2:
            # Tdb and HR
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            Pv = Pvs*HR/100
            w = 0.62198*Pv/(P-Pv)
            mu = w/ws
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
            twb = _twb(tdb, w, P)
        
        elif self.mode == 3:
            # Tdb and Twb
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            w = _W_twb(tdb, twb, ws)
            Pv = w*P/(0.62198+w)
            HR = Pv/Pvs*100
            mu = w/ws
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            tdp = _tdp(Pv)
        
        elif self.mode == 4:
            # Tdb and Tdp
            Pv = _Psat(tdp)
            w = 0.62198*Pv/(P-Pv)
            Pvs = _Psat(tdb)
            ws = 0.62198*Pvs/(P-Pvs)
            mu = w/ws
            v = _v(P, tdb, w)
            h = _h(tdb, w)
            twb = _twb(tdb, w, P)
            
        elif self.mode == 5:
            # Tdp and HR
            pass
            
            
        elif self.mode == 6:
            pass

        elif self.mode == 7:
            pass

#        elif tdp and HR is not None:
#            self.Tdp = tdp
#            self.HR = HR
#            if HR:
#                self.w = self.AireHumedo.Humedad_Absoluta(tdp)
#                self.Hs = self.w/self.HR*100
#            else:
#                self.w = 0
#            self.Tdb = self.AireHumedo.Tdp(self.Hs)
#            self.Twb = self.AireHumedo.Tw(self.Tdb, self.w)
#            self.h = self.AireHumedo.Entalpia(self.Tdb, self.w)
#            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
#            self.Xw = 1-self.Xa
#            self.V = self.AireHumedo.Volumen(self.Tdb, self.Xa)
#            self.rho = unidades.Density(1/self.V)

            
        self.Tdp = unidades.Temperature(tdp)
        self.Tdb = unidades.Temperature(tdb)
        self.Twb = unidades.Temperature(twb)
        self.P = unidades.Pressure(P)
        self.Pvs = unidades.Pressure(Pvs)
        self.Pv = unidades.Pressure(Pv)
        self.ws = unidades.Dimensionless(ws, txt="kgw/kgda")
        self.w = unidades.Dimensionless(w, txt="kgw/kgda")
        self.HR = unidades.Dimensionless(HR, txt="%")
        self.mu = unidades.Dimensionless(mu)
        self.v = unidades.SpecificVolume(v)
        self.h = unidades.Enthalpy(h, "kJkg")
        self.Xa = 1/(1+self.w/0.62198)
        self.Xw = 1-self.Xa

    @property
    def corriente(self, caudal):
        corriente = Corriente(T=self.Twb, P=self.P, caudalMasico=caudal, ids=[62, 475], fraccionMolar=[self.Xw, self.Xa])
        return corriente

            
            
#        elif tdp and HR is not None:
#            self.Tdp = tdp
#            self.HR = HR
#            if HR:
#                self.w = self.AireHumedo.Humedad_Absoluta(tdp)
#                self.Hs = self.w/self.HR*100
#            else:
#                self.w = 0
#            self.Tdb = self.AireHumedo.Tdp(self.Hs)
#            self.Twb = self.AireHumedo.Tw(self.Tdb, self.w)
#            self.h = self.AireHumedo.Entalpia(self.Tdb, self.w)
#            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
#            self.Xw = 1-self.Xa
#            self.V = self.AireHumedo.Volumen(self.Tdb, self.Xa)
#            self.rho = unidades.Density(1/self.V)
#        elif tdb and h:
#            self.Tdb = tdb
#            self.h = unidades.Enthalpy(h, "kJkg")
#            f = lambda h: self.AireHumedo.Entalpia(self.Tdb, h).kJkg-h
#            self.w = fsolve(f, 0.001)
#            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
#            self.Twb = self.AireHumedo.Tw(tdb, self.w)
#            self.HR = self.w/self.Hs*100
#            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
#            self.Xw = 1-self.Xa
#            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
#            self.rho = unidades.Density(1/self.V)

    def _Volume(self, T, Xa, eos=0):
        if eos:
            return _V_Virial(T, Xa)
        else:
            return _V_Ideal(T, Xa)
        
    def _V_Ideal(self, T, Xa):
        """volumen por unidad de masa de aire seco"""
        return unidades.SpecificVolume(R_atml*T/self.P.atm/self.aire.M/Xa)
        
    def _V_Virial(self, T, Xa):
        """volumen por unidad de masa de aire seco"""
        # FIXME: Don't work, for now use ideal gas equation
        Bm, Cm=self.Virial(T, Xa)
        vm=roots([1, -R_atml*T/self.P.atm, -R_atml*T*Bm/self.P.atm,-R_atml*T*Cm/self.P.atm])
        if vm[0].imag==0.0:
            v=vm[0].real
        else:
            v=vm[2].real
        return unidades.SpecificVolume(v/self.aire.M/Xa)
    
    def _h(self, Td, w):
        """Enthalpy calculation procedure"""
        cp_air = self.Air.Cp_Gas_DIPPR(Td)
        cp_water = self.Water.Cp_Gas_DIPPR(Td)
        h = (Td-273.15)*(cp_air+cp_water*w)
        return unidades.Enthalpy(h, "kJkg")
    
    def _HS(self, Td):
        """Saturation humidity calculation procedure"""
        pv = self._Ps(Td)
        return self.Water.M*pv/(self.Air.M*(self.P-pv))

    def _Ps(self, Td):
        """Saturation pressure calculation procedure"""
        if Td < 273.15:
            return _Sublimation_Pressure(Td)
        else:
            return _PSat_T(Td)

    def _Tw(self, Td, w):
        Hv = self.agua.Hv_DIPPR(Td)
        Cs = self.Calor_Especifico_Humedo(Td, w)

        def f(Tw):
            return self.Humedad_Absoluta(Tw)-w-Cs/Hv*(Td-Tw)
        Tw = fsolve(f, Td)
        return unidades.Temperature(Tw)

    
class Punto_Psicrometrico(object):
    def __init__(self, P=1, caudal=None, modo=None, tdb=None, twb=None,
                 tdp=None, w=None, HR=None, V=None, h=None):
        self.AireHumedo = Psychrometry(P)
        self.P = unidades.Pressure(P, "atm")
        if caudal:
            self.caudal = unidades.MassFlow(caudal, "kgh")
        else:
            self.caudal = caudal
        if not modo:
            if tdb and w is not None:
                self.definirPunto(0, tdb=unidades.Temperature(tdb), w=w)
            elif tdb and HR is not None:
                self.definirPunto(1, tdb=unidades.Temperature(tdb), HR=HR)
            elif tdb and tdp:
                self.definirPunto(2, tdb=unidades.Temperature(tdb), tdp=unidades.Temperature(tdp))
            elif tdb and twb:
                self.definirPunto(3, tdb=unidades.Temperature(tdb), twb=unidades.Temperature(twb))
            elif tdp and HR is not None:
                self.definirPunto(4, tdp=unidades.Temperature(tdp), HR=HR)
            elif tdb and h:
                self.definirPunto(5, tdb=unidades.Temperature(tdb), h=h)
            else:
                self.P = 0
                self.w = 0
                self.HR = 0
                self.Tdb = 0
                self.Twb = 0
                self.Tdp = 0
                self.V = 0
                self.rho = 0
                self.h = 0
                self.Xa = 0
                self.Xw = 0
                self.modo = -1
                self.Hs = 0

    def definirPunto(self, modo, tdb=None, twb=None, tdp=None, w=None, HR=None,
                     V=None, h=None):
        """
        modo: datos conocidos
            0   -   tdb, w
            1   -   tdb, HR
            2   -   tdb, tdp
            3   -   tdb, twb
            4   -   tdp, HR
            5   -   tdb, h
            """
        self.modo = modo
        if modo == 0:
            self.Tdb = tdb
            self.w = w
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb = self.AireHumedo.Tw(tdb, w)
            self.HR = w/self.Hs*100
            self.Tdp = self.AireHumedo.Tdp(w)
            self.h = self.AireHumedo.Entalpia(tdb, w)
            self.Xa = 1/(1+w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif modo == 1:
            self.Tdb = tdb
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            self.HR = HR
            self.w = self.Hs*HR/100
            self.Twb = self.AireHumedo.Tw(tdb, self.w)
            self.Tdp = self.AireHumedo.Tdp(self.w)
            self.h = self.AireHumedo.Entalpia(tdb, self.w)
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif modo == 2:
            self.Tdb = tdb
            self.Tdp = tdp
            self.w = self.AireHumedo.Humedad_Absoluta(tdp)
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb = self.AireHumedo.Tw(tdb, self.w)
            self.HR = self.w/self.Hs*100
            self.h = self.AireHumedo.Entalpia(tdb, self.w)
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif modo == 3:
            self.Tdb = tdb
            self.Twb = twb
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            Hv = self.AireHumedo.agua.Hv_DIPPR(tdb)

            def f(w):
                return self.AireHumedo.Humedad_Absoluta(twb)-w - \
                    self.AireHumedo.Calor_Especifico_Humedo(tdb, w)/Hv*(tdb-twb)
            self.w = fsolve(f, 0.001)
            self.Tdp = self.AireHumedo.Tdp(self.w)
            self.HR = self.w/self.Hs*100
            self.h = self.AireHumedo.Entalpia(tdb, self.w)
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif modo == 4:
            self.Tdp = tdp
            self.HR = HR
            if HR:
                self.w = self.AireHumedo.Humedad_Absoluta(tdp)
                self.Hs = self.w/self.HR*100
            else:
                self.w = 0
            self.Tdb = self.AireHumedo.Tdp(self.Hs)
            self.Twb = self.AireHumedo.Tw(self.Tdb, self.w)
            self.h = self.AireHumedo.Entalpia(self.Tdb, self.w)
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(self.Tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif modo == 5:
            self.Tdb = tdb
            self.h = unidades.Enthalpy(h, "kJkg")
            f = lambda w: self.AireHumedo.Entalpia(self.Tdb, w).kJkg-h
            self.w = fsolve(f, 0.001)
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb = self.AireHumedo.Tw(tdb, self.w)
            self.HR = self.w/self.Hs*100
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)

    @property
    def corriente(self):
        corriente = Corriente(T=self.Twb, P=self.P, caudalMasico=self.caudal, fraccionMolar=[self.Xw, self.Xa])
        return corriente

if __name__ == '__main__':
#    aireHumedo=Psychrometry(1)
#    aireHumedo.Tw(298, 0.01)

#    aire = Punto_Psicrometrico(caudal=1000, tdb=300, w=0.1)
#    print aire.__dict__
#    aire = Psy_state(tdb=300, w=0.1)
    aire = Psy_state(tdb=40+273.15, HR=10)
    print aire.Tdb.C, aire.Twb.C, aire.Tdp.C, aire.w, aire.v, aire.mu, aire.h.kJkg, aire.Pv.Pa
#    print aire.__dict__
    
