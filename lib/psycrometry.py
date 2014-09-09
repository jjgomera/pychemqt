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
from lib.corriente import Corriente


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

    def altura(self):
        return unidades.Length(1/2.25577e-5*(1-exp(log(self.P.atm)/5.2559)))
#        return unidades.Length((8430.153*log(101325/self.P))/(1+0.095*log(101325/self.P)))

    def Presion(self, Z):
        return unidades.Pressure((1-2.25577e-5*Z)**5.256, "atm")

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
    """Class to model a psychrometric state with properties"""
    def __init__(self, eos=0, P=101325, **kwargs):
        """
        eos: Index of equation of state to use
            0 - Ideal gas
            1 - Real gas
        P: Pressure
        kwargs definition parameters:
            tdp: dew-point temperature
            tdb: dry-bulb temperature
            twb: web-bulb temperature
            w: Humidity Ratio [kg water/kg dry air]
            HR: Relative humidity
            h: Mixture enthalpy
            v: Mixture specified volume
        """
        self.Water = Componente(62)
        self.Air = Componente(475)
        self.P = unidades.Pressure(P)
        
        tdp = kwargs.get("tdp", 0)
        tdb = kwargs.get("tdb", 0)
        twb = kwargs.get("twb", 0)
        w = kwargs.get("w", None)
        HR = kwargs.get("HR", None)
        h = kwargs.get("h", None)
        v = kwargs.get("v", 0)
        
        if tdb and w is not None:
            self.Tdb = tdb
            self.w = w
            self.Xa = 1/(1+w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.h = self._h(tdb, w)
            self.Hs = self._HS(tdb)
            self.HR = w/self.Hs*100
            self.Twb = self.AireHumedo.Tw(tdb, w)
            self.Tdp = self.AireHumedo.Tdp(w)
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)
        elif tdb and HR is not None:
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
        elif tdb and tdp:
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
        elif tdb and twb:
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
        elif tdp and HR is not None:
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
        elif tdb and h:
            self.Tdb = tdb
            self.h = unidades.Enthalpy(h, "kJkg")
            f = lambda h: self.AireHumedo.Entalpia(self.Tdb, h).kJkg-h
            self.w = fsolve(f, 0.001)
            self.Hs = self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb = self.AireHumedo.Tw(tdb, self.w)
            self.HR = self.w/self.Hs*100
            self.Xa = 1/(1+self.w*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw = 1-self.Xa
            self.V = self.AireHumedo.Volumen(tdb, self.Xa)
            self.rho = unidades.Density(1/self.V)

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
    aire = Psy_state(tdb=300, w=0.1)
    print aire.__dict__
    
