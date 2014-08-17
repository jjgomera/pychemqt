#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from scipy.optimize import fsolve, leastsq, brentq
from scipy.special import erf
from scipy.linalg import det
from scipy import roots, power, log, sqrt, log10, exp, sin, r_, real, imag, zeros, transpose
from pylab import triu, plot, grid, show, figtext, xlabel, ylabel
from PyQt4.QtGui import QApplication

from compuestos import Componente
from bip import SRK, PR, BWRS
from physics import R_atml, R_cal, R, factor_acentrico_octano
from lib import unidades, config, eos, mEoS
from lib.corriente import Corriente


class Psychrometry(object):
    def __init__(self, P=1):
        self.agua=Componente(62)
        self.aire=Componente(475)
        self.P=unidades.Pressure(P, "atm")
        
        self.Volumen(300, 0.999)
    
    def Virial(self, T, Xa):
        """Método que devuelve los coeficientes de la ecuación del virial
        Temperatura en kelvin
        X fracción molar de aire en la mezcla"""
        Baa=0.349568e2-0.668772e4/T-0.210141e7/T**2+0.924746e8/T**3
        Caaa=0.125975e4-0.190905e6/T+0.632467e8/T**2
        Bww=R_atml*T*(0.7e-8-0.147184e-8*exp(1734.29/T))
        Cwww=R_atml**2*T**2*(0.104e-14-0.335297e-17*exp(3645.09/T))**2*Bww**2
        
        Xw=1-Xa
        Baw=0.32366097e2-0.141138e5/T-0.1244535e7/T**2-0.2348789e10/T**4
        Caaw=0.482737e3+0.105678e6/T-0.656394e8/T**2+0.299444e11/T**3-0.319317e13/T**4
        Caww=-1e-6*exp(-0.10728876e2+0.347802e4/T-0.383383e6/T**2+0.33406e8/T**3)
        Bm=Xa**2*Baa+2*Xa*Xw*Baw+Xw**2*Bww
        Cm=Xa**3*Caaa+3*Xa**2*Xw*Caaw+3*Xa*Xw**2*Caww+Xw**3*Cwww
        return Bm, Cm
        
    def Volumen(self, T, Xa):
        """volumen por unidad de masa de aire seco"""
        #FIXME: Usaré de momento gas ideal ya que la ecuación del virial no da valores coherentes
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

    def Presion(self, altura):
        return unidades.Pressure((1-2.25577e-5*altura)**5.2559, "atm")

    def Humedad_Absoluta(self, T):
        pv=self.agua.Pv(T)
        return self.agua.M*pv/(self.aire.M*(self.P-pv))
        
    def Tdp(self, H):
        Pv=H*self.aire.M*self.P/(H*self.aire.M+self.agua.M)
        if Pv:
            return unidades.Temperature(19.5322+13.6626*log(0.00145*Pv)+1.17678*log(0.00145*Pv)**2-0.189693*log(0.00145*Pv)**3+0.087453*log(0.00145*Pv)**4-0.0174053*log(0.00145*Pv)**5+0.00214768*log(0.00145*Pv)**6-0.138343e-3*log(0.00145*Pv)**7+0.38e-5*log(0.00145*Pv)**8+255.38)
        else:
            return unidades.Temperature(19.5322+255.38)
        
    def Calor_Especifico_Humedo(self, T, H):
        return unidades.SpecificHeat(self.aire.Cp_Gas_DIPPR(T)+self.agua.Cp_Gas_DIPPR(T)*H)
        
    def Isoentalpica(self, T, H):
        Hv=self.agua.Hv_DIPPR(T)
        Hs=self.Humedad_Absoluta(T)
        Cs=self.Calor_Especifico_Humedo(T, H)
        return unidades.Temperature(T+(Hs-H)*Hv/Cs)
        
    def Tw(self, Td, H):
        Hv=self.agua.Hv_DIPPR(Td)
        Cs=self.Calor_Especifico_Humedo(Td, H)
        def f(Tw):
            return self.Humedad_Absoluta(Tw)-H-Cs/Hv*(Td-Tw)
        Tw=fsolve(f, Td)
        return unidades.Temperature(Tw)

    def Entalpia(self, Td, H):
        h=1.006*Td.C+H*(2501+1.805*Td.C)
        return unidades.Enthalpy(h, "kJkg")

    def Isocora_H(self, T, v):
        V=unidades.SpecificVolume(v)
        return (self.P.atm*v*self.aire.M/R_atml/T-1)*self.agua.M/self.aire.M
        
    def Isocora_T(self, H, v):
        V=unidades.SpecificVolume(v)
        return unidades.Temperature(self.P.atm*V.lg*self.aire.M/R_atml/(1+H*self.aire.M/self.agua.M))
        
    def definirPunto(self, modo, tdb=None, twb=None, tdp=None, H=None, HR=None, volumen=None, entalpia=None):
        """
        modo: datos conocidos
            0   -   tdb, H
            1   -   tdb, HR
            2   -   tdb, tdp
            3   -   tdb, twb
            4   -   tdp, HR"""
        punto=Punto_Psicrometrico()
        if modo==0:
            punto.P=self.P
            punto.Tdb=tdb
            punto.H=H
            punto.Hs=self.Humedad_Absoluta(tdb)
            punto.Twb=self.Tw(tdb, H)
            punto.HR=H/punto.Hs*100
            punto.Tdp=self.Tdp(H)
            punto.entalpia=self.Entalpia(tdb, H)
            punto.Xa=1/(1+H*self.aire.M/self.agua.M)
            punto.Xw=1-punto.Xa
            punto.volumen=self.Volumen(tdb, punto.Xa)
            punto.densidad=unidades.Density(1/punto.volumen)
        elif modo==1:
            punto.P=self.P
            punto.Tdb=tdb
            punto.Hs=self.Humedad_Absoluta(tdb)
            punto.HR=HR
            punto.H=punto.Hs*HR/100
            punto.Twb=self.Tw(tdb, punto.H)
            punto.Tdp=self.Tdp(punto.H)
            punto.entalpia=self.Entalpia(tdb, punto.H)
            punto.Xa=1/(1+punto.H*self.aire.M/self.agua.M)
            punto.Xw=1-punto.Xa
            punto.volumen=self.Volumen(tdb, punto.Xa)
            punto.densidad=unidades.Density(1/punto.volumen)
        elif modo==2:
            punto.P=self.P
            punto.Tdb=tdb
            punto.Tdp=tdp
            punto.H=self.Humedad_Absoluta(tdp)
            punto.Hs=self.Humedad_Absoluta(tdb)
            punto.Twb=self.Tw(tdb, punto.H)
            punto.HR=punto.H/punto.Hs*100
            punto.entalpia=self.Entalpia(tdb, punto.H)
            punto.Xa=1/(1+punto.H*self.aire.M/self.agua.M)
            punto.Xw=1-punto.Xa
            punto.volumen=self.Volumen(tdb, punto.Xa)
            punto.densidad=unidades.Density(1/punto.volumen)
        elif modo==3:
            punto.P=self.P
            punto.Tdb=tdb
            punto.Twb=twb
            punto.Hs=self.Humedad_Absoluta(tdb)
            Hv=self.agua.Hv_DIPPR(tdb)
            def f(H):
                return self.Humedad_Absoluta(twb)-H-self.Calor_Especifico_Humedo(tdb, H)/Hv*(tdb-twb)
            punto.H= fsolve(f, 0.001)
            punto.Tdp=self.Tdp(punto.H)
            punto.HR=punto.H/punto.Hs*100
            punto.entalpia=self.Entalpia(tdb, punto.H)
            punto.Xa=1/(1+punto.H*self.aire.M/self.agua.M)
            punto.Xw=1-punto.Xa
            punto.volumen=self.Volumen(tdb, punto.Xa)
            punto.densidad=unidades.Density(1/punto.volumen)
        elif modo==4:
            punto.P=self.P
            punto.Tdp=tdp
            punto.HR=HR
            punto.H=self.Humedad_Absoluta(tdp)
            punto.Hs=punto.H/punto.HR*100
            punto.Tdb=self.Tdp(punto.Hs)
            punto.Twb=self.Tw(punto.Tdb, punto.H)
            punto.entalpia=self.Entalpia(punto.Tdb, punto.H)
            punto.Xa=1/(1+punto.H*self.aire.M/self.agua.M)
            punto.Xw=1-punto.Xa
            punto.volumen=self.Volumen(punto.Tdb, punto.Xa)
            punto.densidad=unidades.Density(1/punto.volumen)
        elif modo==5:
            self.Tdb=tdb
            self.entalpia=unidades.Enthalpy(entalpia, "kJkg")
            f=lambda h: self.Entalpia(self.Tdb, h).kJkg-entalpia
            self.H=fsolve(f, 0.001)
            self.Hs=self.Humedad_Absoluta(tdb)
            self.Twb=self.Tw(tdb, self.H)
            self.HR=self.H/self.Hs*100
            self.Xa=1/(1+self.H*self.aire.M/self.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        return punto



class Punto_Psicrometrico(object):
    def __init__(self, P=1, caudal=None, modo=None, tdb=None, twb=None, tdp=None, H=None, HR=None, volumen=None, entalpia=None):
        self.AireHumedo=Psychrometry(P)
        self.P=unidades.Pressure(P, "atm")
        if caudal:
            self.caudal=unidades.MassFlow(caudal, "kgh")
        else: self.caudal=caudal
        if not modo:
            if tdb and H!=None:
                self.definirPunto(0, tdb=unidades.Temperature(tdb), H=H)
            elif tdb and HR!=None:
                self.definirPunto(1, tdb=unidades.Temperature(tdb), HR=HR)
            elif tdb and tdp:
                self.definirPunto(2, tdb=unidades.Temperature(tdb), tdp=unidades.Temperature(tdp))
            elif tdb and twb:
                self.definirPunto(3, tdb=unidades.Temperature(tdb), twb=unidades.Temperature(twb))
            elif tdp and HR!=None:
                self.definirPunto(4, tdp=unidades.Temperature(tdp), HR=HR)
            elif tdb and entalpia:
                self.definirPunto(5, tdb=unidades.Temperature(tdb), entalpia=entalpia)
            else:
                self.P=0
                self.H=0
                self.HR=0
                self.Tdb=0
                self.Twb=0
                self.Tdp=0
                self.volumen=0
                self.densidad=0
                self.entalpia=0
                self.Xa=0
                self.Xw=0
                self.modo=-1
                self.Hs=0
        
    def definirPunto(self, modo, tdb=None, twb=None, tdp=None, H=None, HR=None, volumen=None, entalpia=None):
        """
        modo: datos conocidos
            0   -   tdb, H
            1   -   tdb, HR
            2   -   tdb, tdp
            3   -   tdb, twb
            4   -   tdp, HR
            5   -   tdb, entalpia
            """
        self.modo=modo
        if modo==0:
            self.Tdb=tdb
            self.H=H
            self.Hs=self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb=self.AireHumedo.Tw(tdb, H)
            self.HR=H/self.Hs*100
            self.Tdp=self.AireHumedo.Tdp(H)
            self.entalpia=self.AireHumedo.Entalpia(tdb, H)
            self.Xa=1/(1+H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        elif modo==1:
            self.Tdb=tdb
            self.Hs=self.AireHumedo.Humedad_Absoluta(tdb)
            self.HR=HR
            self.H=self.Hs*HR/100
            self.Twb=self.AireHumedo.Tw(tdb, self.H)
            self.Tdp=self.AireHumedo.Tdp(self.H)
            self.entalpia=self.AireHumedo.Entalpia(tdb, self.H)
            self.Xa=1/(1+self.H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        elif modo==2:
            self.Tdb=tdb
            self.Tdp=tdp
            self.H=self.AireHumedo.Humedad_Absoluta(tdp)
            self.Hs=self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb=self.AireHumedo.Tw(tdb, self.H)
            self.HR=self.H/self.Hs*100
            self.entalpia=self.AireHumedo.Entalpia(tdb, self.H)
            self.Xa=1/(1+self.H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        elif modo==3:
            self.Tdb=tdb
            self.Twb=twb
            self.Hs=self.AireHumedo.Humedad_Absoluta(tdb)
            Hv=self.AireHumedo.agua.Hv_DIPPR(tdb)
            def f(H):
                return self.AireHumedo.Humedad_Absoluta(twb)-H-self.AireHumedo.Calor_Especifico_Humedo(tdb, H)/Hv*(tdb-twb)
            self.H= fsolve(f, 0.001)
            self.Tdp=self.AireHumedo.Tdp(self.H)
            self.HR=self.H/self.Hs*100
            self.entalpia=self.AireHumedo.Entalpia(tdb, self.H)
            self.Xa=1/(1+self.H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        elif modo==4:
            self.Tdp=tdp
            self.HR=HR
            if HR:
                self.H=self.AireHumedo.Humedad_Absoluta(tdp)
                self.Hs=self.H/self.HR*100
            else:
                self.H=0
            self.Tdb=self.AireHumedo.Tdp(self.Hs)
            self.Twb=self.AireHumedo.Tw(self.Tdb, self.H)
            self.entalpia=self.AireHumedo.Entalpia(self.Tdb, self.H)
            self.Xa=1/(1+self.H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(self.Tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
        elif modo==5:
            self.Tdb=tdb
            self.entalpia=unidades.Enthalpy(entalpia, "kJkg")
            f=lambda h: self.AireHumedo.Entalpia(self.Tdb, h).kJkg-entalpia
            self.H=fsolve(f, 0.001)
            self.Hs=self.AireHumedo.Humedad_Absoluta(tdb)
            self.Twb=self.AireHumedo.Tw(tdb, self.H)
            self.HR=self.H/self.Hs*100
            self.Xa=1/(1+self.H*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            self.Xw=1-self.Xa
            self.volumen=self.AireHumedo.Volumen(tdb, self.Xa)
            self.densidad=unidades.Density(1/self.volumen)
            
    @property
    def corriente(self):
        corriente=Corriente(T=self.Twb, P=self.P, caudalMasico=self.caudal, fraccionMolar=[self.Xw, self.Xa])
        return corriente
        
        
if __name__ == '__main__':
#    aireHumedo=Psychrometry(1)
#    aireHumedo.Tw(298, 0.01)

    aire=Punto_Psicrometrico(caudal=1000, tdb=300, H=0.1)
    print aire.__dict__

