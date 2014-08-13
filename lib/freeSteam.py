#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to use freesteam library to calculate water properties with IAPWS-IF97
# to meos calculation or in stream calculations.
# http://freesteam.sourceforge.net/
# This library is optional, pychemqt has a python implementation for IAPWS-IF97
# in file lib/iapws.py, but if freesteam is available in your system this in c++
# is faster

#   - Freesteam: Class with properties calculations
###############################################################################

from math import exp

from PyQt4.QtGui import QApplication
from scipy.constants import R
try:
    import freesteam
except:
    pass

from lib import unidades, mEoS, iapws
from config import fluid


class Freesteam(object):
    """Class to define a water stream
    It can be defined by the pair:
        P,T
        P,h
        P,s
        P,v
        T,s
        T,x

    where:
        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -x: Quality, [-]
        -s: entropy, kJ/kgK
        -h: enthalpy, kJ/kg
        -v: specific volume, m³/kg


    """
    kwargs={"T": 0.0,
            "P": 0.0,
            "x": None,
            "h": 0.0,
            "s": 0.0,
            "v": 0.0}

    status=0
    msg="Unknown variables"

    def __init__(self, **kwargs):
        self.kwargs=Freesteam.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status=1
            self.calculo()
            self.msg=""

    @property
    def calculable(self):
        self._thermo=""
        if self.kwargs["T"] and self.kwargs["P"]:
            self._thermo=0
            self.var1=self.kwargs["P"]
            self.var2=self.kwargs["T"]
        elif self.kwargs["P"] and self.kwargs["h"]:
            self._thermo=1
            self.var1=self.kwargs["P"]
            self.var2=self.kwargs["h"]
        elif self.kwargs["P"] and self.kwargs["s"]:
            self._thermo=2
            self.var1=self.kwargs["P"]
            self.var2=self.kwargs["s"]
        elif self.kwargs["P"] and self.kwargs["v"]:
            self._thermo=3
            self.var1=self.kwargs["P"]
            self.var2=self.kwargs["v"]
        elif self.kwargs["T"] and self.kwargs["s"]:
            self._thermo=4
            self.var1=self.kwargs["T"]
            self.var2=self.kwargs["s"]
        elif self.kwargs["T"] and self.kwargs["x"]!=None:
            self._thermo=5
            self.var1=self.kwargs["T"]
            self.var2=self.kwargs["x"]
        return self._thermo+1


    def calculo(self):
        func=[freesteam.steam_pT, freesteam.steam_ph, freesteam.steam_ps,
               freesteam.steam_pv, freesteam.steam_Ts, freesteam.steam_Tx][self._thermo]
        fluido=func(self.var1, self.var2)

        self.M=unidades.Dimensionless(mEoS.H2O.M)
        self.Pc=unidades.Pressure(freesteam.PCRIT)
        self.Tc=unidades.Temperature(freesteam.TCRIT)
        self.rhoc=unidades.Density(freesteam.RHOCRIT*self.M)
        self.Tt=mEoS.H2O.Tt
        self.Tb=mEoS.H2O.Tb
        self.f_accent=unidades.Dimensionless(mEoS.H2O.f_acent)
        self.momentoDipolar=mEoS.H2O.momentoDipolar

        self.phase=self.getphase(fluido)
        self.x=unidades.Dimensionless(fluido.x)
        self.name=mEoS.H2O.name
        self.synonim=mEoS.H2O.synonym
        self.CAS=mEoS.H2O.CASNumber

        self.T=unidades.Temperature(fluido.T)
        self.P=unidades.Pressure(fluido.p)
        self.rho=unidades.Density(fluido.rho)
        self.v=unidades.SpecificVolume(1./self.rho)

        self.Liquido=fluid()
        self.Gas=fluid()
        if self.x<1:        #Fase Liquida
            liquido=freesteam.steam_Tx(fluido.T, 0.)
            self.fill(self.Liquido, liquido)
            self.Liquido.epsilon=unidades.Tension(iapws._Tension(self.T))
        if self.x>0:
            vapor=freesteam.steam_Tx(fluido.T, 1.)
            self.fill(self.Gas, vapor)

        if self.x in (0, 1):
            self.fill(self, fluido)
        else:
            self.h=unidades.Enthalpy(self.x*self.Vapor.h+(1-self.x)*self.Liquido.h)
            self.s=unidades.SpecificHeat(self.x*self.Vapor.s+(1-self.x)*self.Liquido.s)
            self.u=unidades.SpecificHeat(self.x*self.Vapor.u+(1-self.x)*self.Liquido.u)
            self.a=unidades.Enthalpy(self.x*self.Vapor.a+(1-self.x)*self.Liquido.a)
            self.g=unidades.Enthalpy(self.x*self.Vapor.g+(1-self.x)*self.Liquido.g)

            self.cv=unidades.SpecificHeat(None)
            self.cp=unidades.SpecificHeat(None)
            self.cp_cv=unidades.Dimensionless(None)
            self.w=unidades.Speed(None)

    def fill(self, fase, estado):
        fase.M=self.M
        fase.rho=unidades.Density(estado.rho)
        fase.v=unidades.SpecificVolume(estado.v)
        fase.Z=unidades.Dimensionless(self.P*estado.v/R/1000*self.M/self.T)

        fase.h=unidades.Enthalpy(estado.h)
        fase.s=unidades.SpecificHeat(estado.s)
        fase.u=unidades.Enthalpy(estado.u)
        fase.a=unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g=unidades.Enthalpy(fase.h-self.T*fase.s)

        fase.cv=unidades.SpecificHeat(estado.cv)
        fase.cp=unidades.SpecificHeat(estado.cp)
        fase.cp_cv=unidades.Dimensionless(fase.cp/fase.cv)
        fase.w=unidades.Speed(estado.w)

        fase.mu=unidades.Viscosity(estado.mu)
        fase.k=unidades.ThermalConductivity(estado.k)
        fase.nu=unidades.Diffusivity(fase.mu/fase.rho)
        fase.dielec=unidades.Dimensionless(iapws._Dielectric(estado.rho, self.T))
        fase.Prandt=unidades.Dimensionless(estado.mu*estado.cp/estado.k)

#            fase.joule=unidades.TemperaturePressure(self.Liquido["hjt"], "KkPa")
#        fase.xkappa=unidades.InvPressure(DerivTerms("IsothermalCompressibility", self.T, fase.rho, self.name), "kPa")
#        fase.alfav=unidades.InvTemperature(-estado.PFC.drhodT_constp()/estado.rho)

        cp0=iapws.prop0(self.T, self.P)
        fase.v0=unidades.SpecificVolume(cp0.v)
        fase.h0=unidades.Enthalpy(cp0.h)
        fase.u0=unidades.Enthalpy(fase.h0-self.P*1000*fase.v0)
        fase.s0=unidades.SpecificHeat(cp0.s)
        fase.a0=unidades.Enthalpy(fase.u0-self.T*fase.s0)
        fase.g0=unidades.Enthalpy(fase.h0-self.T*fase.s0)
        fase.cp0=unidades.SpecificHeat(cp0.cp)
        fase.cv0=unidades.SpecificHeat(cp0.cv)
        fase.cp0_cv=unidades.Dimensionless(fase.cp0/fase.cv0)
        fase.w0=cp0.w
#        fase.gamma0=cp0.gamma
        fase.f=unidades.Pressure(self.P*exp((fase.g-fase.g0)/R/self.T))


    def getphase(self, fld):
        u'''Return fluid phase"'''
        #check if fld above critical pressure
        if fld.p > self.Pc:
            #check if fld above critical pressure
            if fld.T > self.Tc:
                return QApplication.translate("pychemqt", "Supercritical fluid")
            else:
                return QApplication.translate("pychemqt", "Compressible liquid")
        #check if fld above critical pressure
        elif fld.T > self.Tc:
            return QApplication.translate("pychemqt", "Gas")
        #check quality
        if fld.x >= 1.:
            if self.kwargs["x"]==1.:
                return QApplication.translate("pychemqt", "Saturated vapor")
            else:
                return QApplication.translate("pychemqt", "Vapor")
        elif 0 < fld.x < 1:
            return QApplication.translate("pychemqt", "Two phases")
        elif fld.x <= 0.:
            if self.kwargs["x"]==0.:
                return QApplication.translate("pychemqt", "Saturated liquid")
            else:
                return QApplication.translate("pychemqt", "Liquid")



if __name__ == '__main__':
#    fluido=Freesteam(T=373.15, x=1)
#    print fluido.h.kJkg, fluido.P
    vapor=freesteam.steam_Tx(370, 1)
    for key, value in vapor.__dict__.iteritems():
        print key, value
    print dir(vapor)
    print ord(vapor.region)
