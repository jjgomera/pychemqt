#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library to use freesteam library to calculate water properties with
# IAPWS-IF97 to meos calculation or in stream calculations.
# http://freesteam.sourceforge.net/
# This library is optional, pychemqt has a python implementation for IAPWS-IF97
# in file lib/iapws.py, but if freesteam is available in your system this is
# faster because is implemented in c++

#   - Freesteam: Class with properties calculations
###############################################################################


from math import exp

from PyQt5.QtWidgets import QApplication
from scipy.constants import R

try:
    import freesteam
except:
    pass

from lib import unidades, mEoS, iapws
from lib.thermo import ThermoWater


class Freesteam(ThermoWater):
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
    kwargs = {"T": 0.0,
              "P": 0.0,
              "x": None,
              "h": 0.0,
              "s": 0.0,
              "v": 0.0,
              "l": 0.5893}
    __doi__ = [
        {"autor": "Freesteam",
         "title": "Revised Release on the IAPWS Industrial Formulation 1997 "
                  "for the Thermodynamic Properties of Water and Steam",
         "ref": "",
         "doi": ""}]

    M = unidades.Dimensionless(mEoS.H2O.M)
    Pc = unidades.Pressure(freesteam.PCRIT)
    Tc = unidades.Temperature(freesteam.TCRIT)
    rhoc = unidades.Density(freesteam.RHOCRIT*M)
    Tt = mEoS.H2O.Tt
    Tb = mEoS.H2O.Tb
    f_accent = unidades.Dimensionless(mEoS.H2O.f_acent)
    momentoDipolar = mEoS.H2O.momentoDipolar

    @property
    def calculable(self):
        self._thermo = -1
        if self.kwargs["T"] and self.kwargs["P"]:
            self._thermo = 0
            self.var1 = self.kwargs["P"]
            self.var2 = self.kwargs["T"]
        elif self.kwargs["P"] and self.kwargs["h"]:
            self._thermo = 1
            self.var1 = self.kwargs["P"]
            self.var2 = self.kwargs["h"]
        elif self.kwargs["P"] and self.kwargs["s"]:
            self._thermo = 2
            self.var1 = self.kwargs["P"]
            self.var2 = self.kwargs["s"]
        elif self.kwargs["P"] and self.kwargs["v"]:
            self._thermo = 3
            self.var1 = self.kwargs["P"]
            self.var2 = self.kwargs["v"]
        elif self.kwargs["T"] and self.kwargs["s"]:
            self._thermo = 4
            self.var1 = self.kwargs["T"]
            self.var2 = self.kwargs["s"]
        elif self.kwargs["T"] and self.kwargs["x"] is not None:
            self._thermo = 5
            self.var1 = self.kwargs["T"]
            self.var2 = self.kwargs["x"]
        return self._thermo+1

    def calculo(self):
        method = [freesteam.steam_pT, freesteam.steam_ph, freesteam.steam_ps,
                  freesteam.steam_pv, freesteam.steam_Ts, freesteam.steam_Tx]
        func = method[self._thermo]
        fluido = func(self.var1, self.var2)

        self.phase = self.getphase(fluido)
        self.x = unidades.Dimensionless(fluido.x)
        self.name = mEoS.H2O.name
        self.synonim = mEoS.H2O.synonym
        self.CAS = mEoS.H2O.CASNumber

        self.T = unidades.Temperature(fluido.T)
        self.P = unidades.Pressure(fluido.p)
        self.Tr = unidades.Dimensionless(self.T/self.Tc)
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.rho = unidades.Density(fluido.rho)
        self.v = unidades.SpecificVolume(1./self.rho)

        cp0 = iapws.prop0(self.T, self.P)
        self._cp0(cp0)

        self.Liquido = ThermoWater()
        self.Gas = ThermoWater()
        if self.x == 0:
            # only liquid phase
            self.fill(self, fluido)
            self.fill(self.Liquido, fluido)
            self.Liquido.sigma = unidades.Tension(freesteam.surftens_T(self.T))

            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)
        elif self.x == 1:
            # only vapor phase
            self.fill(self, fluido)
            self.fill(self.Gas, fluido)

            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)
        else:
            # two phases
            liquido = freesteam.steam_Tx(fluido.T, 0.)
            self.fill(self.Liquido, liquido)
            self.Liquido.sigma = unidades.Tension(freesteam.surftens_T(self.T))
            vapor = freesteam.steam_Tx(fluido.T, 1.)
            self.fill(self.Gas, vapor)

            self.h = unidades.Enthalpy(
                self.x*self.Gas.h+(1-self.x)*self.Liquido.h)
            self.s = unidades.SpecificHeat(
                self.x*self.Gas.s+(1-self.x)*self.Liquido.s)
            self.u = unidades.SpecificHeat(
                self.x*self.Gas.u+(1-self.x)*self.Liquido.u)
            self.a = unidades.Enthalpy(
                self.x*self.Gas.a+(1-self.x)*self.Liquido.a)
            self.g = unidades.Enthalpy(
                self.x*self.Gas.g+(1-self.x)*self.Liquido.g)

            self.cv = unidades.SpecificHeat(None)
            self.cp = unidades.SpecificHeat(None)
            self.cp_cv = unidades.Dimensionless(None)
            self.w = unidades.Speed(None)

            self.Hvap = unidades.Enthalpy(vapor["h"]-liquido["h"], "kJkg")
            self.Svap = unidades.SpecificHeat(vapor["s"]-liquido["s"], "kJkgK")

    def fill(self, fase, estado):
        fase._bool = True
        fase.M = self.M
        fase.rho = unidades.Density(estado.rho)
        fase.v = unidades.SpecificVolume(estado.v)
        fase.Z = unidades.Dimensionless(self.P*estado.v/R/1000*self.M/self.T)

        fase.h = unidades.Enthalpy(estado.h)
        fase.s = unidades.SpecificHeat(estado.s)
        fase.u = unidades.Enthalpy(estado.u)
        fase.a = unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g = unidades.Enthalpy(fase.h-self.T*fase.s)
        fase.f = unidades.Pressure(self.P*exp((fase.g-self.g0)/R/self.T))

        fase.cv = unidades.SpecificHeat(estado.cv)
        fase.cp = unidades.SpecificHeat(estado.cp)
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.w = unidades.Speed(estado.w)

        fase.mu = unidades.Viscosity(estado.mu)
        fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        fase.k = unidades.ThermalConductivity(estado.k)
        fase.alfa = unidades.Diffusivity(fase.k/1000/fase.rho/fase.cp)
        fase.epsilon = unidades.Dimensionless(
            iapws._Dielectric(estado.rho, self.T))
        fase.Prandt = unidades.Dimensionless(estado.mu*estado.cp/estado.k)
        fase.n = unidades.Dimensionless(
            iapws._Refractive(fase.rho, self.T, self.kwargs["l"]))

        fase.alfav = unidades.InvTemperature(1/fase.v*estado.deriv("vTp"))
        fase.kappa = unidades.InvPressure(-1/fase.v*estado.deriv("vpT"))

        fase.joule = unidades.TemperaturePressure(estado.deriv("Tph"))
        fase.deltat = unidades.EnthalpyPressure(estado.deriv("hpT"))
        fase.gamma = unidades.Dimensionless(
            -fase.v/self.P/1000*estado.deriv("pvs"))

        fase.alfap = unidades.Density(fase.alfav/self.P/fase.kappa)
        fase.betap = unidades.Density(-1/self.P/1000*estado.deriv("pvT"))

    def getphase(self, fld):
        """Return fluid phase with translation support"""
        # check if fld above critical pressure
        if fld.p > self.Pc:
            # check if fld above critical pressure
            if fld.T > self.Tc:
                return QApplication.translate("pychemqt", "Supercritical fluid")
            else:
                return QApplication.translate("pychemqt", "Compressible liquid")
        # check if fld above critical pressure
        elif fld.T > self.Tc:
            return QApplication.translate("pychemqt", "Gas")
        # check quality
        if fld.x >= 1.:
            if self.kwargs["x"] == 1.:
                return QApplication.translate("pychemqt", "Saturated vapor")
            else:
                return QApplication.translate("pychemqt", "Vapor")
        elif 0 < fld.x < 1:
            return QApplication.translate("pychemqt", "Two phases")
        elif fld.x <= 0.:
            if self.kwargs["x"] == 0.:
                return QApplication.translate("pychemqt", "Saturated liquid")
            else:
                return QApplication.translate("pychemqt", "Liquid")


if __name__ == '__main__':
    fluido = Freesteam(T=373.15, x=1)
    print(fluido.h.kJkg, fluido.P)
