#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



###############################################################################
# Library to multiparameter equation of state calculation using coolprop
# http://coolprop.sourceforge.net/index.html
# optional method to meos tools calculations and to unicomponent streams
###############################################################################

from PyQt5.QtWidgets import QApplication
from scipy.constants import R

try:
    import CoolProp as CP
except ImportError as e:
    raise e

from lib import unidades
from lib.thermo import Fluid, Thermo


class CoolProp(Thermo):
    """Stream class using coolProp external library
    Parameters needed to define it are:

        -fluido: index of fluid
        -fraccion: molar fraction

        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -rho: Density, kg/m3
        -H: Enthalpy, J/kg
        -S: Entropy, J/kgK
        -x: Quality, -
    """
    kwargs = {"fluido": None,
              "fraccion": [1],

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "rho": 0.0,
              "H": 0.0,
              "S": 0.0}

    @property
    def calculable(self):
        if self.kwargs["fluido"] is not None:
            self._definition = True
        else:
            self._definition = False

        self._thermo = ""
        if self.kwargs["T"] and self.kwargs["P"]:
            self._thermo = "PT"
        elif self.kwargs["T"] and self.kwargs["x"] is not None:
            self._thermo = "TQ"
        elif self.kwargs["P"] and self.kwargs["x"] is not None:
            self._thermo = "PQ"
        elif self.kwargs["T"] and self.kwargs["rho"]:
            self._thermo = "TD"
        elif self.kwargs["P"] and self.kwargs["rho"]:
            self._thermo = "PD"
        elif self.kwargs["P"] and self.kwargs["H"]:
            self._thermo = "PH"
        elif self.kwargs["P"] and self.kwargs["S"]:
            self._thermo = "PS"
        elif self.kwargs["H"] and self.kwargs["S"]:
            self._thermo = "HS"
        return self._definition and self._thermo

    def args(self):
        if "Q" in self._thermo:
            self.kwargs["Q"] = self.kwargs["x"]
        if "D" in self._thermo:
            self.kwargs["D"] = self.kwargs["rho"]

        var1 = self.kwargs[self._thermo[0]]
        var2 = self.kwargs[self._thermo[1]]

        # units conversion to coolprop expected unit:
        # P in kPa, H in kJ/kg, S in kJ/kgK
        if self._thermo[0] == "P":
            var1 /= 1000.
        if self._thermo == "PH":
            var2 /= 1000.
        if self._thermo == "HS":
            var1 /= 1000.
        if "S" in self._thermo:
            var2 /= 1000.

        args = {self._thermo[0]: var1, self._thermo[1]: var2}
        return args

    def calculo(self):
        fluido = CP.FluidsList()[__all__[self.kwargs["fluido"]]]
        args = self.args()
        estado = CP.State(fluido, args)

        self.M = unidades.Dimensionless(estado.Props(CP.iMM))
        self.Pc = unidades.Pressure(estado.Props(CP.iPcrit), "kPa")
        self.Tc = unidades.Temperature(estado.Props(CP.iTcrit))
        self.rhoc = unidades.Density(estado.Props(CP.iRhocrit))

        self.name = fluido
        alias = CP.get_aliases(fluido)
        if len(alias) >= 2:
            self.synonim = alias[1]
        else:
            self.synonim = ""
        self.CAS = CP.get_CAS_code(fluido)

        self.Tt = unidades.Temperature(estado.Props(CP.iPtriple))
        self.Tb = unidades.Temperature(None)
        self.f_accent = unidades.Dimensionless(estado.Props(CP.iAccentric))
        # self.momentoDipolar = unidades.DipoleMoment(estado.Props(CP.iDipole))
        self.momentoDipolar = unidades.DipoleMoment(None)
        self.Rgas = unidades.SpecificHeat(R/self.M)

        self.T = unidades.Temperature(estado.T)
        self.P = unidades.Pressure(estado.p, "kPa")
        self.rho = unidades.Density(estado.rho)
        self.v = unidades.SpecificVolume(1./self.rho)

        self.phase, x = self.getphase(estado)
        self.x = unidades.Dimensionless(x)

        self.Liquido = Fluid()
        self.Gas = Fluid()
        if self.x == 0:
            # liquid phase
            self.fill(self.Liquido, estado)
            self.Liquido.sigma = unidades.Tension(estado.Props(CP.iI))
            self.fill(self, estado)
        elif self.x == 1:
            # vapor phase
            self.fill(self.Gas, estado)
            self.fill(self, estado)
        else:
            # Two phase
            estado.update_Trho(estado.T, estado.PFC.rhoL())
            self.fill(self.Liquido, estado)
            self.Liquido.sigma = unidades.Tension(estado.Props(CP.iI))
            estado.update_Trho(estado.T, estado.PFC.rhoV())
            self.fill(self.Gas, estado)

            self.h = unidades.Enthalpy(self.x*self.Gas.h+(1-self.x)*self.Liquido.h)
            self.s = unidades.SpecificHeat(self.x*self.Gas.s+(1-self.x)*self.Liquido.s)
            self.u = unidades.SpecificHeat(self.x*self.Gas.u+(1-self.x)*self.Liquido.u)
            self.a = unidades.Enthalpy(self.x*self.Gas.a+(1-self.x)*self.Liquido.a)
            self.g = unidades.Enthalpy(self.x*self.Gas.g+(1-self.x)*self.Liquido.g)

            self.cv = unidades.SpecificHeat(None)
            self.cp = unidades.SpecificHeat(None)
            self.cp_cv = unidades.Dimensionless(None)
            self.w = unidades.Speed(None)

    def fill(self, fase, estado):
        fase.M = self.M
        fase.rho = unidades.Density(estado.rho)
        fase.v = unidades.SpecificVolume(1./fase.rho)
        fase.Z = unidades.Dimensionless(CP.DerivTerms("Z", self.T, fase.rho, self.name))

        fase.h = unidades.Enthalpy(estado.h, "kJkg")
        fase.s = unidades.SpecificHeat(estado.s, "kJkgK")
        fase.u = unidades.Enthalpy(estado.u, "kJkg")
        fase.a = unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g = unidades.Enthalpy(estado.Props(CP.iG), "kJkg")

        fase.cv = unidades.SpecificHeat(estado.cv, "kJkgK")
        fase.cp = unidades.SpecificHeat(estado.cp, "kJkgK")
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.w = unidades.Speed(estado.Props(CP.iA))

        fase.mu = unidades.Viscosity(estado.visc)
        fase.k = unidades.ThermalConductivity(estado.k, "kWmK")
        fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        fase.dielec = unidades.Dimensionless(None)
        fase.Prandt = unidades.Dimensionless(estado.Prandtl)

        # fase.joule = unidades.TemperaturePressure(self.Liquido["hjt"], "KkPa")
        fase.kappa = unidades.InvPressure(CP.DerivTerms("IsothermalCompressibility", self.T, fase.rho, self.name), "kPa")
        fase.alfav = unidades.InvTemperature(-estado.PFC.drhodT_constp()/estado.rho)

        fase.cp0 = unidades.SpecificHeat(estado.Props(CP.iC0), "kJkgK")
        fase.cp0_cv = unidades.Dimensionless(fase.cp0/fase.cv)

    def getphase(self, estado):
        """Return fluid phase with translation support"""
        phase = estado.Phase()
        if phase == CP.iphase_supercritical:
            if self.T > self.Tc:
                return QApplication.translate("pychemqt", "Supercritical fluid"), 1.
            else:
                return QApplication.translate("pychemqt", "Compressible liquid"), 1.
        elif phase == CP.iphase_gas:
            if estado.superheat > 0:
                return QApplication.translate("pychemqt", "Vapor"), 1.
            else:
                return QApplication.translate("pychemqt", "Saturated vapor"), 1.
        elif phase == CP.iphase_liquid:
            if estado.subcooling > 0:
                return QApplication.translate("pychemqt", "Liquid"), 0.
            else:
                return QApplication.translate("pychemqt", "Saturated liquid"), 0.
        elif phase == CP.iphase_twophase:
                return QApplication.translate("pychemqt", "Two phases"), estado.Q


noIds = {"ParaHydrogen": 5,
         "OrthoHydrogen": 6,
         "R1234yf": 12,
         "R1234ze(E)": 13,
         "SES36": 16,
         "R236FA": 27,
         "R365MFC": 30,
         "HFE143m": 32,
         "R1234ze(Z)": 49,
         "MDM": 83,
         "MD2M": 84,
         "MD3M": 85,
         "MethylPalmitate": 95,
         "MethylStearate": 96,
         "MethylLinoleate": 98,
         "MethylLinolenate": 99,
         "Deuterium": 104,
         "ParaDeuterium": 105,
         "OrthoDeuterium": 106,
         "R404A": 108,
         "R410A": 109,
         "R407C": 110,
         "R507A": 111}


__all__ = {1: 4,
           258: 36,
           643: 23,
           4: 10,
           645: 14,
           134: 19,
           7: 57,
           8: 75,
           9: 58,
           10: 76,
           11: 77,
           140: 65,
           642: 43,
           14: 55,
           50: 56,
           971: 60,
           994: 63,
           15: 34,
           16: 24,
           146: 112,
           1633: 68,
           1631: 81,
           3: 72,
           51: 67,
           22: 17,
           23: 25,
           24: 91,
           25: 93,
           26: 94,
           5: 74,
           133: 20,
           36: 26,
           6: 73,
           38: 79,
           40: 33,
           41: 62,
           42: 100,
           43: 101,
           44: 102,
           45: 103,
           46: 9,
           48: 53,
           692: 42,
           52: 59,
           671: 70,
           245: 80,
           66: 39,
           68: 92,
           12: 78,
           2: 71,
           13: 61,
           208: 40,
           212: 2,
           215: 45,
           216: 47,
           217: 82,
           218: 46,
           219: 54,
           220: 15,
           225: 52,
           98: 7,
           231: 44,
           232: 48,
           107: 37,
           236: 64,
           110: 66,
           241: 69,
           243: 22,
           117: 41,
           247: 31,
           62: 0,
           47: 3,
           63: 11,
           475: 107,
           919: 97,
           1671: 90,
           1430: 89,
           1835: 88,
           1376: 87,
           1674: 86,
           953: 18,
           1798: 21,
           1873: 28,
           1872: 29,
           1231: 35,
           1629: 38,
           1817: 51}


if __name__ == '__main__':
#    from CoolProp.CoolProp import IProps, get_Fluid_index
#    iPropane = get_Fluid_index('R134a')
#    fluido2=CoolProp(fluido=get_Fluid_index('R134a'), T=273.15, x=1.)
#    print fluido2.h.kJkg
##    In [9]: Props('H','T',273.15,'Q',1,'R134a')
##Out[9]: 398.60345362765497
#    fluido=IAPWS97_PT(101325, 300)
#    print fluido.cp
    fluido = CoolProp(fluido=0, T=300, P=101325)
    print(fluido.M)

