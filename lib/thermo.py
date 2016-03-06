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



#################################################################################
# Module with common thermal utilities:
#   - Fluid: Common class for thermodynamic method
#      Implement properites save/load to stream. Null properties state
#################################################################################


from . import unidades


class Fluid(object):
    """Custom object with null parameter to model a fluid with properties"""
    _bool = False
    h = 0
    s = 0
    cp = 0
    cv = 0
    cp_cv = 0
    cp0 = 0

    def writeStatetoJSON(self, state, fase):
        fluid = {}
        if self._bool:
            fluid["M"] = self.M
            fluid["v"] = self.v

            fluid["h"] = self.h
            fluid["s"] = self.s
            fluid["u"] = self.u
            fluid["a"] = self.a
            fluid["g"] = self.g

            fluid["cv"] = self.cv
            fluid["cp"] = self.cp
            fluid["cp/cv"] = self.cp_cv
            fluid["w"] = self.w

            fluid["Z"] = self.Z
            fluid["alfav"] = self.alfav
            fluid["xkappa"] = self.xkappa

            fluid["mu"] = self.mu
            fluid["k"] = self.k
            fluid["nu"] = self.nu
            fluid["epsilon"] = self.epsilon
            fluid["Prandt"] = self.Prandt
            fluid["n"] = self.n

            fluid["alfa"] = self.alfa
            fluid["joule"] = self.joule
            fluid["deltat"] = self.deltat
            fluid["gamma"] = self.gamma

            fluid["alfap"] = self.alfap
            fluid["betap"] = self.betap

            fluid["vo"] = self.v0
            fluid["ho"] = self.h0
            fluid["uo"] = self.u0
            fluid["so"] = self.s0
            fluid["ao"] = self.a0
            fluid["go"] = self.g0

            fluid["cpo"] = self.cp0
            fluid["cvo"] = self.cv0
            fluid["cpo/cvo"] = self.cp0_cv
            fluid["wo"] = self.w0
            fluid["gammao"] = self.gamma0
            fluid["f"] = self.f

            fluid["volFlow"] = self.Q
            fluid["massFlow"] = self.caudalmasico
            fluid["molarFlow"] = self.caudalmolar
            fluid["fraction"] = self.fraccion
            fluid["massFraction"] = self.fraccion_masica
            fluid["massUnitFlow"] = self.caudalunitariomasico
            fluid["molarUnitFlow"] = self.caudalunitariomolar
        state[fase] = fluid

    def readStatefromJSON(self, fluid):
        if fluid:
            self._bool = True

            self.M = unidades.Dimensionless(fluid["M"])
            self.v = unidades.SpecificVolume(fluid["v"])
            self.rho = unidades.Density(1/self.v)

            self.h = unidades.Enthalpy(fluid["h"])
            self.s = unidades.SpecificHeat(fluid["s"])
            self.u = unidades.Enthalpy(fluid["u"])
            self.a = unidades.Enthalpy(fluid["a"])
            self.g = unidades.Enthalpy(fluid["g"])

            self.cv = unidades.SpecificHeat(fluid["cv"])
            self.cp = unidades.SpecificHeat(fluid["cp"])
            self.cp_cv = unidades.Dimensionless(fluid["cp/cv"])
            self.w = unidades.Speed(fluid["w"])

            self.Z = unidades.Dimensionless(fluid["Z"])
            self.alfav = unidades.InvTemperature(fluid["alfav"])
            self.xkappa = unidades.InvPressure(fluid["xkappa"])

            self.mu = unidades.Viscosity(fluid["mu"])
            self.k = unidades.ThermalConductivity(fluid["k"])
            self.nu = unidades.Diffusivity(fluid["nu"])
            self.epsilon = unidades.Dimensionless(fluid["epsilon"])
            self.Prandt = unidades.Dimensionless(fluid["Prandt"])
            self.n = unidades.Dimensionless(fluid["n"])

            self.alfa = unidades.Diffusivity(fluid["alfa"])
            self.joule = unidades.TemperaturePressure(fluid["joule"])
            self.deltat = unidades.EnthalpyPressure(fluid["deltat"])
            self.gamma = unidades.Dimensionless(fluid["gamma"])

            self.alfap=fluid["alfap"]
            self.betap=fluid["betap"]

            self.v0 = unidades.SpecificVolume(fluid["vo"])
            self.h0 = unidades.Enthalpy(fluid["ho"])
            self.u0 = unidades.Enthalpy(fluid["uo"])
            self.s0 = unidades.SpecificHeat(fluid["so"])
            self.a0 = unidades.Enthalpy(fluid["ao"])
            self.g0 = unidades.Enthalpy(fluid["go"])

            self.cp0 = unidades.SpecificHeat(fluid["cpo"])
            self.cv0 = unidades.SpecificHeat(fluid["cvo"])
            self.cp0_cv = unidades.Dimensionless(fluid["cpo/cvo"])
            self.w0 = unidades.Speed(fluid["wo"])
            self.gamma0 = unidades.Dimensionless(fluid["gammao"])
            self.f = unidades.Pressure(fluid["f"])

            self.Q = unidades.VolFlow(fluid["volFlow"])
            self.caudalmasico = unidades.MassFlow(fluid["massFlow"])
            self.caudalmolar = unidades.MolarFlow(fluid["molarFlow"])
            self.fraccion = fluid["fraction"]
            self.fraccion_masica = fluid["massFraction"]
            self.caudalunitariomasico = fluid["massUnitFlow"]
            self.caudalunitariomolar = fluid["molarUnitFlow"]
