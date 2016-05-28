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
# Module with common thermal utilities:
#   - Fluid: Common class for thermodynamic method
#      Implement properites save/load to stream. Null properties state
###############################################################################


from PyQt5.QtWidgets import QApplication

from lib import unidades


data = [(QApplication.translate("pychemqt", "Temperature"), "T", unidades.Temperature),
        (QApplication.translate("pychemqt", "Reduced temperature"), "Tr", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Pressure"), "P", unidades.Pressure),
        (QApplication.translate("pychemqt", "Reduced Pressure"), "Pr", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Quality"), "x", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Density"), "rho", unidades.Density),
        (QApplication.translate("pychemqt", "Volume"), "v", unidades.SpecificVolume),
        (QApplication.translate("pychemqt", "Enthalpy"), "h", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Entropy"), "s", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Internal Energy"), "u", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Gibbs Free Energy"), "g", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Helmholtz Free Energy"), "a", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Specific isochoric heat capacity"), "cv", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Specific isobaric heat capacity"), "cp", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Heat capacities ratio"), "cp_cv", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Speed sound"), "w", unidades.Speed),
        (QApplication.translate("pychemqt", "Compresibility"), "Z", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Fugacity"), "f", unidades.Pressure),
        (QApplication.translate("pychemqt", "Isoentropic exponent"), "gamma", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Volumetric Expansitivy"), "alfav", unidades.InvTemperature),
        (QApplication.translate("pychemqt", "Isotermic compresibility"), "kappa", unidades.InvPressure),
        (QApplication.translate("pychemqt", "Relative pressure"), "alfap", unidades.InvTemperature),
        (QApplication.translate("pychemqt", "Isothermal stress"), "betap", unidades.Density),
        (QApplication.translate("pychemqt", "Joule-Thomson coefficient"), "joule", unidades.TemperaturePressure),
        (QApplication.translate("pychemqt", "Isothermal throttling coefficient"), "deltat", unidades.EnthalpyPressure),
        (QApplication.translate("pychemqt", "Vaporization heat"), "Hvap", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Vaporization entropy"), "Svap", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Viscosity"), "mu", unidades.Viscosity),
        (QApplication.translate("pychemqt", "Thermal conductivity"), "k", unidades.ThermalConductivity),
        (QApplication.translate("pychemqt", "Kinematic viscosity"), "nu", unidades.Diffusivity),
        (QApplication.translate("pychemqt", "Thermal diffusivity"), "alfa", unidades.Diffusivity),
        (QApplication.translate("pychemqt", "Surface tension"), "sigma", unidades.Tension),
        (QApplication.translate("pychemqt", "Prandtl number"), "Prandt", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Refractive index"), "n", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Ideal gas Specific volume"), "v0", unidades.SpecificVolume),
        (QApplication.translate("pychemqt", "Ideal gas Density"), "rho0", unidades.Density),
        (QApplication.translate("pychemqt", "Ideal gas Specific enthalpy"), "h0", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Ideal gas Specific internal energy"), "u0", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Ideal gas Specific entropy"), "s0", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Ideal gas Specific Helmholtz free energy"), "a0", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Ideal gas Specific Gibbs free energy"), "g0", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Ideal gas Specific isobaric heat capacity"), "cp0", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Ideal gas Specific isochoric heat capacity"), "cv0", unidades.SpecificHeat),
        (QApplication.translate("pychemqt", "Ideal gas heat capacities ratio"), "cp0_cv", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Ideal gas Isoentropic exponent"), "gamma0", unidades.Dimensionless)]


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
            fluid["kappa"] = self.kappa

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
            self.kappa = unidades.InvPressure(fluid["kappa"])

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

            self.alfap = unidades.Dimensionless(fluid["alfap"])
            self.betap = unidades.Dimensionless(fluid["betap"])
            self.f = unidades.Pressure(fluid["f"])

            self.Q = unidades.VolFlow(fluid["volFlow"])
            self.caudalmasico = unidades.MassFlow(fluid["massFlow"])
            self.caudalmolar = unidades.MolarFlow(fluid["molarFlow"])
            self.fraccion = [unidades.Dimensionless(x) for x in fluid["fraction"]]
            self.fraccion_masica = [unidades.Dimensionless(x) for x in fluid["massFraction"]]
            self.caudalunitariomasico = [unidades.MassFlow(x) for x in fluid["massUnitFlow"]]
            self.caudalunitariomolar = [unidades.MolarFlow(x) for x in fluid["molarUnitFlow"]]


class Fluid_MEOS(Fluid):
    """Extended custom object for meos"""
    v = None
    rho = None

    h = None
    s = None
    u = None
    a = None
    g = None

    cp = None
    cv = None
    cp_cv = None
    w = None
    Z = None
    fi = None
    f = None

    rhoM = None
    hM = None
    sM = None
    uM = None
    aM = None
    gM = None
    cvM = None
    cpM = None

    mu = None
    k = None
    nu = None
    Prandt = None
    epsilon = None
    alfa = None
    n = None

    alfap = None
    betap = None
    joule = None
    Gruneisen = None
    alfav = None
    kappa = None
    betas = None
    gamma = None
    Kt = None
    kt = None
    Ks = None
    ks = None
    dpdT_rho = None
    dpdrho_T = None
    drhodT_P = None
    drhodP_T = None
    dhdT_rho = None
    dhdT_P = None
    dhdrho_T = None
    dhdrho_P = None
    dhdP_T = None
    dhdP_rho = None

    Z_rho = None
    IntP = None
    hInput = None

    def writeStatetoJSON(self, state, fase):
        Fluid.writeStatetoJSON(self, state, fase)
        fluid = state[fase]
        if self._bool:
            fluid["rhoM"] = self.rhoM
            fluid["hM"] = self.hM
            fluid["sM"] = self.sM
            fluid["uM"] = self.uM
            fluid["aM"] = self.aM
            fluid["gM"] = self.gM
            fluid["cvM"] = self.cvM
            fluid["cpM"] = self.cpM

            fluid["fi"] = self.fi
            fluid["betas"] = self.betas
            fluid["Gruneisen"] = self.Gruneisen
            fluid["dpdT_rho"] = self.dpdT_rho
            fluid["dpdrho_T"] = self.dpdrho_T
            fluid["drhodT_P"] = self.drhodT_P
            fluid["drhodP_T"] = self.drhodP_T
            fluid["dhdT_rho"] = self.dhdT_rho
            fluid["dhdP_T"] = self.dhdP_T
            fluid["dhdT_P"] = self.dhdT_P
            fluid["dhdrho_T"] = self.dhdrho_T
            fluid["dhdP_rho"] = self.dhdP_rho
            fluid["kt"] = self.kt
            fluid["ks"] = self.ks
            fluid["Ks"] = self.Ks
            fluid["Kt"] = self.Kt
            fluid["IntP"] = self.IntP
            fluid["hInput"] = self.hInput
            fluid["alfa"] = self.alfa

        state[fase] = fluid

    def readStatefromJSON(self, fluid):
        Fluid.readStatefromJSON(self, fluid)
        if self._bool:
            self.rhoM = unidades.MolarDensity(fluid["rhoM"])
            self.hM = unidades.MolarEnthalpy(fluid["hM"])
            self.sM = unidades.MolarSpecificHeat(fluid["sM"])
            self.uM = unidades.MolarEnthalpy(fluid["uM"])
            self.aM = unidades.MolarEnthalpy(fluid["aM"])
            self.gM = unidades.MolarEnthalpy(fluid["gM"])
            self.cvM = unidades.MolarSpecificHeat(fluid["cvM"])
            self.cpM = unidades.MolarSpecificHeat(fluid["cpM"])

            self.fi = unidades.Dimensionless(fluid["fi"])
            self.betas = unidades.TemperaturePressure(fluid["betas"])
            self.Gruneisen = unidades.Dimensionless(fluid["Gruneisen"])
            self.dpdT_rho = unidades.PressureTemperature(fluid["dpdT_rho"])
            self.dpdrho_T = unidades.PressureDensity(fluid["dpdrho_T"])
            self.drhodT_P = unidades.DensityTemperature(fluid["drhodT_P"])
            self.drhodP_T = unidades.DensityPressure(fluid["drhodP_T"])
            self.dhdT_rho = unidades.SpecificHeat(fluid["dhdT_rho"])
            self.dhdP_T = unidades.EnthalpyPressure(fluid["dhdP_T"])
            self.dhdT_P = unidades.SpecificHeat(fluid["dhdT_P"])
            self.dhdrho_T = unidades.EnthalpyDensity(fluid["dhdrho_T"])
            self.dhdP_rho = unidades.EnthalpyPressure(fluid["dhdP_rho"])
            self.kt = unidades.Dimensionless(fluid["kt"])
            self.ks = unidades.InvPressure(fluid["ks"])
            self.Ks = unidades.Pressure(fluid["Ks"])
            self.Kt = unidades.Pressure(fluid["Kt"])
            self.IntP = unidades.Pressure(fluid["IntP"])
            self.hInput = unidades.Enthalpy(fluid["hInput"])
            self.alfa = unidades.Diffusivity(fluid["alfa"])


class Thermo(object):
    """Class with common functionality for special thermo model, children class
    are iapws, coolprop, refprop"""
    # TODO: Add here the common functionality of model
    pass
