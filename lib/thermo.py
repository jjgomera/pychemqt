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
# Module with common thermal utilities:
#   - Fluid: Common class for thermodynamic method
#      Implement properites save/load to stream. Null properties state
###############################################################################


from PyQt5.QtWidgets import QApplication

from lib import unidades


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

    _bool = False
    status = 0
    msg = "Unknown variables"
    kwargs = {}

    h = 0
    s = 0
    cp = 0
    cv = 0
    cp_cv = 0
    cp0 = 0

    def __init__(self, **kwargs):
        self.kwargs = self.__class__.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()
            self.msg = "Solved"

    def calculable(self):
        pass

    def calculo(self):
        pass

    def _cp0(self, cp0):
        "Set ideal properties to state"""
        self.v0 = unidades.SpecificVolume(cp0["v"])
        self.rho0 = unidades.Density(1./cp0["v"])
        self.h0 = unidades.Enthalpy(cp0["h"])
        self.u0 = unidades.Enthalpy(self.h0-self.P*self.v0)
        self.s0 = unidades.SpecificHeat(cp0["s"])
        self.a0 = unidades.Enthalpy(self.u0-self.T*self.s0)
        self.g0 = unidades.Enthalpy(self.h0-self.T*self.s0)

        self.cp0 = unidades.SpecificHeat(cp0["cp"])
        self.cv0 = unidades.SpecificHeat(cp0["cv"])
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv0)
        self.w0 = unidades.Speed(cp0["w"])
        self.gamma0 = self.cp0_cv

    def derivative(self, z, x, y, fase):
        """Calculate generic partial derivative: (δz/δx)y
        where x, y, z can be: P, T, v, u, h, s, g, a"""
        dT = {"P": 0,
              "T": 1,
              "v": fase.v*fase.alfav,
              "u": fase.cp-self.P*fase.v*fase.alfav,
              "h": fase.cp,
              "s": fase.cp/self.T,
              "g": -fase.s,
              "a": -self.P*fase.v*fase.alfav-fase.s}
        dP = {"P": 1,
              "T": 0,
              "v": -fase.v*fase.kappa,
              "u": fase.v*(self.P*fase.kappa-self.T*fase.alfav),
              "h": fase.v*(1-self.T*fase.alfav),
              "s": -fase.v*fase.alfav,
              "g": fase.v,
              "a": self.P*fase.v*fase.kappa}
        return (dP[z]*dT[y]-dT[z]*dP[y])/(dP[x]*dT[y]-dT[x]*dP[y])

    @classmethod
    def properties(cls):
        l = [
            (QApplication.translate("pychemqt", "Temperature"), "T", unidades.Temperature),
            (QApplication.translate("pychemqt", "Reduced temperature"), "Tr", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Pressure"), "P", unidades.Pressure),
            (QApplication.translate("pychemqt", "Reduced Pressure"), "Pr", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Quality"), "x", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Density"), "rho", unidades.Density),
            (QApplication.translate("pychemqt", "Molar Density"), "rhoM", unidades.MolarDensity),
            (QApplication.translate("pychemqt", "Volume"), "v", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "Enthalpy"), "h", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Molar Enthalpy"), "hM", unidades.MolarEnthalpy),
            (QApplication.translate("pychemqt", "Entropy"), "s", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Molar Entropy"), "sM", unidades.MolarSpecificHeat),
            (QApplication.translate("pychemqt", "Internal Energy"), "u", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Molar Internal Energy"), "uM", unidades.MolarEnthalpy),
            (QApplication.translate("pychemqt", "Helmholtz Free Energy"), "a", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Molar Helmholtz Free Energy"), "aM", unidades.MolarEnthalpy),
            (QApplication.translate("pychemqt", "Gibbs Free Energy"), "g", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Molar Gibbs Free Energy"), "gM", unidades.MolarEnthalpy),
            (QApplication.translate("pychemqt", "Specific isochoric heat capacity"), "cv", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Molar Specific isochoric heat capacity"), "cvM", unidades.MolarSpecificHeat),
            (QApplication.translate("pychemqt", "Specific isobaric heat capacity"), "cp", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Molar Specific isobaric heat capacity"), "cpM", unidades.MolarSpecificHeat),
            (QApplication.translate("pychemqt", "Heat  capacities ratio"), "cp_cv", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Speed sound"), "w", unidades.Speed),
            (QApplication.translate("pychemqt", "Compresibility"), "Z", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Fugacity coefficient"), "fi", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Fugacity"), "f", unidades.Pressure),
            (QApplication.translate("pychemqt", "Isoentropic exponent"), "gamma", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Volume Expansivity"), "alfav", unidades.InvTemperature),  # 1/V dV/dt = -1/D dD/dt
            (QApplication.translate("pychemqt", "Isothermal compresibility"), "kappa", unidades.InvPressure),  # -1/V (dV/dP)T = 1/D (dD/dP)T
            (QApplication.translate("pychemqt", "Adiabatic compresibility"), "kappas", unidades.InvPressure),  # -1/V (dV/dP)s = 1/D (dD/dP)s
            (QApplication.translate("pychemqt", "Relative pressure coerricient"), "alfap", unidades.InvTemperature),  # 1/P (dP/dT)v
            (QApplication.translate("pychemqt", "Isothermal stress coefficient"), "betap", unidades.Density),  # -1/P (dP/dv)T = 1/P (dP/dD)T
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
        return l

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
            fluid["Prandt"] = self.Prandt

            fluid["alfa"] = self.alfa
            fluid["joule"] = self.joule
            fluid["deltat"] = self.deltat
            fluid["gamma"] = self.gamma

            fluid["alfap"] = self.alfap
            fluid["betap"] = self.betap
            fluid["fi"] = self.fi
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
            self.Prandt = unidades.Dimensionless(fluid["Prandt"])

            self.alfa = unidades.Diffusivity(fluid["alfa"])
            self.joule = unidades.TemperaturePressure(fluid["joule"])
            self.deltat = unidades.EnthalpyPressure(fluid["deltat"])
            self.gamma = unidades.Dimensionless(fluid["gamma"])

            self.alfap = unidades.Dimensionless(fluid["alfap"])
            self.betap = unidades.Dimensionless(fluid["betap"])
            self.fi = [unidades.Dimensionless(f) for f in fluid["fi"]]
            self.f = [unidades.Pressure(f) for f in fluid["f"]]

            self.Q = unidades.VolFlow(fluid["volFlow"])
            self.caudalmasico = unidades.MassFlow(fluid["massFlow"])
            self.caudalmolar = unidades.MolarFlow(fluid["molarFlow"])
            self.fraccion = [unidades.Dimensionless(x) for x in fluid["fraction"]]
            self.fraccion_masica = [unidades.Dimensionless(x) for x in fluid["massFraction"]]
            self.caudalunitariomasico = [unidades.MassFlow(x) for x in fluid["massUnitFlow"]]
            self.caudalunitariomolar = [unidades.MolarFlow(x) for x in fluid["molarUnitFlow"]]

            self.rhoM = unidades.MolarDensity(self.rho/self.M)
            self.hM = unidades.MolarEnthalpy(self.h/self.M)
            self.sM = unidades.MolarSpecificHeat(self.s/self.M)
            self.uM = unidades.MolarEnthalpy(self.u/self.M)
            self.aM = unidades.MolarEnthalpy(self.a/self.M)
            self.gM = unidades.MolarEnthalpy(self.g/self.M)
            self.cvM = unidades.MolarSpecificHeat(self.cv/self.M)
            self.cpM = unidades.MolarSpecificHeat(self.cp/self.M)


class ThermoWater(Thermo):
    """Custom specified thermo instance to add special properties for water"""
    @classmethod
    def properties(cls):
        prop = Thermo.properties()[:]
        l = [
           (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", unidades.Dimensionless),
           (QApplication.translate("pychemqt", "Refractive index"), "n", unidades.Dimensionless)]
        for p in l:
            prop.insert(-11, p)
        return prop

    def writeStatetoJSON(self, state, fase):
        Thermo.writeStatetoJSON(self, state, fase)
        if self._bool:
            state[fase]["n"] = self.n
            state[fase]["epsilon"] = self.epsilon

    def readStatefromJSON(self, fluid):
        Thermo.readStatefromJSON(self, fluid)
        if fluid:
            self.epsilon = unidades.Dimensionless(fluid["epsilon"])
            self.n = unidades.Dimensionless(fluid["n"])


class ThermoAdvanced(Thermo):
    """Custom specified thermo instance to add special properties for advanced
    model as coolprop, refprop and meos"""
    @classmethod
    def properties(cls):
        prop = Thermo.properties()[:]
        l = [
            (QApplication.translate("pychemqt", "Isentropic temperature-pressure"), "betas", unidades.TemperaturePressure),
            (QApplication.translate("pychemqt", "Gruneisen parameter"), "Gruneisen", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "2nd virial coefficient"), "virialB", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "3er virial coefficient"), "virialC", unidades.SpecificVolume_square),
            ("(dp/dT)_rho", "dpdT_rho", unidades.PressureTemperature),
            ("(dp/drho)_T", "dpdrho_T", unidades.PressureDensity),
            ("(drho/dT)_P", "drhodT_P", unidades.DensityTemperature),
            ("(drho/dP)_T", "drhodP_T", unidades.DensityPressure),
            ("(dh/dT)_rho", "dhdT_rho", unidades.SpecificHeat),
            ("(dh/dP)_T", "dhdP_T", unidades.EnthalpyPressure),
            ("(dh/dT)_P", "dhdT_P", unidades.SpecificHeat),
            ("(dh/drho)_T", "dhdrho_T", unidades.EnthalpyDensity),
            ("(dh/dP)_rho", "dhdP_rho", unidades.EnthalpyPressure),
            (QApplication.translate("pychemqt", "Isothermal expansion coefficient"), "kt", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Isentropic expansion coefficient"), "ks", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Adiabatic bulk modulus"), "Ks", unidades.Pressure),
            (QApplication.translate("pychemqt", "Isothermal bulk modulus"), "Kt", unidades.Pressure),
            #        Z_rho     -   (Z-1) over the density, m³/kg
            (QApplication.translate("pychemqt", "Internal pressure"), "IntP", unidades.Pressure),
            (QApplication.translate("pychemqt", "Negative reciprocal temperature"), "invT", unidades.InvTemperature),
            (QApplication.translate("pychemqt", "Specific heat input"), "hInput", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", unidades.Dimensionless)]
        for p in l:
            prop.insert(34, p)
        return prop

    def writeStatetoJSON(self, state, fase):
        Thermo.writeStatetoJSON(self, state, fase)
        if self._bool:
            state[fase]["betas"] = self.betas
            state[fase]["Gruneisen"] = self.Gruneisen
            state[fase]["virialB"] = self.virialB
            state[fase]["virialC"] = self.virialC
            state[fase]["dpdT_rho"] = self.dpdT_rho
            state[fase]["dpdrho_T"] = self.dpdrho_T
            state[fase]["drhodT_P"] = self.drhodT_P
            state[fase]["drhodP_T"] = self.drhodP_T
            state[fase]["dhdT_rho"] = self.dhdT_rho
            state[fase]["dhdP_T"] = self.dhdP_T
            state[fase]["dhdT_P"] = self.dhdT_P
            state[fase]["dhdrho_T"] = self.dhdrho_T
            state[fase]["dhdP_rho"] = self.dhdP_rho
            state[fase]["kt"] = self.kt
            state[fase]["ks"] = self.ks
            state[fase]["Ks"] = self.Ks
            state[fase]["Kt"] = self.Kt
            state[fase]["IntP"] = self.IntP
            state[fase]["invT"] = self.invT
            state[fase]["hInput"] = self.hInput
            state[fase]["epsilon"] = self.epsilon

    def readStatefromJSON(self, fluid):
        Thermo.readStatefromJSON(self, fluid)
        if fluid:
            self.betas = unidades.TemperaturePressure(fluid["betas"])
            self.Gruneisen = unidades.Dimensionless(fluid["Gruneisen"])
            self.virialB = unidades.SpecificVolume(fluid["virialB"])
            self.virialC = unidades.SpecificVolume_square(fluid["virialC"])
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
            self.invT = unidades.InvTemperature(fluid["invT"])
            self.hInput = unidades.Enthalpy(fluid["hInput"])
            self.epsilon = unidades.Dimensionless(fluid["epsilon"])


class ThermoRefProp(ThermoAdvanced):
    """Custom specified thermo instance to add special properties for advanced
    model as coolprop, refprop and meos"""
    @classmethod
    def properties(cls):
        prop = ThermoAdvanced.properties()[:]
        l = [
            (QApplication.translate("pychemqt", "Ideal Pressure"), "P0", unidades.Pressure),
            (QApplication.translate("pychemqt", "Residual Pressure"), "P_Pideal", unidades.Pressure),
            (QApplication.translate("pychemqt", "K value"), "K", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Heat Capacity along the saturation line"), "csat", unidades.SpecificHeat),
            ("dP/dT [sat]", "dpdt_sat", unidades.PressureTemperature),
            (QApplication.translate("pychemqt", "Cv two phases"), "cv2p", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Excess volume"), "vE", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "Excess internal energy"), "eE", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Excess enthalpy"), "hE", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Excess entropy"), "sE", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Excess Helmholtz energy"), "aE", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Excess Gibbs energy"), "gE", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Residual pressure"), "pr", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "Residual internal energy"), "er", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Residual enthalpy"), "hr", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Residual entropy"), "sr", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Residual Helmholtz energy"), "ar", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Residual Gibbs energy"), "gr", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Residual isobaric heat capacity"), "cpr", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Residual isochoric heat capacity"), "cvr", unidades.SpecificHeat),
            (QApplication.translate("pychemqt", "Supercompressibility factor"), "fpv", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Chemical potential"), "chempot", unidades.Enthalpy),
            (QApplication.translate("pychemqt", "Fourth virial coefficient"), "virialD", unidades.Dimensionless),
            (QApplication.translate("pychemqt", "Second acoustic virial coefficient"), "virialBa", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "Third acoustic virial coefficient"), "virialCa", unidades.SpecificVolume_square),
            ("dC/dT", "dCdt", unidades.Dimensionless),
            ("d²C/dT²", "dCdt2", unidades.Dimensionless),
            ("dB/dT", "dBdt", unidades.Dimensionless),
            ("b12", "b12", unidades.SpecificVolume),
            (QApplication.translate("pychemqt", "Critical flow factor"), "cstar", unidades.Dimensionless)]

        for p in l:
            prop.append(p)
        return prop

    def writeStatetoJSON(self, state, fase):
        Thermo.writeStatetoJSON(self, state, fase)
        if self._bool:
            state[fase]["K"] = self.K
            state[fase]["P0"] = self.P0
            state[fase]["P_Pideal"] = self.P_Pideal
            state[fase]["csat"] = self.csat
            state[fase]["dpdt_sat"] = self.dpdt_sat
            state[fase]["cv2p"] = self.cv2p
            state[fase]["vE"] = self.vE
            state[fase]["eE"] = self.eE
            state[fase]["hE"] = self.hE
            state[fase]["sE"] = self.sE
            state[fase]["aE"] = self.aE
            state[fase]["gE"] = self.gE
            state[fase]["pr"] = self.pr
            state[fase]["er"] = self.er
            state[fase]["hr"] = self.hr
            state[fase]["sr"] = self.sr
            state[fase]["ar"] = self.ar
            state[fase]["gr"] = self.gr
            state[fase]["cpr"] = self.cpr
            state[fase]["cvr"] = self.cvr
            state[fase]["fpv"] = self.fpv
            state[fase]["chempot"] = self.chempot
            state[fase]["viriald"] = self.viriald
            state[fase]["virialba"] = self.virialba
            state[fase]["virialca"] = self.virialca
            state[fase]["dcdt"] = self.dcdt
            state[fase]["dcdt2"] = self.dcdt2
            state[fase]["dbdt"] = self.dbdt
            state[fase]["b12"] = self.b12
            state[fase]["cstar"] = self.cstar

    def readStatefromJSON(self, fluid):
        Thermo.readStatefromJSON(self, fluid)
        if fluid:
            self.K = [unidades.Dimensionless(k) for k in fluid["K"]]
            self.P0 = unidades.Pressure(fluid["P0"])
            self.P_Pideal = unidades.Pressure(fluid["P_Pideal"])
            self.csat = unidades.Dimensionless(fluid["csat"])
            self.dpdt_sat = unidades.Dimensionless(fluid["dpdt_sat"])
            self.cv2p = unidades.Dimensionless(fluid["cv2p"])
            self.vE = unidades.SpecificVolume(fluid["vE"])
            self.eE = unidades.Enthapy(fluid["eE"])
            self.hE = unidades.Enthapy(fluid["hE"])
            self.sE = unidades.SpecificHeat(fluid["sE"])
            self.aE = unidades.Enthapy(fluid["aE"])
            self.gE = unidades.Enthapy(fluid["gE"])
            self.pr = unidades.Pressure(fluid["pr"])
            self.er = unidades.Enthapy(fluid["er"])
            self.hr = unidades.Enthapy(fluid["hr"])
            self.sr = unidades.SpecificHeat(fluid["sr"])
            self.ar = unidades.Enthapy(fluid["ar"])
            self.gr = unidades.Enthapy(fluid["gr"])
            self.cpr = unidades.SpecificHeat(fluid["cpr"])
            self.cvr = unidades.SpecificHeat(fluid["cvr"])
            self.fpv = unidades.Dimensionless(fluid["fpv"])
            self.chempot = [unidades.Dimensionless(m) for m in fluid["chempot"]]
            self.viriald = unidades.Dimensionless(fluid["viriald"])
            self.virialba = unidades.SpecificVolume(fluid["virialba"])
            self.virialca = unidades.SpecificVolume_square(fluid["virialca"])
            self.dcdt = unidades.Dimensionless(fluid["dcdt"])
            self.dcdt2 = unidades.Dimensionless(fluid["dcdt2"])
            self.dbdt = unidades.Dimensionless(fluid["dbdt"])
            self.b12 = unidades.SpecificVolume(fluid["b12"])
            self.cstar = unidades.Dimensionless(fluid["cstar"])



