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
# Library with a wrapper class of iapws97 for IAPWS formulation
###############################################################################


from iapws.iapws97 import IAPWS97 as IAPWS
from iapws._iapws import M, Tc, Pc, rhoc, Tt, Tb, Dipole, f_acent

from lib import unidades
from lib.thermo import ThermoWater


class IAPWS97(ThermoWater):
    """Class to model a state for liquid water or steam with the IAPWS-IF97

    Incoming properties:
    T   -   Temperature, K
    P   -   Pressure, Pa
    h   -   Specific enthalpy, kJ/kg
    s   -   Specific entropy, kJ/kg·K
    x   -   Quality

    Optional:
    l   -   Wavelength of light, for refractive index

    Definitions options:
    T, P    Not valid for two-phases region
    P, h
    P, s
    h, s
    T, x    Only for two-phases region
    P, x    Only for two-phases region


    Properties:
    P        -   Pressure, MPa
    T        -   Temperature, K
    g        -   Specific Gibbs free energy, kJ/kg
    a        -   Specific Helmholtz free energy, kJ/kg
    v        -   Specific volume, m³/kg
    r        -   Density, kg/m³
    h        -   Specific enthalpy, kJ/kg
    u        -   Specific internal energy, kJ/kg
    s        -   Specific entropy, kJ/kg·K
    c        -   Specific isobaric heat capacity, kJ/kg·K
    c        -   Specific isochoric heat capacity, kJ/kg·K
    Z        -   Compression factor
    f        -   Fugacity, MPa
    gamma    -   Isoentropic exponent
    alfav    -   Isobaric cubic expansion coefficient, 1/K
    kappa   -   Isothermal compressibility, 1/MPa
    alfap    -   Relative pressure coefficient, 1/K
    betap    -   Isothermal stress coefficient, kg/m³
    joule    -   Joule-Thomson coefficient, K/MPa
    deltat   -   Isothermal throttling coefficient, kJ/kg·MPa
    region   -   Region

    v0       -   Ideal specific volume, m³/kg
    u0       -   Ideal specific internal energy, kJ/kg
    h0       -   Ideal specific enthalpy, kJ/kg
    s0       -   Ideal specific entropy, kJ/kg·K
    a0       -   Ideal specific Helmholtz free energy, kJ/kg
    g0       -   Ideal specific Gibbs free energy, kJ/kg
    cp0      -   Ideal specific isobaric heat capacity, kJ/kg·K
    cv0      -   Ideal specific isochoric heat capacity, kJ/kg·K
    w0       -   Ideal speed of sound, m/s
    gamma0   -   Ideal isoentropic exponent

    w        -   Speed of sound, m/s
    mu       -   Dynamic viscosity, Pa·s
    nu       -   Kinematic viscosity, m²/s
    k        -   Thermal conductivity, W/m·K
    alfa     -   Thermal diffusivity, m²/s
    sigma    -   Surface tension, N/m
    epsilon  -   Dielectric constant
    n        -   Refractive index
    Prandt   -   Prandtl number
    Pr       -   Reduced Pressure
    Tr       -   Reduced Temperature

    Usage:
    >>> water=IAPWS97(T=170+273.15,x=0.5)
    >>> "%0.4f %0.4f %0.1f %0.2f" % (water.Liquido.cp.kJkgK, \
        water.Gas.cp.kJkgK, water.Liquido.w, water.Gas.w)
    '4.3695 2.5985 1418.3 498.78'
    >>> water=IAPWS97(T=325+273.15,x=0.5)
    >>> "%0.4f %0.8f %0.7f %0.2f %0.2f" % (water.P.MPa, water.Liquido.v, \
        water.Gas.v, water.Liquido.h.kJkg, water.Gas.h.kJkg)
    '12.0505 0.00152830 0.0141887 1493.37 2684.48'
    >>> water=IAPWS97(T=50+273.15,P=611.2127)
    >>> "%0.4f %0.4f %0.2f %0.3f %0.2f" % (water.cp0.kJkgK, water.cv0.kJkgK, \
        water.h0.kJkg, water.s0.kJkgK, water.w0)
    '1.8714 1.4098 2594.66 9.471 444.93'
    """

    __doi__ = [
        {"autor": "IAPWS",
         "title": "Revised Release on the IAPWS Industrial Formulation 1997"
                  "for the Thermodynamic Properties of Water and Steam",
         "ref": "",
         "doi": ""},
        {"autor": "Wagner, W., Cooper, J. R., Dittmann, A., Kijima, J., "
                  "Kretzschmar, H.-J., Kruse, A., Mareš, R., Oguchi, K., Sato,"
                  " H., Stöcker, I., Šifner, O., Takaishi, Y., Tanishita, I., "
                  "Trübenbach, J., and Willkommen, Th.",
         "title": "The IAPWS Industrial Formulation 1997 for the Thermodynamic"
                  " Properties of Water and Steam",
         "ref": "J. Eng. Gas Turbines & Power 122, 150-182 (2000)",
         "doi": "10.1115/1.483186"}]

    M = unidades.Dimensionless(M)
    Pc = unidades.Pressure(Pc, "MPa")
    Tc = unidades.Temperature(Tc)
    rhoc = unidades.Density(rhoc)
    Tt = unidades.Temperature(Tt)
    Tb = unidades.Temperature(Tb)
    f_accent = unidades.Dimensionless(f_acent)
    momentoDipolar = unidades.DipoleMoment(Dipole, "Debye")

    def __init__(self, **kwargs):
        if "P" in kwargs:
            kwargs["P"] /= 1e6
        elif "h" in kwargs:
            kwargs["h"] /= 1e3
        elif "s" in kwargs:
            kwargs["s"] /= 1e3

        st = IAPWS(**kwargs)
        self.status = st.status
        self.msg = st.msg
        if self.status:
            self.calculo(st)

    def calculo(self, st):
        self.x = unidades.Dimensionless(st.x)
        self.region = st.region
        self.phase = self.getphase(phase=st.phase)
        self.name = st.name
        self.synonim = st.synonim
        self.CAS = st.CAS

        self.T = unidades.Temperature(st.T)
        self.P = unidades.Pressure(st.P, "MPa")
        self.Tr = unidades.Dimensionless(st.Tr)
        self.Pr = unidades.Dimensionless(st.Pr)
        self.v = unidades.SpecificVolume(st.v)
        self.rho = unidades.Density(st.rho)

        cp0 = {}
        cp0["v"] = st.v0
        cp0["h"] = st.h0*1000
        cp0["s"] = st.s0*1000
        cp0["cp"] = st.cp0*1000
        cp0["cv"] = st.cv0*1000
        cp0["w"] = st.w0
        self._cp0(cp0)

        self.Liquido = ThermoWater()
        self.Gas = ThermoWater()
        if self.x == 0:
            # only liquid phase
            self.fill(self, st.Liquid)
            self.fill(self.Liquido, st.Liquid)
            self.sigma = unidades.Tension(st.sigma)
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)

        elif self.x == 1:
            # only vapor phase
            self.fill(self, st.Vapor)
            self.fill(self.Gas, st.Vapor)
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)

        else:
            # two phases
            self.fill(self.Liquido, st.Liquid)
            self.sigma = unidades.Tension(st.sigma)
            self.fill(self.Gas, st.Vapor)

            self.h = unidades.Enthalpy(st.h)
            self.s = unidades.SpecificHeat(st.s)
            self.u = unidades.SpecificHeat(st.u)
            self.a = unidades.Enthalpy(st.a)
            self.g = unidades.Enthalpy(st.g)

            self.cv = unidades.SpecificHeat(None)
            self.cp = unidades.SpecificHeat(None)
            self.cp_cv = unidades.Dimensionless(None)
            self.w = unidades.Speed(None)

            self.Hvap = unidades.Enthalpy(st.Hvap, "kJkg")
            self.Svap = unidades.SpecificHeat(st.Svap, "kJkgK")

    def fill(self, fase, st):
        """Fill phase properties"""
        fase._bool = True
        fase.M = self.M
        fase.v = unidades.SpecificVolume(st.v)
        fase.rho = unidades.Density(st.rho)
        fase.Z = unidades.Dimensionless(st.Z)

        fase.h = unidades.Enthalpy(st.h, "kJkg")
        fase.s = unidades.SpecificHeat(st.s, "kJkgK")
        fase.u = unidades.Enthalpy(st.u, "kJkg")
        fase.a = unidades.Enthalpy(st.a, "kJkg")
        fase.g = unidades.Enthalpy(st.g, "kJkg")
        fase.fi = [unidades.Dimensionless(st.fi)]
        fase.f = [unidades.Pressure(st.f, "MPa")]

        fase.cv = unidades.SpecificHeat(st.cv, "kJkgK")
        fase.cp = unidades.SpecificHeat(st.cp, "kJkgK")
        fase.cp_cv = unidades.Dimensionless(st.cp_cv)
        fase.gamma = fase.cp_cv
        fase.w = unidades.Speed(st.w)

        fase.rhoM = unidades.MolarDensity(fase.rho/self.M)
        fase.hM = unidades.MolarEnthalpy(fase.h*self.M)
        fase.sM = unidades.MolarSpecificHeat(fase.s*self.M)
        fase.uM = unidades.MolarEnthalpy(fase.u*self.M)
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(fase.cv*self.M)
        fase.cpM = unidades.MolarSpecificHeat(fase.cp*self.M)

        fase.alfav = unidades.InvTemperature(st.alfav)
        fase.kappa = unidades.InvPressure(st.xkappa, "MPa")
        fase.kappas = unidades.InvPressure(st.kappas, "MPa")

        fase.mu = unidades.Viscosity(st.mu)
        fase.nu = unidades.Diffusivity(st.nu)
        fase.k = unidades.ThermalConductivity(st.k)
        fase.alfa = unidades.Diffusivity(st.alfa)
        fase.epsilon = unidades.Dimensionless(st.epsilon)
        fase.Prandt = unidades.Dimensionless(st.Prandt)
        fase.n = unidades.Dimensionless(st.n)

        fase.joule = unidades.TemperaturePressure(st.joule)
        fase.deltat = unidades.EnthalpyPressure(st.deltat)

        fase.betap = unidades.Density(st.betap)
        fase.alfap = unidades.Density(st.alfap)
        fase.fraccion = [1]
        fase.fraccion_masica = [1]


class IAPWS97_PT(IAPWS97):
    """Derivated class for direct P and T input"""
    def __init__(self, P, T):
        IAPWS97.__init__(self, T=T, P=P)


class IAPWS97_Ph(IAPWS97):
    """Derivated class for direct P and h input"""
    def __init__(self, P, h):
        IAPWS97.__init__(self, P=P, h=h)


class IAPWS97_Ps(IAPWS97):
    """Derivated class for direct P and s input"""
    def __init__(self, P, s):
        IAPWS97.__init__(self, P=P, s=s)


class IAPWS97_Px(IAPWS97):
    """Derivated class for direct P and x input"""
    def __init__(self, P, x):
        IAPWS97.__init__(self, P=P, x=x)


class IAPWS97_Tx(IAPWS97):
    """Derivated class for direct T and x input"""
    def __init__(self, T, x):
        IAPWS97.__init__(self, T=T, x=x)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
