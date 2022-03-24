#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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

try:
    from iapws import __doi__
except ImportError:
    __doi__ = {}

__doi__["1"] = {
    "autor": "Wagner, W., Cooper, J.R., Dittmann, A., Kijima, J., "
             "Kretzschmar, H.-J., Kruse, A., Mareš, R., Oguchi, K., Sato,"
             " H., Stöcker, I., Šifner, O., Takaishi, Y., Tanishita, I., "
             "Trübenbach, J., Willkommen, T.",
    "title": "The IAPWS Industrial Formulation 1997 for the Thermodynamic "
             "Properties of Water and Steam",
    "ref": "J. Eng. Gas Turbines & Power 122 (2000) 150-182",
    "doi": "10.1115/1.483186"}


class IAPWS97(ThermoWater):
    """Class to model a state for liquid water or steam with the IAPWS-IF97

    Parameters
    ----------
    T : float
        Temperature [K]
    P : float
        Pressure [MPa]
    h : float
        Specific enthalpy [kJ/kg]
    s : float
        Specific entropy [kJ/kgK]
    x : float
        Vapor quality [-]
    l : float, optional
        Wavelength of light, for refractive index [nm]

    Notes
    -----
    Definitions options:
        * T, P: Not valid for two-phases region
        * P, h
        * P, s
        * h, s
        * T, x: Only for two-phases region
        * P, x: Only for two-phases region

    Returns
    -------
    The calculated instance has the following properties:
        * P: Pressure [MPa]
        * T: Temperature [K]
        * g: Specific Gibbs free energy [kJ/kg]
        * a: Specific Helmholtz free energy [kJ/kg]
        * v: Specific volume [m³/kg]
        * rho: Density [kg/m³]
        * h: Specific enthalpy [kJ/kg]
        * u: Specific internal energy [kJ/kg]
        * s: Specific entropy [kJ/kg·K]
        * cp: Specific isobaric heat capacity [kJ/kg·K]
        * cv: Specific isochoric heat capacity [kJ/kg·K]
        * Z: Compression factor [-]
        * fi: Fugacity coefficient [-]
        * f: Fugacity [MPa]

        * gamma: Isoentropic exponent [-]
        * alfav: Isobaric cubic expansion coefficient [1/K]
        * xkappa: Isothermal compressibility [1/MPa]
        * kappas: Adiabatic compresibility [1/MPa]
        * alfap: Relative pressure coefficient [1/K]
        * betap: Isothermal stress coefficient [kg/m³]
        * joule: Joule-Thomson coefficient [K/MPa]
        * deltat: Isothermal throttling coefficient [kJ/kg·MPa]
        * region: Region

        * v0: Ideal specific volume [m³/kg]
        * u0: Ideal specific internal energy [kJ/kg]
        * h0: Ideal specific enthalpy [kJ/kg]
        * s0: Ideal specific entropy [kJ/kg·K]
        * a0: Ideal specific Helmholtz free energy [kJ/kg]
        * g0: Ideal specific Gibbs free energy [kJ/kg]
        * cp0: Ideal specific isobaric heat capacity [kJ/kg·K]
        * cv0: Ideal specific isochoric heat capacity [kJ/kg·K]
        * w0: Ideal speed of sound [m/s]
        * gamma0: Ideal isoentropic exponent [-]

        * w: Speed of sound [m/s]
        * mu: Dynamic viscosity [Pa·s]
        * nu: Kinematic viscosity [m²/s]
        * k: Thermal conductivity [W/m·K]
        * alfa: Thermal diffusivity [m²/s]
        * sigma: Surface tension [N/m]
        * epsilon: Dielectric constant [-]
        * n: Refractive index [-]
        * Prandt: Prandtl number [-]
        * Pr: Reduced Pressure [-]
        * Tr: Reduced Temperature [-]

    Examples
    --------
    Region 1, Table A3
    >>> st = IAPWS97(T=300, P=3e6)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '0.00100215168 115.331273 112.324818'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '0.392294792 4.17301218 1507.73921'
    >>> st = IAPWS97(T=300, P=8e7)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '0.000971180894 184.142828 106.448356'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '0.368563852 4.01008987 1634.69054'
    >>> st = IAPWS97(T=500, P=3e6)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '0.001202418 975.542239 971.934985'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '2.58041912 4.65580682 1240.71337'

    Region 2, Table A6
    >>> st = IAPWS97(T=300, P=3.5e3)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '39.4913866 2549.91145 2411.6916'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '8.52238967 1.91300162 427.920172'
    >>> st = IAPWS97(T=700, P=3.5e3)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '92.3015898 3335.68375 3012.62819'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '10.1749996 2.08141274 644.289068'
    >>> st = IAPWS97(T=700, P=3e7)
    >>> "%0.9g %0.9g %0.9g" % (st.v, st.h.kJkg, st.u.kJkg)
    '0.00542946619 2631.49474 2468.61076'
    >>> "%0.9g %0.9g %0.9g" % (st.s.kJkgK, st.cp.kJkgK, st.w)
    '5.17540298 10.3505092 480.386523'

    Region 3, Table A10
    >>> from iapws.iapws97 import _Region3
    >>> st = _Region3(500, 650)
    >>> "%0.9g %0.9g %0.9g" % (st["P"], st["h"], st["h"]-st["P"]*1000*st["v"])
    '25.5837018 1863.43019 1812.26279'
    >>> "%0.9g %0.9g %0.9g" % (st["s"], st["cp"], st["w"])
    '4.05427273 13.8935717 502.005554'
    >>> st = _Region3(200, 650)
    >>> "%0.9g %0.9g %0.9g" % (st["P"], st["h"], st["h"]-st["P"]*1000*st["v"])
    '22.2930643 2375.12401 2263.65868'
    >>> "%0.9g %0.9g %0.9g" % (st["s"], st["cp"], st["w"])
    '4.85438792 44.6579342 383.444594'
    >>> st = _Region3(500, 750)
    >>> "%0.9g %0.9g %0.9g" % (st["P"], st["h"], st["h"]-st["P"]*1000*st["v"])
    '78.3095639 2258.68845 2102.06932'
    >>> "%0.9g %0.9g %0.9g" % (st["s"], st["cp"], st["w"])
    '4.46971906 6.34165359 760.696041'

    Region 4, Table A12
    >>> st1 = IAPWS97(T=300, x=0.5)
    >>> st2 = IAPWS97(T=500, x=0.5)
    >>> st3 = IAPWS97(T=600, x=0.5)
    >>> "%0.9g %0.9g %0.9g" % (st1.P.MPa, st2.P.MPa, st3.P.MPa)
    '0.00353658941 2.63889776 12.3443146'

    Region 4, Table A29
    >>> st1 = IAPWS97(P=1e5, x=0.5)
    >>> st2 = IAPWS97(P=1e6, x=0.5)
    >>> st3 = IAPWS97(P=1e7, x=0.5)
    >>> "%0.9g %0.9g %0.9g" % (st1.T, st2.T, st3.T)
    '372.755919 453.035632 584.149488'

    Other test in paper have been upgraded by new equation, i.e. for region 5
    or for ancillary equation, the testing is done in iapws testing layer

    References
    ----------
    [1]_ Wagner, W., Cooper, J.R., Dittmann, A., Kijima, J., Kretzschmar,
        H.-J., Kruse, A., Mareš, R., Oguchi, K., Sato, H., Stöcker, I.,
        Šifner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., Willkommen, T.
        The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties
        of Water and Steam. J. Eng. Gas Turbines & Power 122 (2000) 150-182.
    """

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

        self._cp0()

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

    def _cp0(self):
        "Set ideal properties to state"""
        self.v0 = unidades.SpecificVolume(None)
        self.rho0 = unidades.Density(None)
        self.h0 = unidades.Enthalpy(None)
        self.u0 = unidades.Enthalpy(None)
        self.s0 = unidades.SpecificHeat(None)
        self.a0 = unidades.Enthalpy(None)
        self.g0 = unidades.Enthalpy(None)

        self.cp0 = unidades.SpecificHeat(None)
        self.cv0 = unidades.SpecificHeat(None)
        self.cp0_cv = unidades.Dimensionless(None)
        self.w0 = unidades.Speed(None)
        self.gamma0 = self.cp0_cv

        self.rhoM0 = unidades.MolarDensity(None)
        self.hM0 = unidades.MolarEnthalpy(None)
        self.uM0 = unidades.MolarEnthalpy(None)
        self.sM0 = unidades.MolarSpecificHeat(None)
        self.aM0 = unidades.MolarEnthalpy(None)
        self.gM0 = unidades.MolarEnthalpy(None)
        self.cpM0 = unidades.MolarSpecificHeat(None)
        self.cvM0 = unidades.MolarSpecificHeat(None)


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
