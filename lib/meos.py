#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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


from itertools import product

from PyQt5.QtWidgets import QApplication
from scipy import exp, log, sinh, cosh, tanh, arctan
from scipy.constants import Boltzmann, pi, Avogadro, R, u
from scipy.optimize import fsolve

from lib import unidades
from lib.utilities import SimpleEq
from lib.physics import R_atml, Collision_Neufeld
from lib.thermo import ThermoAdvanced
from lib.compuestos import RhoL_Costald, Pv_Lee_Kesler, MuG_Chung, MuG_P_Chung
from lib.compuestos import ThG_Chung, ThG_P_Chung, Tension_Pitzer


__doi__ = {
    1:
        {"autor": "Mulero, A., Cachadiña, I., Parra, M.I.",
         "title": "Recommended Correlations for the Surface Tension of "
                  "Common Fluids",
         "ref": "J. Phys. Chem. Ref. Data 41(4) (2012) 043105",
         "doi": "10.1063/1.4768782"},
    2:
        {"autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
         "title": "Generalization of the Friction Theory for Viscosity "
                  "Modeling",
         "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
         "doi": "10.1021/jp0618577"},
    3:
        {"autor": "Quiñones-Cisneros, S.E., Huber, M.L., Deiters, U.K.",
         "title": "Correlation for the Viscosity of Sulfur Hexafluoride (SF6) "
                  "from the Triple Point to 1000 K and Pressures to 50 MPa",
         "ref": "J. Phys. Chem. Ref. Data 41(2) (2012) 023102",
         "doi": "10.1063/1.3702441"},
    4:
        {"autor": "Younglove, B.A.",
         "title": "Thermophysical Properties of Fluids. I. Argon, Ethylene, "
                  "Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen",
         "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
         "doi": ""},
    5:
        {"autor": "Younglove, B.A., Ely, J.F.",
         "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, "
                  "Propane, Isobutane, and Normal Butane",
         "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
         "doi": "10.1063/1.555785"},
    6:
        {"autor": "Lemmon, E.W., Jacobsen, R.T.",
         "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, "
                  "Oxygen, Argon, and Air",
         "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
         "doi": "10.1023/B:IJOT.0000022327.04529.f3"},
    7:
        {"autor": "Tariq, U., Jusoh, A.R.B., Riesco, N., Vesovic, V.",
         "title": "Reference Correlation of the Viscosity of Cyclohexane from "
                  " the Triple Point to 700K and up to 110 MPa",
         "ref": "J. Phys. Chem. Ref. Data 43(3) (2014) 033101",
         "doi": "10.1063/1.4891103"},
    8:
        {"autor": "Neufeld, P.D., Janzen, A.R., Aziz, R.A.",
         "title": "Empirical Equations to Calculate 16 of the Transport "
                  "Collision Integrals Ω for the Lennard-Jones Potential",
         "ref": "J. Chem. Phys. 57(3) (1972) 1100-1102",
         "doi": "10.1063/1.1678363"},

    9:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
    # 1:
        # {"autor": "Reeves, L.E., Scott, G.J., Babb, S.E., Jr.",
         # "title": "Melting Curves of Pressure‐Transmitting Fluids ",
         # "ref": "J. Chem. Phys. 40, 3662 (1964)",
         # "doi": "10.1063/1.1725068"},
    # 2:
        # {"autor": "Harvey, A.H. and Lemmon, E.W.",
         # "title": "Method for Estimating the Dielectric Constant of "
                  # "Natural Gas Mixtures ",
         # "ref": "Int. J. Thermophys., 26(1):31-46, 2005.",
         # "doi": "10.1007/s10765-005-2351-5"},

    10:
        {"autor": "Olchowy, G.A., Sengers, J.V.",
         "title": "A Simplified Representation for the Thermal Conductivity "
                  "of Fluids in the Critical Region",
         "ref": "Int. J. Thermophys. 10(2) (1989) 417-426",
         "doi": "10.1007/bf01133538"},

    49:
        {"autor": "Chung, T.H., Ajlan, M., Lee, L.L., Starling, K.E.",
         "title": "Generalized Multiparameter Correlation for Nonpolar and "
                  "Polar Fluid Transport Properties",
         "ref": "Ind. Eng. Chem. Res. 27(4) (1988) 671-679",
         "doi": "10.1021/ie00076a024"},
    50:
        {"autor": "Chung, T.H., Lee, L.L., Starling, K.E.",
         "title": "Applications of Kinetic Gas Theories and Multiparameter "
                  "Correlation for Prediction of Dilute Gas Viscosity and "
                  "Thermal Conductivity",
         "ref": "Ind. Eng. Chem. Fundam. 23(1) (1984) 8-13",
         "doi": "10.1021/i100013a002"},

        }


class MEoS(ThermoAdvanced):
    r"""General class for implement multiparameter equation of state
    Each child class must define the parameters for the calculations

    Compound definition:
        * name: Name of component
        * CASNumber: CAS Number of component
        * formula: Empiric formula
        * synonym: Alternate formula (Refrigerant name)
        * id: index of component in pychemqt database
        * _refPropName = Codename of compound in RefProp
        * _coolPropName = Coedename of compound in coolProp

    Physical properties:
        * rhoc: Critical density, [kg/m³]
        * Tc: Critical temperature, [K]
        * Pc: Critical pressure, [Pa]
        * M: Molecular weigth, [g/mol]
        * Tt: Temperature of triple point, [K]
        * Tb: Normal boiling point temperature, [K]
        * f_acent: Acentric factor, [-]
        * momentoDipolar: Depole moment, [C·m]

    Parameters of mEoS and transport correlation:
        * eq: List of pointer for mEoS correlations
        * _viscosity: List of pointer for viscosity correlations
        * _thermal: List of pointer for thermal conductivity correlations

    Other properties correlations:
        * _vapor_Pressure: Parameters for vapor pressure ancillary equation
        * _liquid_Density: Parameters for liquid density ancillary equation
        * _vapor_Density: Parameters for vapor density ancillary equation
        * _dielectric: Parameters for dielectric constant calculation
        * _melting: Parameters for melting line calculation
        * _sublimation: Parameters for sublimation line calculation
        * _surface: Parameters for surface tension calculation

    Parameters necessary only for special EoS:
        * _PR: Peneloux volume correction for Peng-Robinson equation of state
        * _Tr: Temperature parameter for generalized equation
        * _rhor: Density parameter for generalized equation
        * _w: Acentric factor for generalized equation
    """

    id = None
    _refPropName = ""
    _coolPropName = ""
    _Tr = None
    _rhor = None
    _w = None

    eq = ()
    _PR = 0.

    _dielectric = None
    _melting = None
    _sublimation = None
    _surface = None
    _vapor_Pressure = None
    _liquid_Density = None
    _vapor_Density = None

    _omega = None
    _viscosity = None
    _thermal = None
    _critical = None

    _test = []

    kwargs = {"T": 0.0,
              "P": 0.0,
              "rho": None,
              "h": None,
              "s": None,
              "u": None,
              "x": None,
              "v": 0.0,

              "eq": 0,
              "visco": 0,
              "thermal": 0,
              "ref": None,
              "refvalues": None,
              "rho0": 0,
              "T0": 0}
    status = 0
    msg = QApplication.translate("pychemqt", "Unknown Variables")

    def __init__(self, **kwargs):
        """Incoming properties:
        T   -   Temperature, K
        P   -   Pressure, MPa
        rho -   Density, kg/m3
        v   -   Specific volume, m3/kg
        h   -   Specific enthalpy, kJ/kg
        s   -   Specific entropy, kJ/kg·K
        u   -   Specific internal energy, kJ/kg·K
        x   -   Quality
        It needs two incoming properties

        eq: Equation to use:
            number: Indice
            "PR": Peng Robinson Cubic equation
            "Generalised": Generalised Multiparameter equation Span and Wagner
            "GERG": Model for mixtures

        visco: Equation for viscosity calculation
        thermal: Equation for thermal conductivity
        ref: Código con el estado de referencia
            OTO:  h,s=0 at 25ºC and 1 atm
            NBP:  h,s=0 saturated liquid at Tb
            IIR:  h=200,s=1 saturated liquid 0ºC
            ASHRAE:  h,s=0 saturated liquid at -40ºC
            CUSTOM
        refvalues: array with custom values of reference state
            [Tref, Pref, ho, so]
        rho0: Initial value for iteration over density
        T0: Initial value for iteration over temperature

    Calculated properties:
        P         -   Pressure, MPa
        Pr        -   Reduce pressure
        T         -   Temperature, K
        Tr        -   Reduced temperature
        x         -   Quality
        v         -   Specific volume, m³/kg
        rho       -   Density, kg/m³
        h         -   Specific enthalpy, kJ/kg
        s         -   Specific entropy, kJ/kg·K
        u         -   Specific internal energy, kJ/kg
        g         -   Specific Gibbs free energy, kJ/kg
        a         -   Specific Helmholtz free energy, kJ/kg
        cp        -   Specific isobaric heat capacity, kJ/kg·K
        cv        -   Specific isochoric heat capacity, kJ/kg·K
        cp_cv     -   Heat capacity ratio
        w         -   Speed of sound, m/s
        Z         -   Compression factor
        fi        -   Fugacity coefficient
        f         -   Fugacity, MPa
        gamma     -   Isoentropic exponent
        Hvap      -   Vaporization heat, kJ/kg
        alfav     -   Thermal expansion coefficient (Volume expansivity), 1/K
        kappa     -   Isothermal compressibility, 1/MPa
        alfap     -   Relative pressure coefficient, 1/K
        betap     -   Isothermal stress coefficient, kg/m³
        betas     -   Isoentropic temperature-pressure coefficient
        joule     -   Joule-Thomson coefficient, K/MPa
        Gruneisen -   Gruneisen parameter
        virialB   -   Second virial coefficient, m³/kg
        virialC   -   Third virial coefficient, m⁶/kg²
        dpdT_rho  -   Derivatives, dp/dT at constant rho, MPa/K
        dpdrho_T  -   Derivatives, dp/drho at constant T, MPa·m³/kg
        drhodT_P  -   Derivatives, drho/dT at constant P, kg/m³·K
        drhodP_T  -   Derivatives, drho/dP at constant T, kg/m³·MPa
        dhdT_rho  -   Derivatives, dh/dT at constant rho, kJ/kg·K
        dhdP_T    -   Isothermal throttling coefficient, kJ/kg·MPa
        dhdT_P    -   Derivatives, dh/dT at constant P, kJ/kg·K
        dhdrho_T  -   Derivatives, dh/drho at constant T, kJ·m³/kg²
        dhdrho_P  -   Derivatives, dh/drho at constant P, kJ·m³/kg²
        dhdP_rho  -   Derivatives, dh/dP at constant rho, kJ/kg·MPa
        kt        -   Isothermal Expansion Coefficient
        ks        -   Adiabatic Expansion Coefficient
        # ks        -   Adiabatic Compressibility, 1/MPa
        Ks        -   Adiabatic bulk modulus, MPa
        Kt        -   Isothermal bulk modulus, MPa

        Z_rho     -   (Z-1) over the density, m³/kg
        IntP      -   Internal pressure
        invT      -   Negative reciprocal temperature
        hInput    -   Specific heat input, kJ/kg
        mu        -   Dynamic viscosity, Pa·s
        nu        -   Kinematic viscosity, m²/s
        k         -   Thermal conductivity, W/m·K
        sigma     -   Surface tension, N/m
        alfa      -   Thermal diffusivity, m²/s
        Pramdt    -   Prandtl number
        epsilon   -   Dielectric constant

        v0        -   Ideal gas Specific volume, m³/kg
        rho0      -   Ideal gas Density, kg/m³
        h0        -   Ideal gas Specific enthalpy, kJ/kg
        u0        -   Ideal gas Specific internal energy, kJ/kg
        s0        -   Ideal gas Specific entropy, kJ/kg·K
        a0        -   Ideal gas Specific Helmholtz free energy, kJ/kg
        g0        -   Ideal gas Specific Gibbs free energy, kJ/kg
        cp0       -   Ideal gas Specific isobaric heat capacity, kJ/kg·K
        cv0       -   Ideal gas Specific isochoric heat capacity, kJ/kg·K
        cp0_cv    -   Ideal gas Heat capacity ratio
        gamma0    -   Ideal gas Isoentropic exponent
        """

        self.kwargs = MEoS.kwargs.copy()
        self.__call__(**kwargs)

        # Define general documentation
        if self._surface and "__doi__" not in self._surface:
            self._surface["__doi__"] = __doi__[1]

    def __call__(self, **kwargs):
        """Make instance callable to let definition one parameters for one"""
        # Let user refer to equation using the internal name of equation
        eq = kwargs.get("eq", 0)
        if isinstance(eq, str) and eq in self.__class__.__dict__:
            eq = self.eq.index(self.__class__.__dict__[eq])
            kwargs["eq"] = eq
        elif isinstance(eq, str) and eq not in self.__class__.__dict__ and \
                eq not in ["PR", "GERG", "Generalised"]:
            raise ValueError("Unknown equation")

        self.kwargs.update(kwargs)

        if eq == "PR":
            self._eq = self._PengRobinson
            self._constants = self.eq[0]
        elif eq == "Generalised":
            self._eq = self._Helmholtz
            self._Generalised()
        elif eq == "GERG":
            try:
                self._constants = self.GERG
            except:
                self._constants = self.eq[0]
            if self._constants["__type__"] == "Helmholtz":
                self._eq = self._Helmholtz
            else:
                self._eq = self._MBWR
        elif self.eq[eq]["__type__"] == "Helmholtz":
            self._eq = self._Helmholtz
            self._constants = self.eq[eq]
        elif self.eq[eq]["__type__"] == "MBWR":
            self._eq = self._MBWR
            self._constants = self.eq[eq]
        elif self.eq[eq]["__type__"] == "ECS":
            self._eq = self._ECS
            self._constants = self.eq[eq]

        visco = self.kwargs["visco"]
        thermal = self.kwargs["thermal"]
        if self._viscosity:
            self._viscosity = self._viscosity[visco]
        if self._thermal:
            self._thermal = self._thermal[thermal]

        # Configure custom parameter from eq
        if "M" in self._constants:
            self.M = self._constants["M"]
        if "Tc" in self._constants:
            self.Tc = unidades.Temperature(self._constants["Tc"])
        if "Pc" in self._constants:
            self.Pc = unidades.Pressure(self._constants["Pc"], "kPa")
        if "rhoc" in self._constants:
            self.rhoc = unidades.Density(self._constants["rhoc"]*self.M)
        if "Tt" in self._constants:
            self.Tt = unidades.Temperature(self._constants["Tt"])

        self.R = unidades.SpecificHeat(self._constants["R"]/self.M, "kJkgK")
        self.Zc = self.Pc/self.rhoc/self.R/self.Tc
        self.cleanOldValues(**kwargs)

        if self.calculable:
            self.calculo()
            if self.status in (1, 3):
                converge = True
                for input in self._mode.split("-"):
                    inp = self.kwargs[input]
                    out = self.__getattribute__(input)._data

                    if abs(inp-out) > 1e-9:
                        converge = False
                        break

                if not converge:
                    self.status = 5
                    self.msg = QApplication.translate(
                            "pychemqt", "Solution don´t converge")
                    err = self.kwargs[input]-self.__getattribute__(input)._data
                    print("State don't converge for %s by %g" % (input, err))

    def cleanOldValues(self, **kwargs):
        """Convert alternative input parameters"""
        # Alternative rho input
        if "rhom" in kwargs:
            kwargs["rho"] = kwargs["rhom"]*self.M
            del kwargs["rhom"]
        elif kwargs.get("v", 0):
            kwargs["rho"] = 1./kwargs["v"]
            del kwargs["v"]
        elif kwargs.get("vm", 0):
            kwargs["rho"] = self.M/kwargs["vm"]
            del kwargs["vm"]
        self.kwargs.update(kwargs)

    @property
    def calculable(self):
        """Check if instance has enough input to be calculated"""

        # Define _mode variable to know in calculo procedure the input pair
        self._mode = ""
        if self.kwargs["T"] and self.kwargs["P"]:
            self._mode = "T-P"
        elif self.kwargs["T"] and self.kwargs["rho"] is not None:
            self._mode = "T-rho"
        elif self.kwargs["T"] and self.kwargs["h"] is not None:
            self._mode = "T-h"
        elif self.kwargs["T"] and self.kwargs["s"] is not None:
            self._mode = "T-s"
        elif self.kwargs["T"] and self.kwargs["u"] is not None:
            self._mode = "T-u"
        elif self.kwargs["P"] and self.kwargs["rho"] is not None:
            self._mode = "P-rho"
        elif self.kwargs["P"] and self.kwargs["h"] is not None:
            self._mode = "P-h"
        elif self.kwargs["P"] and self.kwargs["s"] is not None:
            self._mode = "P-s"
        elif self.kwargs["P"] and self.kwargs["u"] is not None:
            self._mode = "P-u"
        elif self.kwargs["rho"] is not None and self.kwargs["h"] is not None:
            self._mode = "rho-h"
        elif self.kwargs["rho"] is not None and self.kwargs["s"] is not None:
            self._mode = "rho-s"
        elif self.kwargs["rho"] is not None and self.kwargs["u"] is not None:
            self._mode = "rho-u"
        elif self.kwargs["h"] is not None and self.kwargs["s"] is not None:
            self._mode = "h-s"
        elif self.kwargs["h"] is not None and self.kwargs["u"] is not None:
            self._mode = "h-u"
        elif self.kwargs["s"] is not None and self.kwargs["u"] is not None:
            self._mode = "s-u"

        elif self.kwargs["T"] and self.kwargs["x"] is not None:
            self._mode = "T-x"
        elif self.kwargs["P"] and self.kwargs["x"] is not None:
            self._mode = "P-x"

        return bool(self._mode)

    def calculo(self):
        """Calculate procedure"""

        # Get input parameters
        T = self.kwargs["T"]
        rho = self.kwargs["rho"]
        P = self.kwargs["P"]
        s = self.kwargs["s"]
        h = self.kwargs["h"]
        u = self.kwargs["u"]
        x = self.kwargs["x"]

        # Define reference state
        ref = self.kwargs["ref"]
        refvalues = self.kwargs["refvalues"]
        self._ref(ref, refvalues)

        propiedades = None

        if x is None:
            # Methods with iteration necessary to get x
            if self._mode == "T-P":
                rhoo = []
                if self.kwargs["rho0"]:
                    rhoo.append(self.kwargs["rho0"])
                if T < self.Tc:
                    Pv = self._Vapor_Pressure(T)
                    rhov = self._Vapor_Density(T)
                    rhol = self._Liquid_Density(T)
                    if P > Pv:
                        rhoo.append(rhol)
                    else:
                        rhoo.append(rhov)
                else:
                    rhoo.append(self._Liquid_Density(self.Tc))

                rhoo.append(self._constants["rhomax"]*self.M)
                rhoo.append(self.rhoc)
                rhoo.append(P/T/self.R)

                def f(rho):
                    return self._eq(rho, T)["P"]-P
                rho = self.fsolve(f, False, **{"P": P, "T": T, "rho0": rhoo})

            elif self._mode == "T-h":
                def f(rho):
                    return self._eq(rho, T)["h"]*1000-h

                def f2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*1000.*x+liquido["h"]*1000.*(1-x)-h

                rho = self.fsolve(f, True, f2, **{"h": h, "T": T})

            elif self._mode == "T-s":
                def f(rho):
                    return self._eq(rho, T)["s"]*1000.-s

                def f2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["s"]*1000.*x+liquido["s"]*1000.*(1-x)-s

                rho = self.fsolve(f, True, f2, **{"s": s, "T": T})

            elif self._mode == "T-u":
                def f(rho):
                    par = self._eq(rho, T)
                    return par["h"]-par["P"]/1000.*par["v"]-u/1000.

                def f2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vu*x+lu*(1-x)-u

                rho = self.fsolve(f, True, f2, **{"u": u, "T": T})

            elif self._mode == "P-rho":
                def f(T):
                    return self._eq(rho, T)["P"]-P

                def f2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    return Ps-P

                T = self.fsolve(f, True, f2, **{"P": P, "rho": rho})

            elif self._mode == "P-h":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["P"]-P, par["h"]*1000-h

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return Ps-P, vapor["h"]*1000*x+liquido["h"]*1000*(1-x)-h

                rho, T = self.fsolve(f, True, f2, **{"P": P, "h": h})

            elif self._mode == "P-s":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["P"]*1000-P, par["s"]*1000-s

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1./rhol)/(1./rhov-1./rhol)
                    return Ps-P, vapor["s"]*1000*x+liquido["s"]*1000*(1-x)-s

                rho, T = self.fsolve(f, True, f2, **{"P": P, "s": s})

            elif self._mode == "P-u":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]*par["v"]-u/1000., par["P"]-P

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return Ps-P, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(f, True, f2, **{"P": P, "u": u})

            elif self._mode == "rho-h":
                def f(T):
                    return self._eq(rho, T)["h"]*1000-h

                def f2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*x+liquido["h"]*(1-x)-h/1000

                T = self.fsolve(f, True, f2, **{"h": h, "rho": rho})

            elif self._mode == "rho-s":
                def f(T):
                    return self._eq(rho, T)["s"]*1000.-s

                def f2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["s"]*1000.*x+liquido["s"]*1000.*(1-x)-s

                T = self.fsolve(f, True, f2, **{"s": s, "rho": rho})

            elif self._mode == "rho-u":
                def f(T):
                    par = self._eq(rho, T)
                    return par["h"]-par["P"]/1000*par["v"]-u/1000

                def f2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vu*x+lu*(1-x)-u/1000

                T = self.fsolve(f, True, f2, **{"u": u, "rho": rho})

            elif self._mode == "h-s":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]*1000-h, par["s"]*1000.-s

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return (vapor["h"]*1000*x+liquido["h"]*1000*(1-x)-h,
                            vapor["s"]*1000*x+liquido["s"]*1000*(1-x)-s)

                rho, T = self.fsolve(f, True, f2, **{"s": s, "h": h})

            elif self._mode == "h-u":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]/1000*par["v"]-u, par["h"]-h

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*x+liquido["h"]*(1-x)-h, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(f, True, f2, **{"u": u, "h": h})

            elif self._mode == "s-u":
                def f(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]*par["v"]-u, par["s"]*1000.-s

                def f2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    vs = vapor["s"]*1000.
                    ls = liquido["s"]*1000.
                    return vs*x+ls*(1-x)-s, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(f, True, f2, **{"u": u, "s": s})

            if self._mode == "T-rho" and self.kwargs["rho"] == 0:
                self.status = 3
                self.msg = QApplication.translate(
                    "pychemqt", "Ideal condition at zero pressure")
            elif self._constants["Tmin"] <= T <= self._constants["Tmax"] and \
                    0 < rho:  # <= self._constants["rhomax"]*self.M:
                self.status = 1
                self.msg = ""
            else:
                self.status = 5
                self.msg = QApplication.translate(
                    "pychemqt", "input out of range")
                return

            rho = float(rho)
            T = float(T)
            propiedades = self._eq(rho, T)
            if self._mode != "T-P" and T <= self.Tc and rho:
                rhol = self._Liquid_Density(T)
                rhov = self._Vapor_Density(T)
                if rhol > rho > rhov:
                    rhol, rhov, Ps = self._saturation(T)
                    x = (1/rho-1/rhol)/(1/rhov-1/rhol)
                    if x < 0:
                        x = 0
                    elif x > 1:
                        x = 1
                    if rhol > rho > rhov:
                        P = Ps
                elif rho <= rhov:
                    x = 1
                elif rho >= rhol:
                    x = 0

                vapor = self._eq(rhov, T)
                liquido = self._eq(rhol, T)

            elif self._mode == "T-P":
                if T > self.Tc:
                    x = 1
                else:
                    rhol, rhov, Ps = self._saturation(T)
                    if Ps > P:
                        x = 1
                    else:
                        x = 0

            elif T > self.Tc:
                x = 1
            elif rho == 0:
                x = 1
            else:
                raise NotImplementedError("Incoming out of bound")

            if self.kwargs["P"]:
                P = self.kwargs["P"]/1000
            elif not P:
                P = propiedades["P"]/1000.
            else:
                P = P/1000.

        elif self._mode == "T-x":
            # Check input T in saturation range
            if self.Tt > T or self.Tc < T:
                raise ValueError("Wrong input values")

            rhol, rhov, Ps = self._saturation(T)
            rho = 1/(1/rhov*x+1/rhol*(1-x))
            vapor = self._eq(rhov, T)
            liquido = self._eq(rhol, T)
            if x == 0:
                propiedades = liquido
            elif x == 1:
                propiedades = vapor
            P = Ps/1000.
            self.status = 1

        elif self._mode == "P-x":
            # Iterate over saturation routine to get T
            def funcion(T):
                T = float(T)
                rhol, rhov, Ps = self._saturation(T)
                return Ps-P
            T = fsolve(funcion, 0.99*self.Tc)[0]
            rhol, rhov, Ps = self._saturation(T)
            rho = 1/(1/rhov*x+1/rhol*(1-x))
            vapor = self._eq(rhov, T)
            liquido = self._eq(rhol, T)
            if x == 0:
                propiedades = liquido
            elif x == 1:
                propiedades = vapor
            P = P/1000.
            self.status = 1

        self.T = unidades.Temperature(T)
        self.Tr = unidades.Dimensionless(T/self.Tc)
        self.P = unidades.Pressure(P, "kPa")
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.x = unidades.Dimensionless(x)

        # Ideal properties
        cp0 = self._prop0(rho, self.T)
        self.v0 = unidades.SpecificVolume(cp0["v"])
        self.rho0 = unidades.Density(1./self.v0)
        self.rhoM0 = unidades.MolarDensity(self.rho0/self.M)
        self.h0 = unidades.Enthalpy(cp0["h"])
        self.hM0 = unidades.MolarEnthalpy(self.h0/self.M)
        self.u0 = unidades.Enthalpy(self.h0-self.P*self.v0)
        self.uM0 = unidades.MolarEnthalpy(self.u0/self.M)
        self.s0 = unidades.SpecificHeat(cp0["s"])
        self.sM0 = unidades.MolarSpecificHeat(self.s0/self.M)
        self.a0 = unidades.Enthalpy(self.u0-self.T*self.s0)
        self.aM0 = unidades.MolarEnthalpy(self.a0/self.M)
        self.g0 = unidades.Enthalpy(self.h0-self.T*self.s0)
        self.gM0 = unidades.MolarEnthalpy(self.g0/self.M)
        self.cp0 = unidades.SpecificHeat(cp0["cp"])
        self.cpM0 = unidades.MolarSpecificHeat(self.cp0/self.M)
        self.cv0 = unidades.SpecificHeat(cp0["cv"])
        self.cvM0 = unidades.MolarSpecificHeat(self.cv0/self.M)
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv0)
        self.gamma0 = self.cp0_cv

        self.Liquido = ThermoAdvanced()
        self.Gas = ThermoAdvanced()
        if x == 0:
            # liquid phase
            self.fill(self.Liquido, propiedades)
            self.fill(self, propiedades)
            self.fillNone(self.Gas)
        elif x == 1:
            # vapor phase
            self.fill(self.Gas, propiedades)
            self.fill(self, propiedades)
            self.fillNone(self.Liquido)
        else:
            self.fillNone(self)
            self.fill(self.Liquido, liquido)
            self.fill(self.Gas, vapor)

            self.v = unidades.SpecificVolume(x*self.Gas.v+(1-x)*self.Liquido.v)
            self.rho = unidades.Density(1./self.v)

            self.h = unidades.Enthalpy(x*self.Gas.h+(1-x)*self.Liquido.h)
            self.s = unidades.SpecificHeat(x*self.Gas.s+(1-x)*self.Liquido.s)
            self.u = unidades.Enthalpy(x*self.Gas.u+(1-x)*self.Liquido.u)
            self.a = unidades.Enthalpy(x*self.Gas.a+(1-x)*self.Liquido.a)
            self.g = unidades.Enthalpy(x*self.Gas.g+(1-x)*self.Liquido.g)

            self.rhoM = unidades.MolarDensity(self.rho/self.M)
            self.hM = unidades.MolarEnthalpy(self.h*self.M)
            self.sM = unidades.MolarSpecificHeat(self.s*self.M)
            self.uM = unidades.MolarEnthalpy(self.u*self.M)
            self.aM = unidades.MolarEnthalpy(self.a*self.M)
            self.gM = unidades.MolarEnthalpy(self.g*self.M)

        # Calculate special properties useful only for one phase
        if x < 1 and self.Tt <= T <= self.Tc:
            self.sigma = unidades.Tension(self._Surface(T))
        else:
            self.sigma = unidades.Tension(None)

        if 0 < x < 1:
            self.virialB = unidades.SpecificVolume(vapor["B"]/self.rhoc)
            self.virialC = unidades.SpecificVolume_square(
                vapor["C"]/self.rhoc**2)
        else:
            self.virialB = unidades.SpecificVolume(propiedades["B"]/self.rhoc)
            self.virialC = unidades.SpecificVolume_square(
                propiedades["C"]/self.rhoc**2)

        if 0 < self.x < 1:
            self.Hvap = unidades.Enthalpy(vapor["h"]-liquido["h"], "kJkg")
            self.Svap = unidades.SpecificHeat(vapor["s"]-liquido["s"], "kJkgK")
        else:
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)
        self.invT = unidades.InvTemperature(-1/self.T)

    def fsolve(self, function, phases=True, function2phase=None, **kwargs):
        """Iterate to calculate T and rho
        function: function to iterate
        phases: calculate two phases region
        funtion2phase: function to iterate in two phase region"""
        if "T" not in kwargs:
            to = [self.Tc, self._constants["Tmin"], self._constants["Tmax"]]
            if self.kwargs["T0"]:
                to.insert(0, self.kwargs["T0"])
        if "rho" not in kwargs:
            rhov = self._Vapor_Density(self.Tt)
            ro = [rhov, self.rhoc, self._constants["rhomax"]*self.M, 1, 1e-3]
            if kwargs["rho0"]:
                if isinstance(kwargs["rho0"], list):
                    for rho in kwargs["rho0"][-1::-1]:
                        ro.insert(0, rho)
                else:
                    ro.insert(0, kwargs["rho0"])

        rinput = None
        rho, T = 0, 0
        if "T" in kwargs:
            T = kwargs["T"]
            if T >= self.Tc:
                phases = False
            for r in ro:
                try:
                    rinput = fsolve(function, r, full_output=True)
                    rho = rinput[0][0]
                except:
                    pass
                else:
                    f1 = rinput[1]["fvec"]
                    idx = rinput[2]
                    if 0 < rho and abs(f1) < 1e-3 and idx == 1:
                        break
        elif "rho" in kwargs:
            rho = kwargs["rho"]
            for t in to:
                try:
                    rinput = fsolve(function, t, full_output=True)
                    T = rinput[0][0]
                except:
                    pass
                else:
                    f1 = function(T)
                    if T != t and abs(f1) < 1e-6:
                        break
        else:
            for r, t in product(ro, to):
                try:
                    rinput = fsolve(function, [r, t], full_output=True)
                    rho, T = rinput[0]
                except:
                    pass
                else:
                    f1, f2 = function([rho, T])
                    if (rho != r or T != t) and \
                            0 < rho < self._constants["rhomax"]*self.M and \
                            abs(f1) < 1e-6 and abs(f2) < 1e-6:
                        break

        if phases:
            if self.Tt <= T < self.Tc:
                rhol = self._Liquid_Density(T)
                rhov = self._Vapor_Density(T)
                if rinput is None or rinput[2] != 1 or rhov <= rho <= rhol:
                    rhol = self._Liquid_Density(T)
                    rhov = self._Vapor_Density(T)
                    rhom = 1/(0.5/rhol+0.5/rhov)
                    t = (self.Tc+self.Tc)/2
                    rhol2 = self._Liquid_Density(t)
                    rhov2 = self._Vapor_Density(t)
                    rhom2 = 1/(0.5/rhol2+0.5/rhov2)
                    ro = [rhov, rhol, rhom, rhov2, rhol2, rhom2]
                    to = [T, T, T, t, t, t]

                    if "T" in kwargs:
                        for r in ro:
                            try:
                                rinput = fsolve(
                                    function2phase, r, full_output=True)
                                rho = rinput[0][0]
                            except:
                                pass
                            else:
                                f1 = function2phase(rho)
                                if rho != r and abs(f1) < 1e-3:
                                    break
                    elif "rho" in kwargs:
                        for t in to:
                            try:
                                rinput = fsolve(
                                    function2phase, t, full_output=True)
                                T = rinput[0][0]
                            except:
                                pass
                            else:
                                f1 = function2phase(T)
                                if 0 < rho < self._constants["rhomax"]*self.M \
                                        and T != t and abs(f1) < 1e-3:
                                    break
                    else:
                        for r, t in zip(ro, to):
                            try:
                                rinput = fsolve(
                                    function2phase, [r, t], full_output=True)
                                rho, T = rinput[0]
                            except:
                                pass
                            else:
                                f1, f2 = function2phase([rho, T])
                                if 0 < rho < self._constants["rhomax"]*self.M \
                                        and (rho != r or T != t) and \
                                        abs(f1) < 1e-3 and abs(f2) < 1e-3:
                                    break

        if "T" in kwargs:
            return rho
        elif "rho" in kwargs:
            return T
        else:
            return rho, T

    def fill(self, fase, estado):
        """Fill phase properties"""
        fase._bool = True
        fase.M = unidades.Dimensionless(self.M)
        fase.v = unidades.SpecificVolume(estado["v"])
        fase.rho = unidades.Density(1/fase.v)
        fase.Z = unidades.Dimensionless(self.P*fase.v/self.T/self.R)

        fase.h = unidades.Enthalpy(estado["h"], "kJkg")
        fase.s = unidades.SpecificHeat(estado["s"], "kJkgK")
        fase.u = unidades.Enthalpy(fase.h-self.P*fase.v)
        fase.a = unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g = unidades.Enthalpy(fase.h-self.T*fase.s)
        fase.fi = [unidades.Dimensionless(estado["fugacity"])]
        fase.f = [unidades.Pressure(f*self.P) for f in fase.fi]

        fase.cp = unidades.SpecificHeat(estado["cp"], "kJkgK")
        fase.cv = unidades.SpecificHeat(estado["cv"], "kJkgK")
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
#        fase.cps = estado["cps"]
        fase.w = unidades.Speed(estado["w"])

        fase.rhoM = unidades.MolarDensity(fase.rho/self.M)
        fase.hM = unidades.MolarEnthalpy(fase.h*self.M)
        fase.sM = unidades.MolarSpecificHeat(fase.s*self.M)
        fase.uM = unidades.MolarEnthalpy(fase.u*self.M)
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(fase.cv*self.M)
        fase.cpM = unidades.MolarSpecificHeat(fase.cp*self.M)

        fase.alfap = unidades.InvTemperature(estado["alfap"])
        fase.betap = unidades.Density(estado["betap"])

        # if fase.rho:
        #    d2Pdvdt = self.derivative("P", "v", "T", fase))

        if fase.rho:
            fase.gamma = unidades.Dimensionless(
                -fase.v/self.P*self.derivative("P", "v", "s", fase))
        else:
            fase.gamma = self.gamma0

        fase.joule = unidades.TemperaturePressure(
            self.derivative("T", "P", "h", fase))
        fase.Gruneisen = unidades.Dimensionless(
            fase.v/fase.cv*self.derivative("P", "T", "v", fase))

        if fase.rho:
            fase.alfav = unidades.InvTemperature(
                self.derivative("v", "T", "P", fase)/fase.v)
            fase.kappa = unidades.InvPressure(
                -self.derivative("v", "P", "T", fase)/fase.v)
            fase.kappas = unidades.InvPressure(
                -1/fase.v*self.derivative("v", "P", "s", fase))
            fase.betas = unidades.TemperaturePressure(
                self.derivative("T", "P", "s", fase))
            fase.kt = unidades.Dimensionless(
                -fase.v/self.P*self.derivative("P", "v", "T", fase))
            fase.ks = unidades.Dimensionless(
                -fase.v/self.P*self.derivative("P", "v", "s", fase))
            fase.Ks = unidades.Pressure(
                -fase.v*self.derivative("P", "v", "s", fase))
            fase.Kt = unidades.Pressure(
                -fase.v*self.derivative("P", "v", "T", fase))
            fase.dhdT_rho = unidades.SpecificHeat(
                self.derivative("h", "T", "rho", fase))
            fase.dhdT_P = unidades.SpecificHeat(
                self.derivative("h", "T", "P", fase))
            fase.dhdP_T = unidades.EnthalpyPressure(
                self.derivative("h", "P", "T", fase))  # deltat
            fase.deltat = fase.dhdP_T
            fase.dhdP_rho = unidades.EnthalpyPressure(
                self.derivative("h", "P", "rho", fase))
            fase.dhdrho_T = unidades.EnthalpyDensity(estado["dhdrho"])
            fase.dhdrho_P = unidades.EnthalpyDensity(
                estado["dhdrho"]+fase.dhdT_rho/estado["drhodt"])
            fase.dpdrho_T = unidades.PressureDensity(estado["dpdrho"])
            fase.drhodP_T = unidades.DensityPressure(1/estado["dpdrho"])
            fase.drhodT_P = unidades.DensityTemperature(estado["drhodt"])

            fase.Z_rho = unidades.SpecificVolume((fase.Z-1)/fase.rho)
            fase.hInput = unidades.Enthalpy(
                fase.v*self.derivative("h", "v", "P", fase))

        fase.dpdT_rho = unidades.PressureTemperature(
            self.derivative("P", "T", "rho", fase))
        fase.IntP = unidades.Pressure(self.derivative("u", "v", "T", fase))

        fase.virialB = unidades.SpecificVolume(estado["B"]/self.rhoc)
        fase.virialC = unidades.SpecificVolume_square(
            estado["C"]/self.rhoc**2)
        fase.invT = unidades.InvTemperature(-1/self.T)

        fase.mu = self._Viscosity(fase.rho, self.T, fase)
        fase.k = self._ThCond(fase.rho, self.T, fase)
        if fase.mu and fase.rho:
            fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        else:
            fase.nu = unidades.Diffusivity(None)
        if fase.k and fase.rho:
            fase.alfa = unidades.Diffusivity(fase.k/fase.rho/fase.cp)
        else:
            fase.alfa = unidades.Diffusivity(None)
        if fase.mu and fase.k:
            fase.Prandt = unidades.Dimensionless(fase.mu*fase.cp/fase.k)
        else:
            fase.Prandt = unidades.Dimensionless(None)
        fase.epsilon = unidades.Dimensionless(
            self._Dielectric(fase.rho, self.T))
        fase.fraccion = [1]
        fase.fraccion_masica = [1]

    def _saturation(self, T=None):
        """Saturation calculation for two phase search"""
        if not T:
            T = self.T
        if T > self.Tc:
            T = self.Tc
        T = float(T)
        rhoLo = self._Liquid_Density(T)
        rhoGo = self._Vapor_Density(T)

        def f(parr):
            rhol, rhog = parr
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc
            liquido = self._eq(rhol, T)
            vapor = self._eq(rhog, T)
            Jl = deltaL*(1+deltaL*liquido["fird"])
            Jv = deltaG*(1+deltaG*vapor["fird"])
            Kl = deltaL*liquido["fird"]+liquido["fir"]+log(deltaL)
            Kv = deltaG*vapor["fird"]+vapor["fir"]+log(deltaG)
            return Kv-Kl, Jv-Jl

        rhoL, rhoG = fsolve(f, [rhoLo, rhoGo])

        if rhoL == rhoG:
            Ps = self.Pc
        else:
            liquido = self._eq(rhoL, T)
            vapor = self._eq(rhoG, T)
            deltaL = rhoL/self.rhoc
            deltaG = rhoG/self.rhoc
            Ps = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
                liquido["fir"] - vapor["fir"] + log(deltaL/deltaG))
        return rhoL, rhoG, Ps

    def _Helmholtz(self, rho, T):
        """General implementation of Helmholtz free energy multiparameter
        equation of state"""
        delta = rho/self.rhoc
        tau = self.Tc/T

        ideal = self._phi0(self._constants["cp"], tau, delta)
        fio = ideal["fio"]
        fiot = ideal["fiot"]
        fiott = ideal["fiott"]
        fiodt = ideal["fiodt"]

        res = self._phir(tau, delta)
        fir = res["fir"]
        firt = res["firt"]
        firtt = res["firtt"]
        fird = res["fird"]
        firdd = res["firdd"]
        firdt = res["firdt"]
        B = res["B"]
        C = res["C"]

        propiedades = {}
        propiedades["fir"] = fir
        propiedades["fird"] = fird
        propiedades["firdd"] = firdd

        propiedades["T"] = T
        propiedades["P"] = (1+delta*fird)*self.R*T*rho
        if rho:
            propiedades["v"] = 1./rho
        else:
            propiedades["v"] = float("inf")

        propiedades["h"] = self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"] = self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cv"] = -self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["cp"] = self.R.kJkgK*(
            -tau**2*(fiott+firtt) +
            (1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["w"] = (self.R*T*(
            1 + 2*delta*fird+delta**2*firdd -
            (1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"] = (1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"] = rho*(
            1 + (delta*fird+delta**2*firdd)/(1 + delta*fird))
        propiedades["fugacity"] = exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"] = B
        propiedades["C"] = C
        propiedades["dpdrho"] = self.R*T*(1+2*delta*fird+delta**2*firdd)
        propiedades["drhodt"] = -rho*(1+delta*fird-delta*tau*firdt) / \
            (T*(1+2*delta*fird+delta**2*firdd))
        if rho:
            propiedades["dhdrho"] = self.R*T/rho * \
                (tau*delta*(fiodt+firdt)+delta*fird+delta**2*firdd)
        else:
            propiedades["dhdrho"] = 0
#        dbt=-phi11/rho/t
#        propiedades["cps"] = propiedades["cv"]-self.R*(1+delta*fird-delta*tau*firdt)*T/rho*propiedades["drhodt"]
#        propiedades["cps"] = self.R*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)/(1+2*delta*fird+delta**2*firdd)*
#                                    ((1+delta*fird-delta*tau*firdt)-self.rhoc/self.R/delta*self.derivative("P", "T", "rho", propiedades)))
#        propiedades["cps"] = propiedades["cv"] Add cps from Argon pag.27
        return propiedades


    def _ECS(self,  rho, T):
        delta = rho/self.rhoc
        tau = self.Tc/T

        ideal = self._phi0(self._constants["cp"], tau, delta)
        fio = ideal["fio"]
        fiot = ideal["fiot"]
        fiott = ideal["fiott"]

        ref = self._constants["ref"](eq=self._constants["eq"])
        Tr = T/self.Tc
        rhor = rho/self.rhoc

        psi = 1+(self.f_acent-ref.f_acent)*(self._constants["ft"][0]+self._constants["ft"][1]*log(Tr))
        for n, m in zip(self._constants["ft_add"], self._constants["ft_add_exp"]):
            psi += n*Tr**m
        for n, m in zip(self._constants["fd"], self._constants["fd_exp"]):
            psi += n*rhor**m
        T0 = T*ref.Tc/self.Tc/psi

        phi = ref.Zc/self.Zc*(1+(self.f_acent-ref.f_acent)*(self._constants["ht"][0]+self._constants["ht"][1]*log(Tr)))
        for n, m in zip(self._constants["ht_add"], self._constants["ht_add_exp"]):
            phi += n*Tr**m
        for n, m in zip(self._constants["hd"], self._constants["hd_exp"]):
            phi += n*rhor**m
        rho0 = rho*ref.rhoc/self.rhoc*phi

        deltaref = rho0/ref.rhoc
        tauref = ref.Tc/T0
        res = self._phir(tauref, deltaref)
        fir = res["fir"]
        firt = res["firt"]
        firtt = res["firtt"]
        fird = res["fird"]
        firdd = res["firdd"]
        firdt = res["firdt"]
        B = res["B"]
        C = res["C"]

        propiedades = {}
        propiedades["fir"]=fir
        propiedades["fird"]=fird
        propiedades["firdd"]=firdd

        propiedades["T"]=T
        propiedades["P"]=(1+delta*fird)*self.R.JkgK*T*rho
        propiedades["v"]=1/rho
        propiedades["h"]=self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"]=self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cv"]=-self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["cp"]=self.R.kJkgK*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["w"]=(self.R*T*(1+2*delta*fird+delta**2*firdd-(1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"]=(1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"]=rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
        propiedades["fugacity"]=exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"]=B
        propiedades["C"]=C
        propiedades["dpdrho"]=self.R*T*(1+2*delta*fird+delta**2*firdd)
        propiedades["dpdT"]=self.R.kJkgK*rho*(1+delta*fird+delta*tau*firdt)
#        propiedades["cps"]=propiedades["cv"] Add cps from Argon pag.27
        return propiedades

    def _MBWR(self, rho, T):
        """Multiparameter equation of state of Benedict-Webb-Rubin"""
        rho = rho/self.M
        delta = rho/self.rhoc
        tau = self.Tc/T
        b = self._constants["b"]

        ideal = self._phi0(self._constants["cp"], tau, delta)
        fio = ideal["fio"]
        fiot = ideal["fiot"]
        fiott = ideal["fiott"]

        a = [None]
        a.append(self.R*T)
        a.append(b[1]*T+b[2]*T**0.5+b[3]+b[4]/T+b[5]/T**2)
        a.append(b[6]*T+b[7]+b[8]/T+b[9]/T**2)
        a.append(b[10]*T+b[11]+b[12]/T)
        a.append(b[13])
        a.append(b[14]/T+b[15]/T**2)
        a.append(b[16]/T)
        a.append(b[17]/T+b[18]/T**2)
        a.append(b[19]/T**2)
        a.append(b[20]/T**2+b[21]/T**3)
        a.append(b[22]/T**2+b[23]/T**4)
        a.append(b[24]/T**2+b[25]/T**3)
        a.append(b[26]/T**2+b[27]/T**4)
        a.append(b[28]/T**2+b[29]/T**3)
        a.append(b[30]/T**2+b[31]/T**3+b[32]/T**4)

        P = sum([a[n]*rho**n for n in range(1, 10)])
        P += exp(-delta**2)*sum([a[n]*rho**(2*n-17) for n in range(10, 16)])
        P = P*10

        dPdrho = sum([a[n]*n*rho**(n-1) for n in range(1, 10)])
        dPdrho += exp(-delta**2)*sum([(2*n-17-2*delta**2)*a[n]*rho**(2*n-18) for n in range(10, 16)])
        dPdrho = dPdrho*100

        d2Prho = sum([a[n]*n*(n-1)*rho**(n-2) for n in range(1, 10)])
        d2Prho += exp(-delta**2)*sum([(-35*n+2*n**2+153+33*delta**2+2*delta**4-4*n*delta**2)*2*a[n]*rho**(2*n-19) for n in range(10, 16)])
        d2Prho = d2Prho*100

        A = 0
        for n in range(2, 10):
            A += a[n]/(n-1.)*rho**(n-1)

        A -= 0.5*a[10]*self.rhoc**2*(exp(-delta**2)-1)
        A -= 0.5*a[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        A -= 0.5*a[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        A -= 0.5*a[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        A -= 0.5*a[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        A -= 0.5*a[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        A = A*100

        adT = [None, self.R]
        adT.append(b[1]+b[2]/2/T**0.5-b[4]/T**2-2*b[5]/T**3)
        adT.append(b[6]-b[8]/T**2-2*b[9]/T**3)
        adT.append(b[10]-b[12]/T**2)
        adT.append(0)
        adT.append(-b[14]/T**2-2*b[15]/T**3)
        adT.append(-b[16]/T**2)
        adT.append(-b[17]/T**2-2*b[18]/T**3)
        adT.append(-2*b[19]/T**3)
        adT.append(-2*b[20]/T**3-3*b[21]/T**4)
        adT.append(-2*b[22]/T**3-4*b[23]/T**5)
        adT.append(-2*b[24]/T**3-3*b[25]/T**4)
        adT.append(-2*b[26]/T**3-4*b[27]/T**5)
        adT.append(-2*b[28]/T**3-3*b[29]/T**4)
        adT.append(-2*b[30]/T**3-3*b[31]/T**4-4*b[32]/T**5)

        dA = 0
        for n in range(2, 10):
            dA += adT[n]/(n-1.)*rho**(n-1)

        dA -= 0.5*adT[10]*self.rhoc**2*(exp(-delta**2)-1)
        dA -= 0.5*adT[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        dA -= 0.5*adT[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        dA -= 0.5*adT[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        dA -= 0.5*adT[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        dA -= 0.5*adT[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        dA = dA*100

        dPdT = sum([adT[n]*rho**n for n in range(1, 10)])
        dPdT += exp(-delta**2)*sum([adT[n]*rho**(2*n-17) for n in range(10, 16)])
        dPdT = dPdT*100

        adTT = [None, 0]
        adTT.append(-b[2]/4/T**1.5+2*b[4]/T**3+6*b[5]/T**4)
        adTT.append(2*b[8]/T**3+6*b[9]/T**4)
        adTT.append(2*b[12]/T**3)
        adTT.append(0)
        adTT.append(2*b[14]/T**3+6*b[15]/T**4)
        adTT.append(2*b[16]/T**3)
        adTT.append(2*b[17]/T**3+6*b[18]/T**4)
        adTT.append(6*b[19]/T**4)
        adTT.append(6*b[20]/T**4+12*b[21]/T**5)
        adTT.append(6*b[22]/T**4+20*b[23]/T**6)
        adTT.append(6*b[24]/T**4+12*b[25]/T**5)
        adTT.append(6*b[26]/T**4+20*b[27]/T**6)
        adTT.append(6*b[28]/T**4+12*b[29]/T**5)
        adTT.append(6*b[30]/T**4+12*b[31]/T**5+20*b[32]/T**6)

        d2A = 0
        for n in range(2, 10):
            d2A += adTT[n]*rho**(n-1)/(n-1)
        d2A -= 0.5*adTT[10]*self.rhoc**2*(exp(-delta**2)-1)
        d2A -= 0.5*adTT[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        d2A -= 0.5*adTT[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        d2A -= 0.5*adTT[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        d2A -= 0.5*adTT[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        d2A -= 0.5*adTT[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        d2A = d2A*100

        # TODO:
        B, C = 0, 0

        fir = A/self.R/T
        firdt = P/T-dPdT/self.R/rho
        firt = A/T-dA/self.R
        firtt = d2A*T/self.R
        fird = P/self.R/T/rho-1.
        firdd = (dPdrho-2*P/rho)/self.R/T+1.
#        firddd = (d2Prho*rho-4*dPdrho+6*P/rho)/self.R/T-2.

        propiedades = {}
        propiedades["T"] = T
        propiedades["P"] = P*0.1  # converted from bar to MPa
        propiedades["v"] = 1/rho
        propiedades["h"] = self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"] = self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cp"] = self.R.kJkgK*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["cv"] = -self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["w"] = (self.R*T*(1+2*delta*fird+delta**2*firdd-(1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"] = (1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"] = rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
        propiedades["fugacity"] = exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"] = B
        propiedades["C"] = C
        propiedades["fir"] = fir
        propiedades["fird"] = fird
        return propiedades

    def _PengRobinson(self, rho, T):
        """Peng, D.-Y.; Robinson, D.B. A New Two-Constant Equation of State. I&EC Fundam. 1976, 15(1), 59
        Peneloux, A.; Rauzy, E.; Freze, R. A consistent correction for Redlich-Kwong-Soave volumes. Fluid Phase Eq. 1982, 8, 7.
        http://dx.doi.org/10.1021/i160057a011
        http://dx.doi.org/10.1016/0378-3812(82)80002-2"""
        # FIXME: no sale por poco
        delta = rho/self.rhoc
        tau = self.Tc/T
        delta_0 = 1e-50
        rho = rho/self.M

        Tr = 1./tau
        m = 0.37464+1.54226*self.f_acent-0.26992*self.f_acent**2
        alfa = (1+m*(1-Tr**0.5))**2
        a = 0.457235*R_atml**2*self.Tc**2/self.Pc.atm
        b = 0.077796*R_atml*self.Tc/self.Pc.atm

        daT = -a*m/T**0.5/self.Tc**0.5*alfa**0.5
        d2aT = a*m*(1+m)/2/T/(T*self.Tc)**0.5

        ideal = self._phi0(self._constants["cp"], tau, delta)
        fio = ideal["fio"]
        fiot = ideal["fiot"]
        fiott = ideal["fiott"]

        v = 1./rho
        vb = v-b
        q = b*8**0.5
        v1 = 2*v+2*b+q
        v2 = 2*v+2*b-q
        v1n = v1+2*self._PR
        v2n = v2+2*self._PR
        fir = log(v)-log(vb+self._PR)+a/self.R.kJkgK/T*log(v2n/v1n)/q

        phipart = (1./self.R.kJkgK/q)*(daT/T-a/T**2)
        dtdtau = -T**2/self.Tc
        dphidtau = phipart*dtdtau
        term1 = 2+2*b*rho-rho*q+2*rho*self._PR
        term2 = 2+2*b*rho+rho*q+2*rho*self._PR
        phidpart = -4*self.rhoc/self.M*q/term1/term2
        firdt = dphidtau*phidpart

        fird = -1+1./(1-b*rho+rho*self._PR)+a/self.R.kJkgK/T/q*(2/v1n-2/v2n)/rho
        bdt = 1-b*rho+rho*self._PR
        firdd = 1-1./bdt-(self._PR*self.rhoc/self.M-b*self.rhoc/self.M) * \
            delta/bdt**2 + 4*a/self.R.kJkgK/T/q/rho * \
            ((1./v2n-1./v1n)+(1./v1n**2-1./v2n**2)/rho)
        dphidt = 1./self.R.kJkgK/q*log(v2n/v1n)*(daT/T-a/T**2)
        firt = dphidt*dtdtau
        d2phidt2 = (1./self.R.kJkgK/q)*log(v2n/v1n)*(d2aT/T-2*daT/T**2+2*a/T**3)
        d2phid2tau = 1./tau**2*(T**2*d2phidt2+2*T*dphidt)
        firtt = d2phid2tau

#        fir=-log(1-b*rho)+a/(self.R*T*2**1.5*b)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        fird=b/(1-b*rho)+a/(self.R*T*2**1.5*b)*((1-2**0.5)*b/(1+(1-2**0.5)*b*rho)-(1+2**0.5)*b/(1+(1+2**0.5)*b*rho))
#        firdd=(b/(1-b*rho))**2+a/(rho**2*self.R*T*2**1.5*b)*(-((1-2**0.5)*b*rho/(1+(1-2**0.5)*b*rho))**2+((1+2**0.5)*b*rho/(1+(1+2**0.5)*b*rho))**2)
#        firt=1/(self.R*T*2**1.5*b)*(daT-a/T)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        firtt=1/(self.R*2**1.5*b*T)*(d2aT-2/T*daT+2*a/T**2)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        firdt=1/(rho*T*self.R*2**1.5*b)*(daT-a/T)*((1-2**0.5)*b*rho/(1+(1-2**0.5)*b*rho)-(1+2**0.5)*b*rho/(1+(1+2**0.5)*b*rho))
        firdtt = B = C = 0

        propiedades = {}
        propiedades["fir"] = fir
        propiedades["fird"] = fird
        propiedades["firdd"] = firdd

        propiedades["T"] = T
        propiedades["P"] = (1+delta*fird)*self.R.JkgK*T*rho
        propiedades["v"] = 1/rho
        propiedades["h"] = self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"] = self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cp"] = self.R.kJkgK*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["cv"] = -self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["w"] = (self.R*T*(1+2*delta*fird+delta**2*firdd-(1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"] = (1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"] = rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
        propiedades["fugacity"] = exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"] = B
        propiedades["C"] = C
        propiedades["dpdrho"] = self.R*T*(1+2*delta*fird+delta**2*firdd)
        propiedades["dpdT"] = self.R*rho*(1+delta*fird+delta*tau*firdt)
        return propiedades

    def _Generalised(self):
        """Span R., Wagner W., "An accurate empirical three parameter equation of state for nonpolar fluids". To be submitted to Fluid Phase Equilibria. 2000 """
        if self._Tr:
            Tref = self._Tr
        else:
            Tref = self.Tc
        if self._rhor:
            rhoref = self._rhor
        else:
            rhoref = self.rhoc
        if self._w:
            w = self._w
        else:
            w = self.f_acent

        helmholtz = {
            "R": 8.31451,
            "Tc": Tref,
            "rhoc": rhoref,
            "cp": self.eq[0]["cp"],

            "d1": [1, 1, 2, 3, 8],
            "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

            "d2": [2, 3, 1, 4, 3],
            "t2": [0.625, 2, 4.125, 4.125, 17],
            "c2": [1, 1, 2, 2, 3],
            "gamma2": [1]*5}

        c1 = [0.636479524, -0.174667493e1, -0.144442644e-1, 0.6799731e-1,
              0.767320032e-4, 0.218194143, 0.810318494e-1, -0.907368899e-1,
              0.25312225e-1, -0.209937023e-1]
        c2 = [0.82247342, -0.954932692, -0.745462328, 0.182685593,
              0.547120142e-4, 0.761697913, 0.415691324, -0.825206373,
              -0.240558288, -0.643818403e-1]
        c3 = [-0.186193063e1, 0.105083555e2, 0.16403233e1, -0.613747797,
              -0.69318829e-3, -0.705727791e1, -0.290006245e1, -0.232497527,
              -0.282346515, 0.254250643e1]
        nr = [c1[i]+c2[i]*w+c3[i]*w**4 for i in range(10)]
        helmholtz["nr1"] = nr[:5]
        helmholtz["nr2"] = nr[5:]
        self._constants = helmholtz

    def _ref(self, ref, refvalues=None):
        """Define reference state
        Input:
            ref: name of standard
                OTO | NBP | IIR | ASHRAE | CUSTOM
            refvalues: array with custom refvalues
                Tref
                Pref in kPa
                ho in J/mol
                so in J/mol·K
        """
        if ref is None:
            rf = self._constants["ref"]
            if isinstance(rf, str):
                ref = rf
            elif isinstance(rf, dict):
                ref = "CUSTOM"
                refvalues = [rf["Tref"], rf["Pref"], rf["ho"], rf["so"]]
            else:
                ref = "OTO"

        if ref == "OTO":
            self.Tref = 298.15
            self.Pref = 101325.
            self.ho = 0
            self.so = 0
        elif ref == "NBP":
            self.Tref = self.Tb
            self.Pref = 101325.
            self.ho = 0
            self.so = 0
        elif ref == "IIR":
            self.Tref = 273.15
            self.Pref = self._Vapor_Pressure(273.15)
            self.ho = 200
            self.so = 1
        elif ref == "ASHRAE":
            self.Tref = 233.15
            self.Pref = self._Vapor_Pressure(233.15)
            self.ho = 0
            self.so = 0
        elif ref == "CUSTOM":
            if refvalues is None:
                refvalues = [298.15, 101325., 0., 0.]
            self.Tref = unidades.Temperature(refvalues[0])
            self.Pref = unidades.Pressure(refvalues[1], "kPa")
            self.ho = unidades.Enthalpy(refvalues[2]/self.M, "Jg")
            self.so = unidades.SpecificHeat(refvalues[3]/self.M, "JgK")

    def _prop0(self, rho, T):
        """Ideal gas properties"""
        delta = rho/self.rhoc
        tau = self.Tc/T
        ideal = self._phi0(self._constants["cp"], tau, delta)
        fio = ideal["fio"]
        fiot = ideal["fiot"]
        fiott = ideal["fiott"]

        propiedades = {}
        if rho:
            propiedades["v"] = self.R*T/self.P
        else:
            propiedades["v"] = float("inf")
        propiedades["h"] = self.R*T*(1+tau*fiot)
        propiedades["s"] = self.R*(tau*fiot-fio)
        propiedades["cv"] = -self.R*tau**2*fiott
        propiedades["cp"] = self.R*(-tau**2*fiott+1)
        return propiedades

    def _PHIO(self, cp):
        """Convert cp dict in phi0 dict when the cp expression isn't in
        Helmholtz free energy terms"""
        # TODO: Add support for 1/τ terms, give tau*logtau terms and more
        R = cp.get("R", self._constants["R"])/self.M*1000
        rho0 = self.Pref/R/self.Tref
        tau0 = self.Tc/self.Tref
        delta0 = rho0/self.rhoc
        co = cp["ao"]-1
        ti = [-x for x in cp["pow"]]
        ci = [-n/(t*(t+1))*self.Tc**t for n, t in zip(cp["an"], cp["pow"])]
        titao = [fi/self.Tc for fi in cp["exp"]]
        hyp = [fi/self.Tc for fi in cp["hyp"]]
        cI = self.ho/R/self.Tc-cp["ao"]/tau0
        cII = -1+log(tau0/delta0)+cp["ao"]-cp["ao"]*log(tau0)-self.so/R
        for c, t in zip(ci, ti):
            cI -= c*t*tau0**(t-1)
            cII += c*(t-1)*tau0**t
        for ao, tita in zip(cp["ao_exp"], titao):
            cI -= ao*tita/(exp(tita*tau0)-1)
            cII -= ao*(tita*tau0/(1-exp(tita*tau0))+log(1-exp(-tita*tau0)))
        if cp["ao_hyp"]:
            for i in [0, 2]:
                cI -= cp["ao_hyp"][i]*hyp[i]/(tanh(hyp[i]*tau0))
                cII += cp["ao_hyp"][i]*(hyp[i]*tau0/tanh(hyp[i]*tau0) -
                                        log(abs(sinh(hyp[i]*tau0))))
            for i in [1, 3]:
                cI += cp["ao_hyp"][i]*hyp[i]*tanh(hyp[i]*tau0)
                cII -= cp["ao_hyp"][i]*(hyp[i]*tau0*tanh(hyp[i]*tau0) -
                                        log(abs(cosh(hyp[i]*tau0))))

        Fi0 = {"ao_log": [1,  co],
               "pow": [0, 1] + ti,
               "ao_pow": [cII, cI] + ci,
               "ao_exp": cp["ao_exp"],
               "titao": titao,
               "ao_hyp": cp["ao_hyp"],
               "hyp": hyp}
        return Fi0

    def _phi0(self, cp, tau, delta):
        """Ideal gas Helmholtz free energy and derivatives

        Parameters
        ----------
        tau : float
            Inverse reduced temperature, Tc/T [-]
        delta : float
            Reduced density, rho/rhoc [-]

        Returns
        -------
        prop : dictionary with ideal adimensional helmholtz energy and deriv
            fio  [-]
            fiot: [∂fio/∂τ]δ  [-]
            fiod: [∂fio/∂δ]τ  [-]
            fiott: [∂²fio/∂τ²]δ  [-]
            fiodt: [∂²fio/∂τ∂δ]  [-]
            fiodd: [∂²fio/∂δ²]τ  [-]
        """

        if "ao_log" in cp:
            Fi0 = cp
        else:
            Fi0 = self._PHIO(cp)

#        T = self._constants.get("Tref", self.Tc)/tau
#        rho = delta*self.rhoc
#        rho0 = self.Pref/self.R/self.Tref
#        delta0=rho0/self.rhoc
#        tau0=self._constants.get("Tref", self.Tc)/self.Tref
#
#        cp0sav = self._Cp0(cp, T)
#        cpisav = self._dCp(cp, T, self.Tref)
#        cptsav = self._dCpT(cp, T, self.Tref)
#        cpt2sav = self._dCpT2(cp, T, self.Tref)
#
#        fio = self.ho*tau/self.R/self.Tc-self.so/self.R-1+\
#            log(delta*tau0/delta0/tau)-tau*cpt2sav/self.R*1000+cptsav/self.R
#        fiot = (cpisav/self.R*1000/T-1)/tau+self.ho/self.R/self.Tc
#        fiott = (1-cp0sav/self.R*1000)/tau**2

        # FIXME: Reference estate
        fio = Fi0["ao_log"][1]*log(tau)
        fiot = Fi0["ao_log"][1]/tau
        fiott = -Fi0["ao_log"][1]/tau**2

        if delta:
            fiod = 1/delta
            fiodd = -1/delta**2
        else:
            fiod, fiodd = 0, 0
        fiodt = 0

        R_ = cp.get("R", self._constants["R"])
        for n, t in zip(Fi0["ao_pow"], Fi0["pow"]):
            fio += n*tau**t
            if t != 0:
                fiot += t*n*tau**(t-1)
            if t not in [0, 1]:
                fiott += n*t*(t-1)*tau**(t-2)

        for n, g in zip(Fi0["ao_exp"], Fi0["titao"]):
            fio += n*log(1-exp(-g*tau))
            fiot += n*g*((1-exp(-g*tau))**-1-1)
            fiott -= n*g**2*exp(-g*tau)*(1-exp(-g*tau))**-2

        # Special case for τ·ln(τ) terms (i.e. undecane, D2O)
        if "tau*logtau" in Fi0:
            fio += Fi0["tau*logtau"]*tau*log(tau)
            fiot += Fi0["tau*logtau"]*(log(tau)+1)
            fiott += Fi0["tau*logtau"]/tau

        if "tau*logdelta" in Fi0 and delta:
            fio += Fi0["tau*logdelta"]*tau*log(delta)
            fiot += Fi0["tau*logdelta"]*log(delta)
            fiod += Fi0["tau*logdelta"]*tau/delta
            fiodd -= Fi0["tau*logdelta"]*tau/delta**2
            fiodt += Fi0["tau*logdelta"]/delta

        # Special case for Lemmon-Jacobsen mEoS for air
        if "ao_exp2" in Fi0:
            for n, g, C in zip(Fi0["ao_exp2"], Fi0["titao2"], Fi0["sum2"]):
                fio += n*log(C+exp(g*tau))
                fiot += n*g/(C*exp(-g*tau)+1)
                fiott += C*n*g**2*exp(-g*tau)/(C*exp(-g*tau)+1)**2

        if "ao_hyp" in Fi0 and Fi0["ao_hyp"]:
            for i in [0, 2]:
                fio += Fi0["ao_hyp"][i]*log(abs(sinh(Fi0["hyp"][i]*tau)))
                fiot += Fi0["ao_hyp"][i]*Fi0["hyp"][i]/tanh(Fi0["hyp"][i]*tau)
                fiott -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/sinh(
                    Fi0["hyp"][i]*tau)**2

            for i in [1, 3]:
                fio -= Fi0["ao_hyp"][i]*log(abs(cosh(Fi0["hyp"][i]*tau)))
                fiot -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]*tanh(Fi0["hyp"][i]*tau)
                fiott -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/cosh(
                    Fi0["hyp"][i]*tau)**2

        R_ = cp.get("R", self._constants["R"])
        factor = R_/self._constants["R"]
        if delta:
            fio = Fi0["ao_log"][0]*log(delta)+factor*fio
        else:
            fio *= factor

        prop = {}
        prop["fio"] = fio
        prop["fiot"] = fiot*factor
        prop["fiott"] = fiott*factor
        prop["fiod"] = fiod
        prop["fiodd"] = fiodd
        prop["fiodt"] = fiodt
        return prop

    def _Cp0(self, T=False):
        Tc = self._constants.get("Tref", self.Tc)
        if not T:
            T = self.T
        cp = self._constants["cp"]

        if "ao_log" in cp:
            tau = Tc/T
            ideal = self._phi0(self._constants["cp"], tau, 0)
            fiott = ideal["fiott"]
            cpo = (-tau**2*fiott+1)*self.R
            return unidades.SpecificHeat(cpo)
        else:
            tau = Tc/T
            cpo = cp["ao"]
            for a, t in zip(cp["an"], cp["pow"]):
                cpo += a*T**t
            for m, tita in zip(cp["ao_exp"], cp["exp"]):
                cpo += m*(tita/T)**2*exp(tita/T)/(1-exp(tita/T))**2
            if cp["ao_hyp"]:
                for i in [0, 2]:
                    cpo += cp["ao_hyp"][i]*(cp["hyp"][i]/T/(sinh(cp["hyp"][i]/T)))**2
                for i in [1, 3]:
                    cpo += cp["ao_hyp"][i]*(cp["hyp"][i]/T/(cosh(cp["hyp"][i]/T)))**2
            return unidades.SpecificHeat(cpo*self.R/1000)

    def _dCp(self, cp, T, Tref):
        """Calcula la integral de Cp0 entre T y Tref, necesario para calcular la entalpia usando estados de referencia"""
        cpsum = cp.get("ao", 0)*(T-Tref)
        for n, t in zip(cp["an"], cp["pow"]):
            t += 1
            if t == 0:
                cpsum += n*log(T/Tref)
            else:
                cpsum += n*(T**(t+1)-Tref**(t+1))/((t+1)*self.Tc**t)
#        for tita, ao in zip(cp["ao_exp"], cp["exp"]):
#            cpsum += -ao*0.5*tita*((1.+exp(tita/T))/(1.-exp(tita/T))-(1.+exp(tita/Tref))/(1.-exp(tita/Tref)))
        return unidades.Enthalpy(cpsum/self.M/1000)

    def _dCpT(self, cp, T, Tref):
        """Calcula la integral de Cp0/T entre T y Tref, necesario para calcular la entropia usando estados de referencia"""
        cpsum = cp.get("ao", 0)*log(T/Tref)
#        for n, t in zip(cp["an"], cp["pow"]):
#            cpsum += n*(T**t-Tref**t)/(t*self.Tc**t)
#        for tita, ao in zip(cp["ao_exp"], cp["exp"]):
#            cpsum += ao*(log((1.-exp(tita/Tref))/(1.-exp(tita/T)))+tita/T*exp(tita/T)/(exp(tita/T)-1.)-tita/Tref*exp(tita/Tref)/(exp(tita/Tref)-1.))
        return unidades.SpecificHeat(cpsum*self.R*1000)

    def _dCpT2(self, cp, T, Tref):
        """Calcula la integral de Cp0/T² entre T y Tref, necesario para calcular la entropia usando estados de referencia"""
        cpsum = -cp.get("ao", 0)*(1/T-1/Tref)
        return unidades.SpecificHeat(cpsum*self.M*1000)

    def _phir(self, tau, delta):
        """Residual contribution to the free Helmholtz energy

        Parameters
        ----------
        tau : float
            Inverse reduced temperature, Tc/T [-]
        delta : float
            Reduced density, rho/rhoc [-]

        Returns
        -------
        prop : dictionary with residual adimensional helmholtz energy and deriv
            fir  [-]
            firt: [∂fir/∂τ]δ,x  [-]
            fird: [∂fir/∂δ]τ,x  [-]
            firtt: [∂²fir/∂τ²]δ,x  [-]
            firdt: [∂²fir/∂τ∂δ]x  [-]
            firdd: [∂²fir/∂δ²]τ,x  [-]
        """

        delta_0 = 1e-100
        fir = fird = firdd = firt = firtt = firdt = firdtt = B = C = 0

        if delta:
            # Polinomial terms
            nr1 = self._constants.get("nr1", [])
            d1 = self._constants.get("d1", [])
            t1 = self._constants.get("t1", [])
            for n, d, t in zip(nr1, d1, t1):
                fir += n*delta**d*tau**t
                fird += n*d*delta**(d-1)*tau**t
                firdd += n*d*(d-1)*delta**(d-2)*tau**t
                firt += n*t*delta**d*tau**(t-1)
                firtt += n*t*(t-1)*delta**d*tau**(t-2)
                firdt += n*t*d*delta**(d-1)*tau**(t-1)
                firdtt += n*t*d*(t-1)*delta**(d-1)*tau**(t-2)
                B += n*d*delta_0**(d-1)*tau**t
                C += n*d*(d-1)*delta_0**(d-2)*tau**t

            # Exponential terms
            nr2 = self._constants.get("nr2", [])
            d2 = self._constants.get("d2", [])
            g2 = self._constants.get("gamma2", [])
            t2 = self._constants.get("t2", [])
            c2 = self._constants.get("c2", [])
            for n, d, g, t, c in zip(nr2, d2, g2, t2, c2):
                fir += n*delta**d*tau**t*exp(-g*delta**c)
                fird += n*exp(-g*delta**c)*delta**(d-1)*tau**t*(d-g*c*delta**c)
                firdd += n*exp(-g*delta**c)*delta**(d-2)*tau**t * \
                    ((d-g*c*delta**c)*(d-1-g*c*delta**c)-g**2*c**2*delta**c)
                firt += n*t*delta**d*tau**(t-1)*exp(-g*delta**c)
                firtt += n*t*(t-1)*delta**d*tau**(t-2)*exp(-g*delta**c)
                firdt += n*t*delta**(d-1)*tau**(t-1)*(d-g*c*delta**c) * \
                    exp(-g*delta**c)
                firdtt += n*t*(t-1)*delta**(d-1)*tau**(t-2) * \
                    (d-g*c*delta**c) * exp(-g*delta**c)
                B += n*exp(-g*delta_0**c)*delta_0**(d-1)*tau**t * \
                    (d-g*c*delta_0**c)
                C += n*exp(-g*delta_0**c) * \
                    (delta_0**(d-2)*tau**t*((d-g*c*delta_0**c)*(
                        d-1-g*c*delta_0**c)-g**2*c**2 * delta_0**c))

            # Gaussian terms
            nr3 = self._constants.get("nr3", [])
            d3 = self._constants.get("d3", [])
            t3 = self._constants.get("t3", [])
            a3 = self._constants.get("alfa3", [])
            e3 = self._constants.get("epsilon3", [])
            b3 = self._constants.get("beta3", [])
            g3 = self._constants.get("gamma3", [])
            exp1 = self._constants.get("exp1", [2]*len(nr3))
            exp2 = self._constants.get("exp2", [2]*len(nr3))
            for n, d, t, a, e, b, g, ex1, ex2 in zip(
                    nr3, d3, t3, a3, e3, b3, g3, exp1, exp2):
                expr = exp(-a*(delta-e)**ex1-b*(tau-g)**ex2)

                fir += n*delta**d*tau**t * expr
                fird += expr * (n*d*delta**(d-1)*tau**t -
                                n*a*delta**d*(delta-e)**(ex1-1)*ex1*tau**t)
                firdd += expr * (
                    n*a**2*delta**d*(delta-e)**(2*ex1-2)*ex1**2*tau**t -
                    n*a*delta**d*(delta-e)**(ex1-2)*(ex1-1)*ex1*tau**t -
                    2*n*a*d*delta**(d-1)*(delta-e)**(ex1-1)*ex1*tau**t +
                    n*(d-1)*d*delta**(d-2)*tau**t)
                firt += expr * (n*delta**d*t*tau**(t-1) -
                                n*b*delta**d*ex2*tau**t*(tau-g)**(ex2-1))
                firtt += expr * (
                    n*b**2*delta**d*ex2**2*tau**t*(tau-g)**(2*ex2-2) -
                    2*n*b*delta**d*ex2*t*tau**(t-1)*(tau-g)**(ex2-1) -
                    n*b*delta**d*(ex2-1)*ex2*tau**t*(tau-g)**(ex2-2) +
                    n*delta**d*(t-1)*t*tau**(t-2))
                firdt += expr * (
                    n*a*b*delta**d*(delta-e)**(ex1-1)*ex1*ex2*tau**t *
                    (tau-g)**(ex2-1) -
                    n*b*d*delta**(d-1)*ex2*tau**t*(tau-g)**(ex2-1) -
                    n*a*delta**d*(delta-e)**(ex1-1)*ex1*t*tau**(t-1) +
                    n*d*delta**(d-1)*t*tau**(t-1))
                firdtt += expr * (
                    -n*a*b**2*delta**d*(delta-e)**(ex1-1)*ex1*ex2**2*tau**t *
                    (tau-g)**(2*ex2-2) +
                    n*b**2*d*delta**(d-1)*ex2**2*tau**t*(tau-g)**(2*ex2-2) +
                    2*n*a*b*delta**d*(delta-e)**(ex1-1)*ex1*ex2*t*tau**(t-1) *
                    (tau-g)**(ex2-1) -
                    2*n*b*d*delta**(d-1)*ex2*t*tau**(t-1)*(tau-g)**(ex2-1) +
                    n*a*b*delta**d*(delta-e)**(ex1-1)*ex1*(ex2-1)*ex2*tau**t *
                    (tau-g)**(ex2-2) -
                    n*b*d*delta**(d-1)*(ex2-1)*ex2*tau**t*(tau-g)**(ex2-2) -
                    n*a*delta**d*(delta-e)**(ex1-1)*ex1*(t-1)*t*tau**(t-2) +
                    n*d*delta**(d-1)*(t-1)*t*tau**(t-2))

                expr_ = exp(-a*(delta_0-e)**ex1-b*(tau-g)**ex2)
                B += expr_ * (n*d*delta_0**(d-1)*tau**t -
                              n*a*delta_0**d*(delta_0-e)**(ex1-1)*ex1*tau**t)
                C += expr_ * (
                    n*a**2*delta_0**d*(delta_0-e)**(2*ex1-2)*ex1**2*tau**t -
                    n*a*delta_0**d*(delta_0-e)**(ex1-2)*(ex1-1)*ex1*tau**t -
                    2*n*a*d*delta_0**(d-1)*(delta_0-e)**(ex1-1)*ex1*tau**t +
                    n*(d-1)*d*delta_0**(d-2)*tau**t)

            # Non analitic terms
            # FIXME: Invalid value in critical point with this term
            ni = self._constants.get("nr4", [])
            ai = self._constants.get("a4", [])
            bi = self._constants.get("b4", [])
            Ai = self._constants.get("A", [])
            Bi = self._constants.get("B", [])
            Ci = self._constants.get("C", [])
            Di = self._constants.get("D", [])
            b_ = self._constants.get("beta4", [])

            for n, a, b, A, B, C, D, bt in zip(ni, ai, bi, Ai, Bi, Ci, Di, b_):
                Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
                F = exp(-C*(delta-1)**2-D*(tau-1)**2)
                Fd = -2*C*F*(delta-1)
                Fdd = 2*C*F*(2*C*(delta-1)**2-1)
                Ft = -2*D*F*(tau-1)
                Ftt = 2*D*F*(2*D*(tau-1)**2-1)
                Fdt = 4*C*D*F*(delta-1)*(tau-1)
                Fdtt = 4*C*D*F*(delta-1)*(2*D*(tau-1)**2-1)

                Delta = Tita**2+B*((delta-1)**2)**a
                Deltad = (delta-1)*(A*Tita*2/bt*((delta-1)**2)**(0.5/bt-1) +
                                    2*B*a*((delta-1)**2)**(a-1))
                if delta == 1:
                    Deltadd = 0
                else:
                    Deltadd = Deltad/(delta-1)+(delta-1)**2*(
                        4*B*a*(a-1)*((delta-1)**2)**(a-2) +
                        2*A**2/bt**2*(((delta-1)**2)**(0.5/bt-1))**2 +
                        A*Tita*4/bt*(0.5/bt-1)*((delta-1)**2)**(0.5/bt-2))

                DeltaBd = b*Delta**(b-1)*Deltad
                DeltaBdd = b*(Delta**(b-1)*Deltadd +
                              (b-1)*Delta**(b-2)*Deltad**2)
                DeltaBt = -2*Tita*b*Delta**(b-1)
                DeltaBtt = 2*b*Delta**(b-1)+4*Tita**2*b*(b-1)*Delta**(b-2)
                DeltaBdt = -A*b*2/bt*Delta**(b-1)*(delta-1)*((delta-1)**2)**(
                    0.5/bt-1)-2*Tita*b*(b-1)*Delta**(b-2)*Deltad
                DeltaBdtt = 2*b*(b-1)*Delta**(b-2) * \
                    (Deltad*(1+2*Tita**2*(b-2)/Delta)+4*Tita*A *
                     (delta-1)/bt*((delta-1)**2)**(0.5/bt-1))

                fir += n*Delta**b*delta*F
                fird += n*(Delta**b*(F+delta*Fd)+DeltaBd*delta*F)
                firdd += n*(Delta**b*(2*Fd+delta*Fdd)+2*DeltaBd*(F+delta*Fd) +
                            DeltaBdd*delta*F)
                firt += n*delta*(DeltaBt*F+Delta**b*Ft)
                firtt += n*delta*(DeltaBtt*F+2*DeltaBt*Ft+Delta**b*Ftt)
                firdt += n*(Delta**b*(Ft+delta*Fdt)+delta*DeltaBd*Ft +
                            DeltaBt*(F+delta*Fd)+DeltaBdt*delta*F)
                firdtt += n*(
                    (DeltaBtt*F+2*DeltaBt*Ft+Delta**b*Ftt) +
                    delta*(DeltaBdtt*F+DeltaBtt*Fd+2*DeltaBdt*Ft +
                           2*DeltaBt*Fdt + DeltaBt*Ftt+Delta**b*Fdtt))

                Tita_ = (1-tau)+A*((delta_0-1)**2)**(0.5/bt)
                Delta_ = Tita_**2+B*((delta_0-1)**2)**a
                Deltad_ = (delta_0-1)*(A*Tita_*2/bt*((delta_0-1)**2)**(
                    0.5/bt-1)+2*B*a*((delta_0-1)**2)**(a-1))
                Deltadd_ = Deltad_/(delta_0-1) + (delta_0-1)**2*(
                    4*B*a*(a-1)*((delta_0-1)**2)**(a-2) +
                    2*A**2/bt**2*(((delta_0-1)**2)**(0.5/bt-1))**2 +
                    A*Tita_*4/bt*(0.5/bt-1)*((delta_0-1)**2)**(0.5/bt-2))
                DeltaBd_ = b*Delta_**(b-1)*Deltad_
                DeltaBdd_ = b*(Delta_**(b-1)*Deltadd_ +
                               (b-1)*Delta_**(b-2)*Deltad_**2)
                F_ = exp(-C*(delta_0-1)**2-D*(tau-1)**2)
                Fd_ = -2*C*F_*(delta_0-1)
                Fdd_ = 2*C*F_*(2*C*(delta_0-1)**2-1)

                B += n*(Delta_**b*(F_+delta_0*Fd_)+DeltaBd_*delta_0*F_)
                C += n*(Delta_**b*(2*Fd_+delta_0*Fdd_)+2*DeltaBd_ *
                        (F_+delta*Fd_) + DeltaBdd_*delta_0*F_)

            # Hard sphere term
            if self._constants.get("Fi", None):
                f = self._constants["Fi"]
                n = 0.1617
                a = 0.689
                g = 0.3674
                X = n*delta/(a+(1-a)/tau**g)
                Xd = n/(a+(1-a)/tau**g)
                Xt = n*delta*(1-a)*g/tau**(g+1)/(a+(1-a)/tau**g)**2
                Xdt = n*(1-a)*g/tau**(g+1)/(a+(1-a)/tau**g)**2
                Xtt = -n*delta*((1-a)*g/tau**(g+2)*((g+1)*(a+(1-a)/tau**g) -
                                2*g*(1-a)/tau**g))/(a+(1-a)/tau**g)**3
                Xdtt = -n*((1-a)*g/tau**(g+2)*((g+1)*(a+(1-a)/tau**g) -
                           2*g*(1-a)/tau**g))/(a+(1-a)/tau**g)**3

                ahdX = -(f**2-1)/(1-X) + (f**2+3*f+X*(f**2-3*f))/(1-X)**3
                ahdXX = -(f**2-1)/(1-X)**2 + \
                    (3*(f**2+3*f)+(f**2-3*f)*(1+2*X))/(1-X)**4
                ahdXXX = -2*(f**2-1)/(1-X)**3 + \
                    6*(2*(f**2+3*f)+(f**2-3*f)*(1+X))/(1-X)**5

                fir += (f**2-1)*log(1-X)+((f**2+3*f)*X-3*f*X**2)/(1-X)**2
                fird += ahdX*Xd
                firdd += ahdXX*Xd**2
                firt += ahdX*Xt
                firtt += ahdXX*Xt**2+ahdX*Xtt
                firdt += ahdXX*Xt*Xd+ahdX*Xdt
                firdtt += ahdXXX*Xt**2*Xd+ahdXX*(Xtt*Xd+2*Xdt*Xt)*ahdX*Xdtt

                X_virial = n*delta_0/(a+(1-a)/tau**g)
                ahdX_virial = -(f**2-1)/(1-X_virial) + \
                    (f**2+3*f+X_virial*(f**2-3*f))/(1-X_virial)**3
                ahdXX_virial = -(f**2-1)/(1-X_virial)**2 + \
                    (3*(f**2+3*f)+(f**2-3*f)*(1+2*X_virial))/(1-X_virial)**4
                B += ahdX_virial*Xd
                C += ahdXX_virial*Xd**2

            # Special form from Saul-Wagner Water 58 coefficient equation
            if "nr5" in self._constants:
                delta_0 = 1e-200
                if delta < 0.2:
                    factor = 1.6*delta**6*(1-1.2*delta**6)
                else:
                    factor = exp(-0.4*delta**6)-exp(-2*delta**6)

                nr5 = self._constants.get("nr5", [])
                d5 = self._constants.get("d5", [])
                t5 = self._constants.get("t5", [])
                fr, frt, frtt, frdtt1, frdtt2 = 0, 0, 0, 0, 0
                frd1, frd2 = 0, 0
                frdd1, frdd2, frdd3 = 0, 0, 0
                frdt1, frdt2 = 0, 0
                Bsum1, Bsum2, Csum1, Csum2, Csum3 = 0, 0, 0, 0, 0
                for n, d, t in zip(nr5, d5, t5):
                    fr += n*delta**d*tau**t
                    frd1 += n*delta**(d+5)*tau**t
                    frd2 += n*d*delta**(d-1)*tau**t
                    frdd1 += n*delta**(d+10)*tau**t
                    frdd2 += n*(2*d+5)*delta**(d+4)*tau**t
                    frdd3 += n*d*(d-1)*delta**(d-2)*tau**t
                    frt += n*delta**d*t*tau**(t-1)
                    frtt += n*delta**d*t*(t-1)*tau**(t-2)
                    frdt1 += n*delta**(d+5)*t*tau**(t-1)
                    frdt2 += n*d*delta**(d-1)*t*tau**(t-1)
                    frdtt1 += n*delta**(d+5)*t*(t-1)*tau**(t-2)
                    frdtt2 += n*d*delta**(d-1)*t*(t-1)*tau**(t-2)
                    Bsum1 += n*delta_0**(d+5)*tau**t
                    Bsum2 += n*d*delta_0**(d-1)*tau**t
                    Csum1 += n*delta_0**(d+10)*tau**t
                    Csum2 += n*(2*d+5)*delta_0**(d+4)*tau**t
                    Csum3 += n*d*(d-1)*delta_0**(d-2)*tau**t

                fir += factor*fr
                fird += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6))*frd1 + \
                    factor*frd2
                firdd += (5.76*exp(-0.4*delta**6)-144*exp(-2*delta**6)) * \
                    frdd1 + (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6)) * \
                    frdd2 + factor*frdd3
                firt += factor*frt
                firtt += factor*frtt
                firdt += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6)) * \
                    frdt1 + factor*frdt2
                firdtt += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6)) * \
                    frdtt1 + factor*frdtt2
                B += (-2.4*exp(-0.4*delta_0**6)+12*exp(-2*delta_0**6)) * \
                    Bsum1 + (exp(0.4*delta_0**6)-exp(-2*delta_0**6))*Bsum2
                C += (5.76*exp(-0.4*delta_0**6)-144*exp(-2*delta_0**6)) * \
                    Csum1+(-2.4*exp(-0.4*delta_0**6)+12*exp(-2*delta_0**6)) * \
                    Csum2 + (exp(0.4*delta_0**6)-exp(-2*delta_0**6))*Csum3

        prop = {}
        prop["fir"] = fir
        prop["firt"] = firt
        prop["firtt"] = firtt
        prop["fird"] = fird
        prop["firdd"] = firdd
        prop["firdt"] = firdt
        prop["firdtt"] = firdtt
        prop["B"] = B
        prop["C"] = C
        return prop

    def derivative(self, z, x, y, fase):
        """Calculate generic partial derivative: (δz/δx)y
        where x, y, z can be: P, T, v, u, h, s, g, a"""
        dT = {"P": self.P*fase.alfap,
              "T": 1,
              "v": 0,
              "rho": 0,
              "u": fase.cv,
              "h": fase.cv+self.P*fase.v*fase.alfap,
              "s": fase.cv/self.T,
              "g": self.P*fase.v*fase.alfap-fase.s,
              "a": -fase.s}
        dv = {"P": -self.P*fase.betap,
              "T": 0,
              "v": 1,
              "rho": -1,
              "u": self.P*(self.T*fase.alfap-1),
              "h": self.P*(self.T*fase.alfap-fase.v*fase.betap),
              "s": self.P*fase.alfap,
              "g": -self.P*fase.v*fase.betap,
              "a": -self.P}
        return (dv[z]*dT[y]-dT[z]*dv[y])/(dv[x]*dT[y]-dT[x]*dv[y])

    def _Vapor_Pressure(self, T):
        """Vapor Pressure ancillary equation"""
        if self._vapor_Pressure:
            Pr = SimpleEq(self.Tc, T, self._vapor_Pressure)
            Pv = unidades.Pressure(Pr*self.Pc)
        else:
            # If no available data use the Lee-Kesler correlation
            Pv = Pv_Lee_Kesler(T, self.Tc, self.Pc, self.f_acent)
        return Pv

    def _Liquid_Density(self, T):
        """Saturated liquid density ancillary equation"""
        if self._liquid_Density:
            if T < self.Tt:
                T = self.Tt
            if T > self.Tc:
                T = self.Tc

            Pr = SimpleEq(self.Tc, T, self._liquid_Density)
            rho = unidades.Density(Pr*self.rhoc)
        else:
            # Use the Costald method to calculate saturated liquid density
            rho = RhoL_Costald(T, self.Tc, self.f_acent, 1/self.rhoc)
        return rho

    def _Vapor_Density(self, T):
        """Saturated vapor density ancillary equation"""
        if self._vapor_Density:
            if T < self.Tt:
                T = self.Tt
            if T > self.Tc:
                T = self.Tc

            Pr = SimpleEq(self.Tc, T, self._vapor_Density)
            rho = unidades.Density(Pr*self.rhoc)
        else:
            # If no available data use the ideal gas correlation
            P = self._Vapor_Pressure(T)
            rho = P*1e-3/T/R*self.M
        return rho

    def _Surface(self, T):
        r"""Equation for the surface tension

        .. math::
            \sigma(T) = \sum_i \sigma_i\left(1-\frac{T}{T_c}\right)^{n_i}

        The subclass must define the parameters of correlation in _surface::
            * sigma: Coefficient of polynomial term
            * exp: Exponential of polynomial term

        References
        ----------
        [1]_ Mulero, A., Cachadiña, I., Parra, M.I. Recommended Correlations
            for the Surface Tension of Common Fluids. J. Phys. Chem. Ref. Data
            41(4) (2012) 043105.
        """
        if self.Tt <= self.T <= self.Tc:
            if self._surface:
                tau = 1-T/self.Tc
                tension = 0
                for sigma, n in zip(self._surface["sigma"],
                                    self._surface["exp"]):
                    tension += sigma*tau**n
                sigma = unidades.Tension(tension)
            else:
                # If no available data use the Pitzer correlation
                sigma = Tension_Pitzer(T, self.Tc, self.Pc, self.f_acent)
        else:
            # Undefined out of range of two phase region
            sigma = None
        return sigma

    def _Dielectric(self, rho, T):
        """Calculate the dielectric constant"""
        if self._dielectric:
            if rho < 1e-12:
                e = 1.
            else:
                delta = rho/self.M/self._dielectric["rhoref"]
                tau = T/self._dielectric["Tref"]
                a0 = self._dielectric["a0"]
                expt0 = self._dielectric["expt0"]
                expd0 = self._dielectric["expd0"]
                a1 = self._dielectric["a1"]
                expt1 = self._dielectric["expt1"]
                expd1 = self._dielectric["expd1"]
                a2 = self._dielectric["a2"]
                expt2 = self._dielectric["expt2"]
                expd2 = self._dielectric["expd2"]
                c = 0
                for a, expt, expd in zip(a0, expt0, expd0):
                    c += a * tau**expt * delta**expd
                for a, expt, expd in zip(a1, expt1, expd1):
                    c += a * (tau-1)**expt * delta**expd
                for a, expt, expd in zip(a2, expt2, expd2):
                    c += a * (1./tau-1)**expt * delta**expd
                if self._dielectric["eq"] == 3:
                    e = (1+2*c)/(1-c)
                else:
                    e = 0.25*(1+9*c+3*(9*c**2+2*c+1)**0.5)
        else:
            e = None
        return unidades.Dimensionless(e)

    @classmethod
    def _Melting_Pressure(cls, T):
        """Calculate the melting pressure"""
        if cls._melting:
            Tref = cls._melting["Tref"]
            Pref = cls._melting["Pref"]
            Tita = T/Tref
            suma = 0
            for a, t in zip(cls._melting["a1"], cls._melting["exp1"]):
                suma += a*Tita**t
            for a, t in zip(cls._melting["a2"], cls._melting["exp2"]):
                suma += a*(Tita-1)**t
            for a, t in zip(cls._melting["a3"], cls._melting["exp3"]):
                suma += a*log(Tita)**t

            if cls._melting["eq"] == 1:
                P = suma*Pref
            elif cls._melting["eq"] == 2:
                P = exp(suma)*Pref
            return unidades.Pressure(P, "kPa")
        else:
            return None

    @classmethod
    def _Sublimation_Pressure(cls, T):
        """Calculate the sublimation pressure"""
        if cls._sublimation:
            Tref = cls._sublimation["Tref"]
            Pref = cls._sublimation["Pref"]
            Tita = T/Tref
            suma = 0
            for a, t in zip(cls._sublimation["a1"], cls._sublimation["exp1"]):
                suma += a*Tita**t
            for a, t in zip(cls._sublimation["a2"], cls._sublimation["exp2"]):
                suma += a*(1-Tita)**t
            for a, t in zip(cls._sublimation["a3"], cls._sublimation["exp3"]):
                suma += a*log(Tita)**t

            if cls._sublimation["eq"] == 1:
                P = suma*Pref
            elif cls._sublimation["eq"] == 2:
                P = exp(suma)*Pref
            elif cls._sublimation["eq"] == 3:
                P = exp(Tref/T*suma)*Pref
            return unidades.Pressure(P, "kPa")
        else:
            return None

    # Viscosity calculation methods
    def _Viscosity(self, rho, T, fase):
        r"""Viscosity calculation procedure, implement several general method

        The derived class must define a dict object with the parameters for the
        method.
        The key eq define the procedure to use:

            * 0 - Hardcoded for special procedures (i.e.: R23, H2O, ...)
            * 1 - General formulation with different contribution
            * 2 - Younglove formulation
            * 4 - Quiñones-Cisneros friction theory model

        **Hardcoded procedures**

            Special procedures easier to implement as harcoded in subclasses,
            used in component as: He, D2O, H2, H20, R23, Ethylene

                * method: Name of procedure with code

        **General formulation**

            Normal formulation with several contribution, so the more flexible
            for implement correlation

            .. math::
                \eta = \eta^o(T)+\eta^1(T)\rho+\eta^r(\tau,\delta)+\eta^{CP}

            *Initial density terms, second virial coefficient*

                .. math::
                    \eta^1(T) = B_\eta \eta^0

                .. math::
                    B_\eta = \eta_r B^*_\eta

                .. math::
                    B^*_\eta = \sum_in_i\tau^t

                * Tref_virial: Initial density reference temperature
                * muref_virial: Initial density reference viscosity
                * n_virial: Initial density parameter
                * t_virial: Initial density tau exponent

            *Residual fluid contribution*

                .. math::
                    \eta^r(\tau,\delta) = \sum_{i=1}^nN_i\tau^{t_i}\delta^{d_i}
                    \exp\left(-\gamma_i\delta^{c_i}\right) + \frac{\sum_{i}
                    ^nN_i\tau^{t_i}\delta^{d_i}\exp\left(-\gamma_i\delta^{c_i}
                    \right)}{\sum_{i}^nN_i\tau^{t_i}\delta^{d_i}\exp\left(
                    -\gamma_i\delta^{c_i}\right)}

                * Tref_res: Residual viscosity reference temperature
                * rhoref_res: Residual viscosity reference density
                * muref_res: Residual viscosity reference viscosity
                * nr: Residual viscosity parameter
                * tr: Residual viscosity tau exponent
                * dr: Residual viscosity delta exponent
                * gr: Residual viscosity exponential parameter
                * cr: Reisidual viscosity delta exponent in exponential term

                * nr_num: Fractional numerator coefficient
                * tr_num: Fractional numerator temperature exponent
                * dr_num: Fractional numerator density exponent
                * gr_num: Fractional numerator exponential coefficient
                * cr_num: Fractional numerator exponential density exponent

                * nr_den: Fractional denominator coefficient
                * tr_den: Fractional denominator temperature exponent
                * dr_den: Fractional denominator density exponent
                * gr_den: Fractional denominator exponential coefficient
                * cr_den: Fractional denominator exponential density exponent

            *Modified Batschinkski-Hildebrand contribution*

                .. math::
                    \eta^{CP} = f\left(\frac{\delta}{\delta_0(\tau)-\delta}-
                    \frac{\delta}{\delta_0(\tau)}\right)

                .. math::
                    \delta_0(\tau) = g_1\left(1+\sum_i g_i\tau^{t_i}\right)

                * CPf: f parameter for closed packed term
                * CPg1: g1 parameter for closed packed term
                * CPgi: g parameter for aditional term of closed packed term
                * CPti: tau exponent for aditional term of closed packed term

            *Special terms*
                A hardcoded term for any non stndard formulation term

                * special: Name of procedure with hardcoded method


        **Younglove formulation**
            Formulation as explain in [4]_ and [5]_

            .. math::
                \eta = \eta^o(T)+\eta^1(T)\rho+\eta^r(\tau,\delta)

            .. math::
                \eta^1 = F_{[1]}+F_{[2]}\left(F_{[3]}-\ln\left(\frac{T}
                {F_{[4]}}\right)\right)^2

            .. math::
                \eta^2 = \exp(F)-\exp(G)

            .. math::
                G = E_{[1]}+\frac{E_{[2]}}{T}

            .. math::
                H = \frac{\rho^{0.5}\left(\rho-\rho_c\right)}{\rho_c}

            .. math::
                F = G + \left(E_{[3]}+\frac{E_{[4]}}{T^{1.5}}\right)\rho^{0.1}+
                H \left(E_{[5]}+\frac{E_{[6]}}{T}+\frac{E_{[7]}}{T^2}\right)

            The derived class must define the variables:
                * F: η1 contribution parameters (4 terms)
                * E: η2 contribution parameters (7 terms)
                * rhoc: Reducing density, [mol/l]

        **Quiñones-Cisneros formulation**
            Formulation as explain in [2]_ and [3]_

            .. math::
                \eta = \eta^o + \kappa_iP_{id} + \kappa_r\Delta P_r +
                \kappa_aP_a + \kappa_{ii}P_{id}^2 + \kappa_{rr}\Delta P_r^2 +
                \kappa_{aa}P_a^2 + \kappa_{rrr}P_r^3 + \kappa_{aaa}P_a^3

            .. math::
                \kappa_a = \Gamma\left(a_0+a_1\psi_1+a_2\psi_2\right)

            .. math::
                \kappa_{aa} = \Gamma^3\left(A_0+A_1\psi_1+A_2\psi_2\right)

            .. math::
                \kappa_r = \Gamma\left(b_0+b_1\psi_1+b_2\psi_2\right)

            .. math::
                \kappa_{rr} = \Gamma^3\left(B_0+B_1\psi_1+B_2\psi_2\right)

            .. math::
                \kappa_i = \Gamma\left(c_0+c_1\psi_1+c_2\psi_2\right)

            .. math::
                \kappa_{ii} = \Gamma^3\left(C_0+C_1\psi_1+C_2\psi_2\right)

            .. math::
                \kappa_{rrr} = \Gamma\left(D_0+D_1\psi_1+D_2\psi_2\right)

            .. math::
                \kappa_{aaa} = \Gamma\left(E_0+E_1\psi_1+E_2\psi_2\right)

            .. math::
                \Gamma = \frac{T_c}{T}

            .. math::
                \psi_1 = \exp(\Gamma)-1

            .. math::
                \psi_2 = \exp(\Gamma^2)-1

            The derived class must define the variables:
                * a: Ka correlation coefficients
                * A: Kaa correlation coefficients
                * b: Kr correlation coefficients
                * B: Krr correlation coefficients
                * c: Ki correlation coefficients
                * C: Kii correlation coefficients
                * D: Krrr correlation coefficients
                * E: Kaaa correlation coefficients

        In all cases the dilute-gas limit density is calculate in a separate
        procedure _Visco0, see its docs for its parameters

        Parameters
        ----------
        rho : float
            Density [kg/m³]
        T : float
            Temperature [K]
        fase: dict
            phase properties

        References
        ----------
        [2]_ Quiñones-Cisneros, S.E., Deiters, U.K. Generalization of the
            Friction Theory for Viscosity Modeling. J. Phys. Chem. B, 110(25)
            (2006) 12820-12834
        [3]_ Quiñones-Cisneros, S.E., Huber, M.L., Deiters, U.K. Correlation
            for the Viscosity of Sulfur Hexafluoride (SF6) from the Triple
            Point to 1000 K and Pressures to 50 MPa. J. Phys. Chem. Ref. Data
            41(2) (2012) 023102
        [4]_ Younglove, B.A. Thermophysical Properties of Fluids. I. Argon,
            Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen.
            J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)
        [5]_ Younglove, B.A., Ely, J.F. Thermophysical Properties of Fluids.
            II. Methane, Ethane, Propane, Isobutane, and Normal Butane. J.
            Phys. Chem. Ref. Data 16(4) (1987) 577-798
        """
        coef = self._viscosity
        if coef:
            M = coef.get("M", self.M)
            if coef["eq"] == 0:
                # Hardcoded method
                method = self.__getattribute__(coef["method"])
                mu = method(rho, T, fase).muPas

            elif coef["eq"] == 1:
                muo = self._Visco0()

                # Initial-density term, second virial coefficient
                mud = 0
                if rho and "n_virial" in coef:
                    Tc = coef.get("Tref_virial", self.Tc)
                    if "muref_virial" in coef:
                        mur = coef["muref_virial"]
                    else:
                        mur = Avogadro*(coef["sigma"]*1e-9)**3

                    tau = T/Tc
                    B_ = 0
                    for n, t in zip(coef["n_virial"], coef["t_virial"]):
                        B_ += n*tau**t

                    mud = mur*B_*1e3*rho/M*muo

                # Residual term
                Tr = coef.get("Tref_res", self.Tc)
                tau = Tr/T
                rhor = coef.get("rhoref_res", self.rhoc)
                delta = rho/rhor
                mured = coef.get("muref_res", 1)

                mur = 0
                # polynomial term
                if rho and "nr" in coef:
                    if "cr" in coef:
                        for n, t, d, c, g in zip(
                                coef["nr"], coef["tr"], coef["dr"],
                                coef["cr"], coef["gr"]):
                            mur += n*tau**t*delta**d*exp(-g*delta**c)
                    else:
                        for n, t, d in zip(coef["nr"], coef["tr"], coef["dr"]):
                            mur += n*tau**t*delta**d

                # numerator of rational poly; denominator of rat. poly;
                if rho and "nr_num" in coef:
                    num = 0
                    den = 0
                    if "cr_num" in coef:
                        for n, t, d, c, g in zip(
                                coef["nr_num"], coef["tr_num"], coef["dr_num"],
                                coef["cr_num"], coef["gr_num"]):
                            num += n*tau**t*delta**d*exp(-g*delta**c)
                    else:
                        for n, t, d in zip(coef["nr_num"], coef["tr_num"],
                                           coef["dr_num"]):
                            num += n*tau**t*delta**d

                    if "nr_den" in coef:
                        if "cr_den" in coef:
                            for n, t, d, c, g in zip(
                                    coef["nr_den"], coef["tr_den"],
                                    coef["dr_den"], coef["cr_den"],
                                    coef["gr_den"]):
                                den += n*tau**t*delta**d*exp(-g*delta**c)
                        else:
                            for n, t, d in zip(coef["nr_den"], coef["tr_den"],
                                               coef["dr_den"]):
                                den += n*tau**t*delta**d

                    else:
                        den = 1.

                    mur += num/den
                mur *= mured

                # Modified Batschinkski-Hildebrand terms with reduced
                # close-packed density
                muCP = 0
                if "CPf" in coef:
                    f = coef["CPf"]
                    g1 = coef["CPg1"]
                    gi = coef["CPgi"]
                    ti = coef["CPti"]
                    rest = 0
                    for g, t in zip(gi, ti):
                        rest += g*tau**t
                    delta0 = g1*(1+rest)
                    muCP = f*(delta/(delta0-delta)-delta/delta0)

                # Gaussian terms
                # Reference: Propane, ethane correlation from Vogel
                if "nr_gaus" in coef:
                    for n, b, e in zip(coef["nr_gaus"], coef["br_gaus"],
                                       coef["er_gaus"]):
                        mur += n*tau*delta*exp(-b*(delta-1)**2-e*abs(tau-1))

                # Special term hardcoded in component class
                mue = 0
                if "special" in coef:
                    method = self.__getattribute__(coef["special"])
                    mue = method(rho, T, fase)

                # print(muo, mud, mur, mue)
                mu = muo+mud+mur+muCP+mue

            elif coef["eq"] == 2:
                # Younglove form
                rhom = rho/M
                muo = self._Visco0()

                f = coef["F"]
                e = coef["E"]
                mod = coef.get("mod", False)

                mu1 = f[0]+f[1]*(f[2]-log(T/f[3]))**2                   # Eq 21

                G = e[0]+e[1]/T                                         # Eq 23
                H = rhom**0.5*(rhom-coef["rhoc"])/coef["rhoc"]          # Eq 25
                if mod:
                    # Eq 23 in ref 2
                    F = e[0] + e[1]*H + e[2]*rhom**0.1 + e[3]*H/T**2 + \
                        e[4]*rhom**0.1/T**1.5 + e[5]/T + e[6]*H/T
                else:
                    # Eq 24 in ref 1
                    F = G + (e[2]+e[3]/T**1.5)*rhom**0.1 + \
                        H*(e[4]+e[5]/T+e[6]/T**2)
                mu2 = exp(F)-exp(G)                                     # Eq 22

                mu = muo+mu1*rhom+mu2                                   # Eq 18

            elif self._viscosity["eq"] == 3 or self._viscosity["eq"] == "ecs":
                # Huber ECS model

                Tr = T/self.Tc
                rhor = rho/self.rhoc

                fint = 0
                for n, m in zip(self._viscosity["fint"], self._viscosity["fint_t"]):
                    fint += n*T**m

                psi = 0
                for n, t, d in zip(self._viscosity["psi"], self._viscosity["psi_t"], self._viscosity["psi_d"]):
                    psi += n*Tr**t*rhor**d

                mu = None


            elif coef["eq"] == 4:
                # Quiñones-Cisneros correlations
                muo = self._Visco0()

                Gamma = self.Tc/T                                      # Eq 35
                psi1 = exp(Gamma)-1.0                                  # Eq 33
                psi2 = exp(Gamma**2)-1.0                               # Eq 34

                a = coef["a"]
                b = coef["b"]
                c = coef["c"]
                A = coef["A"]
                B = coef["B"]
                C = coef["C"]

                ka = (a[0] + a[1]*psi1 + a[2]*psi2) * Gamma            # Eq 27
                kaa = (A[0] + A[1]*psi1 + A[2]*psi2) * Gamma**3        # Eq 28
                kr = (b[0] + b[1]*psi1 + b[2]*psi2) * Gamma            # Eq 29
                krr = (B[0] + B[1]*psi1 + B[2]*psi2) * Gamma**3        # Eq 30
                ki = (c[0] + c[1]*psi1 + c[2]*psi2) * Gamma            # Eq 31
                kii = (C[0] + C[1]*psi1 + C[2]*psi2) * Gamma**3        # Eq 32

                # All parameteres has pressure units of bar
                Patt = -fase.IntP.bar
                Prep = T*fase.dpdT_rho.barK
                Pid = rho*self.R*self.T/1e5
                delPr = Prep-Pid

                # Eq 26
                mur = ki*Pid + kr*delPr + ka*Patt + kii*Pid**2 + \
                    krr*delPr**2 + kaa*Patt**2

                # 3rd terms only necessary en SF6 Quiñones corelation
                if "D" in coef:
                    D = coef["D"]
                    E = coef["E"]
                    krrr = (D[0] + D[1]*psi1 + D[2]*psi2) * Gamma       # Eq 17
                    kaaa = (E[0] + E[1]*psi1 + E[2]*psi2) * Gamma       # Eq 18
                    mur += krrr*Prep**3 + kaaa*Patt**3

                mu = (muo+mur*1e3)

        else:
            # Use the Chung method as default method for no high quality method
            # muo = MuL_LetsouStiel(T, self.M, self.Tc, self.Pc, self.f_acent)
            muo = MuG_Chung(T, self.Tc, 1/self.rhoc, self.M, self.f_acent,
                            self.momentoDipolar.Debye, 0)
            mu = MuG_P_Chung(T, self.Tc, 1/self.rhoc, self.M, self.f_acent,
                             self.momentoDipolar.Debye, 0, rho, muo).muPas
        return unidades.Viscosity(mu, "muPas")

    def _Visco0(self, coef=None):
        r"""Dilute gas viscosity calculation

        .. math::
            \eta^o(T) = \frac{N_{chapman}\left(MT\right)^{t_{chapman}}}
            {\sigma^2\Omega(T^*)} + \sum_{i}^nN_i\tau^{t_i}
            + \frac{\sum_{i}^nN_i\tau^{t_i}}{\sum_{i}^nN_i\tau^{t_i}}
            + \eta_̣{hc}^o

        .. math::
            T^* = \frac{T}{ε/k}

        .. math::
            \tau = \frac{T}{T_r}

        The collision integral is calculated in separate procedure _Omega, see
        its docs for parameters

        Notes
        -----
        The derived class must define a dict object with the parameters for the
        method

        Champman-Enskop:
            * ek: Lennard-Jones energy parameter, [K]
            * sigma: Lennard-Jones size parameter, [nm]
            * Tref: Dilute gas viscosity reference temperature, default=1
            * rhoref: Dilute gas viscosity reference density, default=1
            * n_chapman: Dilute gas viscosity parameter, default=0.0266958
            * t_chapman: Dilute gas viscosity temperature exponent, default=0.5

        Polynomial terms:
            * Toref: Polynomial terms reference temperature, default=1
            * no: Polynomial terms coefficient
            * to: Polynomial terms tau exponent

        Fractional terms:
            * no_num: Fractional numerator coefficient for dilute gas viscosity
            * to_num: Fractional numerator temperature exponent
            * no_den: Fraction denominator coefficient for dilute gas viscosity
            * to_den: Fraction denominator temperature exponent

        Custom Hardcoded terms:
            * special0: Name of procedure with hardcoded method
        """
        if coef is None:
            coef = self._viscosity

        muo = 0
        M = coef.get("M", self.M)

        # Collision integral calculation
        if coef["omega"]:
            omega = self._Omega(coef)
            N_chap = coef.get("n_chapman", 0.0266958)
            t_chap = coef.get("t_chapman", 0.5)
            Tr = coef.get("Tref", 1.)

            muo += N_chap*(M*self.T/Tr)**t_chap/(coef["sigma"]**2*omega)

        tau = self.T/coef.get("Toref", 1.)
        # other aditional polynomial terms
        if "no" in coef:
            for n, t in zip(coef["no"], coef["to"]):
                muo += n*tau**t

        # Other fractional terms
        if "no_num" in coef:
            num = 0
            den = 0
            for n, t in zip(coef["no_num"], coef["to_num"]):
                num += n*tau**t
            if "no_den" in coef:
                for n, t, in zip(coef["no_den"], coef["to_den"]):
                    den += n*tau**t
            else:
                den = 1
            muo += num/den

        # Special hardcoded method:
        if "special0" in coef:
            method = self.__getattribute__(coef["special0"])
            muo += method()

        return muo

    def _Omega(self, coef):
        r"""Collision integral calculations

        The derived class must define a dict object with the parameters for the
        method:
            * omega: Collision integral procedure calculation
                * 0 - None
                * 1 - Lemmon expression from []: Air, N2, O2...
                * 2 - Younglove expression from []: C1, C2, C3, iC4, nC4
                * 3 - Inverse temperature polinomial form exponential (cC6)
                * 4 - Inverse temperature polinomial form (H2S)
                * 5 - Neufeld correlation for Lennard-Jones 6-12 potential

            * collision: Alternate contributions values for collision integral

        .. math::
            \Omega_1 = \exp\left(\sum_i^nb_i\ln\left(T/εk\right)^i\right)

        .. math::
            \Omega_2 = \frac{1}{\sum_i^n\left(b_i\frac{εk}{T}\right)^{(n+2)/3}}

        .. math::
            \Omega_3 = \exp\left(\sum_{i=0}^n \frac{b_i}{T_r^i}\right)

        .. math::
            \Omega_4 = \sum_{i=0}^n \frac{b_i}{T_r^i}

        .. math::
            \Omega_5 = \frac{1.16145}{T_r^{0.14874}} + \frac{0.52487}
            {\exp(0.7732T_r)} + \frac{2.16178}{\exp(2.4378T_r)} - 6.435e^{-4}
            T_r^{0.14874}\sin\left(\frac{18.0323}{T_r^{0.76830}}-7.27371\right)

        References
        ----------
        [6]_ Lemmon, E.W., Jacobsen, R.T. Viscosity and Thermal Conductivity
            Equations for Nitrogen, Oxygen, Argon, and Air. Int. J.
            Thermophys., 25(1) (2004) 21-69
        [5]_ Younglove, B.A., Ely, J.F. Thermophysical Properties of Fluids.
            II. Methane, Ethane, Propane, Isobutane, and Normal Butane. J.
            Phys. Chem. Ref. Data 16(4) (1987) 577-798
        [7]_ Tariq, U., Jusoh, A.R.B., Rittliesco, N., Vesovic, V. Reference
            Correlation of the Viscosity of Cyclohexane from the Triple Point
            to 700K and up to 110 MPa. J. Phys. Chem. Ref. Data 43(3) (2014)
            033101
        [3]_ Quiñones-Cisneros, S.E., Huber, M.L., Deiters, U.K. Correlation
            for the Viscosity of Sulfur Hexafluoride (SF6) from the Triple
            Point to 1000 K and Pressures to 50 MPa. J. Phys. Chem. Ref. Data
            41(2) (2012) 023102
        [8]_ Neufeld, P.D., Janzen, A.R., Aziz, R.A. Empirical Equations to
            Calculate 16 of the Transport Collision Integrals Ω for the
            Lennard-Jones (12-6) Potential. J. Chem. Phys. 57(3) (1972)
            1100-1102
        """
        b = coef.get("collision", None)
        if coef["omega"] == 1:
            if not b:
                b = [0.431, -0.4623, 0.08406, 0.005341, -0.00331]
            T_ = log(self.T/coef["ek"])
            suma = 0
            for i, bi in enumerate(b):
                suma += bi*T_**i
            omega = exp(suma)
        elif coef["omega"] == 2:
            if not b:
                # Pag 580, Table 2 from [2]
                b = [-3.0328138281, 16.918880086, -37.189364917, 41.288861858,
                     -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390,
                     -0.0096570437074]
            T_ = coef["ek"]/self.T
            suma = 0
            for i, bi in enumerate(b):
                suma += bi*T_**((3.-i)/3.)
            omega = 1./suma

        elif self._viscosity["omega"] == 3:
            Tr = self.T/coef.get("ek", 1)
            suma = 0
            for i, bi in enumerate(b):
                suma += bi/Tr**i
            omega = exp(suma)

        elif self._viscosity["omega"] == 4:
            Tr = self.T/coef.get("ek", 1)
            omega = 0
            for i, bi in enumerate(b):
                omega += bi/Tr**i

        elif self._viscosity["omega"] == 5:
            T_ = self.T/coef["ek"]
            omega = Collision_Neufeld(T_)

        return omega

    # Thermal conductivity methods
    def _ThCond(self, rho, T, fase):
        r"""Thermal conductivity calculation procedure

        The derived class must define a dict object with the parameters for the
        method.
        The key eq define the procedure to use:

            * 0 - Hardcoded for special procedures (i.e.: R23, H2O, ...)
            * 1 - General formulation with different contribution
            * 2 - Younglove formulation

        **Hardcoded procedures**

            Special procedures easier to implement as harcoded in subclasses,
            used in component as:

                * method: Name of procedure with code

        **General formulation**

            Normal formulation with several contribution, so the more flexible
            for implement correlations

            .. math::
                \eta = \eta^o(T)+\eta^1(T)\rho+\eta^r(\tau,\delta)+\eta^{CP}

            *Initial density terms, second virial coefficient*

                .. math::
                    \eta^1(T) = B_\eta \eta^0

                .. math::
                    B_\eta = \eta_r B^*_\eta

                .. math::
                    B^*_\eta = \sum_in_i\tau^t

                * Tref_virial: Initial density reference temperature
                * muref_virial: Initial density reference viscosity
                * n_virial: Initial density parameter
                * t_virial: Initial density tau exponent

            *Residual fluid contribution*

                .. math::
                    \eta^r(\tau,\delta) = \sum_{i=1}^nN_i\tau^{t_i}\delta^{d_i}
                    \exp\left(-\gamma_i\delta^{c_i}\right) + \frac{\sum_{i}
                    ^nN_i\tau^{t_i}\delta^{d_i}\exp\left(-\gamma_i\delta^{c_i}
                    \right)}{\sum_{i}^nN_i\tau^{t_i}\delta^{d_i}\exp\left(
                    -\gamma_i\delta^{c_i}\right)}

                * Tref_res: Residual viscosity reference temperature
                * rhoref_res: Residual viscosity reference density
                * muref_res: Residual viscosity reference viscosity
                * nr: Residual viscosity parameter
                * tr: Residual viscosity tau exponent
                * dr: Residual viscosity delta exponent
                * gr: Residual viscosity exponential parameter
                * cr: Reisidual viscosity delta exponent in exponential term

                * nr_num: Fractional numerator coefficient
                * tr_num: Fractional numerator temperature exponent
                * dr_num: Fractional numerator density exponent
                * gr_num: Fractional numerator exponential coefficient
                * cr_num: Fractional numerator exponential density exponent

                * nr_den: Fractional denominator coefficient
                * tr_den: Fractional denominator temperature exponent
                * dr_den: Fractional denominator density exponent
                * gr_den: Fractional denominator exponential coefficient
                * cr_den: Fractional denominator exponential density exponent

            *Modified Batschinkski-Hildebrand contribution*

                .. math::
                    \eta^{CP} = f\left(\frac{\delta}{\delta_0(\tau)-\delta}-
                    \frac{\delta}{\delta_0(\tau)}\right)

                .. math::
                    \delta_0(\tau) = g_1\left(1+\sum_i g_i\tau^{t_i}\right)

                * CPf: f parameter for closed packed term
                * CPg1: g1 parameter for closed packed term
                * CPgi: g parameter for aditional term of closed packed term
                * CPti: tau exponent for aditional term of closed packed term

            *Special terms*
                A hardcoded term for any non stndard formulation term

                * special: Name of procedure with hardcoded method


        **Younglove formulation**
            Formulation as explain in [4]_ and [5]_

            .. math::
                \lambda = \lambda_0 + \frac{\left(F_0+F_1\rho\right)\rho}{1-F_2\rho}+\lambda_c

            .. math::
                \lambda_0 = 1000\eta_0\left(C_p^0-2.5R\right)\left(G_1+G_2\epsilon/kT\right]

            .. math::
                F_0 = \sum_{n=1}^3E_nT^{1-n}

            .. math::
                F_1 = \sum_{n=4}^6E_nT^{1-n}

            .. math::
                F_2 = \sum_{n=7}^8E_nT^{1-n}

            The derived class must define the variables:
                * G: λ0 contribution parameters (2 terms)
                * E: λ1 contribution parameters (8 terms)

        Parameters
        ----------
        rho : float
            Density [kg/m³]
        T : float
            Temperature [K]
        fase: dict
            phase properties

        References
        ----------
        [4]_ Younglove, B.A. Thermophysical Properties of Fluids. I. Argon,
            Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen.
            J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)
        [5]_ Younglove, B.A., Ely, J.F. Thermophysical Properties of Fluids.
            II. Methane, Ethane, Propane, Isobutane, and Normal Butane. J.
            Phys. Chem. Ref. Data 16(4) (1987) 577-798
        """
        coef = self._thermal

        if coef:
            # Let define viscosity coefficient to use in thermal conductivity
            if "visco" in coef:
                visco = coef["visco"]
            else:
                visco = None

            if coef["eq"] == 0:
                # Hardcoded method
                method = self.__getattribute__(coef["method"])
                k = method(rho, T, fase)

            elif coef["eq"] == 1:
                Tr = T/coef.get("Toref", 1)

                # Dilute gas terms
                kg = 0

                # Special term with dilute-gas limit viscosity value
                # Reference from Lemmon for Air, N2, O2...
                if "no_visco" in coef:
                    muo = self._Visco0(visco)
                    kg += coef["no_visco"]*muo

                # Special term with dilute-gas viscosity and ideal gas isobaric
                # heat capacity
                # From Frield echane meos
                if "no_viscoCp" in coef:
                    # Contribution, Eq 14
                    f = 0
                    for n, t in zip(coef["no_viscoCp"], coef["to_viscoCp"]):
                        f += n*Tr**t

                    muo = self._Visco0(visco)
                    # 1e-6 factor because viscosity is in μPa·s
                    n = 1e-6*self.R/u/Avogadro
                    kg += muo*n*(3.75+f*(self.cp0/self.R-2.5))

                # Polynomial terms
                if "to" in coef:
                    for n, t in zip(coef["no"], coef["to"]):
                        kg += n*Tr**t

                # Fractional terms
                if "no_num" in coef:
                    num = 0
                    for n, t in zip(coef["no_num"], coef["to_num"]):
                        num += n*Tr**t

                    den = 0
                    for n, t in zip(coef["no_den"], coef["to_den"]):
                        den += n*Tr**t
                    kg += num/den

                        # if c == -99:
                            # cpi = 1.+n*(self.cp0.kJkgK-2.5*self.R.kJkgK)
                            # kg *= cpi
                        # elif c == -98:
                            # muo = self._Visco0()
                            # cpi = self.cp0-R
                            # kg += n*cpi*muo
                        # elif c == -96:
                            # cpi = self.cp0/self.R-2.5
                            # muo = self._Visco0()
                            # kg = (kg*cpi+15./4.)*self.R.kJkgK*muo/self.M
                        # else:

                kg *= coef.get("koref", 1)

                # Backgraund terms
                kr = 0
                kc = 0
                if rho > 0:
                    if "nr" in coef:
                        Tr = coef.get("Tref_res", self.Tc)
                        tau = Tr/T
                        rhor = coef.get("rhoref_res", self.rhoc)
                        delta = rho/rhor

                        if "cr" in coef:
                            for n, t, d, c, g in zip(
                                    coef["nr"], coef["tr"], coef["dr"],
                                    coef["cr"], coef["gr"]):
                                kr += n*tau**t*delta**d*exp(-g*delta**c)
                        else:
                            for n, t, d, in zip(
                                    coef["nr"], coef["tr"], coef["dr"]):
                                kr += n*tau**t*delta**d

                    # Special term with density of saturated vapor factor
                    # Rererence from Friend correlation for methane
                    if "nr_s" in coef:
                        if T < self.Tc and rho < self.rhoc:
                            delta_s = self._Vapor_Density(T)
                        else:
                            delta_s = 11

                        for n, t, d, in zip(
                                coef["nr_s"], coef["tr_s"], coef["dr_s"]):
                            kr += n*tau**t*delta**d/delta_s


                        # if "nbden" in self._thermal:
                            # den = 0
                            # for n, t, d in zip(self._thermal["nbden"], self._thermal["tbden"], self._thermal["dbden"]):
                                # den += n*tau**t*delta**d
                            # kb /= den

                    kr *= coef.get("kref_res", 1)

                    # Critical enhancement
                    kc = self._KCritical(rho, T, fase)

                # Special term hardcoded in component class
                ke = 0
                if "special" in coef:
                    method = self.__getattribute__(coef["special"])
                    ke = method(rho, T, fase)

                # print(kg, kr, kc, ke, kg+kr)
                k = kg+kr+kc+ke

            elif self._thermal["eq"] == 2:
                # Younglove form
                muo = self._Visco0(coef["visco"])
                print(muo, self.cp0, self.R)

                # Eq 27
                kg = 1e-6*muo*(self.cp0-2.5*self.R) * \
                    (coef["G"][0]+coef["G"][1]*coef["visco"]["ek"]/T)

                if rho:
                    E = self._thermal["E"]
                    F0 = E[0] + E[1]/T + E[2]/T**2                      # Eq 28
                    F1 = E[3]/T**3 + E[4]/T**4 + E[5]/T**5              # Eq 29
                    F2 = E[6]/T**6 + E[7]/T**7                          # Eq 30

                    rhom = rho/self.M
                    kb = (F0+F1*rhom)*rhom/(1-F2*rhom)

                    # Critical enhancement
                    kc = self._KCritical(rho, T, fase)

                else:
                    kb = 0
                    kc = 0

                print(kg, kb, kc)
                k = kg+kb+kc

            elif self._thermal["eq"] == 3:
                T_ = self._thermal["ek"]/T
                rhom = rho/self.M
                suma = 0
                for i, bi in enumerate(self._thermal["b"]):
                    suma += bi*T_**((3.-i)/3.)
                ko = self._thermal["Nchapman"]*T**self._thermal["tchapman"]/(self._thermal["sigma"]**2/suma)

                f = self._thermal["F"]
                e = self._thermal["E"]
                k1 = f[0]+f[1]*(f[2]-log(T/f[3]))**2
                kg = ko+k1*rhom

                G = e[0]+e[1]/T
                H = rhom**0.5*(rhom-self._thermal["rhoc"])/self._thermal["rhoc"]
                F = G+(e[2]+e[3]*T**-1.5)*rhom**0.1+H*(e[4]+e[5]/T+e[6]/T**2)
                kb = exp(F)-exp(G)

                bl = self._thermal["ff"]*(self._thermal["rm"]**5*rho/1000.*Avogadro/self.M*self._thermal["Nchapman"]/T)**0.5
                y = 6.0*pi*fase.mu/100000.*bl*(Boltzmann*T*rho/1000.*Avogadro/self.M)**0.5
                deltaL = 0.0
                der = self.derivative("P", "rho", "T", fase)
                if der > 0:
                    deltaL = Boltzmann*(T*fase.dpdT_rho)**2/(rho/1000.*der)**0.5/y
                else:
                    deltaL = 0.0
                kc = deltaL*exp(-18.66*((rho-self.rhoc)/self.rhoc)**4-4.25*((T-self.Tc)/self.Tc)**2)

                k = kg+kb+kc

            elif self._thermal["eq"] == "ecs":
                # TODO
                k = 0

        else:
            # Use the Chung method as default method for no high quality method
            # ko = ThL_Pachaiyappan(T, self.Tc, self.M, rho, self.branched)
            # muo = MuL_LetsouStiel(T, self.M, self.Tc, self.Pc, self.f_acent)
            # muo = MuG_Chung(T, self.Tc, 1/self.rhoc, self.M, self.f_acent,
                            # self.momentoDipolar.Debye, 0)
            # ko = ThG_Chung(T, self.Tc, self.M, self.f_acent, fase.cv, muo)
            # k = ThG_P_Chung(T, self.Tc, 1/self.rhoc, self.M, self.f_acent,
                            # self.momentoDipolar.Debye, 0, rho, ko)
            k = 0
        return unidades.ThermalConductivity(k)

    def _KCritical(self, rho, T, fase):
        r"""Enchancement thermal conductivity calculation for critical region
        The model is defined in _thermal["critical"]. The defined models are:
            * 0 : No critical enhancement
            * 1 : Gaussian term, used in Laesecke correlation for R123
            * 3 : Olchowy-Sengers

        **Olchowy-Sengers**

            .. math::
                \lambda_c = \rho c_p \frac{R_ok_BT}{6\pi\eta\xi}
                \left(\bar{\Omega}-\bar{\Omega}_0\right)

            .. math::
                \bar{\Omega} = \frac{2}{\pi}\left[\left(\frac{c_p-c_v}{c_p}
                \right)\arctan\left(q_D\xi\right)+\frac{c_v}{c_p}q_D\xi\right]

            .. math::
                \bar{\Omega}_0 = \frac{2}{\pi}\left[1-\exp\left(-\frac{1}
                {\frac{1}{q_D\xi}+\frac{(q_D\xi\rho_c)^2}{\rho^23}}\right)
                \right]

            .. math::
                \xi = \xi_0\left(\frac{\Delta\tilde{\chi}}{\Gamma}\right)
                ^{v/\gamma}

            .. math::
                \Delta\tilde{\chi} = \tilde{\chi}(T,\rho) - \tilde{\chi}
                (T_R,\rho)\frac{T_R}{T}

            .. math::
                \tilde{\chi}(T,\rho) = \frac{P_c\rho}{\rho_c^2}\left(\frac
                {\partial\rho}{\partial P}\right)_T

            The derived class must define the variables:
                * Tcref: Reference temperature far above the Tc, [K]
                * gnu: critical exponent parameter, [-]
                * gam0: critical exponent parameter, [-]
                * Xio: Amplitude of the asymptotic power laws for ξ, [m]
                * gamma: Amplitude of the asymptotic power laws for χ, [-]
                * R0: Amplitude parameter, [-]
                * qd : finite upper cutoff, [m]

        ** Gaussian terms**

            .. math::
                \ln\frac{\lambda_c}{\lambda_r} = \sum n\left(\tau-\alfa\right)
                ^t \left(\delta-\beta\right)^d

            The derived class must define the variables:
               * Trefc: Reference temperature, [K]
               * rhorefc: Reference density, [kg/m³]
               * krefc: Reference thermal conductivity, [W/m·K]
               * nc: Polynomial coefficient
               * alfac: τ term translation
               * tc: τ exponent
               * betac: δ term translation
               * dc: δ exponent

        References
        ----------
                   "autor": "Laesecke, A., Perkins, R.A., Howley, J.B.",
                   "title": "An improved correlation for the thermal "
                            "conductivity of HCFC123 (2,2-dichloro-1,1,1-"
                            "trifluoroethane)",
                   "ref": "Int. J. Refrigeration 19(4) (1996) 231-238",
        [] Olchowy, G.A., Sengers, J.V. A Simplified Representation for the
            Thermal Conductivity of Fluids in the Critical Region. Int. J.
            Thermophys. 10(2) (1989) 417-426

        """
        if self._thermal["critical"] == 0:
            # No critical enhancement
            tc = 0

        elif self._thermal["critical"] == 1:
            # Empirical formulation, referenced in Laesecke for R123
            tau = self._thermal["Trefc"]/T
            delta = rho/self._thermal["rhorefc"]

            expo = 0
            for n, alfa, t, beta, d in zip(
                    self._thermal["nc"], self._thermal["alfac"],
                    self._thermal["tc"], self._thermal["betac"],
                    self._thermal["dc"]):
                expo += n*(tau+alfa)**t*(delta+beta)**d
            tc = self._thermal["krefc"]*exp(expo)

        elif self._thermal["critical"] == 2:
            # Younglove form for Propane
            rhom = rho/self.M
            Tc = self._thermal.get("Tcref", self.Tc)
            rhoc = self._thermal.get("rhocref", self.rhoc)

            # Eq D4
            delT = abs(T-Tc)/Tc
            delrho = abs(rho-rhoc)/rhoc

            X = self._thermal["X"]
            xi = self.Pc*rho/self.rhoc**2/fase.drhodP_T                 # Eq D3

            # Eq D2
            DLc = X[3]*Boltzmann/self.Pc * \
                (T/Tc*fase.drhodT_P*self.rhoc/rho)**2*xi**X[2]
            tc = DLc*exp(-X[0]*delT**4 - X[1]*delrho**4)
            #/(6*pi*fase.mu*self._thermal["Z"])

        elif self._thermal["critical"] == 3:
            # fase.mu = unidades.Viscosity(213.0207196, "muPas")
            # Olchowy-Sengers
            qd = self._thermal["qd"]
            Tref = self._thermal["Tcref"]
            Xio = self._thermal["Xio"]
            gam0 = self._thermal["gam0"]
            gnu = self._thermal["gnu"]
            gamma = self._thermal["gamma"]
            R0 = self._thermal["R0"]

            # Optional critical parameters
            Pc = self._thermal.get("Pc", self.Pc)
            rhoc = self._thermal.get("rhoc", self.rhoc)

            Xi = Pc*rho/rhoc**2*fase.drhodP_T

            st = self._Helmholtz(rho, Tref)
            drho = 1/st["dpdrho"]
            Xi_Tr = Pc*rho/rhoc**2*drho

            delchi = Xi-Xi_Tr*Tref/T                                    # Eq 10
            if delchi <= 0:
                tc = 0
            else:
                Xi = Xio*(delchi/gam0)**(gnu/gamma)                      # Eq 9
                Xq = Xi/qd

                # Eq 7
                omega = 2/pi*((fase.cp-fase.cv)/fase.cp*arctan(Xq) +
                              fase.cv/fase.cp*Xq)

                # Eq 8
                omega0 = 2/pi*(1-exp(-1/(1/Xq+Xq**2/3*(rhoc/rho)**2)))

                # Eq 6
                Kb = 1.380658e-23
                tc = rho*fase.cp*Kb*R0*T/(6*pi*Xi*fase.mu) * \
                    (omega-omega0)

        elif self._thermal["critical"] == 4:
            # Younglove form
            rhom = rho/self.M
            Tc = self._thermal["Tcref"]
            rhoc = self._thermal["rhocref"]

            # Eq D4
            delT = abs(T-Tc)/Tc
            delrho = abs(rhom-rhoc)/rhoc

            Xt = (rhom*self._thermal["Pcref"]/self._thermal["rhocref"]**2/self.derivative("P", "rho", "T", fase))**self._thermal["expo"]
            parterm = self._thermal["alfa"]*Boltzmann/self._thermal["Pcref"]*(T*fase.dpdT*self._thermal["rhocref"]/rhom)**2*Xt*1e21
            expterm = exp(-(self._thermal["alfa"]*delT**2+self._thermal["beta"]*delrho**4))
            tc = parterm*expterm/(6*pi*self._thermal["Xio"]*fase.mu.muPas)*self._thermal["kcref"]

               # "critical": 4,
               # "Tcref": 150.86, "Pcref": 4905.8, "rhocref": 13.41, "kcref": 1e-3,
               # "gamma": 1.02,
               # "expo": 0.46807, "alfa": 39.8, "beta": 5.45, "Xio": 6.0795e-1}

        elif self._thermal["critical"] == "CH4":
            tau = self.Tc/T
            delta = rho/self.rhoc
            ts = (self.Tc-T)/self.Tc
            ds = (self.rhoc-rho)/self.rhoc
            xt = 0.28631*delta*tau/self.derivative("P", "rho", "T", fase)
            ftd = exp(-2.646*abs(ts)**0.5+2.678*ds**2-0.637*ds)
            tc = 91.855/fase.mu/tau**2*fase.dpdT_rho**2*xt**0.4681*ftd*1e-3

        elif isinstance(self._thermal["critical"], str):
            # Hardcoded method
            method = self.__getattribute__(self._thermal["critical"])
            tc = method(rho, T, fase)
        return tc


class MEoSBlend(MEoS):
    """Special meos class to implement pseudocomponent blend and defining its
    ancillary dew and bubble point"""
    @classmethod
    def _dewP(cls, T, eq=0):
        """Using ancillary equation return the pressure of dew point"""
        c = cls.eq[eq]["dew"]
        Tj = cls.eq[eq]["Tj"]
        Pj = cls.eq[eq]["Pj"]
        Tita = 1-T/Tj

        suma = 0
        for i, n in zip(c["i"], c["n"]):
            suma += n*Tita**(i/2.)
        P = Pj*exp(Tj/T*suma)
        return unidades.Pressure(P, "MPa")

    @classmethod
    def _bubbleP(cls, T, eq=0):
        """Using ancillary equation return the pressure of bubble point"""
        c = cls.eq[eq]["bubble"]
        Tj = cls.eq[eq]["Tj"]
        Pj = cls.eq[eq]["Pj"]
        Tita = 1-T/Tj

        suma = 0
        for i, n in zip(c["i"], c["n"]):
            suma += n*Tita**(i/2.)
        P = Pj*exp(Tj/T*suma)
        return unidades.Pressure(P, "MPa")

    @classmethod
    def _Vapor_Pressure(cls, T, eq=0):
        # Use a mean value between dew and bubble pressure to get a correct
        # initial phase checking
        dew = cls._dewP(T, eq)
        bubble = cls._bubbleP(T, eq)
        Pv = unidades.Pressure((dew+bubble)/2)
        return Pv


data = MEoS.properties()
propiedades = [p[0] for p in data]
keys = [p[1] for p in data]
units = [p[2] for p in data]
properties = dict(list(zip(keys, propiedades)))
inputData = [data[0], data[2], data[4], data[5], data[6], data[7], data[8],
             data[9]]


if __name__ == "__main__":
#Calculate estado estandard  OTO (h,s=0 at 25ºC and 1 atm)
#    import pickle
#    std={}
#    for compuesto in __all__:
#        method={}
#        for metodo in range(len(compuesto.eq)):
#            cmp=compuesto(T=298.15, P=0.101325, eq=metodo, ref=None)
#            method[str(metodo)]=(cmp.h.kJkg, cmp.s.kJkgK)
#        std[compuesto.__name__]=method
#    cPickle.dump(std, open("/home/jjgomera/Programacion/pychemqt/oto.pkl", "w"))

#Calculate estado estandard  NBP (h,s=0 saturated liquid at Tb)
#    std={}
#    for compuesto in __all__:
#        method={}
#        for metodo in range(len(compuesto.eq)):
#            try:
#                cmp=compuesto(T=compuesto.Tb, x=0.,ref=None)
#                method[str(metodo)]=(cmp.h.kJkg, cmp.s.kJkgK)
#            except:
#                print compuesto.__name__, metodo
#                method[str(metodo)]=(0, 0)
#        std[compuesto.__name__]=method
#    cPickle.dump(std, open("/home/jjgomera/Programacion/pychemqt/oto.pkl", "w"))

#Calculate estado estandard  IIR (h=200,s=1 saturated liquid 0ºC)
#    std={}
#    for compuesto in __all__:
#        cmp=compuesto(T=273.15, x=0., ref=None)
#        std[compuesto.__name__]=(cmp.h.kJkg+200, cmp.s.kJkgK+1)
#    cPickle.dump(std, open("/home/jjgomera/Programacion/pychemqt/IIR.pkl", "w"))

#Calculate estado estandard  ASHRAE (h,s=0 saturated liquid at -40ºC)
#    std={}
#    for compuesto in __all__:
#        try:
#            cmp=compuesto(T=233.15, x=0., ref=None)
#            std[compuesto.__name__]=(cmp.h.kJkg, cmp.s.kJkgK)
#        except:
#            print compuesto.__name__
#            std[compuesto.__name__]=(0, 0)
#    cPickle.dump(std, open("/home/jjgomera/Programacion/pychemqt/ASHRAE.pkl", "w"))

    from lib.mEoS import Ethylene
    st = Ethylene(T=105, P=1e5, eq="jahangiri")
    print(st.x, st.h.kJkg, st.s.kJkgK, st.msg)
    print(st._Cp0()*st.M)
