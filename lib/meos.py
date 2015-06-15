#!/usr/bin/python
# -*- coding: utf-8 -*-

#############################################################################
# Implemented multiparameter equation of state
#   o   Ecuación de estado Setzmann-Wagner, basada en la energía de Helmholtz
#   o   Ecuación MBWR
#   o   Ecuación Peng-Robinson con translación de Peneloux
#############################################################################

import cPickle
import os
from itertools import product
from PyQt4.QtGui import QApplication
from scipy import exp, log, log10, sin, sinh, cosh, tanh, arctan, __version__
if int(__version__.split(".")[1]) < 10:
    from scipy.constants import Bolzmann as Boltzmann
else:
    from scipy.constants import Boltzmann
from scipy.constants import pi, Avogadro, R
from scipy.optimize import fsolve

from lib import unidades, compuestos
from physics import R_atml
from config import Fluid

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
        (QApplication.translate("pychemqt", "Fugacity coef."), "fi", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Fugacity"), "f", unidades.Pressure),
        (QApplication.translate("pychemqt", "Isoentropic exponent"), "gamma", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Vaporization heat"), "Hvap", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Volumetric Expansitivy"), "alfav", unidades.InvTemperature),
        (QApplication.translate("pychemqt", "Isotermic compresibility"), "kappa", unidades.InvPressure),
        (QApplication.translate("pychemqt", "Relative pressure"), "alfap", unidades.InvTemperature),
        (QApplication.translate("pychemqt", "Isothermal stress"), "betap", unidades.Density),
        (QApplication.translate("pychemqt", "Isentropic temperature-pressure"), "betas", unidades.TemperaturePressure),
        (QApplication.translate("pychemqt", "Joule-Thomson coefficient"), "joule", unidades.TemperaturePressure),
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
        (QApplication.translate("pychemqt", "Isothermal expansion"), "kt", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Adiabatic compresibility"), "ks", unidades.InvPressure),
        (QApplication.translate("pychemqt", "Isothermal modulus"), "Ks", unidades.Pressure),
        (QApplication.translate("pychemqt", "Adiabatic modulus"), "Kt", unidades.Pressure),
#        Z_rho     -   (Z-1) over the density, m³/kg
        (QApplication.translate("pychemqt", "Internal pressure"), "IntP", unidades.Pressure),
        (QApplication.translate("pychemqt", "Negative reciprocal temperature"), "invT", unidades.InvTemperature),
        (QApplication.translate("pychemqt", "Specific heat input"), "hInput", unidades.Enthalpy),
        (QApplication.translate("pychemqt", "Viscosity"), "mu", unidades.Viscosity),
        (QApplication.translate("pychemqt", "Thermal conductivity"), "k", unidades.ThermalConductivity),
        (QApplication.translate("pychemqt", "Kinematic viscosity"), "nu", unidades.Diffusivity),
        (QApplication.translate("pychemqt", "Surface tension"), "sigma", unidades.Tension),
        (QApplication.translate("pychemqt", "Thermal diffusivity"), "alfa", unidades.Diffusivity),
        (QApplication.translate("pychemqt", "Prandtl number"), "Prandt", unidades.Dimensionless),
        (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", unidades.Dimensionless),
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

propiedades = [p[0] for p in data]
keys = [p[1] for p in data]
units = [p[2] for p in data]
properties = dict(zip(keys, propiedades))
inputData = [data[0], data[2], data[4], data[5], data[6], data[7], data[8], data[9]]

class _fase(object):
    """Class to implement a null phase"""
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


class MEoS(_fase):
    """General class for implement multiparameter equation of state
    Each child class must define parameters for do calculations:
        name: Name of component
        CASNumber: CAS Number of component
        formula: Empiric formula
        synonym: Alternate formula (Refrigerant name or sigles

        rhoc: Critical density, instance of unidades.Density
        Tc: Critical temperature, instance of unidades.Temperature
        Pc: Critical pressure, instance of unidades.Pressure
        M: Molecular weigth, g/mol
        Tt: Temperature of triple point, instance of unidades.Temperature
        Tb: Normal boiling point temperature, instance of unidades.Temperature
        f_acent: Acentric factor
        momentoDipolar: Depole moment, instance of unidades.DipoleMoment
        id: index of component in database if exist

        _Tr: Temperature parameter for generalized equation
        _rhor: Density parameter for generalized equation
        _w: Acentric factor for generalized equation

        eq: array of pointer for equations parameters
        _PR: Peneloux volume correction for Peng-Robinson equation of state

        _dielectric: Data for dielectric constant calculation
        _melting: Data for melting line calculation
        _sublimation: Data for sublimation line calculation
        _surface: Data for surface tension calculation
        _vapor_Pressure: Data for vapor pressure ancillary equation
        _liquid_Density: Data for liquid density ancillary equation
        _vapor_Density: Data for vapor density ancillary equation

        _viscosity: Tuple with viscosity equation
        _thermal: Tuple with thermal conductivity equation
    """
    id = None
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
              "T0": 0, 
              "recursion": True}
    status = 0
    msg = QApplication.translate("pychemqt", "Unknown Variables")
    __doi__ = {"surface":
        {"autor": "Mulero, A., Cachadiña, I., and Parra, M.I.",
         "title": "Recommended Correlations for the Surface Tension of Common Fluids", 
         "ref": "J. Phys. Chem. Ref. Data 41, 043105 (2012)",
         "doi": "10.1063/1.4768782"}, 
        }

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
        ks        -   Adiabatic Compressibility, 1/MPa
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
            self._surface["__doi__"] = self.__doi__["surface"]

    def __call__(self, **kwargs):
        self.cleanOldValues(**kwargs)
        
        self._constants = self.eq[self.kwargs["eq"]]
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

        if self.calculable:
            self.calculo()
            if self.status in (1, 3):
                converge = True
                for input in self._mode.split("-"):
                    if abs(self.kwargs[input]-self.__getattribute__(input)._data) > 1e-9:
                        converge = False
                        break
                if not converge:
                    self.status = 5
                    self.msg = QApplication.translate("pychemqt", "Solution don´t converge", None, QApplication.UnicodeUTF8)
                    print "dont converge for %s by %g" %(input, self.kwargs[input]-self.__getattribute__(input)._data)            
            
    def cleanOldValues(self, **kwargs):
        """Convert alternative rho input to correct rho value"""
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
        T = self.kwargs["T"]
        rho = self.kwargs["rho"]
        P = self.kwargs["P"]
        s = self.kwargs["s"]
        h = self.kwargs["h"]
        u = self.kwargs["u"]
        x = self.kwargs["x"]
        eq = self.kwargs["eq"]
        visco = self.kwargs["visco"]
        thermal = self.kwargs["thermal"]
        ref = self.kwargs["ref"]
        refvalues = self.kwargs["refvalues"]
        
        self._ref(ref, refvalues)

        if self.id:
            self.componente = compuestos.Componente(self.id)

        # Opcion de aceptar el nombre interno de la ecuacion
        if isinstance(eq, str) and eq in self.__class__.__dict__:
            eq = self.eq.index(self.__class__.__dict__[eq])

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

        if self._viscosity:
            self._viscosity = self._viscosity[visco]
        if self._thermal:
            self._thermal = self._thermal[thermal]

        propiedades = None

        if x is None:
            # Method with iteration necessary to get x
            if self._mode == "T-P":
                
                if self.kwargs["rho0"]:
                    rhoo = self.kwargs["rho0"]
                elif T < 0.99*self.Tc and \
                        self._Vapor_Pressure(T) < P:
                    rhoo = self._Liquid_Density(T)
                elif T < 0.99*self.Tc:
                    rhoo = self._Vapor_Density(T)
                elif T > 2*self.Tc or P > 2*self.Pc:
                    rhoo = self.eq[eq]["rhomax"]*self.M
                elif 0.99*self.Tc <= T < self.Tc and self.Pc*0.9 < P < self.Pc:
                    rhoo = self.rhoc
                else:
                    rhoo = P/T/self.R
                rinput = fsolve(lambda rho: self._eq(rho, T)["P"]-P, rhoo, full_output=True)
            
                if rinput[2] != 1:
                    self.status = 0
                    return
                rho = rinput[0][0]

            elif self._mode == "T-h":
                funcion = lambda rho: self._eq(rho, T)["h"]*1000-h
                def funcion2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*1000.*x+liquido["h"]*1000.*(1-x)-h

                rho = self.fsolve(funcion, True, funcion2, **{"h": h, "T": T})


            elif self._mode == "T-s":
                funcion = lambda rho: self._eq(rho, T)["s"]*1000.-s
                def funcion2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["s"]*1000.*x+liquido["s"]*1000.*(1-x)-s

                rho = self.fsolve(funcion, True, funcion2, **{"s": s, "T": T})

            elif self._mode == "T-u":
                def funcion(rho):
                    par = self._eq(rho, T)
                    return par["h"]-par["P"]/1000.*par["v"]-u/1000.
                def funcion2(rho):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vu*x+lu*(1-x)-u

                rho = self.fsolve(funcion, True, funcion2, **{"u": u, "T": T})

            elif self._mode == "P-rho":
                funcion = lambda T: self._eq(rho, T)["P"]-P
                def funcion2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    return Ps-P

                T = self.fsolve(funcion, True, funcion2, **{"P": P, "rho": rho})

            elif self._mode == "P-h":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["P"]-P, par["h"]*1000-h
                def funcion2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return Ps-P, vapor["h"]*1000*x+liquido["h"]*1000*(1-x)-h

                rho, T = self.fsolve(funcion, True, funcion2, **{"P": P, "h": h})

            elif self._mode == "P-s":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["P"]*1000-P, par["s"]*1000-s
                def funcion2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1./rhol)/(1./rhov-1./rhol)
                    return Ps-P, vapor["s"]*1000*x+liquido["s"]*1000*(1-x)-s

                rho, T = self.fsolve(funcion, True, funcion2, **{"P": P, "s": s})

            elif self._mode == "P-u":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]*par["v"]-u/1000., par["P"]-P
                def funcion(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return Ps-P, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(funcion, True, funcion2, **{"P": P, "u": u})

            elif self._mode == "rho-h":
                funcion = lambda T: self._eq(rho, T)["h"]*1000-h
                def funcion2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*x+liquido["h"]*(1-x)-h/1000

                T = self.fsolve(funcion, True, funcion2, **{"h": h, "rho": rho})

            elif self._mode == "rho-s":
                funcion = lambda T: self._eq(rho, T)["s"]*1000.-s
                def funcion2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["s"]*1000.*x+liquido["s"]*1000.*(1-x)-s
                
                T = self.fsolve(funcion, True, funcion2, **{"s": s, "rho": rho})
                
            elif self._mode == "rho-u":
                def funcion(T):
                    par = self._eq(rho, T)
                    return par["h"]-par["P"]/1000*par["v"]-u/1000
                def funcion2(T):
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vu*x+lu*(1-x)-u/1000

                T = self.fsolve(funcion, True, funcion2, **{"u": u, "rho": rho})

            elif self._mode == "h-s":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]*1000-h, par["s"]*1000.-s
                def funcion2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return (vapor["h"]*1000*x+liquido["h"]*1000*(1-x)-h, 
                            vapor["s"]*1000*x+liquido["s"]*1000*(1-x)-s)
                            
                rho, T = self.fsolve(funcion, True, funcion2, **{"s": s, "h": h})

            elif self._mode == "h-u":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]/1000*par["v"]-u, par["h"]-h
                def funcion(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["h"]*x+liquido["h"]*(1-x)-h, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(funcion, True, funcion2, **{"u": u, "h": h})

            elif self._mode == "s-u":
                def funcion(parr):
                    par = self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]*par["v"]-u, par["s"]*1000.-s
                def funcion2(parr):
                    rho, T = parr
                    rhol, rhov, Ps = self._saturation(T)
                    vapor = self._eq(rhov, T)
                    liquido = self._eq(rhol, T)
                    vu = vapor["h"]-Ps/rhov
                    lu = liquido["h"]-Ps/rhol
                    x = (1./rho-1/rhol)/(1/rhov-1/rhol)
                    return vapor["s"]*1000.*x+liquido["s"]*1000.*(1-x)-s, vu*x+lu*(1-x)-u

                rho, T = self.fsolve(funcion, True, funcion2, **{"u": u, "s": s})

            if self._mode == "T-rho" and self.kwargs["rho"] == 0:
                self.status = 3
                self.msg = QApplication.translate("pychemqt", "Ideal condition at zero pressure")
            elif self._constants["Tmin"]<=T<=self._constants["Tmax"] and \
                    0 < rho:# <= self._constants["rhomax"]*self.M:
                self.status = 1
                self.msg = ""
            else:
                self.status = 5
                self.msg = QApplication.translate("pychemqt", "input out of range")
                return 

            rho = float(rho)
            T = float(T)
            propiedades = self._eq(rho, T)
            if T <= self.Tc:
                rhol = self._Liquid_Density(T)
                rhov = self._Vapor_Density(T)
                if rhol > rho > rhov:
                    rhol, rhov, Ps = self._saturation(T)
                    x = (1/rho-1/rhol)/(1/rhov-1/rhol)
                    if x < 0:
                        x = 0
                    elif x > 1:
                        x = 1
                    P = Ps
                elif rho <= rhov:
                    x = 1
                elif rho >= rhol:
                    x = 0

                vapor = self._eq(rhov, T)
                liquido = self._eq(rhol, T)

            elif T > self.Tc:
                x = 1
            else:
                raise NotImplementedError("Incoming out of bound")

            if not P:
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
            T = fsolve(funcion, 0.9*self.Tc)[0]
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
        self.v0 = unidades.SpecificVolume(cp0.v)
        self.rho0 = unidades.Density(1./self.v0)
        self.rhoM0 = unidades.MolarDensity(self.rho0/self.M)
        self.h0 = unidades.Enthalpy(cp0.h)
        self.hM0 = unidades.MolarEnthalpy(self.h0/self.M)
        self.u0 = unidades.Enthalpy(self.h0-self.P*self.v0)
        self.uM0 = unidades.MolarEnthalpy(self.u0/self.M)
        self.s0 = unidades.SpecificHeat(cp0.s)
        self.sM0 = unidades.MolarSpecificHeat(self.s0/self.M)
        self.a0 = unidades.Enthalpy(self.u0-self.T*self.s0)
        self.aM0 = unidades.MolarEnthalpy(self.a0/self.M)
        self.g0 = unidades.Enthalpy(self.h0-self.T*self.s0)
        self.gM0 = unidades.MolarEnthalpy(self.g0/self.M)
        self.cp0 = unidades.SpecificHeat(cp0.cp)
        self.cpM0 = unidades.MolarSpecificHeat(self.cp0/self.M)
        self.cv0 = unidades.SpecificHeat(cp0.cv)
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv0)
        if self.rho0:
            self.gamma0 = unidades.Dimensionless(-self.v0/self.P/1000*self.derivative("P", "v", "s", cp0))
        else:
            self.gamma0 = 0

        self.Liquido = _fase()
        self.Gas = _fase()
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

            self.Z = unidades.Dimensionless(x*self.Gas.Z+(1-x)*self.Liquido.Z)
            self.f = unidades.Pressure(x*self.Gas.f+(1-x)*self.Liquido.f)
            self.Z_rho = unidades.SpecificVolume(x*self.Gas.Z_rho+(1-x)*self.Liquido.Z_rho)
            self.IntP = unidades.Pressure(x*self.Gas.IntP+(1-x)*self.Liquido.IntP)

        # Calculate special properties useful only for one phase
        if x < 1 and self.Tt <= T <= self.Tc:
            self.sigma = unidades.Tension(self._Surface())
        else:
            self.sigma = unidades.Tension(None)

        if 0 < x < 1:
            self.virialB = unidades.SpecificVolume(vapor["B"]/self.rhoc)
            self.virialC = unidades.SpecificVolume_square(vapor["C"]/self.rhoc**2)
        else:
            self.virialB = unidades.SpecificVolume(propiedades["B"]/self.rhoc)
            self.virialC = unidades.SpecificVolume_square(propiedades["C"]/self.rhoc**2)
            
        if self.Tt <= T <= self.Tc:
            self.Hvap = unidades.Enthalpy(vapor["h"]-liquido["h"], "kJkg")
        else:
            self.Hvap = unidades.Enthalpy(None)
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
            ro = [322, rhov, self.rhoc, self._constants["rhomax"]*self.M, 1, 1e-3]
            if self.kwargs["rho0"]:
                ro.insert(0, self.kwargs["rho0"])
        
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
                    f1 = function(rho)
                    if rho != r and 0 < rho < self._constants["rhomax"]*self.M and abs(f1) < 1e-6:
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
                    if (rho != r or T != t) and 0 < rho < self._constants["rhomax"]*self.M and abs(f1) < 1e-6 and abs(f2) < 1e-6:
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
                                rinput = fsolve(function2phase, r, full_output=True)
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
                                rinput = fsolve(function2phase, t, full_output=True)
                                T = rinput[0][0]
                            except:
                                pass
                            else:
                                f1 = function2phase(T)
                                if T != t and 0 < rho < self._constants["rhomax"]*self.M and abs(f1) < 1e-3:
                                    break
                    else:
                        for r, t in zip(ro, to):
                            try:
                                rinput = fsolve(function2phase, [r, t], full_output=True)
                                rho, T = rinput[0]
                            except:
                                pass
                            else:
                                f1, f2 = function2phase([rho, T])
                                print rho, T, f1, f2
                                if (rho != r or T != t) and 0 < rho < self._constants["rhomax"]*self.M and abs(f1) < 1e-3 and abs(f2) < 1e-3:
                                    break
                                    
        if "T" in kwargs:
            return rho
        elif "rho" in kwargs:
            return T
        else:
            return rho, T

    def fillNone(self, fase):
        """Fill properties in null phase with a explicative msg"""
        if self.x == 0:
            txt = QApplication.translate("pychemqt", "Subcooled")
        elif self.Tr < 1 and self.Pr < 1: 
            txt = QApplication.translate("pychemqt", "Superheated")
        elif self.Tr == 1 and self.Pr == 1: 
            txt = QApplication.translate("pychemqt", "Critic point")
        else:
            txt = QApplication.translate("pychemqt", "Supercritical")
        for key in _fase.__dict__:
            if key[0] != "_":
                fase.__setattr__(key, txt)
            
    def fill(self, fase, estado):
        """Fill phase properties"""
        fase.M = unidades.Dimensionless(self.M)
        fase.v = unidades.SpecificVolume(estado["v"])
        fase.rho = unidades.Density(1/fase.v)

        fase.h = unidades.Enthalpy(estado["h"], "kJkg")
        fase.s = unidades.SpecificHeat(estado["s"], "kJkgK")
        fase.u = unidades.Enthalpy(fase.h-self.P*fase.v)
        fase.a = unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g = unidades.Enthalpy(fase.h-self.T*fase.s)

        fase.Z = unidades.Dimensionless(self.P*fase.v/self.T/self.R)
        fase.fi = unidades.Dimensionless(estado["fugacity"])
        fase.f = unidades.Pressure(fase.fi*self.P)
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

        fase.joule = unidades.TemperaturePressure(self.derivative("T", "P", "h", fase))
        fase.Gruneisen = unidades.Dimensionless(fase.v/fase.cv*self.derivative("P", "T", "v", fase))
        
        if fase.rho:
            fase.alfav = unidades.InvTemperature(self.derivative("v", "T", "P", fase)/fase.v)
            fase.kappa = unidades.InvPressure(-self.derivative("v", "P", "T", fase)/fase.v)
            fase.betas = unidades.TemperaturePressure(self.derivative("T", "P", "s", fase))
            fase.gamma = unidades.Dimensionless(-fase.v/self.P*self.derivative("P", "v", "s", fase))

            fase.kt = unidades.Dimensionless(-fase.v/self.P*self.derivative("P", "v", "T", fase))
            fase.ks = unidades.InvPressure(-self.derivative("v", "P", "s", fase)/fase.v)
            fase.Kt = unidades.Pressure(-fase.v*self.derivative("P", "v", "s", fase))
            fase.Ks = unidades.Pressure(-fase.v*self.derivative("P", "v", "T", fase))
            fase.dhdT_rho = unidades.SpecificHeat(self.derivative("h", "T", "rho", fase))
            fase.dhdT_P = unidades.SpecificHeat(self.derivative("h", "T", "P", fase))
            fase.dhdP_T = unidades.EnthalpyPressure(self.derivative("h", "P", "T", fase)) #deltat
            fase.dhdP_rho = unidades.EnthalpyPressure(self.derivative("h", "P", "rho", fase))
            fase.dhdrho_T = unidades.EnthalpyDensity(estado["dhdrho"])
            fase.dhdrho_P = unidades.EnthalpyDensity(estado["dhdrho"]+fase.dhdT_rho/estado["drhodt"])
            fase.dpdT_rho = unidades.PressureTemperature(self.derivative("P", "T", "rho", fase))
            fase.dpdrho_T = unidades.PressureDensity(estado["dpdrho"])
            fase.drhodP_T = unidades.DensityPressure(1/estado["dpdrho"])
            fase.drhodT_P = unidades.DensityTemperature(estado["drhodt"])

            fase.Z_rho = unidades.SpecificVolume((fase.Z-1)/fase.rho)
            fase.IntP = unidades.Pressure(self.T*self.derivative("P", "T", "rho", fase)-self.P)
            fase.hInput = unidades.Enthalpy(fase.v*self.derivative("h", "v", "P", fase))
        if self.kwargs["recursion"]:
            fase.mu = self._Viscosity(fase.rho, self.T, fase)
            fase.k = self._ThCond(fase.rho, self.T, fase)
            if fase.mu and fase.rho:
                fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
            else:
                fase.nu = None
            if fase.k and fase.rho:
                fase.alfa = unidades.Diffusivity(fase.k/1000/fase.rho/fase.cp)
            else:
                fase.alfa = None
            if fase.mu and fase.k:
                fase.Prandt = unidades.Dimensionless(fase.mu*fase.cp*1000/fase.k)
            else:
                fase.Prandt = None
            fase.epsilon = unidades.Dimensionless(self._Dielectric(fase.rho, self.T))

    def _saturation(self, T=None):
        """Saturation calculation for two phase search"""
        if not T:
            T = self.T
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
            Ps = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(liquido["fir"]-vapor["fir"]+log(deltaL/deltaG))
        return rhoL, rhoG, Ps

    def _Helmholtz(self, rho, T):
        """Implementación general de la ecuación de estado Setzmann-Wagner, ecuación de estado de multiparámetros basada en la energía libre de Helmholtz"""
        delta = rho/self.rhoc
        tau = self.Tc/T
        
        fio, fiot, fiott, fiod, fiodd, fiodt = self._phi0(self._constants["cp"], tau, delta)
        fir, firt, firtt, fird, firdd, firdt, firdtt, B, C=self._phir(tau, delta)
        
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
        propiedades["cp"] = self.R.kJkgK*(-tau**2*(fiott+firtt) + 
            (1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["w"] = abs(self.R*T*(1+2*delta*fird+delta**2*firdd - 
            (1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"] = (1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"] = rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
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

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta)

        ref = self._constants["ref"](eq=self._constants["eq"])
        Tr = T/Tc
        rhor = rho/rhoc

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
        fir, firt, firtt, fird, firdd, firdt, firdtt, B, C=ref._phir(tauref, deltaref)

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
        """Multimaparameter euation of state of Benedict-Webb-Rubin"""
        rho = rho/self.M
        delta = rho/self.rhoc
        tau = self.Tc/T
        b = self._constants["b"]

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta)

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
        firddd = (d2Prho*rho-4*dPdrho+6*P/rho)/self.R/T-2.

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

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta)

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
            refeq = self._constants["ref"]
            if isinstance(refeq, str):
                ref = refeq
            elif isinstance(refeq, dict):
                ref = "CUSTOM"
                refvalues = [refeq["Tref"], refeq["Pref"], refeq["ho"], refeq["so"]]
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
        fio, fiot, fiott, fiod, fiodd, fiodt = self._phi0(self._constants["cp"], tau, delta)

        propiedades = _fase()
        if rho:
            propiedades.v = self.R*T/self.P
        else:
            propiedades.v = float("inf")
        propiedades.h = self.R*T*(1+tau*fiot)
        propiedades.s = self.R*(tau*fiot-fio)
        propiedades.cv = -self.R*tau**2*fiott
        propiedades.cp = self.R*(-tau**2*fiott+1)
        propiedades.alfap = 1/T
        propiedades.betap = rho
        return propiedades
    
    def _PHIO(self, cp):
        """Convert cp dict in phi0 dict"""
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
            cI-=c*t*tau0**(t-1)
            cII+=c*(t-1)*tau0**t
        for ao, tita in zip(cp["ao_exp"], titao):
            cI-=ao*tita/(exp(tita*tau0)-1)
            cII-=ao*(tita*tau0/(1-exp(tita*tau0))+log(1-exp(-tita*tau0)))
        if cp["ao_hyp"]:
            for i in [0, 2]:
                cI-=cp["ao_hyp"][i]*hyp[i]/(tanh(hyp[i]*tau0))
                cII+=cp["ao_hyp"][i]*(hyp[i]*tau0/tanh(hyp[i]*tau0)-log(abs(sinh(hyp[i]*tau0))))
            for i in [1, 3]:
                cI+=cp["ao_hyp"][i]*hyp[i]*tanh(hyp[i]*tau0)
                cII-=cp["ao_hyp"][i]*(hyp[i]*tau0*tanh(hyp[i]*tau0)-log(abs(cosh(hyp[i]*tau0))))
        
        Fi0 = {"ao_log": [1,  co],
               "pow": [0, 1] + ti,
               "ao_pow": [cII, cI] + ci,
               "ao_exp": cp["ao_exp"],
               "titao": titao,
               "ao_hyp": cp["ao_hyp"],
               "hyp": hyp}
        return Fi0

    def _phi0(self, cp, tau, delta):
        if "ao_log" in cp:
            Fi0 = cp
        else:
            Fi0 = self._PHIO(cp)
#        print Fi0
#        
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
        fio=Fi0["ao_log"][1]*log(tau)
        fiot=Fi0["ao_log"][1]/tau
        fiott=-Fi0["ao_log"][1]/tau**2

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
        if "ao_exp2" in Fi0:
            for n, g, sum in zip(Fi0["ao_exp2"], Fi0["titao2"], Fi0["sum2"]):
                fio += n*log(sum+exp(g*tau))
                fiot += n*g/(sum*exp(-g*tau)+1)
                fiott += sum*n*g**2*exp(-g*tau)/(sum*exp(-g*tau)+1)**2

        if "ao_hyp" in Fi0 and Fi0["ao_hyp"]:
            for i in [0, 2]:
                fio += Fi0["ao_hyp"][i]*log(abs(sinh(Fi0["hyp"][i]*tau)))
                fiot += Fi0["ao_hyp"][i]*Fi0["hyp"][i]/tanh(Fi0["hyp"][i]*tau)
                fiott -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/sinh(Fi0["hyp"][i]*tau)**2

            for i in [1, 3]:
                fio -= Fi0["ao_hyp"][i]*log(abs(cosh(Fi0["hyp"][i]*tau)))
                fiot -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]*tanh(Fi0["hyp"][i]*tau)
                fiott -= Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/cosh(Fi0["hyp"][i]*tau)**2

        if delta:
            fiod = 1/delta
            fiodd = -1/delta**2
        else:
            fiod, fiodd = 0, 0
        fiodt = 0
        R_ = cp.get("R", self._constants["R"])
        factor = R_/self._constants["R"]
        if delta:
            fio = Fi0["ao_log"][0]*log(delta)+factor*fio
        else:
            fio *= factor
        return fio, factor*fiot, factor*fiott, fiod, fiodd, fiodt

    def _Cp0(self, T=False):
        Tc = self._constants.get("Tref", self.Tc)
        if not T:
            T = self.T
        cp = self._constants["cp"]
            
        if "ao_log" in cp:
            tau = Tc/T
            fio, fiot, fiott, fiod, fiodd, fiodt = self._phi0(cp, tau, 0)
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
        delta_0 = 1e-200
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
                firdt += n*t*delta**(d-1)*tau**(t-1)*(d-g*c*delta**c)*exp(-g*delta**c)
                firdtt += n*t*(t-1)*delta**(d-1)*tau**(t-2)*(d-g*c*delta**c) * \
                    exp(-g*delta**c)
                B += n*exp(-g*delta_0**c)*delta_0**(d-1)*tau**t*(d-g*c*delta_0**c)
                C += n*exp(-g*delta_0**c)*(delta_0**(d-2)*tau**t *
                    ((d-g*c*delta_0**c)*(d-1-g*c*delta_0**c)-g**2*c**2*delta_0**c))

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
            for n, d, t, a, e, b, g, ex1, ex2 in zip(nr3, d3, t3, a3, e3, b3, g3, exp1, exp2):
                fir += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2)
                fird += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2) * \
                    (d/delta-2*a*(delta-e))
                firdd += n*tau**t*exp(-a*(delta-e)**ex1-b *
                    (tau-g)**ex2)*(-2*a*delta**d+4*a**2 *
                    delta**d*(delta-e)**ex1-4*d*a*delta**(d-1) *
                    (delta-e)+d*(d-1)*delta**(d-2))
                firt += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2) * \
                    (t/tau-2*b*(tau-g))
                firtt += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2) * \
                    ((t/tau-2*b*(tau-g))**ex2-t/tau**2-2*b)
                firdt += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2) * \
                    (t/tau-2*b*(tau-g))*(d/delta-2*a*(delta-e))
                firdtt += n*delta**d*tau**t * \
                    exp(-a*(delta-e)**ex1-b*(tau-g)**ex2) * \
                    ((t/tau-2*b*(tau-g))**ex2-t/tau**2-2*b) * \
                    (d/delta-2*a*(delta-e))
                B += n*delta_0**d*tau**t * \
                    exp(-a*(delta_0-e)**ex1-b*(tau-g)**ex2) * \
                    (d/delta_0-2*a*(delta_0-e))
                C += n*tau**t * \
                    exp(-a*(delta_0-e)**ex1-b*(tau-g)**ex2) * \
                    (-2*a*delta_0**d+4*a**2*delta_0**d *
                    (delta_0-e)**ex1-4*d*a*delta_0**2*(delta_0-e) +
                    d*2*delta_0)

            # Non analitic terms
            nr4 = self._constants.get("nr4", [])
            a4 = self._constants.get("a4", [])
            b = self._constants.get("b4", [])
            A = self._constants.get("A", [])
            Bi = self._constants.get("B", [])
            Ci = self._constants.get("C", [])
            D = self._constants.get("D", [])
            bt = self._constants.get("beta4", [])
            for i in range(len(nr4)):
                Tita = (1-tau)+A[i]*((delta-1)**2)**(0.5/bt[i])
                F = exp(-Ci[i]*(delta-1)**2-D[i]*(tau-1)**2)
                Fd = -2*Ci[i]*F*(delta-1)
                Fdd = 2*Ci[i]*F*(2*Ci[i]*(delta-1)**2-1)
                Ft = -2*D[i]*F*(tau-1)
                Ftt = 2*D[i]*F*(2*D[i]*(tau-1)**2-1)
                Fdt = 4*Ci[i]*D[i]*F*(delta-1)*(tau-1)
                Fdtt = 4*Ci[i]*D[i]*F*(delta-1)*(2*D[i]*(tau-1)**2-1)

                Delta = Tita**2+Bi[i]*((delta-1)**2)**a4[i]
                Deltad = (delta-1)*(A[i]*Tita*2/bt[i]*((delta-1)**2)**(0.5/bt[i]-1) \
                    + 2*Bi[i]*a4[i]*((delta-1)**2)**(a4[i]-1))
                if delta == 1:
                    Deltadd = 0
                    DeltaBd = 0
                    DeltaBdd = 0
                    DeltaBt = 0
                    DeltaBtt = 0
                    DeltaBdt = 0
                    DeltaBdtt = 0
                else:
                    Deltadd = Deltad/(delta-1)+(delta-1)**2 * \
                        (4*Bi[i]*a4[i]*(a4[i]-1)*((delta-1)**2)**(a4[i]-2) +
                        2*A[i]**2/bt[i]**2*(((delta-1)**2)**(0.5/bt[i]-1))**2 +
                        A[i]*Tita*4/bt[i]*(0.5/bt[i]-1)*((delta-1)**2)**(0.5/bt[i]-2))
                    DeltaBd = b[i]*Delta**(b[i]-1)*Deltad
                    DeltaBdd = b[i]*(Delta**(b[i]-1)*Deltadd+(b[i]-1)*Delta**(b[i]-2)*Deltad**2)
                    DeltaBt = -2*Tita*b[i]*Delta**(b[i]-1)
                    DeltaBtt = 2*b[i]*Delta**(b[i]-1)+4*Tita**2*b[i]*(b[i]-1)*Delta**(b[i]-2)
                    DeltaBdt = -A[i]*b[i]*2/bt[i]*Delta**(b[i]-1)*(delta-1)*((delta-1)**2) ** \
                        (0.5/bt[i]-1)-2*Tita*b[i]*(b[i]-1)*Delta**(b[i]-2)*Deltad
                    DeltaBdtt = 2*b[i]*(b[i]-1)*Delta**(b[i]-2) * \
                        (Deltad*(1+2*Tita**2*(b[i]-2)/Delta)+4*Tita*A[i] *
                        (delta-1)/bt[i]*((delta-1)**2)**(0.5/bt[i]-1))

                fir += nr4[i]*Delta**b[i]*delta*F
                fird += nr4[i]*(Delta**b[i]*(F+delta*Fd)+DeltaBd*delta*F)
                firdd += nr4[i]*(Delta**b[i]*(2*Fd+delta*Fdd)+2*DeltaBd *
                    (F+delta*Fd)+DeltaBdd*delta*F)
                firt += nr4[i]*delta*(DeltaBt*F+Delta**b[i]*Ft)
                firtt += nr4[i]*delta*(DeltaBtt*F+2*DeltaBt*Ft+Delta**b[i]*Ftt)
                firdt += nr4[i]*(Delta**b[i]*(Ft+delta*Fdt)+delta*DeltaBd*Ft +
                    DeltaBt*(F+delta*Fd)+DeltaBdt*delta*F)
                firdtt += nr4[i]*((DeltaBtt*F+2*DeltaBt*Ft+Delta**b[i]*Ftt) +
                    delta*(DeltaBdtt*F+DeltaBtt*Fd+2*DeltaBdt*Ft+2*DeltaBt*Fdt +
                    DeltaBt*Ftt+Delta**b[i]*Fdtt))

                Tita_ = (1-tau)+A[i]*((delta_0-1)**2)**(0.5/bt[i])
                Delta_ = Tita_**2+Bi[i]*((delta_0-1)**2)**a4[i]
                Deltad_ = (delta_0-1)*(A[i]*Tita_*2/bt[i]*((delta_0-1)**2) **
                    (0.5/bt[i]-1)+2*Bi[i]*a4[i]*((delta_0-1)**2)**(a4[i]-1))
                Deltadd_ = Deltad_/(delta_0-1)+(delta_0-1)**2 * \
                    (4*Bi[i]*a4[i]*(a4[i]-1)*((delta_0-1)**2)**(a4[i]-2)+2*A[i]**2 /
                    bt[i]**2*(((delta_0-1)**2)**(0.5/bt[i]-1))**2+A[i]*Tita_ *
                    4/bt[i]*(0.5/bt[i]-1)*((delta_0-1)**2)**(0.5/bt[i]-2))
                DeltaBd_ = b[i]*Delta_**(b[i]-1)*Deltad_
                DeltaBdd_ = b[i]*(Delta_**(b[i]-1)*Deltadd_ +
                    (b[i]-1)*Delta_**(b[i]-2)*Deltad_**2)
                F_ = exp(-Ci[i]*(delta_0-1)**2-D[i]*(tau-1)**2)
                Fd_ = -2*Ci[i]*F_*(delta_0-1)
                Fdd_ = 2*Ci[i]*F_*(2*Ci[i]*(delta_0-1)**2-1)

                B += nr4[i]*(Delta_**b[i]*(F_+delta_0*Fd_) +
                    DeltaBd_*delta_0*F_)
                C += nr4[i]*(Delta_**b[i]*(2*Fd_+delta_0*Fdd_) +
                    2*DeltaBd_*(F_+delta_0*Fd_) +
                    DeltaBdd_*delta_0*F_)

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
                Xtt = -n*delta*((1-a)*g/tau**(g+2)*((g+1)*(a+(1-a)/tau**g)-2*g*(1-a)/tau**g))/(a+(1-a)/tau**g)**3
                Xdtt = -n*((1-a)*g/tau**(g+2)*((g+1)*(a+(1-a)/tau**g)-2*g*(1-a)/tau**g))/(a+(1-a)/tau**g)**3

                ahdX = -(f**2-1)/(1-X)+(f**2+3*f+X*(f**2-3*f))/(1-X)**3
                ahdXX = -(f**2-1)/(1-X)**2+(3*(f**2+3*f)+(f**2-3*f)*(1+2*X))/(1-X)**4
                ahdXXX = -2*(f**2-1)/(1-X)**3+6*(2*(f**2+3*f)+(f**2-3*f)*(1+X))/(1-X)**5

                fir += (f**2-1)*log(1-X)+((f**2+3*f)*X-3*f*X**2)/(1-X)**2
                fird += ahdX*Xd
                firdd += ahdXX*Xd**2
                firt += ahdX*Xt
                firtt += ahdXX*Xt**2+ahdX*Xtt
                firdt += ahdXX*Xt*Xd+ahdX*Xdt
                firdtt += ahdXXX*Xt**2*Xd+ahdXX*(Xtt*Xd+2*Xdt*Xt)*ahdX*Xdtt

                X_virial = n*delta_0/(a+(1-a)/tau**g)
                ahdX_virial = -(f**2-1)/(1-X_virial)+(f**2+3*f+X_virial*(f**2-3*f))/(1-X_virial)**3
                ahdXX_virial = -(f**2-1)/(1-X_virial)**2+(3*(f**2+3*f)+(f**2-3*f)*(1+2*X_virial))/(1-X_virial)**4
                B += ahdX_virial*Xd
                C += ahdXX_virial*Xd**2
            
            # Special form from Saul, A. and Wagner, W. Water 58 coefficient equation
            if "nr5" in self._constants:
                delta_0 = 1e-200
                if delta < 0.2:
                    factor = 1.6*delta**6*(1-1.2*delta**6)
                else:
                    factor = exp(0.4*delta**6)-exp(-2*delta**6)

                nr5 = self._constants.get("nr5", [])
                d5 = self._constants.get("d5", [])
                t5 = self._constants.get("t5", [])
                fr, frt, frtt, frdtt1, frdtt2 = 0, 0, 0, 0, 0
                frd1, frd2 = 0, 0
                frdd1, frdd2, frdd3 = 0, 0, 0
                frdt1, frdt2, frdt3 = 0, 0, 0
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
                fird += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6))*frd1 +\
                    factor*frd2
                firdd += (5.76*exp(-0.4*delta**6)-144*exp(-2*delta**6))*frdd1 +\
                    (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6))*frdd2 +\
                    factor*frdd3
                firt += factor*frt
                firtt += factor*frtt
                firdt += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6))*frdt1 +\
                    factor*frdt2
                firdtt += (-2.4*exp(-0.4*delta**6)+12*exp(-2*delta**6))*frdtt1 +\
                    factor*frdtt2
                B += (-2.4*exp(-0.4*delta_0**6)+12*exp(-2*delta_0**6))*Bsum1 +\
                    (exp(0.4*delta_0**6)-exp(-2*delta_0**6))*Bsum2
                C += (5.76*exp(-0.4*delta_0**6)-144*exp(-2*delta_0**6))*Csum1 +\
                    (-2.4*exp(-0.4*delta_0**6)+12*exp(-2*delta_0**6))*Csum2 +\
                    (exp(0.4*delta_0**6)-exp(-2*delta_0**6))*Csum3
                    
        return fir, firt, firtt, fird, firdd, firdt, firdtt, B, C


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

    def _Dielectric(self, rho, T):
        """Harvey, A.H. and Lemmon, E.W. "Method for Estimating the Dielectric Constant of Natural Gas Mixtures," Int. J. Thermophys., 26(1):31-46, 2005.
        http://dx.doi.org/10.1007/s10765-005-2351-5"""
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
            e = 0
        return unidades.Dimensionless(e)

    @classmethod
    def _Melting_Pressure(cls, T):
        if cls._melting:
            Tref = cls._melting["Tref"]
            Pref = cls._melting["Pref"]
            Tita = T/Tref
            suma = 0
            for ai, expi in zip(cls._melting["a1"], cls._melting["exp1"]):
                suma += ai*Tita**expi
            for ai, expi in zip(cls._melting["a2"], cls._melting["exp2"]):
                suma += ai*(Tita-1)**expi
            for ai, expi in zip(cls._melting["a3"], cls._melting["exp3"]):
                suma += ai*log(Tita)**expi

            if cls._melting["eq"] == 1:
                P = suma*Pref
            elif cls._melting["eq"] == 2:
                P = exp(suma)*Pref
            return unidades.Pressure(P, "kPa")
        else:
            return None

    @classmethod
    def _Sublimation_Pressure(cls, T):
        if cls._sublimation:
            Tref = cls._sublimation["Tref"]
            Pref = cls._sublimation["Pref"]
            Tita = T/Tref
            suma = 0
            for ai, expi in zip(cls._sublimation["a1"], cls._sublimation["exp1"]):
                suma += ai*Tita**expi
            for ai, expi in zip(cls._sublimation["a2"], cls._sublimation["exp2"]):
                suma += ai*(1-Tita)**expi
            for ai, expi in zip(cls._sublimation["a3"], cls._sublimation["exp3"]):
                suma += ai*log(Tita)**expi

            if cls._sublimation["eq"] == 1:
                P = suma*Pref
            elif cls._sublimation["eq"] == 2:
                P = exp(suma)*Pref
            elif cls._sublimation["eq"] == 3:
                P = exp(Tref/T*suma)*Pref
            return unidades.Pressure(P, "kPa")
        else:
            return None

    def _Vapor_Pressure(self, T):
        if self._vapor_Pressure:
            eq = self._vapor_Pressure["eq"]
            Tita = 1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita = Tita**0.5
            suma = sum([n*Tita**x for n, x in zip(
                self._vapor_Pressure["ao"], self._vapor_Pressure["exp"])])
            if eq in [1, 2]:
                Pr = suma+1
            elif eq in [3, 4]:
                Pr = exp(suma)
            else:
                Pr = exp(self.Tc/T*suma)
            Pv = unidades.Pressure(Pr*self.Pc)
        else:
            Pv = self.componente.Pv(T)
        return Pv

    def _Liquid_Density(self, T):
        if self._liquid_Density:
            eq = self._liquid_Density["eq"]
            Tita = 1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita = Tita**(1./3)
            suma = sum([n*Tita**x for n, x in zip(
                self._liquid_Density["ao"], self._liquid_Density["exp"])])
            if eq in [1, 2]:
                Pr = suma+1
            elif eq in [3, 4]:
                Pr = exp(suma)
            else:
                Pr = exp(self.Tc/T*suma)
            rho = unidades.Density(Pr*self.rhoc)
        else:
            rho = self.componente.RhoL_DIPPR(T)
        return rho

    def _Vapor_Density(self, T):
        if self._vapor_Density:
            eq = self._vapor_Density["eq"]
            Tita = 1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita = Tita**(1./3)
            suma = sum([n*Tita**x for n, x in zip(
                self._vapor_Density["ao"], self._vapor_Density["exp"])])
            if eq in [1, 2]:
                Pr = suma+1
            elif eq in [3, 4]:
                Pr = exp(suma)
            else:
                Pr = exp(self.Tc/T*suma)
            rho = unidades.Density(Pr*self.rhoc)
        else:
            rho = self._Vapor_Density_Chouaieb(T)
        return rho

    def _Vapor_Density_Chouaieb(self, T):
        """O. Chouaieb, J. Ghazouani, A. Bellagi, Simple correlations for saturated liquid and vapor densities of pure fluids, Thermochimica Acta, Volume 424, Issues 1–2, 15 December 2004, Pages 43-51, ISSN 0040-6031,
        http://dx.doi.org/10.1016/j.tca.2004.05.017"""
        coeff = {"Ne": (4.4960428, 1.1861523, -0.485720905),
                 "Ar": (3.1601166, 1.1183779, -0.5515808),
                 "Kr": (3.4377347, 1.1459094, -0.5312377),
                 "Xe": (.2336283, 1.1272991, -0.5475391),
                 "CH4": (2.9407241, 1.1012273, -0.5697413),
                 "O2": (2.9401387, 1.1115040, -0.5649451),
                 "N2": (2.8539512, 1.1131452, -0.5673304),
                 "F2": (2.5237557, 1.0739354, -0.6052029),
                 "CO": (2.1804625, 1.0460487, -0.6091577),
                 "Ethylene": (2.6032071, 1.1016984, -0.5971817),
                 "C2": (2.6379996, 1.1048142, -0.6018817),
                 "NF3": (2.7962186, 1.1483585, -0.5712755),
                 "Propylene": (3.1010584, 1.1872850, -0.5517088),
                 "C3": (2.5491601, 1.1169908, -0.6110244),
                 "iC4": (2.5860317, 1.1380881, -0.6069818),
                 "nC4": (2.4703930, 1.1187258, -0.629651003),
                 "R22": (2.4474190, 1.1245385, -0.6346948),
                 "CO2": (2.4686277, 1.1345838, -0.6240188),
                 "nC5": (2.4715148, 1.1326580, -0.6412635),
                 "NH3": (2.9025748, 1.1747326, -0.6213074),
                 "R143a": (2.4390945, 1.1210065, -0.6625750),
                 "R152a": (2.4036109, 1.1186849, -0.6709599),
                 "R32": (2.4973193, 1.1207203, -0.6786065),
                 "R123": (2.3683423, 1.1374290, -0.6444662),
                 "R124": (2.4354428, 1.1506277, -0.6362176),
                 "nC6": (2.5036259, 1.1549903, -0.6410813),
                 "R125": (2.3507937, 1.1406387, -0.6493031),
                 "R134a": (2.4203360, 1.1516298, -0.6546052),
                 "H2O": (2.3609558, 1.0916682, -0.7452828),
                 "nC7": (2.5477720, 1.1770021, -0.6419602)}
        tau = 1-T/self.Tc

        if self.__class__.__name__ in coeff:
            m, n, p = coeff[self.__class__.__name__]
        else:
            Zc = self.rhoc*self.R*self.Tc/self.Pc
            Ni = [-0.1497547, 0.6006565]
            Pi = [-19.348354, -41.060325, 1.1878726]
            m = 3-1.71*self.f_acent  # Simple custom regresion of Fig.8
            p = Zc**2/(Pi[0]+Pi[1]*log10(Zc)*Pi[2])
            n = p+1./(Ni[0]*self.f_acent+Ni[1])

        alfa = exp(tau**(1./3)+tau**0.5+tau+tau**m)
        rhog = self.rhoc*exp(p*(alfa**n-exp(1-alfa)))
        return rhog

    def _Surface(self, T=None):
        """Equation for the surface tension"""
        if self.Tt <= self.T <= self.Tc:
            if self._surface:
                if not T:
                    T = self.T
                tau = 1-T/self.Tc
                tension = 0
                for sigma, n in zip(self._surface["sigma"],
                                    self._surface["exp"]):
                    tension += sigma*tau**n
                sigma = unidades.Tension(tension)
            else:
                sigma = self.componente.Tension(self.T)
        else:
            sigma = None
        return sigma

    def _Omega(self):
        """Collision integral calculations
            0 - None
            1 - Lemmon: nitrogen. oxygen, argon, aire  and custom collison parameter, co2
            2 - Younglove: C1, C2, C3, iC4, nC4
            3 - CI0 from NIST:  Chung
        Ref:
            Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.
            Younglove, B.A. and Ely, J.F. (1987). Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane. J. Phys. Chem. Ref. Data  16: 577-798.
            Fenghour, A., Wakeham, W.A., Vesovic, V., Watson, J.T.R., Millat, J., and Vogel, E., "The viscosity of ammonia," J. Phys. Chem. Ref. Data, 24:1649-1667, 1995.
            Fenghour, A., Wakeham, W.A., Vesovic, V., "The Viscosity of Carbon Dioxide," J. Phys. Chem. Ref. Data, 27:31-44, 1998.
            T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679,
        """
        if self._viscosity["omega"] == 1:
            if self._viscosity.has_key("collision"):
                b = self._viscosity["collision"]
            else:
                b = [0.431, -0.4623, 0.08406, 0.005341, -0.00331]
            T_ = log(self.T/self._viscosity["ek"])
            suma = 0
            for i, bi in enumerate(b):
                suma += bi*T_**i
            omega = exp(suma)
        elif self._viscosity["omega"] == 2:
            if "collision" in self._viscosity:
                b = self._viscosity["collision"]
            else:
                b = [-3.0328138281, 16.918880086, -37.189364917, 41.288861858,
                     -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390,
                     -0.009657043707]
            T_ = self._viscosity["ek"]/self.T
            suma = 0
            for i, bi in enumerate(b):
                suma += bi*T_**((3.-i)/3.)
            omega = 1./suma
        elif self._viscosity["omega"] == 3:
            # FIXME: reference from https://github.com/thorade/HelmholtzMedia/blob/master/Interfaces/PartialHelmholtzMedium/Transport/dynamicViscosity_dilute.mo#L27
            T_ = self.T/self._viscosity.get("ek", self.Tc/1.2593)
            omega = 1.16145/T_**0.14874 + 0.52487/exp(0.77320*T_) + \
                2.16178/exp(2.4378*T_) - \
                6.435e-4*T_**0.14874*sin(18.0323*T_**-0.76830-7.27371)
        return omega

    def _Visco0(self):
        if self._viscosity["eq"] == 3:
            Tc = self._viscosity.get("Tref", 1.)
            rhoc = self._viscosity.get("rhoref", 1.)
            tau = self.T/Tc
            delta = self.rho/self.M/rhoc
            muo = 0
            for n, t in zip(self._viscosity["n_poly"], self._viscosity["t_poly"]):
                muo += n*tau**t
            if "n_polyden" in self._viscosity:
                den = 1
                for n, t, d in zip(self._viscosity["n_polyden"], self._viscosity["t_polyden"], self._viscosity["d_polyden"]):
                    den *= n*tau**t*delta**d
                    muo /= den
        elif self._viscosity["eq"] == 5:
            omega = self._Omega()
            Fc = 1-0.2756*self._viscosity["w"] + \
                0.059035*self._viscosity["mur"]**4+self._viscosity["k"]
            muo = 4.0795e-5*(self.M*self.T)**0.5*self.rhoc**(2./3.)/omega*Fc
        else:
            omega = self._Omega()
            Nchapman = self._viscosity.get("n_chapman", 0.0266958)
            tchapman = self._viscosity.get("t_chapman", 0.5)
            muo = Nchapman*(self.M*self.T)**tchapman/(self._viscosity["sigma"]**2*omega)

            # other adittional empirical terms
            if "n_ideal" in self._viscosity:
                tau = self.T/self._viscosity["Tref"]
                for n, t in zip(self._viscosity["n_ideal"], self._viscosity["t_ideal"]):
                    muo += n*tau**t
        return muo

    def _Viscosity(self, rho, T, fase):
        if self._viscosity:
            if self._viscosity["eq"] == 0:
                # Hardcoded method
                mu = self.__getattribute__(self._viscosity["method"])(rho, T, fase)*1e6

            elif self._viscosity["eq"] == 1:
                muo = self._Visco0()

                if rho > 0:
                    # second virial
                    Tc = self._viscosity.get("Tref_virial", self.Tc)
                    etar = self._viscosity.get("etaref_virial", 1.)
                    tau = T/Tc
                    mud = 0
                    if self._viscosity.has_key("n_virial"):
                        muB = 0
                        for n, t in zip(self._viscosity["n_virial"],
                                        self._viscosity["t_virial"]):
                            muB += n*tau**t
                        mud = etar*muB*rho/self.M*muo

                    Tc = self._viscosity.get("Tref_res", self.Tc)
                    rhoc = self._viscosity.get("rhoref_res", self.rhoc)
                    mured = self._viscosity.get("etaref_res", 1.)
                    tau = Tc/T
                    delta = rho/rhoc
                    if abs(delta-1) <= 0.001:
                        expdel = rho/self.rhoc
                    else:
                        expdel = delta

                    mur = 0
                    # close-packed density;
                    if "n_packed" in self._viscosity:
                        del0 = 0
                        for n, t in zip(self._viscosity["n_packed"],
                                        self._viscosity["t_packed"]):
                            del0 += n*tau**t
                    else:
                        del0 = 1.

                    # polynomial term
                    if "n_poly" in self._viscosity:
                        for n, t, d, c, g in zip(
                            self._viscosity["n_poly"], self._viscosity["t_poly"],
                            self._viscosity["d_poly"], self._viscosity["c_poly"],
                            self._viscosity["g_poly"]):
                            vis = n*tau**t*delta**d*del0**g
                            if c:
                                vis *= exp(-expdel**c)
                            mur += vis
                    # numerator of rational poly; denominator of rat. poly;
                    num = 0
                    if "n_num" in self._viscosity:
                        for n, t, d, c, g in zip(
                            self._viscosity["n_num"], self._viscosity["t_num"],
                            self._viscosity["d_num"], self._viscosity["c_num"],
                            self._viscosity["g_num"]):
                            num += n*tau**t*delta**d*del0**g
                            if c:
                                num *= exp(-expdel**c)
                    if "n_den" in self._viscosity:
                        den = 0
                        for n, t, d, c, g in zip(
                            self._viscosity["n_den"], self._viscosity["t_den"],
                            self._viscosity["d_den"], self._viscosity["c_den"],
                            self._viscosity["g_den"]):
                            den += n*tau**t*delta**d*del0**g
                            if c:
                                den *= exp(-expdel**c)
                    else:
                        den = 1.
                    mur += num/den

                    # numerator of exponential; denominator of exponential
                    num = 0
                    if "n_numexp" in self._viscosity:
                        for n, t, d, c, g in zip(
                            self._viscosity["n_numexp"], self._viscosity["t_numexp"],
                            self._viscosity["d_numexp"], self._viscosity["c_numexp"],
                            self._viscosity["g_numexp"]):
                            num += n*tau**t*delta**d*del0**g
                    if "n_denexp" in self._viscosity:
                        den = 0
                        for n, t, d, c, g in zip(
                            self._viscosity["n_denexp"], self._viscosity["t_denexp"],
                            self._viscosity["d_denexp"], self._viscosity["c_denexp"],
                            self._viscosity["g_denexp"]):
                            den += n*tau**t*delta**d*del0**g
                        else:
                            den = 1.
                    if "n_numexp" in self._viscosity:
                        mur += exp(num/den)

                    mur *= mured
                else:
                    mud, mur = 0, 0
                mu = muo+mud+mur

            elif self._viscosity["eq"] == 2:
                muo = self._Visco0()
                f = self._viscosity["F"]
                e = self._viscosity["E"]
                mu1 = f[0]+f[1]*(f[2]-log(T/f[3]))**2

                rhom = rho/self.M
                G = e[0]+e[1]/T
                H = rhom**0.5*(rhom-self._viscosity["rhoc"])/self._viscosity["rhoc"]
                F = G+(e[2]+e[3]*T**-1.5)*rhom**0.1+H*(e[4]+e[5]/self.T+e[6]/self.T**2)
                mu2 = exp(F)-exp(G)
                mu = muo+mu1*rhom+mu2

            elif self._viscosity["eq"] == 3:
                Tc = self._viscosity.get("Tref", 1.)
                rhoc = self._viscosity.get("rhoref", 1.)
                muref = self._viscosity.get("muref", 1.)
                tau = T/Tc
                delta = rho/self.M/rhoc

                muo = self._Visco0()

                mur = 0
                for n, t, d in zip(self._viscosity["n_num"],
                                   self._viscosity["t_num"], self._viscosity["d_num"]):
                    muo += n*tau**t*delta**d
                if "n_den" in self._viscosity:
                    den = 1
                    for n, t, d in zip(self._viscosity["n_den"],
                                       self._viscosity["t_den"], self._viscosity["d_den"]):
                        den *= n*tau**t*delta**d
                    mur /= den

                mu = (muo+mur)*muref

            elif self._viscosity["eq"] == 4:
                muo = self._Visco0()
                mur = 0
                Gamma = self.Tc/T
                psi1 = exp(Gamma)-1.0
                psi2 = exp(Gamma**2)-1.0
                a = self._viscosity["a"]
                b = self._viscosity["b"]
                c = self._viscosity["c"]
                A = self._viscosity["A"]
                B = self._viscosity["B"]
                C = self._viscosity["C"]
                D = self._viscosity["D"]
                ka = (a[0]+a[1]*psi1+a[2]*psi2)*Gamma
                kaa = (A[0]+A[1]*psi1+A[2]*psi2)*Gamma**3
                kr = (b[0]+b[1]*psi1+b[2]*psi2)*Gamma
                krr = (B[0]+B[1]*psi1+B[2]*psi2)*Gamma**3
                ki = (c[0]+c[1]*psi1+c[2]*psi2)*Gamma
                kii = (C[0]+C[1]*psi1+C[2]*psi2)*Gamma**3

                Prep = T*fase.dpdT_rho.barK
                Patt = self.P.bar-Prep
                Pid = rho*self.R*self.T/1e5
                delPr = Prep-Pid
                mur = kr*delPr + ka*Patt + krr*delPr**2 + kaa*Patt**2 + ki*Pid + kii*Pid**2 + D[0]*Prep**3*Gamma**2

                mu = (muo+mur*1e3)

            elif self._viscosity["eq"] == 5:
                a0 = (6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540,
                      -1.89992, 24.2745, 0.79716, -0.23816, 0.68629e-1)
                a1 = (50.4119, -0.11536e-2, 254.209, 38.0957, 7.63034,
                      -12.5367, 3.44945, 1.11764, 0.67695e-1, 0.34793)
                a2 = (-51.6801, -0.62571e-2, -168.481, -8.46414, -14.3544,
                      4.98529, -11.2913, 0.12348e-1, -0.8163, 0.59256)
                a3 = (1189.02, 0.37283e-1, 3898.27, 31.4178, 31.5267,
                      -18.1507, 69.3466, -4.11661, 4.02528, -0.72663)
                A = [None]
                for i in range(10):
                    A.append(a0[i]+a1[i]*self._viscosity["w"]+a2[i]*self._viscosity["mur"]**4+a3[i]*self._viscosity["k"])

                muo = self._Visco0()
                Y = rho/self.rhoc/6
                T_ = self.T/self._viscosity.get("ek", self.Tc/1.2593)
                G1 = (1-0.5*Y)/(1-Y)**3
                G2 = (A[1]*(1-exp(-A[4]*Y))/Y+A[2]*G1*exp(A[5]*Y)+A[3]*G1)/(A[1]*A[4]+A[2]+A[3])
                muk = muo*(1/G2+A[6]*Y)
                mup = 36.344e-6*(self.M*self.Tc)**0.5*self.rhoc**(2./3.)*A[7]*Y**2*G2*exp(A[8]+A[9]/T_+A[10]/T_**2)
                mu = muk+mup

            elif self._viscosity["eq"] == "ecs":
#                rhoc=self._constants.get("rhoref", self.rhoc)
#                Tc=self._constants.get("Tref", self.Tc)
#                delta=rho/rhoc
#                tau=Tc/T
#
#                fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta)
#
#                ref=self._constants["ref"](eq=self._constants["eq"])
#                Tr=T/Tc
#                rhor=rho/rhoc
#
#                psi=1+(self.f_acent-ref.f_acent)*(self._constants["ft"][0]+self._constants["ft"][1]*log(Tr))
#                for n, m in zip(self._constants["ft_add"], self._constants["ft_add_exp"]):
#                    psi+=n*Tr**m
#                for n, m in zip(self._constants["fd"], self._constants["fd_exp"]):
#                    psi+=n*rhor**m
#                T0=T*ref.Tc/self.Tc/psi
#
#                phi=ref.Zc/self.Zc*(1+(self.f_acent-ref.f_acent)*(self._constants["ht"][0]+self._constants["ht"][1]*log(Tr)))
#                for n, m in zip(self._constants["ht_add"], self._constants["ht_add_exp"]):
#                    phi+=n*Tr**m
#                for n, m in zip(self._constants["hd"], self._constants["hd_exp"]):
#                    phi+=n*rhor**m
#                rho0=rho*ref.rhoc/self.rhoc*phi
#
#                ref_rhoc=ref._constants.get("rhoref", ref.rhoc)
#                ref_Tc=ref._constants.get("Tref", ref.Tc)
#                deltaref=rho0/ref_rhoc
#                tauref=ref_Tc/T0
#                fir, firt, firtt, fird, firdd, firdt, firdtt, B, C=ref._phir(tauref, deltaref)

                Tr = T/self.Tc
                rhor = rho/self.rhoc

                fint = 0
                for n, m in zip(self._viscosity["fint"], self._viscosity["fint_t"]):
                    fint += n*T**m

                psi = 0
                for n, t, d in zip(self._viscosity["psi"], self._viscosity["psi_t"], self._viscosity["psi_d"]):
                    psi += n*Tr**t*rhor**d

                mu = None

        else:
            mu = None
        return unidades.Viscosity(mu, "muPas")

    def _KCritical(self, rho, T, fase):
        """Enchancement thermal conductivity calculation for critical region"""
        if self._thermal["critical"] == 0:
            tc = 0

        elif self._thermal["critical"] == 1:
            tc = 0
            if "crit_num_n" in self._thermal:
                Tref = self._thermal["crit_num_Tref"]
                if Tref < 0:
                    tau = Tref/T
                else:
                    tau = T/Tref
                delta = rho/self.M/self._thermal["crit_num_rhoref"]

                for n, alfa, t, beta, d in zip(self._thermal["crit_num_n"], self._thermal["crit_num_alfa"], self._thermal["crit_num_t"], self._thermal["crit_num_beta"], self._thermal["crit_num_d"]):
                    tc += n*(tau+alfa)**t*(delta+beta)**d

                if "crit_den_n" in self._thermal:
                    den = 0
                    for n, alfa, t, beta, d, c in zip(self._thermal["crit_den_n"], self._thermal["crit_den_alfa"], self._thermal["crit_den_t"], self._thermal["crit_den_beta"], self._thermal["crit_den_d"], self._thermal["crit_den_c"]):
                        if c == 99:
                            tcmax = max(tau, alfa-tau)
                            den += n*tcmax**t*(delta+beta)**d
                        else:
                            den += n*(tau+alfa)**t*(delta+beta)**d
                    tc /= den

            if "crit_exp_n" in self._thermal:
                Tref = self._thermal["crit_exp_Tref"]
                if Tref < 0:
                    tau = Tref/T
                else:
                    tau = T/Tref
                delta = rho/self.M/self._thermal["crit_exp_rhoref"]
                expo = 0
                for n, alfa, t, beta, d, c in zip(self._thermal["crit_exp_n"], self._thermal["crit_exp_alfa"], self._thermal["crit_exp_t"], self._thermal["crit_exp_beta"], self._thermal["crit_exp_d"], self._thermal["crit_exp_c"]):
                    expo += n*(tau+alfa)**t*(delta+beta)**d
                tc *= exp(expo)

            tc *= self._thermal["crit_num_k"]

        elif self._thermal["critical"] == 2:
            X = self._thermal["X"]
            xi = self.Pc*rho/self.rhoc**2/self.derivative("P", "rho", "T", fase)
            normterm = X[3]*Boltzmann/self.Pc*(T*fase.dpdT*self.rhoc/rho)**2*xi**X[2]
            delT = abs(T-self.Tc)/self.Tc
            delrho = abs(rho-self.rhoc)/self.rhoc
            expterm = exp(-(X[0]*delT**4+X[1]*delrho**4))
            mu = self._Viscosity()
            tc = normterm*expterm/(6*pi*mu*self._thermal["Z"])

        elif self._thermal["critical"] == 3:
            qd = self._thermal["qd"]
            Tref = self._thermal["Tcref"]
#            ref = self.__class__(T=Tref, rho=rho, recursion=False)
            x_T = self.Pc*rho/self.rhoc**2*self.derivative("rho", "P", "T", fase)
            x_Tr = self.Pc*rho/self.rhoc**2*self.derivative("rho", "P", "T", fase)*Tref/T
            delchi = x_T-x_Tr
            if delchi <= 0:
                tc = 0
            else:
                Xi = self._thermal["Xio"]*(delchi/self._thermal["gam0"])**(self._thermal["gnu"]/self._thermal["gamma"])
                omega = 2/pi*((fase.cp-fase.cv)/fase.cp*arctan(Xi/qd)+fase.cv/fase.cp*Xi/qd)
                omega0 = 2/pi*(1-exp(-1/(1./qd/Xi+Xi**2*qd**2/3*(self.rhoc/rho)**2)))
                tc = rho/self.M*fase.cp*Boltzmann*self._thermal["R0"]*T/(6*pi*Xi*fase.mu)*(omega-omega0)

        elif self._thermal["critical"] == 4:
            rhom = rho/self.M
            Xt = (rhom*self._thermal["Pcref"]/self._thermal["rhocref"]**2/self.derivative("P", "rho", "T", fase))**self._thermal["expo"]
            parterm = self._thermal["alfa"]*Boltzmann/self._thermal["Pcref"]*(T*fase.dpdT*self._thermal["rhocref"]/rhom)**2*Xt*1e21
            delT = abs(T-self._thermal["Tcref"])/self._thermal["Tcref"]
            delrho = abs(rhom-self._thermal["rhocref"])/self._thermal["rhocref"]
            expterm = exp(-(self._thermal["alfa"]*delT**2+self._thermal["beta"]*delrho**4))
            tc = parterm*expterm/(6*pi*self._thermal["Xio"]*fase.mu.muPas)*self._thermal["kcref"]

        elif self._thermal["critical"] == "NH3":
            tr = abs(T-405.4)/405.4
            trr = tr
            if 404.4 < T < 406.5 and (rho/self.M<9.6 or rho/self.M>18): #to avoid infinite value in critical region
                trr = 0.002
            etab = 1.0e-5*(2.6+1.6*tr)
            dPT = 1.0e5*(2.18-0.12/exp(17.8*tr))
            if trr == 0.0:
                tcrhoc = 1.e20
                dtcid = tcrhoc
                xcon = -1.e20
            else:
                tcrhoc = 1.2*1.38066e-23*T**2*dPT**2*0.423e-8/trr**1.24*\
                    (1.0+1.429*tr**0.5)/(6.0*pi*etab*(1.34e-10/trr**0.63*(1.0+1.0*tr**0.5)))
                dtcid = tcrhoc*exp(-36.0*tr**2)
                xcon = 0.61*235+16.5*log(trr)
            if rho/self.rhoc < 0.6:
                tccsw = dtcid*xcon**2/(xcon**2+(141.0-0.96*235.0)**2)
                tc = tccsw*rho**2/141.0**2
            else:
                tc = dtcid*xcon**2/(xcon**2+(rho-0.96*235.0)**2)

        elif self._thermal["critical"] == "CH4":
            tau = self.Tc/T
            delta = rho/self.rhoc
            ts = (self.Tc-T)/self.Tc
            ds = (self.rhoc-rho)/self.rhoc
            xt = 0.28631*delta*tau/self.derivative("P", "rho", "T", fase)
            ftd = exp(-2.646*abs(ts)**0.5+2.678*ds**2-0.637*ds)
            tc = 91.855/fase.mu/tau**2*fase.dpdT_rho**2*xt**0.4681*ftd*1e-3

        return tc

    def _ThCond(self, rho, T, fase):
        """Thermal conductivity calculation"""
        if self._thermal:
            if self._thermal["eq"] == 0:
                k = self.__getattribute__(self._thermal["method"])(rho, T, fase)*1e3

            elif self._thermal["eq"] == 1:
                # Dilute gas terms
                kg = 0
                if "no" in self._thermal:
                    tau = self._thermal["Tref"]/T
                    for n, c in zip(self._thermal["no"], self._thermal["co"]):
                        if c == -99:
                            cpi = 1.+n*(self.cp0.kJkgK-2.5*self.R.kJkgK)
                            kg *= cpi
                        elif c == -98:
                            muo = self._Visco0()
                            cpi = self.cp0-R
                            kg += n*cpi*muo
                        elif c == -97:
                            muo = self._Visco0()
                            kg += n*muo
                        elif c == -96:
                            cpi = self.cp0/self.R-2.5
                            muo = self._Visco0()
                            kg = (kg*cpi+15./4.)*self.R.kJkgK*muo/self.M
                        else:
                            kg += n*tau**c

                    if "noden" in self._thermal:
                        den = 0
                        for n, t in zip(self._thermal["noden"], self._thermal["toden"]):
                            den += n*tau**t
                        kg /= den

                    kg *= self._thermal["kref"]

                # Backgraund terms
                kb = 0
                kc = 0
                if rho > 0:
                    if "nb" in self._thermal:
                        tau = self._thermal["Trefb"]/T
                        delta = rho/self.M/self._thermal["rhorefb"]
                        for n, t, d, c in zip(self._thermal["nb"], self._thermal["tb"], self._thermal["db"], self._thermal["cb"]):
                            if c == -99:
                                if tau < 1:
                                    th = (1.-tau)**(1./3.)
                                    kb /= exp(-1.880284*th**1.062 -
                                              2.8526531*th**2.5-3.000648*th**4.5 -
                                              5.251169*th**7.5-13.191869*th**12.5 -
                                              37.553961*th**23.5)
                            else:
                                if c != 0:
                                    kb += n*tau**t*delta**d*exp(-delta**c)
                                else:
                                    kb += n*tau**t*delta**d

                        if "nbden" in self._thermal:
                            den = 0
                            for n, t, d in zip(self._thermal["nbden"], self._thermal["tbden"], self._thermal["dbden"]):
                                den += n*tau**t*delta**d
                            kb /= den

                        kb *= self._thermal["krefb"]

                    # Critical enhancement
                    kc = self._KCritical(rho, T, fase)

                k = kg+kb+kc

            elif self._thermal["eq"] == 2:
                self._viscosity = self._thermal["visco"]
                muo = self._Visco0()
                G = self._thermal["G"]
                kg = 1e-3*muo/self.M*(3.75*self.R+(self.cp0.kJkgK-2.5*self.R)*(G[0]+G[1]*self._viscosity["ek"]/T))

                E = self._thermal["E"]
                F0 = E[0]+E[1]/T+E[2]/T**2
                F1 = E[3]+E[4]/T+E[5]/T**2
                F2 = E[6]+E[7]/T
                rhom = rho/self.M
                kb = (F0+F1*rhom)*rhom/(1+F2*rhom)

                # Critical enhancement
                kc = self._KCritical(rho, T, fase)

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
                k = None
                # TODO:

        else:
            k = None
        return unidades.ThermalConductivity(k)

    def txt(self):
        """Return a text repr of class with all properties"""
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Advanced MEoS properties")
        txt += "-------------------#"+os.linesep
        doc = self._constants["__doi__"]["autor"] + "; " + \
              self._constants["__doi__"]["title"] + "; " + \
              self._constants["__doi__"]["ref"]
        txt += os.linesep + doc + os.linesep
        
        if 0 < self.x < 1:
            param = "%-40s\t%20s\t%20s"
        else:
            param = "%-40s\t%s"
        if self.x == 0:
            txtphases = "%60s" % QApplication.translate("pychemqt", "Liquid")+os.linesep
            phases = [self.Liquido]
        elif self.x == 1:
            txtphases = "%60s" % QApplication.translate("pychemqt", "Gas")+os.linesep
            phases = [self.Gas]
        else:
            txtphases = "%60s\t%20s" % (QApplication.translate("pychemqt", "Liquid"), 
                                 QApplication.translate("pychemqt", "Gas"))+os.linesep
            phases = [self.Liquido, self.Gas]
            
        complejos = ""
        for propiedad, key, unit in data:
            if key in _fase.__dict__:
                values = [propiedad]
                for phase in phases:
                    values.append(phase.__getattribute__(key).str)
                complejos += param % tuple(values) +os.linesep
            else:
                txt += os.linesep
                txt += "%-40s\t%s" % (propiedad, self.__getattribute__(key).str)
        txt += os.linesep + os.linesep + txtphases + complejos
        return txt

class MEoSBlend(MEoS):
    """Special meos class to implement pseudocomponente blend and defining its
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


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()


#    water=H2O(T=300., P=0.1)
#    print "%0.1f %0.2f %0.4f %0.9f %0.9f %0.5f %0.4f %0.9f %0.4f %0.4f %0.9f %0.6f" % (water.T, water.P.MPa, water.rho, water.kappa, water.alfav, water.n, water.kt, water.ks, water.Kt.MPa, water.Ks.MPa, water.deltat, water.Gruneisen)
#    print -water.derivative("P", "v", "T")/1e6, water.betap
#        alfap        -   Relative pressure coefficient, 1/K
#        betap       -   Isothermal stress coefficient, kg/m³
#        betas       -   Isoentropic temperature-pressure coefficient

#    water=H2O(T=300., x=0.5)
#    print water.x
#    print water.Liquido.Z
#    print water.Gas.Z
#    print water.Z

#    aire=Air( P=0.1, T=300)
#    print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w), aire.P.MPa
#

#    p=unidades.Pressure(1, "atm")
#    water1=H2O(T=400, P=p.MPa)
#    water2=H2O(T=450, P=p.MPa)
#    print "%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK)
#    0.0992418352 4.13018112 1501.51914 0.393062643
#    print water2.h.MJkg-water1.h.MJkg

#    aire=R507A(T=500, P=0.1)
#    print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w), aire.P.MPa

#    agua=H2O(T=350, x=0.5)
#    print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w), aire.P.MPa

#    kwargs={"T": 110., "x": 0.0}
#    aire=N2(**kwargs)
#    print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.Liquido.cv.kJkgK, aire.Liquido.cp.kJkgK, aire.P.MPa)
#    for key in keys:
#        print key, aire.__getattribute__(key)

##    T=50.
##    P=N2._Sublimation_Pressure(T).MPa
##    print T, P
##    aire=N2(T=T, P=P)
##    print "%0.1f %0.4f %0.3f %0.3f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.P.MPa)
#
##    for t in range(80, 160, 10):
##    methanol=R23(T=300, P=1)
##    print  methanol.mu.muPas
#

#    aire=pH2(T=300, P=1)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas

#    aire=O2(T=300, P=1.)
#    aire2=O2(T=300, P=1., visco=1, thermal=1)
#    aire3=O2(T=300, P=1., visco=2, thermal=2)
#    aire4=O2(T=300, P=1., visco=3, thermal=3)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.T, aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas
#    print  aire4.T, aire4.P.MPa, aire4.rho, aire4.k.mWmK, aire4.mu.muPas

#    aire=CO2(T=300, P=0.1)
#    aire2=CO2(T=300, P=0.1, visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=CO(T=300, P=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=He(T=300, P=1, eq=1)
#    aire2=He(T=300, P=1, eq=1, thermal=1)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas, aire.k.mWmK, aire2.k.mWmK

#    aire=Ne(T=300, P=1.)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas

#    aire=Ar(T=300, P=1., visco=0, thermal=0)
#    aire2=Ar(T=300, P=1., visco=1, thermal=1)
#    aire3=Ar(T=300, P=1., visco=2, thermal=2)
#    aire4=Ar(T=300, P=1., visco=3, thermal=3)
#    print  aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas
#    print  aire4.P.MPa, aire4.rho, aire4.k.mWmK, aire4.mu.muPas
#    print  aire.T, aire.P.MPa, aire.rho, aire.h.kJkg, aire.s.kJkgK

#    aire=Xe(T=300, P=1.)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=N2(T=300, P=1., visco=0)
#    aire2=N2(T=300, P=1., visco=1, thermal=1)
#    aire3=N2(T=300, P=1., visco=2, thermal=2)
#    aire4=N2(T=300, P=1., visco=3, thermal=3)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.T, aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas
#    print  aire4.T, aire4.P.MPa, aire4.rho, aire4.k.mWmK, aire4.mu.muPas

#    aire=nC4(T=300, P=1., visco=0)
#    aire2=nC4(T=300, P=1., visco=1, thermal=1)
#    aire3=nC4(T=300, P=1., visco=2)
#    print  aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas

#    aire=C3(T=300, P=0.1, visco=0)
#    aire2=C3(T=300, P=0.1, visco=1, thermal=1)
#    aire3=C3(T=300, P=0.1, visco=2)
#    print  aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas

#    aire=Air(T=300, P=1., visco=0, thermal=0)
#    aire2=Air(T=300, P=1., visco=1, thermal=1)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas, aire2.mu.muPas, aire.k.mWmK, aire2.k.mWmK

#    aire=NH3(T=300, P=1., visco=0)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas, aire.k.mWmK

#    aire=nC12(T=500, P=0.1)
#    print  aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=Ethylene(T=300, P=1, visco=0, thermal=0)
#    aire1=Ethylene(T=300, P=1, visco=1, thermal=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire1.T, aire1.P.MPa, aire1.rho, aire1.k.mWmK, aire1.mu.muPas

#    aire=Ethanol(T=300, P=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas


#    aire=C2(T=300, P=1., visco=0)
#    aire2=C2(T=300, P=1., visco=1)
#    aire3=C2(T=300, P=1., visco=2)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.T, aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas

#    aire=H2S(T=300, P=1., visco=0)
#    aire2=H2S(T=300, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=nC6(T=300, P=1., visco=0)
#    aire2=nC6(T=300, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=nC7(T=500, P=1., visco=0)
#    aire2=nC7(T=500, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=H2(T=300, P=1, visco=0)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas

#    aire=iC5(T=300, P=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=iC6(T=300, P=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=iC4(T=300, P=1., visco=0)
#    aire2=iC4(T=300, P=1., visco=1, thermal=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=CH4(T=300, P=1, visco=0)
#    aire2=CH4(T=300, P=1., visco=1, thermal=1)
#    aire3=CH4(T=300, P=1., visco=2)
#    aire4=CH4(T=300, P=1., visco=3)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.T, aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas
#    print  aire4.T, aire4.P.MPa, aire4.rho, aire4.k.mWmK, aire4.mu.muPas

#    aire=Methanol(T=300, P=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=nC9(T=300, P=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=nC8(T=300, P=1., visco=0)
#    aire2=nC8(T=300, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=nC5(T=300, P=1., visco=0)
#    aire2=nC5(T=300, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=H2O(T=300, P=1., visco=0)
#    aire2=H2O(T=300, P=1., visco=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas

#    aire=D2O(T=300, P=1., visco=0)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas, aire.k.mWmK

#    aire=R23(T=300, P=1., visco=0)
#    print  aire.P.MPa, aire.rho, aire.mu.muPas, aire.k.mWmK

#    aire=R32(T=300, P=1.)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=R152a(T=300, P=1., eq=1, visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=R134a(T=300, P=0.1, visco=0)
#    aire2=R134a(T=300, P=0.1, visco=1)
#    aire3=R134a(T=300, P=0.1, visco=2)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas
#    print  aire2.T, aire2.P.MPa, aire2.rho, aire2.k.mWmK, aire2.mu.muPas
#    print  aire3.T, aire3.P.MPa, aire3.rho, aire3.k.mWmK, aire3.mu.muPas

#    aire=R125(T=300, P=1., visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=R123(T=500, P=1., eq=1)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=R245fa(T=300, P=1., visco=0)
#    print  aire.T, aire.P.MPa, aire.rho, aire.k.mWmK, aire.mu.muPas

#    aire=iButene(T=298.15, P=0.101325)
#    print  aire.T, aire.P.MPa, aire.rho, aire.h.kJkg, aire.s.kJkgK
#
#    aire=H2O(T=298.15, P=0.101325)
#    print "%0.2f %0.6f %0.10f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.P.MPa, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w)
#
#    aire=H2O(T=500, P=1)
#    print "%0.2f %0.6f %0.10f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.P.MPa, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w)
#    print aire.T, aire.P.MPa, aire.rho, aire.h.kJkg, aire.s.kJkgK

#Calculate estado estandard  OTO (h,s=0 at 25ºC and 1 atm)
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

#    octano=nC8(T=500, P=0.1, eq=2)
#    print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (octano.T, octano.rho, octano.u.kJkg, octano.h.kJkg, octano.s.kJkgK, octano.cv.kJkgK, octano.cp.kJkgK, octano.w)

#    argon=Ar(**{'P': 0.101325, 'T': 490.0})
#    print argon.T, argon.rho

#    agua=nC8(T=300, x=1.)
#    print agua.T, agua.rho, agua._Vapor_Density()

#    helio=He(T=4, x=0.0, eq=1)
#    print "%0.1f %0.6f %0.5f %0.4f %0.4f %0.4f %0.4f" % (helio.T, helio.P.MPa, helio.Liquido.rho, helio.h.kJkg, helio.s.kJkgK, helio.Liquido.cp.kJkgK, helio.Gas.cp.kJkgK)

#    helio=H2O(T=300, P=0.1)
#    print "%0.1f %0.6f %0.5f %0.4f %0.4f" % (helio.T, helio.P.MPa, helio.rho, helio.h.kJkg, helio.s.kJkgK)
#    print helio.cp.kJkgK, helio.x

#    aire=nC10(T=273.15+20, P=0.101325)
#    print  aire.mu.muPas

#    agua=H2O(T=298.15, P=0.101325, visco=0, thermal=0)
#    print  agua.P.MPa, agua.rho

    dme=MEoS(T=400, P=0.1, eq=1)
    print "%0.1f %0.5f %0.3f %0.5f %0.5f %0.5f %0.2f" % (dme.T, dme.rho, dme.h.kJkg, dme.s.kJkgK, dme.cv.kJkgK, dme.cp.kJkgK, dme.w)

#300,00	0,10000	1,6025	0,31238	0,52152	0,52033	3,6153	-0,00038001	0,00000066198	17,837	22,741

