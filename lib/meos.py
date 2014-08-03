#!/usr/bin/python
# -*- coding: utf-8 -*-

#############################################################################
#Implementación de la ecuación de estado de multiparámetros
#   o   Ecuación de estado Setzmann-Wagner, basada en la energía de Helmholtz
#   o   Ecuación MBWR
#   o   Ecuación Peng-Robinson con translación de Peneloux
#############################################################################

import cPickle, os

from PyQt4.QtGui import QApplication
from scipy import exp, log, log10, sin, sinh, cosh, tanh, arctan, __version__
if int(__version__.split(".")[1])<10:
    from scipy.constants import Bolzmann as Boltzmann
else:
    from scipy.constants import  Boltzmann
from scipy.constants import pi, Avogadro, R
from scipy.optimize import fsolve

from lib import  unidades, compuestos
from physics import R_atml
from config import fluid


Tref=298.15
Pref=101325.
#so=0
#ho=0

data=[(QApplication.translate("pychemqt", "Temperature"), "T"),
            (QApplication.translate("pychemqt", "Pressure"), "P"),
            (QApplication.translate("pychemqt", "Density"), "rho"),
            (QApplication.translate("pychemqt", "Volume"), "v"),
            (QApplication.translate("pychemqt", "Internal Energy"), "u"),
            (QApplication.translate("pychemqt", "Enthalpy"), "h"),
            (QApplication.translate("pychemqt", "Entropy"), "s"),
            ("Cv", "cv"),
            ("Cp", "cp"),
            ("Cp0", "cp0"),
            ("Cp/Cv", "cp_cv"),
            ("Cp0/Cv", "cp0_cv"),
            (QApplication.translate("pychemqt", "Gibbs Free Energy"), "g"),
            (QApplication.translate("pychemqt", "Helmholtz Free Energy"), "a"),
            ("Hvapor", "Hvap"),
            (QApplication.translate("pychemqt", "Fugacity"), "f"),
            (QApplication.translate("pychemqt", "Fugacity coef."), "fi"),
            (QApplication.translate("pychemqt", "Speed sound"), "w"),
            (QApplication.translate("pychemqt", "Compresibility"), "Z"),
            (QApplication.translate("pychemqt", "Quality"), "x"),
            ("Tr", "Tr"),
            ("Pr", "Pr"),
            (QApplication.translate("pychemqt", "Joule-Thomson coefficient"), "joule"),
            (QApplication.translate("pychemqt", "2nd virial coefficient"), "virialB"),
            (QApplication.translate("pychemqt", "3er virial coefficient"), "virialC"),

            (QApplication.translate("pychemqt", "Viscosity"), "mu"),
            (QApplication.translate("pychemqt", "Thermal conductivity"), "k"),
            (QApplication.translate("pychemqt", "Kinematic viscosity"), "nu"),
            (QApplication.translate("pychemqt", "Surface tension"), "sigma"),
            (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon"),
            (QApplication.translate("pychemqt", "Prandtl number"), "Pr"),
            (QApplication.translate("pychemqt", "Thermal diffusivity"), "alfa"),

            (QApplication.translate("pychemqt", "Isoentropic exponent"), "gamma"),
            (QApplication.translate("pychemqt", "Volumetric Expansitivy"), "alfav"), #Thermal expansion coefficient
            (QApplication.translate("pychemqt", "Isotermic compresibility"), "kappa"),
            (QApplication.translate("pychemqt", "Relative pressure coefficient"), "alfap"),  #dPdT_rho
            (QApplication.translate("pychemqt", "Isothermal stress"), "betap"),
            (QApplication.translate("pychemqt", "Isentropic coefficient"), "betas"),
            (QApplication.translate("pychemqt", "Isothermal throttle coefficient"), "deltat"),
            (QApplication.translate("pychemqt", "Isentropic expansion"), "n"),
            (QApplication.translate("pychemqt", "Isothermal expansion"), "kt"),
            (QApplication.translate("pychemqt", "Adiabatic compresibility"), "ks"),
            (QApplication.translate("pychemqt", "Isothermal modulus"), "Ks"),
            (QApplication.translate("pychemqt", "Adiabatic modulus"), "Kt"),
            ("Gruneisen", "Gruneisen"),
            (QApplication.translate("pychemqt", "Internal pressure"), "InternalPressure"),
            ("dhdp_rho", "dhdp_rho"),
            ("dhdT_rho", "dhdT_rho"),
            ("dhdT_P", "dhdT_P"),
            ("dhdP_T", "dhdP_T")]    #Isothermal throttle coefficient

propiedades=[p[0] for p in data]
keys=[p[1] for p in data]
properties=dict(zip(keys, propiedades))


class MEoS(object):
    """Implementación general de ecuaciones de estado de multiparámetros
    Superclase de la que derivan todas los componentes que dispongan de parámetros para ella
    """

    _dielectric=None
    _melting=None
    _sublimation=None
    _surface=None
    _vapor_Pressure=None
    _liquid_Density=None
    _vapor_Density=None
    _omega=None
    _viscosity=None,
    _thermal=None,
    _critical=None
    _PR=0.
    id=None

    _Tr=None
    _rhor=None
    _w=None

    kwargs={"T": 0.0,
                    "P": 0.0,
                    "rho": 0.0,
                    "v": 0.0,
                    "h": 0.0,
                    "s": 0.0,
                    "u": 0.0,
                    "x": None,
                    "v": 0.0,

                    "eq": 0,
                    "visco": 0,
                    "thermal": 0,
                    "ref": None,
                    "recursion": True
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
            None: Se usan ho,so=0, válido para calcular los valores de referencia

        Properties:
        P            -   Pressure, MPa
        T            -   Temperature, K
        g            -   Specific Gibbs free energy, kJ/kg
        a            -   Specific Helmholtz free energy, kJ/kg
        v            -   Specific volume, m³/kg
        rho         -   Density, kg/m³
        h            -   Specific enthalpy, kJ/kg
        u            -   Specific internal energy, kJ/kg
        s             -   Specific entropy, kJ/kg·K
        x             -   Quatity
        cp           -   Specific isobaric heat capacity, kJ/kg·K
        cv           -   Specific isochoric heat capacity, kJ/kg·K
        Z             -   Compression factor
        f              -   Fugacity, MPa
        fi             -   Fugacity coefficient
        gamma    -   Isoentropic exponent
        alfav        -   Thermal expansion coefficient, 1/K
        kappa      -   Isothermal compressibility, 1/MPa
        alfap        -   Relative pressure coefficient, 1/K
        betap       -   Isothermal stress coefficient, kg/m³
        betas       -   Isoentropic temperature-pressure coefficient
        joule         -   Joule-Thomson coefficient, K/MPa
        dhdP_T     -   Isothermal throttling coefficient, kJ/kg·MPa
        n              -   Isentropic Expansion Coefficient
        kt             -   Isothermal Expansion Coefficient
        ks             -   Adiabatic Compressibility, 1/MPa
        Ks             -   Adiabatic bulk modulus, MPa
        Kt             -   Isothermal bulk modulus, MPa

        virialB       -   Second virial coefficient
        virialC       -   Third virial coefficient

        mu           -   Dynamic viscosity, Pa·s
        nu            -   Kinematic viscosity, m²/s
        k              -   Thermal conductivity, W/m·K
        sigma       -   Surface tension, N/m
        alfa           -   Thermal diffusivity, m²/s
        Pr             -   Prandtl number
        """
        self._constants=self.eq[self.kwargs["eq"]]
        self.R=unidades.SpecificHeat(self._constants["R"]/self.M, "kJkgK")
        self.Zc=self.Pc/self.rhoc/self.R/self.Tc
        self.kwargs=MEoS.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status=1
            self.calculo()
            self.msg=""

    @property
    def calculable(self):
        self._mode=""
        if self.kwargs["T"] and self.kwargs["P"]:
            self._mode="TP"
        elif self.kwargs["T"] and self.kwargs["rho"]:
            self._mode="Trho"
        elif self.kwargs["P"] and self.kwargs["h"]:
            self._mode="Ph"
        elif self.kwargs["P"] and self.kwargs["s"]:
            self._mode="Ps"
        elif self.kwargs["P"] and self.kwargs["v"]:
            self._mode="Pv"
        elif self.kwargs["T"] and self.kwargs["s"]:
            self._mode="Ts"
        elif self.kwargs["T"] and self.kwargs["x"]!=None:
            self._mode="Tx"
        return bool(self._mode)


    def calculo(self):

        T=self.kwargs["T"]
        rho=self.kwargs["rho"]
        P=self.kwargs["P"]
        v=self.kwargs["v"]
        s=self.kwargs["s"]
        h=self.kwargs["h"]
        u=self.kwargs["u"]
        x=self.kwargs["x"]
        recursion=self.kwargs["recursion"]
        eq=self.kwargs["eq"]
        visco=self.kwargs["visco"]
        thermal=self.kwargs["thermal"]
        ref=self.kwargs["ref"]

        if ref==None:
            ho=so=0
        else:
            if ref=="CUSTOM":
                file=conf_dir+"std.pkl"
            else:
                file=os.environ["pychemqt"]+"dat/%s.pkl" % ref
            with open(file) as archivo:
                std=cPickle.load(archivo)
                ho, so=std[self.__class__.__name__][str(eq)]

        if self.id:
            self.componente=compuestos.Componente(self.id)

        #Opcion de aceptar el nombre interno de la ecuacion
        if isinstance(eq, str) and eq in self.__class__.__dict__:
            eq=self.eq.index(self.__class__.__dict__[eq])

        if eq=="PR":
            self._eq=self._PengRobinson
            self._constants=self.eq[0]
        elif eq=="Generalised":
            self._eq=self._Helmholtz
            self._Generalised()
        elif eq=="GERG":
            try:
                self._constants=self.GERG
            except:
                self._constants=self.eq[0]
            if self._constants["__type__"]=="Helmholtz":
                self._eq=self._Helmholtz
            else:
                self._eq=self._MBWR
        elif self.eq[eq]["__type__"]=="Helmholtz":
            self._eq=self._Helmholtz
            self._constants=self.eq[eq]
        elif self.eq[eq]["__type__"]=="MBWR":
            self._eq=self._MBWR
            self._constants=self.eq[eq]
        elif self.eq[eq]["__type__"]=="ECS":
            self._eq=self._ECS
            self._constants=self.eq[eq]

        self._viscosity=self._viscosity[visco]
        self._thermal=self._thermal[thermal]

        self.R=unidades.SpecificHeat(self._constants["R"]/self.M, "kJkgK")

        propiedades=None
        if v and not rho:
            rho=1./v

        if T and x!=None:
            if self.Tt>T or self.Tc<T:
                raise ValueError ("Wrong input values")

            self.T=unidades.Temperature(T)
            self.x=x
            self.Tr=self.T/self.Tc

            rhol, rhov, Ps=self._saturation()
            self.P=unidades.Pressure(Ps)

            self.Liquido=self.__class__(T=T, rho=rhol, recursion=False)
            self.Gas=self.__class__(T=T, rho=rhov, recursion=False)

            self.v=unidades.SpecificVolume(self.Liquido.v*(1-x)+self.Gas.v*x)
            self.rho=unidades.Density(1./self.v)
            self.s=unidades.SpecificHeat(self.Liquido.s*(1-x)+self.Gas.s*x)
            self.h=unidades.Enthalpy(self.Liquido.h*(1-x)+self.Gas.h*x)
            self.u=unidades.Enthalpy(self.Liquido.u*(1-x)+self.Gas.u*x)
            self.a=unidades.Enthalpy(self.Liquido.a*(1-x)+self.Gas.a*x)
            self.g=unidades.Enthalpy(self.Liquido.g*(1-x)+self.Gas.g*x)
            self.Z=None
            self.cp=None
            self.cv=None
            self.w=None
            self.alfap=None
            self.betap=None
            self.alfav=None
            self.kappa=None
            self.joule=None
            self.dhdP_T=None
            self.betas=None
            self.gamma=None
            self.fi=self.Liquido.fi
            self.f=unidades.Pressure(self.fi*self.P)
            self.virialB=self.Liquido.virialB
            self.virialC=self.Liquido.virialC

        else:
            if T and P:
                if T<self.Tc and P<self.Pc.MPa and self._Vapor_Pressure(T).MPa<P:
                    rhoo=self._Liquid_Density(T)
                elif T<self.Tc and P<self.Pc.MPa:
                    rhoo=self._Vapor_Density(T)
                else:
                    rhoo=self.rhoc
#                    print rhoo
#                    rhoo=2.
#                    rhoo=self.componente.RhoG_Lee_Kesler(T, P/0.101325)
                rho=fsolve(lambda rho: self._eq(rho, T)["P"]-P*1e6, rhoo)
                #propiedades=self._eq(rho[0], T)
                propiedades=self._eq(rho, T)
            elif T and rho:
                propiedades=self._eq(rho, T)
            elif T and h!=None:
                rho=fsolve(lambda rho: self._eq(rho, T)["h"]-h, 200)
                propiedades=self._eq(rho[0], T)
            elif T and s!=None:
                rho=fsolve(lambda rho: self._eq(rho, T)["s"]-s, 200)
                propiedades=self._eq(rho[0], T)
            elif T and u!=None:
                def funcion(rho):
                    par=self._eq(rho, T)
                    return par["h"]-par["P"]/1000*par["v"]-u
                rho=fsolve(funcion, 200)
                propiedades=self._eq(rho[0], T)

            elif P and rho:
                T=fsolve(lambda T: self._eq(rho, T)["P"]-P*1e6, 600)
                propiedades=self._eq(rho, T[0])
            elif P and h!=None:
                rho, T=fsolve(lambda par: (self._eq(par[0], par[1])["P"]-P*1e6, self._eq(par[0], par[1])["h"]-h), [200, 600])
                propiedades=self._eq(rho[0], T[0])
            elif P and s!=None:
                rho, T=fsolve(lambda par: (self._eq(par[0], par[1])["P"]-P*1e6, self._eq(par[0], par[1])["s"]-s), [200, 600])
                propiedades=self._eq(rho[0], T[0])
            elif P and u!=None:
                def funcion(parr):
                    par=self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]/1000*par["v"]-u, par["P"]-P*1e6
                rho, T=fsolve(funcion, [200, 600])
                propiedades=self._eq(rho[0], T[0])

            elif rho and h!=None:
                T=fsolve(lambda T: self._eq(rho, T)["h"]-h, 600)
                propiedades=self._eq(rho, T[0])
            elif rho and s!=None:
                T=fsolve(lambda T: self._eq(rho, T)["s"]-s, 600)
                propiedades=self._eq(rho, T[0])
            elif rho and u!=None:
                def funcion(T):
                    par=self._eq(rho, T)
                    return par["h"]-par["P"]/1000*par["v"]-u
                T=fsolve(funcion, 600)
                propiedades=self._eq(rho, T[0])

            elif h!=None and s!=None:
                rho, T=fsolve(lambda par: (self._eq(par[0], par[1])["h"]-h, self._eq(par[0], par[1])["s"]-s), [200, 600])
                propiedades=self._eq(rho[0], T[0])
            elif h!=None and u!=None:
                def funcion(parr):
                    par=self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]/1000*par["v"]-u, par["h"]-h
                rho, T=fsolve(funcion, [200, 600])
                propiedades=self._eq(rho[0], T[0])

            elif s!=None and u!=None:
                def funcion(parr):
                    par=self._eq(parr[0], parr[1])
                    return par["h"]-par["P"]/1000*par["v"]-u, par["s"]-s
                rho, T=fsolve(funcion, [200, 600])
                propiedades=self._eq(rho[0], T[0])

#            if propiedades["P"]<=0:
#                raise ValueError ("Wrong input values")

            self.T=unidades.Temperature(T)
            self.Tr=unidades.Dimensionless(T/self.Tc)
            self.s=unidades.SpecificHeat(propiedades["s"]-so, "kJkgK")
            self.P=unidades.Pressure(propiedades["P"])
            self.Pr=unidades.Dimensionless(self.P/self.Pc)
            self.v=unidades.SpecificVolume(propiedades["v"])
            self.rho=unidades.Density(1/self.v)
            self.h=unidades.Enthalpy(propiedades["h"]-ho, "kJkg")
            self.u=unidades.Enthalpy(self.h-self.P*self.v)
            self.a=unidades.Enthalpy(self.u-self.T*self.s)
            self.g=unidades.Enthalpy(self.h-self.T*self.s)
            self.Z=unidades.Dimensionless(self.P*self.v/self.T/self.R.kJkgK/1000)
            self.cp=unidades.SpecificHeat(propiedades["cp"], "kJkgK")
            self.cv=unidades.SpecificHeat(propiedades["cv"], "kJkgK")
            self.w=unidades.Speed(propiedades["w"])
            self.alfap=unidades.PressureTemperature(propiedades["alfap"], "kPaK")
            self.betap=unidades.Density(propiedades["betap"])

            self.cp0=self._Cp0(self._constants["cp"])
            self.gamma=unidades.Dimensionless(-self.v/self.P.kPa*self.derivative("P", "v", "s"))
            self.fi=unidades.Dimensionless(propiedades["fugacity"])
            self.f=unidades.Pressure(self.fi*self.P)
            self.virialB=unidades.SpecificVolume(propiedades["B"]/self.rhoc)
            self.virialC=unidades.SpecificVolume_square(propiedades["C"]/self.rhoc**2)

            self.cp_cv=unidades.Dimensionless(self.cp/self.cv)
            self.cp0_cv=unidades.Dimensionless(self.cp0/self.cv)
            self.alfav=unidades.InvTemperature(self.derivative("v", "T", "P")/self.v)
            self.kappa=unidades.InvPressure(-self.derivative("v", "P", "T")/self.v)
            self.joule=unidades.TemperaturePressure(self.derivative("T", "P", "h"))
            self.betas=unidades.TemperaturePressure(self.derivative("T", "P", "s"))
            self.n=unidades.Dimensionless(-self.v/self.P*self.derivative("P", "v", "s"))
            self.kt=unidades.Dimensionless(-self.v/self.P*self.derivative("P", "v", "T"))
            self.ks=unidades.InvPressure(-self.derivative("v", "P", "s")/self.v)
            self.Kt=unidades.Pressure(-self.v*self.derivative("P", "v", "s"))
            self.Ks=unidades.Pressure(-self.v*self.derivative("P", "v", "T"))
            self.Gruneisen=unidades.Dimensionless(self.v/self.cv*self.derivative("P", "T", "v"))
            self.dhdp_rho=unidades.EnthalpyPressure(self.derivative("h", "P", "v"))
            self.dhdT_rho=unidades.SpecificHeat(self.derivative("h", "T", "v"))
            self.dhdT_P=unidades.SpecificHeat(self.derivative("h", "T", "P"))
            self.dhdP_T=unidades.EnthalpyPressure(self.derivative("h", "P", "T"))
#            self.dpdT=unidades.PressureTemperature(propiedades["dpdT"])
            self.dpdT=unidades.PressureTemperature(self.derivative("P", "T", "rho")*1e-3)

            self.invT=unidades.InvTemperature(-1/self.T)
            self.Z_rho=unidades.SpecificVolume((self.Z-1)/self.rho)
            self.InternalPressure=unidades.Pressure(self.T*self.derivative("P", "T", "rho")*1e-3-self.P)

#            print propiedades["B"]/self.rho/self.M, propiedades["C"]/self.rho**2/self.M**2
#            print -self.rho*self.cp*self.derivative("P", "rho", "T")/self.derivative("P", "T", "rho"), self.rho*self.derivative("h", "rho", "P")
#            print self.derivative("h", "rho", "T")

#            print self.gamma, self.cp0_cv
#            print self.virialB, self.virialC
#            print self.gamma, self.cp_cv
#            print propiedades["dpdT"], self.derivative("P", "T", "v")

            if recursion:
                if self.Tt<self.T<self.Tc and 0<self.P.atm<self.Pc.atm:

                    rhol, rhov, Ps=self._saturation()
                    self.Liquido=fluid()
                    self.Gas=fluid()

                    vapor=self._eq(rhov, self.T)
                    liquido=self._eq(rhol, self.T)
                    self.Hvap=unidades.Enthalpy(vapor["h"]-liquido["h"], "kJkg")

                    if self.rho>rhol:
                        self.x=unidades.Dimensionless(0)
                    elif self.rho<rhov:
                        self.x=unidades.Dimensionless(1)
                    else:
                        self.x=unidades.Dimensionless((self.s-liquido["s"])/(vapor["s"]-liquido["s"]))

                    if self.x<1:
                        self.Liquido=self.__class__(T=T, rho=rhol, recursion=False)
                    if self.x>0:
                        self.Gas=self.__class__(T=T, rho=rhov, recursion=False)

                else:
                    self.Gas=self.__class__(T=self.T, rho=self.rho, recursion=False)
                    self.Liquido=fluid()
                    self.x=unidades.Dimensionless(1)
                    self.Hvap=unidades.Enthalpy(None)

            self.mu=unidades.Viscosity(self._Viscosity())
            self.k=unidades.ThermalConductivity(self._ThCond())
            self.sigma=unidades.Tension(self._Surface())
            self.epsilon=self._Dielectric()

            if self.mu:
                self.nu=unidades.Diffusivity(self.mu/self.rho)
            else:
                self.nu=unidades.Diffusivity(0)
            if self.mu and self.k:
                self.Prandt=unidades.Dimensionless(self.mu*self.cp/self.k)
            else:
                self.Prandt=unidades.Dimensionless(0)
            if self.k:
                self.alfa=unidades.Diffusivity(self.k/self.rho/self.cp)
            else:
                self.alfa=unidades.Diffusivity(0)



    def _saturation(self, T=None):
        """Akasaka (2008) "A Reliable and Useful Method to Determine the Saturation State from Helmholtz Energy Equations of State", Journal of Thermal Science and Technology, 3, 442-451
        http://dx.doi.org/10.1299/jtst.3.442"""
        if not T:
            T=self.T

        Ps=self._Vapor_Pressure(T)
        rhoL=self._Liquid_Density(T)
        rhoG=self._Vapor_Density(T, Ps)
        g=1000.
        erroro=1e6
        rholo=rhoL
        rhogo=rhoG
        while True:
            deltaL=rhoL/self.rhoc
            deltaG=rhoG/self.rhoc
            liquido=self._eq(rhoL, T)
            vapor=self._eq(rhoG, T)
            Jl=deltaL*(1+deltaL*liquido["fird"])
            Jv=deltaG*(1+deltaG*vapor["fird"])
            Kl=deltaL*liquido["fird"]+liquido["fir"]+log(deltaL)
            Kv=deltaG*vapor["fird"]+vapor["fir"]+log(deltaG)
            Jdl=1+2*deltaL*liquido["fird"]+deltaL**2*liquido["firdd"]
            Jdv=1+2*deltaG*vapor["fird"]+deltaG**2*vapor["firdd"]
            Kdl=2*liquido["fird"]+deltaL*liquido["firdd"]+1/deltaL
            Kdv=2*vapor["fird"]+deltaG*vapor["firdd"]+1/deltaG
            Delta=Jdv*Kdl-Jdl*Kdv
            error=abs(Kv-Kl)+abs(Jv-Jl)

            if error<1e-12:
                break
            elif error>erroro:
                rhoL=rholo
                rhoG=rhogo
                g=g/2.
            else:
                erroro=error
                rholo=rhoL
                rhogo=rhoG
                rhoL=rhoL+g/Delta*((Kv-Kl)*Jdv-(Jv-Jl)*Kdv)
                rhoG=rhoG+g/Delta*((Kv-Kl)*Jdl-(Jv-Jl)*Kdl)
        Ps=self.R*T*rhoL*rhoG/(rhoL-rhoG)*(liquido["fir"]-vapor["fir"]+log(deltaL/deltaG))

        return rhoL, rhoG, Ps

    def _Helmholtz(self, rho, T):
        """Implementación general de la ecuación de estado Setzmann-Wagner, ecuación de estado de multiparámetros basada en la energía libre de Helmholtz"""
        rhoc=self._constants.get("rhoref", self.rhoc)
        Tc=self._constants.get("Tref", self.Tc)
        delta=rho/rhoc
        tau=Tc/T

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta, Tref, Pref)

        fir, firt, firtt, fird, firdd, firdt, firdtt, B, C=self._phir(tau, delta)

        propiedades={}
        propiedades["fir"]=fir
        propiedades["fird"]=fird
        propiedades["firdd"]=firdd

        propiedades["T"]=T
        propiedades["P"]=(1+delta*fird)*self.R.JkgK*T*rho
        propiedades["v"]=1./rho
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


    def _ECS(self,  rho, T):
        rhoc=self._constants.get("rhoref", self.rhoc)
        Tc=self._constants.get("Tref", self.Tc)
        delta=rho/rhoc
        tau=Tc/T

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta, Tref, Pref)

        ref=self._constants["ref"](eq=self._constants["eq"])
        Tr=T/Tc
        rhor=rho/rhoc

        psi=1+(self.f_acent-ref.f_acent)*(self._constants["ft"][0]+self._constants["ft"][1]*log(Tr))
        for n, m in zip(self._constants["ft_add"], self._constants["ft_add_exp"]):
            psi+=n*Tr**m
        for n, m in zip(self._constants["fd"], self._constants["fd_exp"]):
            psi+=n*rhor**m
        T0=T*ref.Tc/self.Tc/psi

        phi=ref.Zc/self.Zc*(1+(self.f_acent-ref.f_acent)*(self._constants["ht"][0]+self._constants["ht"][1]*log(Tr)))
        for n, m in zip(self._constants["ht_add"], self._constants["ht_add_exp"]):
            phi+=n*Tr**m
        for n, m in zip(self._constants["hd"], self._constants["hd_exp"]):
            phi+=n*rhor**m
        rho0=rho*ref.rhoc/self.rhoc*phi

        ref_rhoc=ref._constants.get("rhoref", ref.rhoc)
        ref_Tc=ref._constants.get("Tref", ref.Tc)
        deltaref=rho0/ref_rhoc
        tauref=ref_Tc/T0
        fir, firt, firtt, fird, firdd, firdt, firdtt, B, C=ref._phir(tauref, deltaref)

        propiedades={}
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
        """Implementación general de la ecuación de estado de multiparámetros de Benedict-Webb-Rubin"""
        rho=rho/self.M
        delta=rho/self.rhoc
        tau=self.Tc/T

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta, Tref, Pref)

        a=[None]
        a.append(self.R*T)
        a.append(self._constants["b"][1]*T+self._constants["b"][2]*T**0.5+self._constants["b"][3]+self._constants["b"][4]/T+self._constants["b"][5]/T**2)
        a.append(self._constants["b"][6]*T+self._constants["b"][7]+self._constants["b"][8]/T+self._constants["b"][9]/T**2)
        a.append(self._constants["b"][10]*T+self._constants["b"][11]+self._constants["b"][12]/T)
        a.append(self._constants["b"][13])
        a.append(self._constants["b"][14]/T+self._constants["b"][15]/T**2)
        a.append(self._constants["b"][16]/T)
        a.append(self._constants["b"][17]/T+self._constants["b"][18]/T**2)
        a.append(self._constants["b"][19]/T**2)
        a.append(self._constants["b"][20]/T**2+self._constants["b"][21]/T**3)
        a.append(self._constants["b"][22]/T**2+self._constants["b"][23]/T**4)
        a.append(self._constants["b"][24]/T**2+self._constants["b"][25]/T**3)
        a.append(self._constants["b"][26]/T**2+self._constants["b"][27]/T**4)
        a.append(self._constants["b"][28]/T**2+self._constants["b"][29]/T**3)
        a.append(self._constants["b"][30]/T**2+self._constants["b"][31]/T**3+self._constants["b"][32]/T**4)

        P=sum([a[n]*rho**n for n in range(1, 10)])
        P+=exp(-delta**2)*sum([a[n]*rho**(2*n-17) for n in range(10, 16)])
        P=P*10

        dPdrho=sum([a[n]*n*rho**(n-1) for n in range(1, 10)])
        dPdrho+=exp(-delta**2)*sum([(2*n-17-2*delta**2)*a[n]*rho**(2*n-18) for n in range(10, 16)])
        dPdrho=dPdrho*100

        d2Prho=sum([a[n]*n*(n-1)*rho**(n-2) for n in range(1, 10)])
        d2Prho+=exp(-delta**2)*sum([(-35*n+2*n**2+153+33*delta**2+2*delta**4-4*n*delta**2)*2*a[n]*rho**(2*n-19) for n in range(10, 16)])
        d2Prho=d2Prho*100

        A=0
        for n in range(2, 10):
            A+=a[n]/(n-1.)*rho**(n-1)

        A-=0.5*a[10]*self.rhoc**2*(exp(-delta**2)-1)
        A-=0.5*a[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        A-=0.5*a[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        A-=0.5*a[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        A-=0.5*a[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        A-=0.5*a[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        A=A*100

        adT=[None, self.R]
        adT.append(self._constants["b"][1]+self._constants["b"][2]/2/T**0.5-self._constants["b"][4]/T**2-2*self._constants["b"][5]/T**3)
        adT.append(self._constants["b"][6]-self._constants["b"][8]/T**2-2*self._constants["b"][9]/T**3)
        adT.append(self._constants["b"][10]-self._constants["b"][12]/T**2)
        adT.append(0)
        adT.append(-self._constants["b"][14]/T**2-2*self._constants["b"][15]/T**3)
        adT.append(-self._constants["b"][16]/T**2)
        adT.append(-self._constants["b"][17]/T**2-2*self._constants["b"][18]/T**3)
        adT.append(-2*self._constants["b"][19]/T**3)
        adT.append(-2*self._constants["b"][20]/T**3-3*self._constants["b"][21]/T**4)
        adT.append(-2*self._constants["b"][22]/T**3-4*self._constants["b"][23]/T**5)
        adT.append(-2*self._constants["b"][24]/T**3-3*self._constants["b"][25]/T**4)
        adT.append(-2*self._constants["b"][26]/T**3-4*self._constants["b"][27]/T**5)
        adT.append(-2*self._constants["b"][28]/T**3-3*self._constants["b"][29]/T**4)
        adT.append(-2*self._constants["b"][30]/T**3-3*self._constants["b"][31]/T**4-4*self._constants["b"][32]/T**5)

        dA=0
        for n in range(2, 10):
            dA+=adT[n]/(n-1.)*rho**(n-1)

        dA-=0.5*adT[10]*self.rhoc**2*(exp(-delta**2)-1)
        dA-=0.5*adT[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        dA-=0.5*adT[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        dA-=0.5*adT[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        dA-=0.5*adT[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        dA-=0.5*adT[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        dA=dA*100

        dPdT=sum([adT[n]*rho**n for n in range(1, 10)])
        dPdT+=exp(-delta**2)*sum([adT[n]*rho**(2*n-17) for n in range(10, 16)])
        dPdT=dPdT*100

        adTT=[None, 0]
        adTT.append(-self._constants["b"][2]/4/T**1.5+2*self._constants["b"][4]/T**3+6*self._constants["b"][5]/T**4)
        adTT.append(2*self._constants["b"][8]/T**3+6*self._constants["b"][9]/T**4)
        adTT.append(2*self._constants["b"][12]/T**3)
        adTT.append(0)
        adTT.append(2*self._constants["b"][14]/T**3+6*self._constants["b"][15]/T**4)
        adTT.append(2*self._constants["b"][16]/T**3)
        adTT.append(2*self._constants["b"][17]/T**3+6*self._constants["b"][18]/T**4)
        adTT.append(6*self._constants["b"][19]/T**4)
        adTT.append(6*self._constants["b"][20]/T**4+12*self._constants["b"][21]/T**5)
        adTT.append(6*self._constants["b"][22]/T**4+20*self._constants["b"][23]/T**6)
        adTT.append(6*self._constants["b"][24]/T**4+12*self._constants["b"][25]/T**5)
        adTT.append(6*self._constants["b"][26]/T**4+20*self._constants["b"][27]/T**6)
        adTT.append(6*self._constants["b"][28]/T**4+12*self._constants["b"][29]/T**5)
        adTT.append(6*self._constants["b"][30]/T**4+12*self._constants["b"][31]/T**5+20*self._constants["b"][32]/T**6)

        d2A=0
        for n in range(2, 10):
            d2A+=adTT[n]*rho**(n-1)/(n-1)

        d2A-=0.5*adTT[10]*self.rhoc**2*(exp(-delta**2)-1)
        d2A-=0.5*adTT[11]*self.rhoc**4*(exp(-delta**2)*(delta**2+1)-1)
        d2A-=0.5*adTT[12]*self.rhoc**6*(exp(-delta**2)*(delta**4+2*delta**2+2)-2)
        d2A-=0.5*adTT[13]*self.rhoc**8*(exp(-delta**2)*(delta**6+3*delta**4+6*delta**2+6)-6)
        d2A-=0.5*adTT[14]*self.rhoc**10*(exp(-delta**2)*(delta**8+4*delta**6+12*delta**4+24*delta**2+24)-24)
        d2A-=0.5*adTT[15]*self.rhoc**12*(exp(-delta**2)*(delta**10+5*delta**8+20*delta**6+60*delta**4+120*delta**2+120)-120)
        d2A=d2A*100

        B, C=0, 0

        fir=A/self.R/T
        firdt=P/T-dPdT/self.R/rho
        firt=A/T-dA/self.R
        firtt=d2A*T/self.R
        fird=P/self.R/T/rho-1.
        firdd=(dPdrho-2*P/rho)/self.R/T+1.
        firddd=(d2Prho*rho-4*dPdrho+6*P/rho)/self.R/T-2.

        propiedades={}
        propiedades["T"]=T
        propiedades["P"]=P*0.1      #convertido de bar a MPa
        propiedades["v"]=1/rho
        propiedades["h"]=self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"]=self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cp"]=self.R.kJkgK*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["cv"]=-self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["w"]=(self.R*T*(1+2*delta*fird+delta**2*firdd-(1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"]=(1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"]=rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
        propiedades["fugacity"]=exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"]=B
        propiedades["C"]=C
        propiedades["fir"]=fir
        propiedades["fird"]=fird
        return propiedades

    def _PengRobinson(self, rho, T):
        """Peng, D.-Y.; Robinson, D.B. A New Two-Constant Equation of State. I&EC Fundam. 1976, 15(1), 59
        Peneloux, A.; Rauzy, E.; Freze, R. A consistent correction for Redlich-Kwong-Soave volumes. Fluid Phase Eq. 1982, 8, 7.
        http://dx.doi.org/10.1021/i160057a011
        http://dx.doi.org/10.1016/0378-3812(82)80002-2"""
        #FIXME: no sale por poco
        delta=rho/self.rhoc
        tau=self.Tc/T
        delta_0=1e-50
        rho=rho/self.M

        Tr=1./tau
        m=0.37464+1.54226*self.f_acent-0.26992*self.f_acent**2
        alfa=(1+m*(1-Tr**0.5))**2
        a=0.457235*R_atml**2*self.Tc**2/self.Pc.atm
        b=0.077796*R_atml*self.Tc/self.Pc.atm

        daT=-a*m/T**0.5/self.Tc**0.5*alfa**0.5
        d2aT=a*m*(1+m)/2/T/(T*self.Tc)**0.5

        fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta, Tref, Pref)

        v=1./rho
        vb=v-b
        q=b*8**0.5
        v1=2*v+2*b+q
        v2=2*v+2*b-q
        v1n=v1+2*self._PR
        v2n=v2+2*self._PR
        fir=log(v)-log(vb+self._PR)+a/self.R.kJkgK/T*log(v2n/v1n)/q

        phipart=(1./self.R.kJkgK/q)*(daT/T-a/T**2)
        dtdtau=-T**2/self.Tc
        dphidtau=phipart*dtdtau
        term1=2+2*b*rho-rho*q+2*rho*self._PR
        term2=2+2*b*rho+rho*q+2*rho*self._PR
        phidpart=-4*self.rhoc/self.M*q/term1/term2
        firdt=dphidtau*phidpart

        fird=-1+1./(1-b*rho+rho*self._PR)+a/self.R.kJkgK/T/q*(2/v1n-2/v2n)/rho
        bdt=1-b*rho+rho*self._PR
        firdd=1-1./bdt-(self._PR*self.rhoc/self.M-b*self.rhoc/self.M)*delta/bdt**2 + 4*a/self.R.kJkgK/T/q/rho*((1./v2n-1./v1n)+(1./v1n**2-1./v2n**2)/rho)
        dphidt=1./self.R.kJkgK/q*log(v2n/v1n)*(daT/T-a/T**2)
        firt=dphidt*dtdtau
        d2phidt2=(1./self.R.kJkgK/q)*log(v2n/v1n)*(d2aT/T-2*daT/T**2+2*a/T**3)
        d2phid2tau=1./tau**2*(T**2*d2phidt2+2*T*dphidt)
        firtt=d2phid2tau

#        fir=-log(1-b*rho)+a/(self.R*T*2**1.5*b)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        fird=b/(1-b*rho)+a/(self.R*T*2**1.5*b)*((1-2**0.5)*b/(1+(1-2**0.5)*b*rho)-(1+2**0.5)*b/(1+(1+2**0.5)*b*rho))
#        firdd=(b/(1-b*rho))**2+a/(rho**2*self.R*T*2**1.5*b)*(-((1-2**0.5)*b*rho/(1+(1-2**0.5)*b*rho))**2+((1+2**0.5)*b*rho/(1+(1+2**0.5)*b*rho))**2)
#        firt=1/(self.R*T*2**1.5*b)*(daT-a/T)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        firtt=1/(self.R*2**1.5*b*T)*(d2aT-2/T*daT+2*a/T**2)*log((1+(1-2**0.5)*b*rho)/(1+(1+2**0.5)*b*rho))
#        firdt=1/(rho*T*self.R*2**1.5*b)*(daT-a/T)*((1-2**0.5)*b*rho/(1+(1-2**0.5)*b*rho)-(1+2**0.5)*b*rho/(1+(1+2**0.5)*b*rho))
        firdtt=B=C=0

        propiedades={}
        propiedades["fir"]=fir
        propiedades["fird"]=fird
        propiedades["firdd"]=firdd

        propiedades["T"]=T
        propiedades["P"]=(1+delta*fird)*self.R.JkgK*T*rho
        propiedades["v"]=1/rho
        propiedades["h"]=self.R.kJkgK*T*(1+tau*(fiot+firt)+delta*fird)
        propiedades["s"]=self.R.kJkgK*(tau*(fiot+firt)-fio-fir)
        propiedades["cp"]=self.R.kJkgK*(-tau**2*(fiott+firtt)+(1+delta*fird-delta*tau*firdt)**2/(1+2*delta*fird+delta**2*firdd))
        propiedades["cv"]=-self.R.kJkgK*tau**2*(fiott+firtt)
        propiedades["w"]=(self.R*T*(1+2*delta*fird+delta**2*firdd-(1+delta*fird-delta*tau*firdt)**2/tau**2/(fiott+firtt)))**0.5
        propiedades["alfap"]=(1-delta*tau*firdt/(1+delta*fird))/T
        propiedades["betap"]=rho*(1+(delta*fird+delta**2*firdd)/(1+delta*fird))
        propiedades["fugacity"]=exp(fir+delta*fird-log(1+delta*fird))
        propiedades["B"]=B
        propiedades["C"]=C
        propiedades["dpdrho"]=self.R*T*(1+2*delta*fird+delta**2*firdd)
        propiedades["dpdT"]=self.R*rho*(1+delta*fird+delta*tau*firdt)

        return propiedades

    def _Generalised(self):
        """Span R., Wagner W., "An accurate empirical three parameter equation of state for nonpolar fluids". To be submitted to Fluid Phase Equilibria. 2000 """
        if self._Tr:
            Tref=self._Tr
        else:
            Tref=self.Tc
        if self._rhor:
            rhoref=self._rhor
        else:
            rhoref=self.rhoc
        if self._w:
            w=self._w
        else:
            w=self.f_acent

        helmholtz={
            "R": 8.31451,
            "Tref": Tref,
            "rhoref": rhoref,
            "cp": self.eq[0]["cp"],

            "d1": [1, 1, 2, 3, 8],
            "t1": [0.125, 1.125, 1.25, 0.25, 0.75],

            "d2": [2, 3, 1, 4, 3],
            "t2": [0.625, 2, 4.125, 4.125, 17],
            "c2": [1, 1, 2, 2, 3],
            "gamma2": [1]*5}

        c1=[0.636479524, -0.174667493e1, -0.144442644e-1, 0.6799731e-1, 0.767320032e-4, 0.218194143, 0.810318494e-1, -0.907368899e-1, 0.25312225e-1, -0.209937023e-1]
        c2=[0.82247342, -0.954932692, -0.745462328, 0.182685593, 0.547120142e-4, 0.761697913, 0.415691324, -0.825206373, -0.240558288, -0.643818403e-1]
        c3=[-0.186193063e1, 0.105083555e2, 0.16403233e1, -0.613747797, -0.69318829e-3, -0.705727791e1, -0.290006245e1, -0.232497527, -0.282346515, 0.254250643e1]
        nr=[c1[i]+c2[i]*w+c3[i]*w**4 for i in range(10)]
        helmholtz["nr1"]=nr[:5]
        helmholtz["nr2"]=nr[5:]
        self._constants=helmholtz


    def _phi0(self, cp, tau, delta, To, Po):
        R=cp.get("R", self._constants["R"])/self.M*1000
        rhoc=self._constants.get("rhoref", self.rhoc)
        Tc=self._constants.get("Tref", self.Tc)
        rho0=Po/R/To
        tau0=Tc/To
        delta0=rho0/rhoc
        co=cp["ao"]-1
        ti=[-x for x in cp["pow"]]
        ci=[-n/(t*(t+1))*Tc**t for n, t in zip(cp["an"], cp["pow"])]
        titao=[fi/Tc for fi in cp["exp"]]
        hyp=[fi/Tc for fi in cp["hyp"]]
        cI=-(1+co)/tau0
        cII=co*(1-log(tau0))-log(delta0)
#        for c, t in zip(ci, ti):
#            cI-=c*t*tau0**(t-1)
#            cII+=c*(t-1)*tau0**t
#        for ao, tita in zip(cp["ao_exp"], titao):
#            cI-=ao*tita*(1/(1-exp(-tita*tau0))-1)
#            cII+=ao*tita*(tau0*(1/(1-exp(-tita*tau0))-1)-log(1-exp(-tita*tau0)))
#        if cp["ao_hyp"]:
#            for i in [0, 2]:
#                cI-=cp["ao_hyp"][i]*hyp[i]/(tanh(hyp[i]*tau0))
#                cII+=cp["ao_hyp"][i]*(hyp[i]*tau0/tanh(hyp[i]*tau0)-log(abs(sinh(hyp[i]*tau0))))
#            for i in [1, 3]:
#                cI+=cp["ao_hyp"][i]*hyp[i]*tanh(hyp[i]*tau0)
#                cII-=cp["ao_hyp"][i]*(hyp[i]*tau0*tanh(hyp[i]*tau0)-log(abs(cosh(hyp[i]*tau0))))

        Fi0={"ao_log": [1,  co],
                "pow": [0, 1] + ti,
                "ao_pow": [cII, cI] + ci,
                "ao_exp": cp["ao_exp"],
                "titao": titao,
                "ao_hyp": cp["ao_hyp"],
                "hyp": hyp}

        #FIXME: Reference estate
        T=self._constants.get("Tref", self.Tc)/tau
        rho=delta*self.rhoc
        cp0sav=self._Cp0(cp, T)
        cpisav=self._dCp(cp, T, To)
        cptsav=self._dCpT(cp, T, To)
        fio=cpisav/R/T-cptsav/R+log(T*rho/rho0/Tref)-1
        fiot=(cpisav-1)/tau
        fiott=(1-cp0sav/R)/tau**2

#        fio=Fi0["ao_log"][0]*log(delta)+Fi0["ao_log"][1]*log(tau)
#        fiot=+Fi0["ao_log"][1]/tau
#        fiott=-Fi0["ao_log"][1]/tau**2

        fiod=1/delta
        fiodd=-1/delta**2
        fiodt=0

        for n, t in zip(Fi0["ao_pow"], Fi0["pow"]):
            fio+=n*tau**t
            if t!=0:
                fiot+=t*n*tau**(t-1)
            if t not in [0, 1]:
                fiott+=n*t*(t-1)*tau**(t-2)

        for i in range(len(Fi0["ao_exp"])):
            fio+=Fi0["ao_exp"][i]*log(1-exp(-tau*Fi0["titao"][i]))
            fiot+=Fi0["ao_exp"][i]*Fi0["titao"][i]*((1-exp(-Fi0["titao"][i]*tau))**-1-1)
            fiott-=Fi0["ao_exp"][i]*Fi0["titao"][i]**2*exp(-Fi0["titao"][i]*tau)*(1-exp(-Fi0["titao"][i]*tau))**-2

        if Fi0["ao_hyp"]:
            for i in [0, 2]:
                fio+=Fi0["ao_hyp"][i]*log(abs(sinh(Fi0["hyp"][i]*tau)))
                fiot+=Fi0["ao_hyp"][i]*Fi0["hyp"][i]/tanh(Fi0["hyp"][i]*tau)
                fiott-=Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/sinh(Fi0["hyp"][i]*tau)**2

            for i in [1, 3]:
                fio-=Fi0["ao_hyp"][i]*log(abs(cosh(Fi0["hyp"][i]*tau)))
                fiot-=Fi0["ao_hyp"][i]*Fi0["hyp"][i]*tanh(Fi0["hyp"][i]*tau)
                fiott-=Fi0["ao_hyp"][i]*Fi0["hyp"][i]**2/cosh(Fi0["hyp"][i]*tau)**2

        R_=cp.get("R", self._constants["R"])
        factor=R_/self._constants["R"]
        return factor*fio, factor*fiot, factor*fiott, factor*fiod, factor*fiodd, factor*fiodt


    def _phir(self, tau, delta):
        delta_0=1e-50

        fir=fird=firdd=firt=firtt=firdt=firdtt=B=C=0
        for i in range(len(self._constants.get("nr1", []))):        #Polinomial terms
            fir+=self._constants["nr1"][i]*delta**self._constants["d1"][i]*tau**self._constants["t1"][i]
            fird+=self._constants["nr1"][i]*self._constants["d1"][i]*delta**(self._constants["d1"][i]-1)*tau**self._constants["t1"][i]
            firdd+=self._constants["nr1"][i]*self._constants["d1"][i]*(self._constants["d1"][i]-1)*delta**(self._constants["d1"][i]-2)*tau**self._constants["t1"][i]
            firt+=self._constants["nr1"][i]*self._constants["t1"][i]*delta**self._constants["d1"][i]*tau**(self._constants["t1"][i]-1)
            firtt+=self._constants["nr1"][i]*self._constants["t1"][i]*(self._constants["t1"][i]-1)*delta**self._constants["d1"][i]*tau**(self._constants["t1"][i]-2)
            firdt+=self._constants["nr1"][i]*self._constants["t1"][i]*self._constants["d1"][i]*delta**(self._constants["d1"][i]-1)*tau**(self._constants["t1"][i]-1)
            firdtt+=self._constants["nr1"][i]*self._constants["t1"][i]*self._constants["d1"][i]*(self._constants["t1"][i]-1)*delta**(self._constants["d1"][i]-1)*tau**(self._constants["t1"][i]-2)
            B+=self._constants["nr1"][i]*self._constants["d1"][i]*delta_0**(self._constants["d1"][i]-1)*tau**self._constants["t1"][i]
            C+=self._constants["nr1"][i]*self._constants["d1"][i]*(self._constants["d1"][i]-1)*delta_0**(self._constants["d1"][i]-2)*tau**self._constants["t1"][i]

        for i in range(len(self._constants.get("nr2", []))):    #Exponential terms
            fir+=self._constants["nr2"][i]*delta**self._constants["d2"][i]*tau**self._constants["t2"][i]*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])
            fird+=self._constants["nr2"][i]*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])*delta**(self._constants["d2"][i]-1)*tau**self._constants["t2"][i]*(self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta**self._constants["c2"][i])
            firdd+=self._constants["nr2"][i]*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])*delta**(self._constants["d2"][i]-2)*tau**self._constants["t2"][i]*((self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta**self._constants["c2"][i])*(self._constants["d2"][i]-1-self._constants["gamma2"][i]*self._constants["c2"][i]*delta**self._constants["c2"][i])-self._constants["gamma2"][i]**2*self._constants["c2"][i]**2*delta**self._constants["c2"][i])
            firt+=self._constants["nr2"][i]*self._constants["t2"][i]*delta**self._constants["d2"][i]*tau**(self._constants["t2"][i]-1)*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])
            firtt+=self._constants["nr2"][i]*self._constants["t2"][i]*(self._constants["t2"][i]-1)*delta**self._constants["d2"][i]*tau**(self._constants["t2"][i]-2)*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])
            firdt+=self._constants["nr2"][i]*self._constants["t2"][i]*delta**(self._constants["d2"][i]-1)*tau**(self._constants["t2"][i]-1)*(self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta**self._constants["c2"][i])*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])
            firdtt+=self._constants["nr2"][i]*self._constants["t2"][i]*(self._constants["t2"][i]-1)*delta**(self._constants["d2"][i]-1)*tau**(self._constants["t2"][i]-2)*(self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta**self._constants["c2"][i])*exp(-self._constants["gamma2"][i]*delta**self._constants["c2"][i])
            B+=self._constants["nr2"][i]*exp(-self._constants["gamma2"][i]*delta_0**self._constants["c2"][i])*delta_0**(self._constants["d2"][i]-1)*tau**self._constants["t2"][i]*(self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta_0**self._constants["c2"][i])
            C+=self._constants["nr2"][i]*exp(-self._constants["gamma2"][i]*delta_0**self._constants["c2"][i])*(delta_0**(self._constants["d2"][i]-2)*tau**self._constants["t2"][i]*((self._constants["d2"][i]-self._constants["gamma2"][i]*self._constants["c2"][i]*delta_0**self._constants["c2"][i])*(self._constants["d2"][i]-1-self._constants["gamma2"][i]*self._constants["c2"][i]*delta_0**self._constants["c2"][i])-self._constants["gamma2"][i]**2*self._constants["c2"][i]**2*delta_0**self._constants["c2"][i]))

        for i in range(len(self._constants.get("nr3", []))):    #Gaussian terms
            exp1=self._constants.get("exp1", [2]*len(self._constants["nr3"]))
            exp2=self._constants.get("exp2", [2]*len(self._constants["nr3"]))
            fir+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])
            fird+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(self._constants["d3"][i]/delta-2*self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i]))
            firdd+=self._constants["nr3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(-2*self._constants["alfa3"][i]*delta**self._constants["d3"][i]+4*self._constants["alfa3"][i]**2*delta**self._constants["d3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-4*self._constants["d3"][i]*self._constants["alfa3"][i]*delta**2*(delta-self._constants["epsilon3"][i])+self._constants["d3"][i]*2*delta)
            firt+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(self._constants["t3"][i]/tau-2*self._constants["beta3"][i]*(tau-self._constants["gamma3"][i]))
            firtt+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*((self._constants["t3"][i]/tau-2*self._constants["beta3"][i]*(tau-self._constants["gamma3"][i]))**exp2[i]-self._constants["t3"][i]/tau**2-2*self._constants["beta3"][i])
            firdt+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(self._constants["t3"][i]/tau-2*self._constants["beta3"][i]*(tau-self._constants["gamma3"][i]))*(self._constants["d3"][i]/delta-2*self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i]))
            firdtt+=self._constants["nr3"][i]*delta**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*((self._constants["t3"][i]/tau-2*self._constants["beta3"][i]*(tau-self._constants["gamma3"][i]))**exp2[i]-self._constants["t3"][i]/tau**2-2*self._constants["beta3"][i])*(self._constants["d3"][i]/delta-2*self._constants["alfa3"][i]*(delta-self._constants["epsilon3"][i]))
            B+=self._constants["nr3"][i]*delta_0**self._constants["d3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta_0-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(self._constants["d3"][i]/delta_0-2*self._constants["alfa3"][i]*(delta_0-self._constants["epsilon3"][i]))
            C+=self._constants["nr3"][i]*tau**self._constants["t3"][i]*exp(-self._constants["alfa3"][i]*(delta_0-self._constants["epsilon3"][i])**exp1[i]-self._constants["beta3"][i]*(tau-self._constants["gamma3"][i])**exp2[i])*(-2*self._constants["alfa3"][i]*delta_0**self._constants["d3"][i]+4*self._constants["alfa3"][i]**2*delta_0**self._constants["d3"][i]*(delta_0-self._constants["epsilon3"][i])**exp1[i]-4*self._constants["d3"][i]*self._constants["alfa3"][i]*delta_0**2*(delta_0-self._constants["epsilon3"][i])+self._constants["d3"][i]*2*delta_0)

        for i in range(len(self._constants.get("nr4", []))):    #Non analitic terms
            Tita=(1-tau)+self._constants["A"][i]*((delta-1)**2)**(1/2/self._constants["beta4"][i])
            F=exp(-self._constants["C"][i]*(delta-1)**2-self._constants["D"][i]*(tau-1)**2)
            Fd=-2*self._constants["C"][i]*F*(delta-1)
            Fdd=2*self._constants["C"][i]*F*(2*self._constants["C"][i]*(delta-1)**2-1)
            Ft=-2*self._constants["D"][i]*F*(tau-1)
            Ftt=2*self._constants["D"][i]*F*(2*self._constants["D"][i]*(tau-1)**2-1)
            Fdt=4*self._constants["C"][i]*self._constants["D"][i]*F*(delta-1)*(tau-1)
            Fdtt=4*self._constants["C"][i]*self._constants["D"][i]*F*(delta-1)*(2*self._constants["D"][i]*(tau-1)**2-1)

            Delta=Tita**2+self._constants["B"][i]*((delta-1)**2)**self._constants["a4"][i]
            Deltad=(delta-1)*(self._constants["A"][i]*Tita*2/self._constants["beta4"][i]*((delta-1)**2)**(1/2/self._constants["beta4"][i]-1)+2*self._constants["B"][i]*self._constants["a4"][i]*((delta-1)**2)**(self._constants["a4"][i]-1))
            Deltadd=Deltad/(delta-1)+(delta-1)**2*(4*self._constants["B"][i]*self._constants["a4"][i]*(self._constants["a4"][i]-1)*((delta-1)**2)**(self._constants["a4"][i]-2)+2*self._constants["A"][i]**2/self._constants["beta4"][i]**2*(((delta-1)**2)**(1/2/self._constants["beta4"][i]-1))**2+self._constants["A"][i]*Tita*4/self._constants["beta4"][i]*(1/2/self._constants["beta4"][i]-1)*((delta-1)**2)**(1/2/self._constants["beta4"][i]-2))

            DeltaBd=self._constants["b"][i]*Delta**(self._constants["b"][i]-1)*Deltad
            DeltaBdd=self._constants["b"][i]*(Delta**(self._constants["b"][i]-1)*Deltadd+(self._constants["b"][i]-1)*Delta**(self._constants["b"][i]-2)*Deltad**2)
            DeltaBt=-2*Tita*self._constants["b"][i]*Delta**(self._constants["b"][i]-1)
            DeltaBtt=2*self._constants["b"][i]*Delta**(self._constants["b"][i]-1)+4*Tita**2*self._constants["b"][i]*(self._constants["b"][i]-1)*Delta**(self._constants["b"][i]-2)
            DeltaBdt=-self._constants["A"][i]*self._constants["b"][i]*2/self._constants["beta4"][i]*Delta**(self._constants["b"][i]-1)*(delta-1)*((delta-1)**2)**(1/2/self._constants["beta4"][i]-1)-2*Tita*self._constants["b"][i]*(self._constants["b"][i]-1)*Delta**(self._constants["b"][i]-2)*Deltad
            DeltaBdtt=2*self._constants["b"][i]*(self._constants["b"][i]-1)*Delta**(self._constants["b"][i]-2)*(Deltad*(1+2*Tita**2*(self._constants["b"][i]-2)/Delta)+4*Tita*self._constants["A"][i]*(delta-1)/self._constants["beta4"][i]*((delta-1)**2)**(1/2/self._constants["beta4"][i]-1))

            fir+=self._constants["nr4"][i]*Delta**self._constants["b"][i]*delta*F
            fird+=self._constants["nr4"][i]*(Delta**self._constants["b"][i]*(F+delta*Fd)+DeltaBd*delta*F)
            firdd+=self._constants["nr4"][i]*(Delta**self._constants["b"][i]*(2*Fd+delta*Fdd)+2*DeltaBd*(F+delta*Fd)+DeltaBdd*delta*F)
            firt+=self._constants["nr4"][i]*delta*(DeltaBt*F+Delta**self._constants["b"][i]*delta*Ft)
            firtt+=self._constants["nr4"][i]*delta*(DeltaBtt*F+2*DeltaBt*Ft+Delta**self._constants["b"][i]*Ftt)
            firdt+=self._constants["nr4"][i]*(Delta**self._constants["b"][i]*(Ft+delta*Fdt)+delta*DeltaBd*Ft+DeltaBt*(F+delta*Fd)+DeltaBdt*delta*F)
            firdtt+=self._constants["nr4"][i]*((DeltaBtt*F+2*DeltaBt*Ft+Delta**self._constants["b"][i]*Ftt)+delta*(DeltaBdtt*F+DeltaBtt*Fd+2*DeltaBdt*Ft+2*DeltaBt*Fdt+DeltaBt*Ftt+Delta**self._constants["b"][i]*Fdtt))

            Tita_virial=(1-tau)+self._constants["A"][i]*((delta_0-1)**2)**(1/2/self._constants["beta4"][i])
            Delta_Virial=Tita_virial**2+self._constants["B"][i]*((delta_0-1)**2)**self._constants["a4"][i]
            Deltad_Virial=(delta_0-1)*(self._constants["A"][i]*Tita_virial*2/self._constants["beta4"][i]*((delta_0-1)**2)**(1/2/self._constants["beta4"][i]-1)+2*self._constants["B"][i]*self._constants["a4"][i]*((delta_0-1)**2)**(self._constants["a4"][i]-1))
            Deltadd_Virial=Deltad_Virial/(delta_0-1)+(delta_0-1)**2*(4*self._constants["B"][i]*self._constants["a4"][i]*(self._constants["a4"][i]-1)*((delta_0-1)**2)**(self._constants["a4"][i]-2)+2*self._constants["A"][i]**2/self._constants["beta4"][i]**2*(((delta_0-1)**2)**(1/2/self._constants["beta4"][i]-1))**2+self._constants["A"][i]*Tita_virial*4/self._constants["beta4"][i]*(1/2/self._constants["beta4"][i]-1)*((delta_0-1)**2)**(1/2/self._constants["beta4"][i]-2))
            DeltaBd_Virial=self._constants["b"][i]*Delta_Virial**(self._constants["b"][i]-1)*Deltad_Virial
            DeltaBdd_Virial=self._constants["b"][i]*(Delta_Virial**(self._constants["b"][i]-1)*Deltadd_Virial+(self._constants["b"][i]-1)*Delta_Virial**(self._constants["b"][i]-2)*Deltad_Virial**2)
            F_virial=exp(-self._constants["C"][i]*(delta_0-1)**2-self._constants["D"][i]*(tau-1)**2)
            Fd_virial=-2*self._constants["C"][i]*F_virial*(delta_0-1)
            Fdd_virial=2*self._constants["C"][i]*F_virial*(2*self._constants["C"][i]*(delta_0-1)**2-1)

            B+=self._constants["nr4"][i]*(Delta_Virial**self._constants["b"][i]*(F_virial+delta_0*Fd_virial)+DeltaBd_Virial*delta_0*F_virial)
            C+=self._constants["nr4"][i]*(Delta_Virial**self._constants["b"][i]*(2*Fd_virial+delta_0*Fdd_virial)+2*DeltaBd_Virial*(F_virial+delta_0*Fd_virial)+DeltaBdd_Virial*delta_0*F_virial)

        if self._constants.get("Fi", None): #Hard sphere term
            f=self._constants["Fi"]
            n=0.1617
            a=0.689
            gamma=0.3674
            X=n*delta/(a+(1-a)/tau**gamma)
            Xd=n/(a+(1-a)/tau**gamma)
            Xt=n*delta*(1-a)*gamma/tau**(gamma+1)/(a+(1-a)/tau**gamma)**2
            Xdt=n*(1-a)*gamma/tau**(gamma+1)/(a+(1-a)/tau**gamma)**2
            Xtt=-n*delta*((1-a)*gamma/tau**(gamma+2)*((gamma+1)*(a+(1-a)/tau**gamma)-2*gamma*(1-a)/tau**gamma))/(a+(1-a)/tau**gamma)**3
            Xdtt=-n*((1-a)*gamma/tau**(gamma+2)*((gamma+1)*(a+(1-a)/tau**gamma)-2*gamma*(1-a)/tau**gamma))/(a+(1-a)/tau**gamma)**3

            ahdX=-(f**2-1)/(1-X)+(f**2+3*f+X*(f**2-3*f))/(1-X)**3
            ahdXX=-(f**2-1)/(1-X)**2+(3*(f**2+3*f)+(f**2-3*f)*(1+2*X))/(1-X)**4
            ahdXXX=-2*(f**2-1)/(1-X)**3+6*(2*(f**2+3*f)+(f**2-3*f)*(1+X))/(1-X)**5

            fir+=(f**2-1)*log(1-X)+((f**2+3*f)*X-3*f*X**2)/(1-X)**2
            fird+=ahdX*Xd
            firdd+=ahdXX*Xd**2
            firt+=ahdX*Xt
            firtt+=ahdXX*Xt**2+ahdX*Xtt
            firdt+=ahdXX*Xt*Xd+ahdX*Xdt
            firdtt+=ahdXXX*Xt**2*Xd+ahdXX*(Xtt*Xd+2*Xdt*Xt)*ahdX*Xdtt

            X_virial=n*delta_0/(a+(1-a)/tau**gamma)
            ahdX_virial=-(f**2-1)/(1-X_virial)+(f**2+3*f+X_virial*(f**2-3*f))/(1-X_virial)**3
            ahdXX_virial=-(f**2-1)/(1-X_virial)**2+(3*(f**2+3*f)+(f**2-3*f)*(1+2*X_virial))/(1-X_virial)**4
            B+=ahdX_virial*Xd
            C+=ahdXX_virial*Xd**2
        return fir, firt, firtt, fird, firdd, firdt, firdtt, B, C


    def _Cp0(self, cp, T=False):
        if not T:
            T=self.T
        tau=self.Tc/T
        cpo=cp["ao"]
        for a, t in zip(cp["an"], cp["pow"]):
            cpo+=a*T**t
        for m, tita in zip(cp["ao_exp"], cp["exp"]):
            cpo+=m*(tita/T)**2*exp(tita/T)/(1-exp(tita/T))**2
        if cp["ao_hyp"]:
            for i in [0, 2]:
                cpo+=cp["ao_hyp"][i]*(cp["hyp"][i]/T/(sinh(cp["hyp"][i]/T)))**2
            for i in [1, 3]:
                cpo+=cp["ao_hyp"][i]*(cp["hyp"][i]/T/(cosh(cp["hyp"][i]/T)))**2

        return unidades.SpecificHeat(cpo/self.M*1000)

    def _dCp(self, cp, T, Tref):
        """Calcula la integral de Cp0 entre T y Tref, necesario para calcular la entalpia usando estados de referencia"""
        cpsum=cp.get("ao", 0)*(T-Tref)
        for n, t in zip(cp["an"], cp["pow"]):
            t+=1
            if t==0:
                cpsum+=n*log(T/Tref)
            else:
                cpsum+=n*(T**t-Tref**t)/t
        for tita, ao in zip(cp["ao_exp"], cp["exp"]):
            cpsum+=-ao*0.5*tita*((1.+exp(tita/T))/(1.-exp(tita/T))-(1.+exp(tita/Tref))/(1.-exp(tita/Tref)))
#            cpsum+=ao*-0.5*tita*((1.+exp(tita/T))/(1.-exp(tita/T))-(1.+exp(tita/Tref))/(1.-exp(tita/Tref)))

        return cpsum*self.R

    def _dCpT(self, cp, T, Tref):
        """Calcula la integral de Cp0/T entre T y Tref, necesario para calcular la entropia usando estados de referencia"""
        cpsum=cp.get("ao", 0)*log(T/Tref)
        for n, t in zip(cp["an"], cp["pow"]):
            cpsum+=n*(T**t-Tref**t)/t
        for tita, ao in zip(cp["ao_exp"], cp["exp"]):
            cpsum+=ao*(log((1.-exp(tita/Tref))/(1.-exp(tita/T)))+tita/T*exp(tita/T)/(exp(tita/T)-1.)-tita/Tref*exp(tita/Tref)/(exp(tita/Tref)-1.))
        return cpsum*self.R


    def derivative(self, z, x, y):
        """Calculate generic partial derivative: (δz/δx)y where x, y, z can be: P, T, v, rho, u, h, s, g, a"""
        dT={"P": self.P*self.alfap,
                "T": 1,
                "v": 0,
                "rho": 0,
                "u": self.cv,
                "h": self.cv+self.P*self.v*self.alfap,
                "s": self.cv/self.T,
                "g": self.P*self.v*self.alfap-self.s,
                "a": -self.s}
        dv={"P": -self.P*self.betap,
                "T": 0,
                "v": 1,
                "rho": -1,
                "u": self.P*(self.T*self.alfap-1),
                "h": self.P*(self.T*self.alfap-self.v*self.betap),
                "s": self.P*self.alfap,
                "g": -self.P*self.v*self.betap,
                "a": -self.P}
        return (dv[z]*dT[y]-dT[z]*dv[y])/(dv[x]*dT[y]-dT[x]*dv[y])


    def _Dielectric(self):
        """Harvey, A.H. and Lemmon, E.W. "Method for Estimating the Dielectric Constant of Natural Gas Mixtures," Int. J. Thermophys., 26(1):31-46, 2005.
        http://dx.doi.org/10.1007/s10765-005-2351-5"""
        if self._dielectric:
            if self.rho<1e-12:
                e=1.
            else:
                delta=self.rho/self.M/self._dielectric["rhoref"]
                tau=self.T/self._dielectric["Tref"]
                c=0
                for a, expt, expd in zip(self._dielectric["a0"], self._dielectric["expt0"], self._dielectric["expd0"]):
                    c+=a * tau**expt * delta**expd
                for a, expt, expd in zip(self._dielectric["a1"], self._dielectric["expt1"], self._dielectric["expd1"]):
                    c+=a * (tau-1)**expt * delta**expd
                for a, expt, expd in zip(self._dielectric["a2"], self._dielectric["expt2"], self._dielectric["expd2"]):
                    c+=a * (1./tau-1)**expt * delta**expd
                if self._dielectric["eq"]==3:
                    e=(1+2*c)/(1-c)
                else:
                    e=0.25*(1+9*c+3*(9*c**2+2*c+1)**0.5)
        else:
            e=0
        return unidades.Dimensionless(e)

    @classmethod
    def _Melting_Pressure(cls, T):
        if cls._melting:
            Tref=cls._melting["Tref"]
            Pref=cls._melting["Pref"]
            Tita=T/Tref
            suma=0
            for ai, expi in zip(cls._melting["a1"], cls._melting["exp1"]):
                suma+=ai*Tita**expi
            for ai, expi in zip(cls._melting["a2"], cls._melting["exp2"]):
                suma+=ai*(Tita-1)**expi
            for ai, expi in zip(cls._melting["a3"], cls._melting["exp3"]):
                suma+=ai*log(Tita)**expi

            if cls._melting["eq"]==1:
                P=suma*Pref
            elif cls._melting["eq"]==2:
                P=exp(suma)*Pref
            return unidades.Pressure(P, "kPa")

        else:
            return None

    @classmethod
    def _Sublimation_Pressure(cls, T):
        if cls._sublimation:
            Tref=cls._sublimation["Tref"]
            Pref=cls._sublimation["Pref"]
            Tita=T/Tref
            suma=0
            for ai, expi in zip(cls._sublimation["a1"], cls._sublimation["exp1"]):
                suma+=ai*Tita**expi
            for ai, expi in zip(cls._sublimation["a2"], cls._sublimation["exp2"]):
                suma+=ai*(1-Tita)**expi
            for ai, expi in zip(cls._sublimation["a3"], cls._sublimation["exp3"]):
                suma+=ai*log(Tita)**expi

            if cls._sublimation["eq"]==1:
                P=suma*Pref
            elif cls._sublimation["eq"]==2:
                P=exp(suma)*Pref
            elif cls._sublimation["eq"]==3:
                P=exp(Tref/T*suma)*Pref
            return unidades.Pressure(P, "kPa")

        else:
            return None

    def _Vapor_Pressure(self, T):
        if self._vapor_Pressure:
            eq=self._vapor_Pressure["eq"]
            Tita=1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita=Tita**0.5
            suma=sum([n*Tita**x for n, x in zip(self._vapor_Pressure["ao"], self._vapor_Pressure["exp"])])
            if eq in [1, 2]:
                Pr=suma+1
            elif eq in [3, 4]:
                Pr=exp(suma)
            else:
                Pr=exp(self.Tc/T*suma)
            Pv=unidades.Pressure(Pr*self.Pc)
        else:
            Pv=self.componente.Pv(T)
        return Pv

    def _Liquid_Density(self, T=None):
        if not T:
            T=self.T

        if self._liquid_Density:
            eq=self._liquid_Density["eq"]
            Tita=1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita=Tita**(1./3)
            suma=sum([n*Tita**x for n, x in zip(self._liquid_Density["ao"], self._liquid_Density["exp"])])
            if eq in [1, 2]:
                Pr=suma+1
            elif eq in [3, 4]:
                Pr=exp(suma)
            else:
                Pr=exp(self.Tc/T*suma)
            rho=unidades.Density(Pr*self.rhoc)
        else:
            rho=self.componente.RhoL_DIPPR(T)
        return rho

    def _Vapor_Density(self, T=None, P=None):
        if not T:
            T=self.T

        if self._vapor_Density:
            eq=self._vapor_Density["eq"]
            Tita=1-T/self.Tc
            if eq in [2, 4, 6]:
                Tita=Tita**(1./3)
            suma=sum([n*Tita**x for n, x in zip(self._vapor_Density["ao"], self._vapor_Density["exp"])])
            if eq in [1, 2]:
                Pr=suma+1
            elif eq in [3, 4]:
                Pr=exp(suma)
            else:
                Pr=exp(self.Tc/T*suma)
            rho=unidades.Density(Pr*self.rhoc)
        else:
#            rho=0.5
            rho=self._Vapor_Density_Chouaieb(T)
#            rho=self.componente.RhoG_Lee_Kesler(T, P)
        return rho


    def _Vapor_Density_Chouaieb(self, T):
        """O. Chouaieb, J. Ghazouani, A. Bellagi, Simple correlations for saturated liquid and vapor densities of pure fluids, Thermochimica Acta, Volume 424, Issues 1–2, 15 December 2004, Pages 43-51, ISSN 0040-6031,
        http://dx.doi.org/10.1016/j.tca.2004.05.017"""
        coeff={"Ne": (4.4960428, 1.1861523, -0.485720905),
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
        tau=1-T/self.Tc

        if self.__class__.__name__ in coeff:
            m, n, p=coeff[self.__class__.__name__]
        else:
            Zc= self.rhoc*self.R*self.Tc/self.Pc
            Ni=[-0.1497547, 0.6006565]
            Pi=[-19.348354, -41.060325, 1.1878726]
            m=3-1.71*self.f_acent # Simple custom regresion of Fig.8
            p=Zc**2/(Pi[0]+Pi[1]*log10(Zc)*Pi[2])
            n=p+1./(Ni[0]*self.f_acent+Ni[1])

        alfa=exp(tau**(1./3)+tau**0.5+tau+tau**m)
        rhog=self.rhoc*exp(p*(alfa**n-exp(1-alfa)))
        return rhog


    def _Surface(self, T=None):
        """Equation for the surface tension"""
        if self.Tt<=self.T<=self.Tc:
            if self._surface:
                if not T:
                    T=self.T
                tau=1-T/self.Tc
                tension=0
                for sigma, n in zip(self._surface["sigma"], self._surface["exp"]):
                    tension+=sigma*tau**n
                sigma=unidades.Tension(tension)
            else:
                sigma=self.componente.Tension(self.T)
        else:
            sigma=None
        return sigma


    def _Omega(self):
        """Collision integral calculations
            0   -   None
            1   -   Lemmon: nitrogen. oxygen, argon, aire  and custom collison parameter, co2
            2   -   Younglove: C1, C2, C3, iC4, nC4
            3   -   CI0 from NIST:  Chung
        Ref:
            Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.
            Younglove, B.A. and Ely, J.F. (1987). Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane. J. Phys. Chem. Ref. Data  16: 577-798.
            Fenghour, A., Wakeham, W.A., Vesovic, V., Watson, J.T.R., Millat, J., and Vogel, E., "The viscosity of ammonia," J. Phys. Chem. Ref. Data, 24:1649-1667, 1995.
            Fenghour, A., Wakeham, W.A., Vesovic, V., "The Viscosity of Carbon Dioxide," J. Phys. Chem. Ref. Data, 27:31-44, 1998.
            T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E. "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties" Ind. Eng. Chem. Res. 1998, 27, 671-679,
        """
        if self._viscosity["omega"]==1:
            if self._viscosity.has_key("collision"):
                b=self._viscosity["collision"]
            else:
                b=[0.431, -0.4623, 0.08406, 0.005341, -0.00331]
            T_=log(self.T/self._viscosity["ek"])
            suma=0
            for i, bi in enumerate(b):
                suma+=bi*T_**i
            omega=exp(suma)
        elif self._viscosity["omega"]==2:
            if self._viscosity.has_key("collision"):
                b=self._viscosity["collision"]
            else:
                b=[-3.0328138281, 16.918880086, -37.189364917, 41.288861858, -24.615921140, 8.9488430959, -1.8739245042, 0.20966101390, -0.009657043707]
            T_=self._viscosity["ek"]/self.T
            suma=0
            for i, bi in enumerate(b):
                suma+=bi*T_**((3.-i)/3.)
            omega=1./suma
        elif self._viscosity["omega"]==3:
            #FIXME: reference from https://github.com/thorade/HelmholtzMedia/blob/master/Interfaces/PartialHelmholtzMedium/Transport/dynamicViscosity_dilute.mo#L27
            T_ =self.T/self._viscosity.get("ek", self.Tc/1.2593)
            omega=1.16145/T_**0.14874 + 0.52487/exp(0.77320*T_) + 2.16178/exp(2.4378*T_) - 6.435e-4*T_**0.14874*sin(18.0323*T_**-0.76830-7.27371)
        return omega


    def _Visco0(self):
        if self._viscosity["eq"]==3:
            Tc=self._viscosity.get("Tref", 1.)
            rhoc=self._viscosity.get("rhoref", 1.)
            tau=self.T/Tc
            delta=self.rho/self.M/rhoc
            muo=0
            for n, t in zip(self._viscosity["n_poly"], self._viscosity["t_poly"]):
                muo+=n*tau**t
            if self._viscosity.has_key("n_polyden"):
                den=1
                for n, t, d in zip(self._viscosity["n_polyden"], self._viscosity["t_polyden"], self._viscosity["d_polyden"]):
                    den*=n*tau**t*delta**d
                    muo/=den
        elif self._viscosity["eq"]==5:
            omega=self._Omega()
            Fc=1-0.2756*self._viscosity["w"]+0.059035*self._viscosity["mur"]**4+self._viscosity["k"]
            muo=4.0795e-5*(self.M*self.T)**0.5*self.rhoc**(2./3.)/omega*Fc
        else:
            omega=self._Omega()
            Nchapman=self._viscosity.get("n_chapman", 0.0266958)
            tchapman=self._viscosity.get("t_chapman", 0.5)
            muo=Nchapman*(self.M*self.T)**tchapman/(self._viscosity["sigma"]**2*omega)

            #other adittional empirical terms
            if self._viscosity.has_key("n_ideal"):
                tau=self.T/self._viscosity["Tref"]
                for n, t in zip(self._viscosity["n_ideal"], self._viscosity["t_ideal"]):
                    muo+=n*tau**t
        return muo


    def _Viscosity(self):
        if self._viscosity:
            if self._viscosity["eq"]==0:
                mu=self.__getattribute__(self._viscosity["method"])().muPas

            elif self._viscosity["eq"]==1:
                muo=self._Visco0()

                #second virial
                Tc=self._viscosity.get("Tref_virial", self.Tc)
                etar=self._viscosity.get("etaref_virial", 1.)
                tau=self.T/Tc
                mud=0
                if self._viscosity.has_key("n_virial"):
                    muB=0
                    for n, t in zip(self._viscosity["n_virial"], self._viscosity["t_virial"]):
                        muB+=n*tau**t
                    mud=etar*muB*self.rho/self.M*muo

                Tc=self._viscosity.get("Tref_res", self.Tc)
                rhoc=self._viscosity.get("rhoref_res", self.rhoc)
                mured=self._viscosity.get("etaref_res", 1.)
                tau=self.T/Tc
                delta=self.rho/rhoc
                if abs(delta-1)<=0.001:
                    expdel=self.rho/self.rhoc
                else:
                    expdel=delta

                mur=0
                #close-packed density;
                if self._viscosity.has_key("n_packed"):
                    del0=0
                    for n, t in zip(self._viscosity["n_packed"], self._viscosity["t_packed"]):
                        del0+=n*tau**t
                else:
                    del0=1.

                #polynomial term
                if self._viscosity.has_key("n_poly"):
                    for n, t, d, c, g in zip(self._viscosity["n_poly"], self._viscosity["t_poly"], self._viscosity["d_poly"], self._viscosity["c_poly"], self._viscosity["g_poly"]):
                        vis=n*tau**t*delta**d*del0**g
                        if c:
                            vis*=exp(-expdel**c)
                        mur+=vis
                #numerator of rational poly; denominator of rat. poly;
                num=0
                if self._viscosity.has_key("n_num"):
                    for n, t, d, c, g in zip(self._viscosity["n_num"], self._viscosity["t_num"],self._viscosity["d_num"], self._viscosity["c_num"], self._viscosity["g_num"]):
                        num+=n*tau**t*delta**d*del0**g
                        if c:
                            num*=exp(-expdel**c)
                if self._viscosity.has_key("n_den"):
                    den=0
                    for n, t, d, c, g in zip(self._viscosity["n_den"], self._viscosity["t_den"],self._viscosity["d_den"], self._viscosity["c_den"], self._viscosity["g_den"]):
                        den+=n*tau**t*delta**d*del0**g
                        if c:
                            den*=exp(-expdel**c)
                else:
                    den=1.
                mur+=num/den

                #numerator of exponential; denominator of exponential
                num=0
                if self._viscosity.has_key("n_numexp"):
                    for n, t, d, c, g in zip(self._viscosity["n_numexp"], self._viscosity["t_numexp"],self._viscosity["d_numexp"], self._viscosity["c_numexp"], self._viscosity["g_numexp"]):
                        num+=n*tau**t*delta**d*del0**g
                if self._viscosity.has_key("n_denexp"):
                    den=0
                    for n, t, d, c, g in zip(self._viscosity["n_denexp"], self._viscosity["t_denexp"],self._viscosity["d_denexp"], self._viscosity["c_denexp"], self._viscosity["g_denexp"]):
                        den+=n*tau**t*delta**d*del0**g
                else:
                    den=1.
                if self._viscosity.has_key("n_numexp"):
                    mur+=exp(num/den)

                mur*=mured
                mu=muo+mud+mur

            elif self._viscosity["eq"]==2:
                muo=self._Visco0()
                f=self._viscosity["F"]
                e=self._viscosity["E"]
                mu1=f[0]+f[1]*(f[2]-log(self.T/f[3]))**2

                rho=self.rho/self.M
                G=e[0]+e[1]/self.T
                H=rho**0.5*(rho-self._viscosity["rhoc"])/self._viscosity["rhoc"]
                F=G+(e[2]+e[3]*self.T**-1.5)*rho**0.1+H*(e[4]+e[5]/self.T+e[6]/self.T**2)
                mu2=exp(F)-exp(G)
                mu=muo+mu1*rho+mu2

            elif self._viscosity["eq"]==3:
                Tc=self._viscosity.get("Tref", 1.)
                rhoc=self._viscosity.get("rhoref", 1.)
                muref=self._viscosity.get("muref", 1.)
                tau=self.T/Tc
                delta=self.rho/self.M/rhoc

                muo=self._Visco0()

                mur=0
                for n, t, d in zip(self._viscosity["n_num"], self._viscosity["t_num"], self._viscosity["d_num"]):
                    muo+=n*tau**t*delta**d
                if self._viscosity.has_key("n_den"):
                    den=1
                    for n, t, d in zip(self._viscosity["n_den"], self._viscosity["t_den"], self._viscosity["d_den"]):
                        den*=n*tau**t*delta**d
                    mur/=den

                mu=(muo+mur)*muref

            elif self._viscosity["eq"]==4:
                muo=self._Visco0()
                mur=0
                Gamma=self.Tc/self.T
                psi1=exp(Gamma)-1.0
                psi2= exp(Gamma**2)-1.0
                a=self._viscosity["a"]
                b=self._viscosity["b"]
                c=self._viscosity["c"]
                A=self._viscosity["A"]
                B=self._viscosity["B"]
                C=self._viscosity["C"]
                D=self._viscosity["D"]
                ka=(a[0]+a[1]*psi1+a[2]*psi2)*Gamma
                kaa=(A[0]+A[1]*psi1+A[2]*psi2)*Gamma**3
                kr=(b[0]+b[1]*psi1+b[2]*psi2)*Gamma
                krr=(B[0]+B[1]*psi1+B[2]*psi2)*Gamma**3
                ki=(c[0]+c[1]*psi1+c[2]*psi2)*Gamma
                kii=(C[0]+C[1]*psi1+C[2]*psi2)*Gamma**3

                Prep=self.T*self.dpdT.barK
                Patt=self.P.bar-Prep
                Pid=self.rho*self.R*self.T/1e5
                delPr=Prep-Pid
                mur = kr*delPr + ka*Patt + krr*delPr**2 + kaa*Patt**2 + ki*Pid + kii*Pid**2 + D[0]*Prep**3*Gamma**2

                mu=(muo+mur*1e3)

            elif self._viscosity["eq"]==5:
                a0=6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540, -1.89992, 24.2745, 0.79716, -0.23816, 0.68629e-1
                a1=50.4119, -0.11536e-2, 254.209, 38.0957, 7.63034, -12.5367, 3.44945, 1.11764, 0.67695e-1, 0.34793
                a2=-51.6801, -0.62571e-2, -168.481, -8.46414, -14.3544, 4.98529, -11.2913, 0.12348e-1, -0.8163, 0.59256
                a3=1189.02, 0.37283e-1, 3898.27, 31.4178, 31.5267, -18.1507, 69.3466, -4.11661, 4.02528, -0.72663
                A=[None]
                for i in range(10):
                    A.append(a0[i]+a1[i]*self._viscosity["w"]+a2[i]*self._viscosity["mur"]**4+a3[i]*self._viscosity["k"])

                muo=self._Visco0()
                Y=self.rho/self.rhoc/6
                T_ =self.T/self._viscosity.get("ek", self.Tc/1.2593)
                G1=(1-0.5*Y)/(1-Y)**3
                G2=(A[1]*(1-exp(-A[4]*Y))/Y+A[2]*G1*exp(A[5]*Y)+A[3]*G1)/(A[1]*A[4]+A[2]+A[3])
                muk=muo*(1/G2+A[6]*Y)
                mup=36.344e-6*(self.M*self.Tc)**0.5*self.rhoc**(2./3.)*A[7]*Y**2*G2*exp(A[8]+A[9]/T_+A[10]/T_**2)
                mu=muk+mup

            elif self._viscosity["eq"]=="ecs":
#                rhoc=self._constants.get("rhoref", self.rhoc)
#                Tc=self._constants.get("Tref", self.Tc)
#                delta=rho/rhoc
#                tau=Tc/T
#
#                fio, fiot, fiott, fiod, fiodd, fiodt=self._phi0(self._constants["cp"], tau, delta, Tref, Pref)
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

                Tr=self.T/self.Tc
                rhor=self.rho/self.rhoc

                fint=0
                for n, m in zip(self._viscosity["fint"], self._viscosity["fint_t"]):
                    fint+=n*self.T**m

                psi=0
                for n, t, d in zip(self._viscosity["psi"], self._viscosity["psi_t"], self._viscosity["psi_d"]):
                    psi+=n*Tr**t*rhor**d

                mu=None

        else:
            mu=None
        return unidades.Viscosity(mu, "muPas")


    def _KCritical(self):
        if self._thermal["critical"]==0:
            tc=0

        elif self._thermal["critical"]==1:
            tc=0
            if self._thermal.has_key("crit_num_n"):
                Tref=self._thermal["crit_num_Tref"]
                if Tref<0:
                    tau=Tref/self.T
                else:
                    tau=self.T/Tref
                delta=self.rho/self.M/self._thermal["crit_num_rhoref"]

                for n, alfa, t, beta, d in zip(self._thermal["crit_num_n"], self._thermal["crit_num_alfa"], self._thermal["crit_num_t"], self._thermal["crit_num_beta"], self._thermal["crit_num_d"]):
                    tc+=n*(tau+alfa)**t*(delta+beta)**d

                if self._thermal.has_key("crit_den_n"):
                    den=0
                    for n, alfa, t, beta, d, c in zip(self._thermal["crit_den_n"], self._thermal["crit_den_alfa"], self._thermal["crit_den_t"], self._thermal["crit_den_beta"], self._thermal["crit_den_d"], self._thermal["crit_den_c"]):
                        if c==99:
                            tcmax=max(tau, alfa-tau)
                            den+=n*tcmax**t*(delta+beta)**d
                        else:
                            den+=n*(tau+alfa)**t*(delta+beta)**d
                    tc/=den

            if self._thermal.has_key("crit_exp_n"):
                Tref=self._thermal["crit_exp_Tref"]
                if Tref<0:
                    tau=Tref/self.T
                else:
                    tau=self.T/Tref
                delta=self.rho/self.M/self._thermal["crit_exp_rhoref"]
                expo=0
                for n, alfa, t, beta, d, c in zip(self._thermal["crit_exp_n"], self._thermal["crit_exp_alfa"], self._thermal["crit_exp_t"], self._thermal["crit_exp_beta"], self._thermal["crit_exp_d"], self._thermal["crit_exp_c"]):
                    expo+=n*(tau+alfa)**t*(delta+beta)**d
                tc*=exp(expo)

            tc*=self._thermal["crit_num_k"]

        elif self._thermal["critical"]==2:
            X=self._thermal["X"]
            xi=self.Pc*self.rho/self.rhoc**2/self.derivative("P", "rho", "T")
            normterm=X[3]*Boltzmann/self.Pc*(self.T*self.dpdT*self.rhoc/self.rho)**2*xi**X[2]
            delT=abs(self.T-self.Tc)/self.Tc
            delrho=abs(self.rho-self.rhoc)/self.rhoc
            expterm=exp(-(X[0]*delT**4+X[1]*delrho**4))
            mu=self._Viscosity()
            tc=normterm*expterm/(6*pi*mu*self._thermal["Z"])

        elif self._thermal["critical"]==3:
            qd=self._thermal["qd"]
            Tref=self._thermal["Tcref"]
            x_T=self.Pc*self.rho*self.M/self.rhoc**2/self.derivative("P", "rho", "T")
            x_Tr=self.Pc*self.rho*self.M/self.rhoc**2/self.derivative("P", "rho", "T")*Tref/self.T
            delchi=x_T-x_Tr
            if delchi <= 0:
                tc=0
            else:
                Xi=self._thermal["Xio"]*(delchi/self._thermal["gam0"])**(self._thermal["gnu"]/self._thermal["gamma"])
                omega=2/pi*((self.cp-self.cv)/self.cp*arctan(Xi*qd)+self.cv/self.cp*Xi*qd)
                omega0=2/pi*(1-exp(-1/(1./qd/Xi+Xi**2*qd**2/3*(self.rhoc/self.rho)**2)))
                tc=self.rho/self.M*1e9*self.cp*Boltzmann*self._thermal["R0"]*self.T/(6*pi*Xi*self.mu.muPas)*(omega-omega0)

        elif self._thermal["critical"]==4:
            rho=self.rho/self.M
            Xt=(rho*self._thermal["Pcref"]/self._thermal["rhocref"]**2/self.derivative("P", "rho", "T"))**self._thermal["expo"]
            parterm=self._thermal["alfa"]*Boltzmann/self._thermal["Pcref"]*(self.T*self.dpdT*self._thermal["rhocref"]/rho)**2*Xt*1e21
            delT=abs(self.T-self._thermal["Tcref"])/self._thermal["Tcref"]
            delrho=abs(rho-self._thermal["rhocref"])/self._thermal["rhocref"]
            expterm=exp(-(self._thermal["alfa"]*delT**2+self._thermal["beta"]*delrho**4))
            tc=parterm*expterm/(6*pi*self._thermal["Xio"]*self.mu.muPas)*self._thermal["kcref"]

        elif self._thermal["critical"]=="NH3":
            tr=abs(self.T-405.4)/405.4
            trr=tr
            if 404.4 < self.T < 406.5 and (self.rho/self.M<9.6 or self.rho/self.M>18): #to avoid infinite value in critical region
                trr=0.002
            etab=1.0e-5*(2.6+1.6*tr)
            dPT=1.0e5*(2.18-0.12/exp(17.8*tr))
            if trr == 0.0:
                tcrhoc=1.e20
                dtcid=tcrhoc
                xcon=-1.e20
            else:
                tcrhoc=1.2*1.38066e-23*self.T**2*dPT**2*0.423e-8/(trr**1.24)*(1.0+1.429*tr**0.5)/(6.0*pi*etab*(1.34e-10/trr**0.63*(1.0+1.0*tr**0.5)))
                dtcid=tcrhoc*exp(-36.0*tr**2)
                xcon=0.61*235+16.5*log(trr)
            if self.rho/self.rhoc < 0.6:
                tccsw=dtcid*xcon**2/(xcon**2+(141.0- 0.96*235.0)**2)
                tc=tccsw*self.rho**2/141.0**2
            else:
                tc=dtcid*xcon**2/(xcon**2+(self.rho-0.96*235.0)**2)

        elif self._thermal["critical"]=="CH4":
            tau=self.Tc/self.T
            delta=self.rho/self.rhoc
            ts=(self.Tc-self.T)/self.Tc
            ds=(self.rhoc-self.rho)/self.rhoc
            xt=0.28631*delta*tau/self.derivative("P", "rho", "T")
            ftd=exp(-2.646*abs(ts)**0.5+2.678*ds**2-0.637*ds)
            tc=91.855/self.mu/tau**2*self.dpdT**2*xt**0.4681*ftd*1e-3

        return tc

    def _ThCond(self):
        if self._thermal:
            if self._thermal["eq"]==0:
                k=self.__getattribute__(self._thermal["method"])().mWmK

            elif self._thermal["eq"]==1:
                #Dilute gas terms
                kg=0
                if self._thermal.has_key("no"):
                    tau=self.T/self._thermal["Tref"]
                    for n, c in zip(self._thermal["no"], self._thermal["co"]):
                        if c == -99:
                            cpi=1.+n*(self.cp0.kJkgK-2.5*self.R.kJkgK)
                            kg*=cpi
                        elif c == -98:
                            muo=self._Visco0()
                            cpi=self.cp0-R
                            kg+=n*cpi*muo
                        elif c == -97:
                            muo=self._Visco0()
                            kg+=n*muo
                        elif c == -96:
                            cpi=self.cp0/self.R-2.5
                            muo=self._Visco0()
                            kg=(kg*cpi+15./4.)*self.R.kJkgK*muo/self.M
                        else:
                            kg+=n*tau**c

                    if self._thermal.has_key("noden"):
                        den=0
                        for n, t in zip(self._thermal["noden"], self._thermal["toden"]):
                            den+=n*tau**t
                        kg/=den

                    kg*=self._thermal["kref"]

                #Backgraund terms
                kb=0
                if self._thermal.has_key("nb"):
                    tau=self.T/self._thermal["Trefb"]
                    delta=self.rho/self.M/self._thermal["rhorefb"]
                    for n, t, d, c in zip(self._thermal["nb"], self._thermal["tb"], self._thermal["db"], self._thermal["cb"]):
                        if c == -99:
                            if tau<1:
                                th=(1.-tau)**(1./3.)
                                kb/=exp(-1.880284*th**1.062-2.8526531*th**2.5-3.000648*th**4.5-5.251169*th**7.5-13.191869*th**12.5-37.553961*th**23.5)
                        else:
                            if c!=0:
                                kb+=n*tau**t*delta**d*exp(-delta**c)
                            else:
                                kb+=n*tau**t*delta**d

                    if self._thermal.has_key("nbden"):
                        den=0
                        for n, t, d in zip(self._thermal["nbden"], self._thermal["tbden"], self._thermal["dbden"]):
                            den+=n*tau**t*delta**d
                        kb/=den

                    kb*=self._thermal["krefb"]

                #Critical enhancement
                kc=self._KCritical()

                k=kg+kb+kc


            elif self._thermal["eq"]==2:
                self._viscosity=self._thermal["visco"]
                muo=self._Visco0()
                G=self._thermal["G"]
                kg=1e-3*muo/self.M*(3.75*self.R+(self.cp0.kJkgK-2.5*self.R)*(G[0]+G[1]*self._viscosity["ek"]/self.T))

                E=self._thermal["E"]
                F0=E[0]+E[1]/self.T+E[2]/self.T**2
                F1=E[3]+E[4]/self.T+E[5]/self.T**2
                F2=E[6]+E[7]/self.T
                rho=self.rho/self.M
                kb=(F0+F1*rho)*rho/(1+F2*rho)

                #Critical enhancement
                kc=self._KCritical()

                k=kg+kb+kc

            elif self._thermal["eq"]==3:
                T_=self._thermal["ek"]/self.T
                rho=self.rho/self.M
                suma=0
                for i, bi in enumerate(self._thermal["b"]):
                    suma+=bi*T_**((3.-i)/3.)
                ko=self._thermal["Nchapman"]*self.T**self._thermal["tchapman"]/(self._thermal["sigma"]**2/suma)

                f=self._thermal["F"]
                e=self._thermal["E"]
                k1=f[0]+f[1]*(f[2]-log(self.T/f[3]))**2
                kg=ko+k1*rho

                G=e[0]+e[1]/self.T
                H=rho**0.5*(rho-self._thermal["rhoc"])/self._thermal["rhoc"]
                F=G+(e[2]+e[3]*self.T**-1.5)*rho**0.1+H*(e[4]+e[5]/self.T+e[6]/self.T**2)
                kb=exp(F)-exp(G)

                bl=self._thermal["ff"]*(self._thermal["rm"]**5*self.rho.gcc*Avogadro/self.M*self._thermal["Nchapman"]/self.T)**0.5
                y=6.0*pi*self.mu/100000.*bl*(Boltzmann*self.T*self.rho.gcc*Avogadro/self.M)**0.5
                deltaL=0.0
                if self.derivative("P", "rho", "T") > 0:
                    deltaL=Boltzmann*(self.T*self.dpdT)**2/(self.rho.gcc*self.derivative("P", "rho", "T") )**0.5/y
                else:
                    deltaL=0.0
                kc=deltaL*exp(-18.66*((self.rho-self.rhoc)/self.rhoc)**4-4.25*((self.T-self.Tc)/self.Tc)**2)

                k=kg+kb+kc

            elif self._thermal["eq"]=="ecs":
                k=None
                #TODO:

        else:
            k=None
        return unidades.ThermalConductivity(k)

    @classmethod
    def __test__(cls):
        pruebas={}
        for i, test in enumerate(cls._test):
            prueba=">>> for value1 in {}:".format(test["value1"])+os.linesep
            prueba+="...  for value2 in {}:".format(test["value2"])+os.linesep
            prueba+="...   fluido={}({}=value1, {}=value2)".format(cls.__name__, test["var1"], test["var2"])+os.linesep
            prueba+="...   print("+"'"+"{: .5g} "*len(test["prop"])+"'"+".format("
            for propiedad, unidad in zip(test["prop"], test["unit"]):
                prueba+="fluido.{}.{}, ".format(propiedad, unidad)
            prueba+="))"+os.linesep
            prueba+=test["result"]
        pruebas[str(i)]=prueba
        return pruebas





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

