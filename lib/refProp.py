#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Library to multiparameter equation of state calculation using refprop
# https://github.com/BenThelen/python-refprop
# refprop dll must be installed from NIST package, license required
# optional method to meos tools calculations
###############################################################################


import sys

try:
    import refprop
except ImportError:
    pass

from lib import unidades, mEoS
# from lib.config import Preferences
from lib.mezcla import mix_unitmassflow, mix_unitmolarflow
from lib.thermo import ThermoRefProp


# Automatic loading of refprop name from meos subclass _refPropName property
__all__ = {}
noIds = []
for cmp in mEoS.__all__:
    if cmp.id and cmp._refPropName:
        __all__[cmp.id] = cmp._refPropName
    elif cmp._refPropName:
        noIds.append(cmp._refPropName)


class RefProp(ThermoRefProp):
    """
    Stream class using refProp external library
    Parameters needed to define it are:

        -ref: reference state
        -ids: index of fluid
        -fraccionMolar: molar fraction

        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -rho: Density, kg/m3
        -H: Enthalpy, J/kg
        -S: Entropy, J/kgK
        -U: Internal energy, J/kg
        -x: Quality, -

    setref parameters
    setref(hrf='DEF', ixflag=1, x0=[1], h0=0, s0=0, t0=273, p0=100):
        hrf--reference state for thermodynamic calculations [character*3]
            'NBP':  h,s = 0 at normal boiling point(s)
            'ASH':  h,s = 0 for sat liquid at -40 C (ASHRAE convention)
            'IIR':  h = 200, s = 1.0 for sat liq at 0 C (IIR convention)
            'DEF':  default reference state as specified in fluid file is
                applied to each component (ixflag = 1 is used)
            'OTH':  other, as specified by h0, s0, t0, p0 (real gas state)
            'OT0':  other, as specified by h0, s0, t0, p0 (ideal gas state)
            '???':  change hrf to the current reference state and exit.
        ixflag--composition flag:
            1 = ref state applied to pure components
            2 = ref state applied to mixture icomp
        following input has meaning only if ixflag = 2
            x0--composition for which h0, s0 apply; list(1:nc) [mol frac]
                this is useful for mixtures of a predefined composition, e.g.
                refrigerant blends such as R410A
        following inputs have meaning only if hrf = 'OTH'
            h0--reference state enthalpy at t0,p0 {icomp} [J/mol]
            s0--reference state entropy at t0,p0 {icomp} [J/mol-K]
            t0--reference state temperature [K]
                t0 = -1 indicates saturated liquid at normal boiling point
                    (bubble point for a mixture)
            p0--reference state pressure [kPa]
                p0 = -1 indicates saturated liquid at t0 {and icomp}
                p0 = -2 indicates saturated vapor at t0 {and icomp}

    setmod parameters:
    setmod(htype='NBS', hmix='NBS', *hcomp):
    inputs 'in string format':
        htype - flag indicating which models are to be set [character*3]:
            'EOS':  equation of state for thermodynamic properties
            'ETA':  viscosity
            'TCX':  thermal conductivity
            'STN':  surface tension
            'NBS':  reset all of the above model types and all subsidiary
                component models to 'NBS'; values of hmix and hcomp are ignored
        hmix--mixture model to use for the property specified in htype
            [character*3]:
            ignored if number of components = 1
            some allowable choices for hmix:
                'NBS':  use NIST recommendation for specified fluid/mixture
                'HMX':  mixture Helmholtz model for thermodynamic properties
                'ECS':  extended corresponding states for viscosity or th cond
                'STX':  surface tension mixture model
        hcomp--component model(s) to use for property specified in htype
            [array (1..nc) of character*3]:
                'NBS':  NIST recommendation for specified fluid/mixture
            some allowable choices for an equation of state:
                'FEQ':  Helmholtz free energy model
                'BWR':  pure fluid modified Benedict-Webb-Rubin (MBWR)
                'ECS':  pure fluid thermo extended corresponding states
            some allowable choices for viscosity:
                'ECS':  extended corresponding states (all fluids)
                'VS1':  the 'composite' model for R134a, R152a, NH3, etc.
                'VS2':  Younglove-Ely model for hydrocarbons
                'VS4':  Generalized friction theory of Quinones-Cisneros and
                    Deiters
                'VS5':  Chung et al. (1988) predictive model
            some allowable choices for thermal conductivity:
                'ECS':  extended corresponding states (all fluids)
                'TC1':  the 'composite' model for R134a, R152a, etc.
                'TC2':  Younglove-Ely model for hydrocarbons
                'TC5':  Chung et al. (1988) predictive model
            some allowable choices for surface tension:
                'ST1':  surface tension as f(tau); tau = 1 - T/Tc

    setktv parameters
    setktv(icomp, jcomp, hmodij, fij=([0] * _nmxpar), hfmix='HMX.BNC'):
        icomp--component
        jcomp--component j
        hmodij--mixing rule for the binary pair i,j [character*3] e.g.:
            'LJ1' (Lemmon-Jacobsen model)
            'LM1' (modified Lemmon-Jacobsen model) or
            'LIN' (linear mixing rules)
            'RST' indicates reset all pairs to values from original call to
                SETUP (i.e. those read from file) [all other inputs are
                ignored]
        fij--binary mixture parameters [array of dimension nmxpar; currently
            nmxpar is set to 6] the parameters will vary depending on hmodij;
            for example, for the Lemmon-Jacobsen model
                (LJ1):
                    fij(1) = zeta
                    fij(2) = xi
                    fij(3) = Fpq
                    fij(4) = beta
                    fij(5) = gamma
                    fij(6) = 'not used'
        hfmix--file name [character*255] containing generalized parameters
            for the binary mixture model; this will usually be the same as the
            corresponding input to SETUP (e.g.,':fluids:HMX.BNC')
    """
    kwargs = {"ids": [],
              "fraccionMolar": [],
              "fraccionMasica": [],
              "caudalUnitarioMolar": [],
              "caudalUnitarioMasico": [],

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "Q": None,
              "rho": 0.0,
              "D": 0.0,
              "H": 0.0,
              "S": 0.0,
              "u": 0.0,
              "E": 0.0,

              # Configuration parameters
              "preos": False,
              "aga": False,
              "gerg": False,

              "hrf": "DEF",
              "ixflag": 1,
              "x0": [1],
              "h0": 0,
              "s0": 0,
              "t0": 273,
              "p0": 1e5,

              "htype": "NBS",
              "hmix": "NBS",
              "hcomp": "",
              # setktv don't implemented
              }

    __doi__ = [
        {"autor": "Lemmon, E.W., Huber, M.L., McLinden, M.O.",
         "title": "NIST Standard Reference Database 23:  Reference Fluid"
                  "Thermodynamic and Transport Properties-REFPROP, Version"
                  "9.1, National Institute of Standards and Technology,"
                  "Standard Reference Data Program, Gaithersburg, 2013.",
         "ref": "",
         "doi": ""}]

    @property
    def calculable(self):
        """Check in the class is fully defined. The correct definition require
        two parts:
        * Thermodynamics definition: whatever pair properties between T, P,
            h, s, v, rho, x
        * Compound definition: Ids for compound identification, and in
            multicomponent defintion any class of mixture definition
        """
        # Check mix state
        self._multicomponent = False
        if len(self.kwargs["ids"]) > 1:
            self._multicomponent = True

        # Check supported fluid
        REFPROP_available = True
        for id in self.kwargs["ids"]:

            if id not in __all__ and id not in noIds:
                REFPROP_available = False
                if not REFPROP_available:
                    raise ValueError

        # Mix definition
        self._mix = 0
        if len(self.kwargs["fraccionMolar"]) == len(self.kwargs["ids"]):
            self._mix = 1
        elif len(self.kwargs["fraccionMasica"]) == len(self.kwargs["ids"]):
            self._mix = 2
        elif len(self.kwargs["caudalUnitarioMolar"]) == \
                len(self.kwargs["ids"]):
            self._mix = 3
        elif len(self.kwargs["caudalUnitarioMasico"]) == \
                len(self.kwargs["ids"]):
            self._mix = 4

        # Check correct fluid definition
        if self._multicomponent:
            if self.kwargs["ids"] and self._mix:
                self._definition = True
            else:
                self._definition = False
        else:
            self._definition = True

        # Update the kwargs with the special refprop namespace
        if self.kwargs["x"] != RefProp.kwargs["x"]:
            self.kwargs["Q"] = self.kwargs["x"]
        if self.kwargs["rho"] != RefProp.kwargs["rho"]:
            self.kwargs["D"] = self.kwargs["rho"]
        if self.kwargs["u"] != RefProp.kwargs["u"]:
            self.kwargs["E"] = self.kwargs["u"]

        # Check thermo definition
        self._thermo = ""
        for def_ in ["TP", "TQ", "PQ", "TD", "PD", "PH", "PS", "HS", "TH",
                     "TS", "TE", "PE", "ES", "DH", "DS", "DE"]:
            if self.kwargs[def_[0]] != RefProp.kwargs[def_[0]] and \
                    self.kwargs[def_[1]] != RefProp.kwargs[def_[1]]:
                self._thermo = def_

        return self._definition and self._thermo

    def args(self):
        x = self._x()
        # Convert to float because unit subclass raise InputError
        var1 = float(self.kwargs[self._thermo[0]])
        var2 = float(self.kwargs[self._thermo[1]])

        # unit conversion to refprop accepted units
        # P in kPa, U,H in kJ/kmol, S in kJ/kmolK
        if self._thermo[0] == "P":
            var1 /= 1000.
        elif self._thermo[0] in ("E", "H", "S"):
            var1 /= 1000. / self.M
        elif self._thermo[0] == "D":
            var1 /= self.M

        if self._thermo[1] == "P":
            var2 /= 1000.
        elif self._thermo[1] in ("E", "H", "S"):
            var2 /= 1000. / self.M
        elif self._thermo[1] == "D":
            var2 /= self.M

        return self._thermo, var1, var2, x

    def _name(self):
        name = []
        for fld in self.kwargs["ids"]:
            if fld in __all__:
                name.append(__all__[fld])
            elif fld in noIds:
                name.append(fld)
        return name

    def _x(self):
        if self._mix == 1:
            x = self.kwargs["fraccionMolar"]
        elif self._mix == 2:
            x = refprop.xmole(self.kwargs["fraccionMasica"])["x"]
        elif self._mix == 3:
            kw = mix_unitmolarflow(self.kwargs["caudalUnitarioMolar"])
            x = kw["fraccionMolar"]
        elif self._mix == 3:
            kw = mix_unitmassflow(self.kwargs["caudalUnitarioMasico"])
            x = kw["fraccionMolar"]
        else:
            x = [1]
        return x

    def _fixed(self):
        """Calculate fixed properties with the chemical compound ids specified,
        valid only for one component"""

        x = self._x()
        fluido = self._name()

        refprop.setup("def", fluido)

        info = refprop.info()
        self.M = unidades.Dimensionless(info["wmm"])
        self.Tt = unidades.Temperature(info["ttrp"])
        self.Tb = unidades.Temperature(info["tnbpt"])
        self.Tc = unidades.Temperature(info["tcrit"])

        self.R = unidades.SpecificHeat(info["Rgas"]/self.M)
        self.f_accent = unidades.Dimensionless(info["acf"])
        self.momentoDipolar = unidades.DipoleMoment(info["dip"], "Debye")
        self.rhoc = unidades.Density(info["Dcrit"]*self.M)
        self.zc = unidades.Dimensionless(info["zcrit"])
        self.Pc = unidades.Pressure(self.R*self.Tc*self.rhoc*self.zc, "kPa")

        name = refprop.name()
        self.name = name["hname"]
        self.CAS = name["hcas"]

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

        if kwargs.get("s", 0):
            kwargs["S"] = kwargs["s"]
            del kwargs["s"]

        if kwargs.get("h", 0):
            kwargs["H"] = kwargs["h"]
            del kwargs["h"]
        self.kwargs.update(kwargs)

    def calculo(self):
        try:
            self._calculo()
        except refprop.RefpropdllError as e:
            self.status = 4
            self.msg = e

    def _calculo(self):
        # TODO: Add configuration section to Preferences
        # preos = Preferences.getboolean("refProp", "preos")
        # aga = Preferences.getboolean("refProp", "aga")
        # gerg = Preferences.getboolean("refProp", "gerg")
        preos = self.kwargs["preos"]
        aga = self.kwargs["aga"]
        gerg = self.kwargs["gerg"]

        x = self._x()
        fluido = self._name()

        kwmod = [self.kwargs[k] for k in ('htype', 'hmix', 'hcomp')]
        refprop.setmod(*kwmod)
        if gerg:
            refprop.gerg04(ixflag=1)
        refprop.setup("def", fluido)
        # refprop.setktv()
        if preos:
            refprop.preos(ixflag=2)
        elif aga:
            refprop.setaga()
        kwref = {k: self.kwargs[k] for k in (
            'hrf', 'ixflag', 'x0', 'h0', 's0', 't0', 'p0')}
        refprop.setref(**kwref)

        m = refprop.wmol(x)["wmix"]
        self.M = unidades.Dimensionless(m)
        crit = refprop.critp(x)
        self.Pc = unidades.Pressure(crit["pcrit"], "kPa")
        self.Tc = unidades.Temperature(crit["tcrit"])
        self.rhoc = unidades.Density(crit["Dcrit"]*self.M)

        args = self.args()
        flash = refprop.flsh(*args)

        # check if ['q'] in fld
        if 'q' in flash.keys():
            x = flash['q']
        elif 'h' in flash.keys():
            x = refprop.flsh('ph', flash['p'], flash['h'], flash['x'])['q']
        elif 's' in flash.keys():
            x = refprop.flsh('ps', flash['p'], flash['s'], flash['x'])['q']
        if 0 < x < 1:
            region = 4
        else:
            region = 1

        if x < 0:
            x = 0
        elif x > 1:
            x = 1
        self.x = unidades.Dimensionless(x)
        self.T = unidades.Temperature(flash["t"])
        self.P = unidades.Pressure(flash["p"], "kPa")
        self.Tr = unidades.Dimensionless(self.T/self.Tc)
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.rho = unidades.Density(flash["D"]*self.M)
        self.v = unidades.SpecificVolume(1./self.rho)
        self.phase = self.getphase(Tc=self.Tc, Pc=self.Pc, T=self.T, P=self.Pc,
                                   x=self.x, region=region)

        if flash["nc"] == 1:
            name = refprop.name(flash["nc"])
            self.name = name["hname"]
            self.synonim = name["hn80"]
            self.CAS = name["hcas"]

            info = refprop.info(flash["nc"])
            self.R = unidades.SpecificHeat(info["Rgas"]/self.M)
            self.Tt = unidades.Temperature(info["ttrp"])
            self.Tb = unidades.Temperature(info["tnbpt"])
            self.f_accent = unidades.Dimensionless(info["acf"])
            self.momentoDipolar = unidades.DipoleMoment(info["dip"], "Debye")
            self._doc = {}
            for htype in ['EOS', 'CP0', 'ETA', 'VSK', 'TCX', 'TKK', 'STN',
                          'DE ', 'MLT', 'SBL', 'PS ', 'DL ', 'DV ']:
                self._doc[htype] = refprop.getmod(flash["nc"], htype)["hcite"]
        else:
            self.name = ""
            self.synonim = ""
            self.CAS = ""

            rmix = refprop.rmix2(flash["x"])
            self.R = unidades.SpecificHeat(rmix["Rgas"]/self.M)
            self.Tt = unidades.Temperature(None)
            self.Tb = unidades.Temperature(None)
            self.f_accent = unidades.Dimensionless(None)
            self.momentoDipolar = unidades.DipoleMoment(None)
            self._doc = {}

        self._cp0(flash)

        self.Liquido = ThermoRefProp()
        self.Gas = ThermoRefProp()
        if self.x == 0:
            # liquid phase
            self.fill(self.Liquido, flash["t"], flash["D"], flash["x"])
            self.fill(self, flash["t"], flash["D"], flash["x"])
            self.fillNone(self.Gas)
        elif self.x == 1:
            # vapor phase
            self.fill(self.Gas, flash["t"], flash["D"], flash["x"])
            self.fill(self, flash["t"], flash["D"], flash["x"])
            self.fillNone(self.Liquido)
        else:
            # Two phase
            self.fillNone(self)
            self.fill(self.Liquido, flash["t"], flash["Dliq"], flash["xliq"])
            self.fill(self.Gas, flash["t"], flash["Dvap"], flash["xvap"])

            self.v = unidades.SpecificVolume(x*self.Gas.v+(1-x)*self.Liquido.v)
            self.rho = unidades.Density(1./self.v)

            self.u = unidades.Enthalpy(flash["e"]/self.M, "Jg")
            self.h = unidades.Enthalpy(flash["h"]/self.M, "Jg")
            self.s = unidades.SpecificHeat(flash["s"]/self.M, "JgK")
            self.a = unidades.Enthalpy(self.u-self.T*self.s)
            self.g = unidades.Enthalpy(self.h-self.T*self.s)

            self.rhoM = unidades.MolarDensity(self.rho/self.M)
            self.hM = unidades.MolarEnthalpy(self.h*self.M)
            self.sM = unidades.MolarSpecificHeat(self.s*self.M)
            self.uM = unidades.MolarEnthalpy(self.u*self.M)
            self.aM = unidades.MolarEnthalpy(self.a*self.M)
            self.gM = unidades.MolarEnthalpy(self.g*self.M)

        if self.x < 1 and self.T <= self.Tc:
            surten = refprop.surten(flash["t"], flash["Dliq"], flash["Dvap"],
                                    flash["xliq"], flash["xvap"])
            self.sigma = unidades.Tension(surten["sigma"])
        else:
            self.sigma = unidades.Tension(None)

        if 0 < self.x < 1:
            self.Hvap = unidades.Enthalpy(self.Gas.h-self.Liquido.h)
            self.Svap = unidades.SpecificHeat(self.Gas.s-self.Liquido.s)
            self.K = []
            for x, y in zip(self.Liquido.fraccion, self.Gas.fraccion):
                self.K.append(unidades.Dimensionless(y/x))
        else:
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)
            self.K = [unidades.Dimensionless(1)]*flash["nc"]
        self.invT = unidades.InvTemperature(-1/self.T)

        # NOT supported on Windows
        if sys.platform != "win32":
            excess = refprop.excess(flash["t"], flash["D"], flash["x"])
            self.vE = unidades.Volume(excess["vE"]/self.M)
            self.uE = unidades.Enthalpy(excess["eE"]/self.M, "Jg")
            self.hE = unidades.Enthalpy(excess["hE"]/self.M, "Jg")
            self.sE = unidades.SpecificHeat(excess["sE"]/self.M, "JgK")
            self.aE = unidades.Enthalpy(excess["aE"]/self.M, "Jg")
            self.gE = unidades.Enthalpy(excess["gE"]/self.M, "Jg")
        else:
            self.vE = unidades.Volume(0)
            self.uE = unidades.Enthalpy(0)
            self.hE = unidades.Enthalpy(0)
            self.sE = unidades.SpecificHeat(0)
            self.aE = unidades.Enthalpy(0)
            self.gE = unidades.Enthalpy(0)

        self.csat = []
        self.dpdt_sat = []
        self.cv2p = []
        if self.Tt <= flash["t"] <= self.Tc:
            for i in range(1, flash["nc"]+1):
                dat = refprop.dptsatk(i, flash["t"], kph=2)
                cs = unidades.SpecificHeat(dat["csat"]/self.M, "JgK")
                self.csat.append(cs)
                self.dpdt_sat.append(
                    unidades.PressureTemperature(dat["dpdt"], "kPaK"))
                cv2 = refprop.cv2pk(i, flash["t"], flash["D"])
                cv = unidades.SpecificHeat(cv2["cv2p"]/self.M, "JgK")
                self.cv2p.append(cv)

    def _cp0(self, flash):
        "Set ideal properties to state"""
        cp0 = refprop.therm0(flash["t"], flash["D"], flash["x"])
        self.v0 = unidades.SpecificVolume(self.R*self.T/self.P.kPa)
        self.rho0 = unidades.Density(1/self.v0)
        self.P0 = unidades.Pressure(cp0["p"], "kPa")
        self.P_Pideal = unidades.Pressure(self.P-self.P0)

        self.h0 = unidades.Enthalpy(cp0["h"]/self.M, "kJkg")
        self.u0 = unidades.Enthalpy(cp0["e"]/self.M, "kJkg")
        self.s0 = unidades.SpecificHeat(cp0["s"]/self.M, "kJkgK")
        self.a0 = unidades.Enthalpy(cp0["A"]/self.M, "kJkg")
        self.g0 = unidades.Enthalpy(cp0["G"]/self.M, "kJkg")
        self.w0 = unidades.Speed(cp0["w"])

        cp0 = refprop.therm0(float(self.T), self.rho0/self.M, flash["x"])
        self.cp0 = unidades.SpecificHeat(cp0["cp"]/self.M, "kJkgK")
        self.cv0 = unidades.SpecificHeat(cp0["cv"]/self.M, "kJkgK")
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv0)
        self.gamma0 = self.cp0_cv

        self.rhoM0 = unidades.MolarDensity(self.rho0/self.M)
        self.hM0 = unidades.MolarEnthalpy(self.h0*self.M)
        self.uM0 = unidades.MolarEnthalpy(self.u0*self.M)
        self.sM0 = unidades.MolarSpecificHeat(self.s0*self.M)
        self.aM0 = unidades.MolarEnthalpy(self.a0*self.M)
        self.gM0 = unidades.MolarEnthalpy(self.g0*self.M)
        self.cpM0 = unidades.MolarSpecificHeat(self.cp0*self.M)
        self.cvM0 = unidades.MolarSpecificHeat(self.cv0*self.M)

    def fill(self, fase, T, rho, x):
        if sum(x) != 1:
            x = [round(xi, 10) for xi in x]
        mol = refprop.wmol(x)
        thermo = refprop.therm2(T, rho, x)
        thermo3 = refprop.therm3(T, rho, x)

        fase._bool = True
        fase.M = unidades.Dimensionless(mol["wmix"])
        fase.rho = unidades.Density(rho*fase.M)
        fase.v = unidades.SpecificVolume(1./fase.rho)
        fase.Z = unidades.Dimensionless(thermo["Z"])

        fase.u = unidades.Enthalpy(thermo["e"]/fase.M, "Jg")
        fase.h = unidades.Enthalpy(thermo["h"]/fase.M, "Jg")
        fase.s = unidades.SpecificHeat(thermo["s"]/fase.M, "JgK")
        fase.a = unidades.Enthalpy(thermo["A"]/fase.M, "Jg")
        fase.g = unidades.Enthalpy(thermo["G"]/fase.M, "Jg")

        fase.cv = unidades.SpecificHeat(thermo["cv"]/fase.M, "JgK")
        fase.cp = unidades.SpecificHeat(thermo["cp"]/fase.M, "JgK")
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.gamma = fase.cp_cv
        fase.w = unidades.Speed(thermo["w"])

        fase.rhoM = unidades.MolarDensity(fase.rho/self.M)
        fase.hM = unidades.MolarEnthalpy(fase.h*self.M)
        fase.sM = unidades.MolarSpecificHeat(fase.s*self.M)
        fase.uM = unidades.MolarEnthalpy(fase.u*self.M)
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(fase.cv*self.M)
        fase.cpM = unidades.MolarSpecificHeat(fase.cp*self.M)

        residual = refprop.residual(T, rho, x)
        fase.pr = unidades.Pressure(residual["pr"], "kPa")
        fase.ur = unidades.Enthalpy(residual["er"]/fase.M, "Jg")
        fase.hr = unidades.Enthalpy(residual["hr"]/fase.M, "Jg")
        fase.sr = unidades.SpecificHeat(residual["sr"]/fase.M, "JgK")
        fase.ar = unidades.Enthalpy(residual["Ar"]/fase.M, "Jg")
        fase.gr = unidades.Enthalpy(residual["Gr"]/fase.M, "Jg")
        fase.cvr = unidades.SpecificHeat(residual["cvr"]/fase.M, "JgK")
        fase.cpr = unidades.SpecificHeat(residual["cpr"]/fase.M, "JgK")

        fase.alfav = unidades.InvTemperature(thermo["beta"])
        fase.kappa = unidades.InvPressure(thermo["xkappa"], "kPa")
        fase.kappas = unidades.InvPressure(thermo3["betas"], "kPa")
        fase.alfap = unidades.Density(fase.alfav/self.P/fase.kappa)
        fase.deltat = unidades.EnthalpyPressure(
            thermo3["thrott"]/fase.M, "kJkgkPa")
        fase.joule = unidades.TemperaturePressure(thermo["hjt"], "KkPa")
        fase.betas = unidades.TemperaturePressure(
            self.derivative("T", "P", "s", fase))
        fase.betap = unidades.Density(
            -1/self.P*self.derivative("P", "v", "T", fase))

        fase.Kt = unidades.Pressure(thermo3["xkkt"], "kPa")
        fase.Ks = unidades.Pressure(thermo3["bs"], "kPa")
        fase.kt = unidades.Dimensionless(thermo3["xkt"])
        fase.ks = unidades.Dimensionless(thermo3["xisenk"])

        dh = refprop.dhd1(T, rho, x)
        fase.dhdT_rho = unidades.SpecificHeat(dh["dhdt_D"]/fase.M, "kJkgK")
        fase.dhdT_P = unidades.SpecificHeat(dh["dhdt_p"]/fase.M, "kJkgK")
        fase.dhdP_T = unidades.EnthalpyPressure(dh["dhdp_t"]/fase.M, "kJkgkPa")
        # dhdtp_D : fix in library
        fase.dhdP_rho = unidades.EnthalpyPressure(
            dh["dhdtp_D"]/fase.M, "kJkgkPa")
        fase.dhdrho_T = unidades.EnthalpyDensity(
            dh["dhdD_t"]/fase.M**2, "kJkgkgm3")
        fase.dhdrho_P = unidades.EnthalpyDensity(
            dh["dhdD_p"]/fase.M**2, "kJkgkgm3")

        fase.dpdT_rho = unidades.PressureTemperature(thermo["dpdt"], "kPaK")
        fase.dpdrho_T = unidades.PressureDensity(
            thermo["dpdD"]/fase.M, "kPakgm3")
        # TODO: Add unit for derivative d^2p/dD^2 [kPa-L^2/mol^2]
        # MPa·m⁶/kg²
        fase.d2pdrho2 = unidades.Dimensionless(thermo["d2pdD2"]/fase.M**2/1000)
        fase.drhodP_T = unidades.DensityPressure(
            thermo["dDdp"]*fase.M, "kgm3kPa")
        fase.drhodT_P = unidades.DensityTemperature(thermo["dDdt"]*fase.M)
        fase.Gruneisen = unidades.Dimensionless(fase.v/fase.cv*fase.dpdT_rho)

        fase.Z_rho = unidades.SpecificVolume((fase.Z-1)/fase.rho)
        fase.IntP = unidades.Pressure(thermo3["pint"], "kPa")
        fase.hInput = unidades.Enthalpy(thermo3["spht"]/fase.M, "kJkg")
        fase.invT = unidades.InvTemperature(-1/self.T)

        fpv = refprop.fpv(T, rho, self.P.kPa, x)["Fpv"]
        fase.fpv = unidades.Dimensionless(fpv)

        chempot = refprop.chempot(T, rho, x)["u"]
        fase.chempot = [unidades.Enthalpy(c/fase.M) for c in chempot]
        fi = refprop.fugcof(T, rho, x)["f"]
        fase.fi = [unidades.Dimensionless(f) for f in fi]
        f = refprop.fgcty(T, rho, x)["f"]
        fase.f = [unidades.Pressure(f_i, "kPa") for f_i in f]

        b = refprop.virb(T, x)["b"]
        fase.virialB = unidades.SpecificVolume(b/self.M)
        c = refprop.virc(T, x)["c"]
        fase.virialC = unidades.SpecificVolume_square(c/self.M**2)

        # viriald don't supported for windows
        if sys.platform != "win32":
            d = refprop.vird(T, x)["d"]
            fase.virialD = unidades.Dimensionless(d/self.M**3)
        else:
            fase.virialD = unidades.Dimensionless(0)

        ba = refprop.virba(T, x)["ba"]
        fase.virialBa = unidades.SpecificVolume(ba/self.M)
        ca = refprop.virca(T, x)["ca"]
        fase.virialCa = unidades.SpecificVolume_square(ca/self.M**2)
        dcdt = refprop.dcdt(T, x)["dct"]
        fase.dCdt = unidades.Dimensionless(dcdt/self.M**2)
        dcdt2 = refprop.dcdt2(T, x)["dct2"]
        fase.dCdt2 = unidades.Dimensionless(dcdt2/self.M**2)
        dbdt = refprop.dbdt(T, x)["dbt"]
        fase.dBdt = unidades.Dimensionless(dbdt/self.M)

        b12 = refprop.b12(T, x)["b"]
        fase.b12 = unidades.SpecificVolume(b12*fase.M)
        try:
            cstar = refprop.cstar(T, self.P.kPa, 0, x)["cs"]
            fase.cstar = unidades.Dimensionless(cstar)
        except refprop.RefpropdllError:
            fase.cstar = unidades.Dimensionless(None)

        fase.fraccion = [unidades.Dimensionless(xi) for xi in x]
        xg = refprop.xmass(x)["xkg"]
        fase.fraccion_masica = [unidades.Dimensionless(xi) for xi in xg]

        try:
            transport = refprop.trnprp(T, rho, x)
            fase.mu = unidades.Viscosity(transport["eta"], "muPas")
            fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
            fase.k = unidades.ThermalConductivity(transport["tcx"])
            fase.alfa = unidades.Diffusivity(fase.k/fase.rho/fase.cp)
            fase.Prandt = unidades.Dimensionless(fase.mu*fase.cp/fase.k)
        except refprop.RefpropdllError:
            fase.mu = unidades.Viscosity(None)
            fase.nu = unidades.Diffusivity(None)
            fase.k = unidades.ThermalConductivity(None)
            fase.alfa = unidades.Diffusivity(None)
            fase.Prandt = unidades.Dimensionless(None)

        dielec = refprop.dielec(T, rho, x)
        fase.epsilon = unidades.Dimensionless(dielec["de"])


if __name__ == '__main__':
    # fluido = RefProp(ids=[140], T=300, P=1e6,)
    # fluido2 = RefProp(ids=[140], T=300, P=1e6, preos=True)
    # print(fluido.rho, fluido2.rho)
    # from lib.mEoS import Acetone
    # fluido = Acetone(T=300, P=1e6)
    # from pprint import pprint
    # pprint(fluido._Helmholtz(fluido.Tc/fluido.T, fluido.rho/fluido.rhoc))
    # pprint(fluido._PengRobinson(fluido.rho, fluido.T))
    # fluido2 = Acetone(T=300, P=1e6, eq="PR")
    # print(fluido.rho, fluido2.rho)

    # fluido = RefProp(ids=[47], T=100, P=6e4, htype="EOS", hcomp="BWR")
    # print(fluido.rho, fluido.cpM.JmolK)

    

    # See fortran/Manual.txt and fortran/Prop_sub.for file for explanation of refprop functionality

    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g"
          # % (fluido.rho, fluido.u.kJkg, fluido.h.kJkg, fluido.s.kJkgK,
             # fluido.cv.kJkgK, fluido.cp.kJkgK, fluido.cp0.kJkgK,
             # fluido.cp_cv, fluido.w, fluido.Z, fluido.a, fluido.g))
    # fluido = RefProp(ids=[62], T=300, P=1e6, htype="EOS", hcomp="FEK")
    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g"
          # % (fluido.rho, fluido.u.kJkg, fluido.h.kJkg, fluido.s.kJkgK,
             # fluido.cv.kJkgK, fluido.cp.kJkgK, fluido.cp0.kJkgK,
             # fluido.cp_cv, fluido.w, fluido.Z, fluido.a, fluido.g))
    # 996.960022694 112.478998007 113.482047254 0.392813902092 4.12721730499 4.17810360562 1.86484159263 1.01232944545 0.00724456495227 -5.36517262056 -4.3621233736

    # print("%0.12g %0.12g %0.12g %0.12g %0.12g"
    #       % (fluido.virialb, fluido.virialc, fluido.dbdt, fluido.virialba, fluido.virialca))
    # -0.0666822898709 -0.0130427067934 0.00115953207464 -0.120787648647 0.012705502088

    # print("%0.12g %0.12g"
    #       % (fluido.kappa.MPa, fluido.alfav))
    # 0.000449474822287 0.000275696573736
    # print("%0.12g %0.12g"
    #       % (fluido.IntP.MPa, fluido.hInput))
    # 183.012469708 15154.7171914
    # print("%0.12g %0.12g" % (fluido.Kt.MPa, fluido.Ks.MPa))
    # 2224.81872268 2252.24950376
    # print("%0.12g" % (fluido.fpv))
    # 21.6356265292
    # print("%0.12g %0.12g %0.12g"
          # % (fluido.chempot[0], fluido.fi[0], fluido.f[0].MPa))
    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g"
          # % (fluido.rho0, fluido.h0.kJkg, fluido.s0.kJkgK, fluido.g0.kJkg, fluido.a0.kJkg, fluido.P_Pideal.MPa))
    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g"
          # % (fluido.dhdT_rho.kJkgK, fluido.dhdT_P.kJkgK, fluido.dhdrho_T.kJkgkgm3, fluido.dhdrho_P.kJkgkgm3, fluido.dhdP_T.kJkgMPa, fluido.dhdP_rho.kJkgMPa))
    # print("%0.12g %0.12g" % (fluido.Gruneisen, fluido.Z_rho))

    # fluido = RefProp(ids=[62], T=300, x=.5)
    # print("%0.12g %0.12g" % (fluido.csat[0].kJkgK, fluido.dpdt_sat[0].MPaK))

#    fluido=RefProp(fluido=[u'hydrogen', u'nitrogen', u'oxygen', u'water'], fraccionMolar=[0.03, 0.95, 0.01, 0.01], T=300, P=1e5)
#    print(fluido.rho, fluido.cp.kJkgK, fluido.cv.kJkgK)

    # from lib.corriente import Corriente
    # kw = {"MEoS": True, "refprop": True, "ids": [5, 6, 7, 8, 10]}
    # entrada = Corriente(T=300, x=0.7, caudalMasico=0.01,
                        # fraccionMolar=[.3, 0.25, 0.05, 0.15, 0.25], **kw)

    # fluido = RefProp(ids=[5, 6, 7, 8, 10], fraccionMolar=[0.3, 0.25, 0.15, 0.05, 0.25], T=300, x=0.5)
    # import pprint
    # pprint.pprint(fluido.__dict__)
    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g" % (fluido.vE, fluido.uE.kJkg, fluido.hE.kJkg, fluido.sE.kJkgK, fluido.aE.kJkg, fluido.gE.kJkg))
    # print(fluido.M, fluido.Tc, fluido.Pc.MPa, fluido.rhoc)

#     # Code for advanced thermal test
    # from lib.mEoS.H2O import H2O
    # from lib.coolProp import CoolProp
    # from lib.iapws97 import IAPWS97
    # from lib.freeSteam import Freesteam

    # P = 2e4
    # T = 300

    # # Fluid definition, r refprop used as reference fluid, choose a m fluid to
    # # test against reference
    # r = RefProp(ids=[62], T=T, P=P)

    # # m = H2O(T=T, P=P)
    # # ierr = 1e-5
    # m = CoolProp(ids=[62], fraccionMolar=[1], T=T, P=P)
    # ierr = 1e-5
    # # m = IAPWS97(T=T, P=P)
    # # ierr = 1
    # # m = Freesteam(T=T, P=P)
    # # ierr = 1

    # for prop, key, unit in m.properties():
        # if key == "sigma":
            # try:
                # error = abs(r.Liquido.__getattribute__(key)-m.Liquido.__getattribute__(key))/r.Liquido.__getattribute__(key)*100
                # if error > ierr:
                    # msg = "ERROR %0.5g%s" % (error, "%")
                    # print("%s: %0.9g %0.9g %s" % (key, r.Liquido.__getattribute__(key), m.Liquido.__getattribute__(key), msg))
            # except ZeroDivisionError:
                # print(key, "no implementado")
        # elif key == "n":
            # pass
        # elif isinstance(r.__getattribute__(key), list):
            # error = 0
            # for v1, v2 in zip(r.__getattribute__(key), m.__getattribute__(key)):
                # error2 = abs(v1-v2)/v1*100
                # if error2 > error:
                    # error = error2
                # if error > ierr:
                    # msg = "ERROR %0.5g%s" % (error, "%")
                    # print("%s: " % key, r.__getattribute__(key), m.__getattribute__(key), msg)
        # else:
            # try:
                # error = abs(r.__getattribute__(key)-m.__getattribute__(key))/r.__getattribute__(key)*100
                # if error > ierr:
                    # msg = "ERROR %0.5g%s" % (error, "%")
                    # print("%s: %0.9g %0.9g %s" % (key, r.__getattribute__(key), m.__getattribute__(key), msg))
            # except ZeroDivisionError:
                # if r.__getattribute__(key) == m.__getattribute__(key):
                    # pass

    # st = RefProp(ids=[62], T=300, P=1e5)
    # st2 = st._new(T=300, P=2e5)
    # print(st.rho, st2.rho)
    # print(st.P, st2.P)


    # st = RefProp(ids=[62], T=286.04443, x=0.5)
    # from pprint import pprint
    # pprint(st.__dict__)
    # print(st.rhoM)

    T = 463
    P = 1e5

    st = RefProp(ids=[62], T=T, P=P)
    # st = RefProp(ids=[62], T=T, x=0.5)
    print(st.T, st.P.bar, st.v, st.u)
    st = RefProp(ids=[62], u=st.u, v=st.v)
    print(st.status, st.msg)
    print(st.T, st.P.bar, st.v, st.u)
