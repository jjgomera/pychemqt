#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Library to multiparameter equation of state calculation using coolprop
# http://coolprop.sourceforge.net/index.html
# optional method to meos tools calculations and to multicomponent streams
###############################################################################


__doi__ = {
    1:
        {"autor": "Bell, I.H., Wronski, J., Quoilin, S., Lemort, V.",
         "title": "Pure and Pseudo-pure Fluid Thermophysical Property"
                  "Evaluation and the Open-Source Thermophysical Property"
                  "Library CoolProp",
         "ref": "Ind. Eng. Chem. Res. 53(6) (2014) 2498-2508",
         "doi": "10.1021/ie4033999"}}


import os

from tools.qt import QtWidgets

try:
    import CoolProp as CP
except ImportError as e:
    pass

from lib import unidades, mEoS
from lib.thermo import ThermoAdvanced
from lib.compuestos import Componente


# Automatic loading of coolProp name from meos subclass _coolPropName property
__all__ = {}
noIds = []
for cmp in mEoS.__all__:
    if cmp.id and cmp._coolPropName:
        __all__[cmp.id] = cmp._coolPropName
    elif cmp._coolPropName:
        noIds.append(cmp._coolPropName)


class CoolProp(ThermoAdvanced):
    """Stream class using coolProp external library
    Parameters needed to define it are:

        -ids: index of fluid
        -fraccionMolar: molar fraction

        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -rho: Density, kg/m3
        -h: Enthalpy, J/kg
        -s: Entropy, J/kgK
        -x: Quality, -
    """
    kwargs = {"ids": [],
              "fraccionMolar": [],

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "Q": None,
              "rho": 0.0,
              "v": 0.0,
              "Dmass": 0.0,
              "h": None,
              "Hmass": None,
              "u": None,
              "Umass": None,
              "s": None,
              "Smass": None}

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)

        if self.calculable:
            # try:
                self.calculo()
            # except ValueError as e:
                # self.msg = e
                # self.status = 0
            # else:
                self.status = 1
                self.msg = "Solved"

        elif self._definition and not self._multicomponent and "ids" in kwargs:
            if os.environ["CoolProp"] == "True":
                fluido = self._name()
                estado = CP.AbstractState(b"HEOS", fluido.encode())
                self.Tc = unidades.Temperature(estado.T_critical())
                self.Pc = unidades.Pressure(estado.p_critical())
                self.rhoc = unidades.Density(estado.rhomass_critical())

                self.M = unidades.Dimensionless(estado.molar_mass()*1000)
                self.R = unidades.SpecificHeat(estado.gas_constant()/self.M)
                self.Tt = unidades.Temperature(estado.Ttriple())
                self.f_accent = unidades.Dimensionless(
                        estado.acentric_factor())

                self.name = fluido
                self.CAS = estado.fluid_param_string(b"CAS")
                self.synonim = estado.fluid_param_string(b"aliases")
                self.formula = estado.fluid_param_string(b"formula")

            self.eq = self._limit(fluido, estado)

    def _limit(self, name, estado):
        eq = {}
        eq["Tmin"] = estado.Tmin()
        eq["Tmax"] = estado.Tmax()
        eq["Pmin"] = CP.CoolProp.PropsSI(b"pmin", name.encode())
        eq["Pmax"] = estado.pmax()
        return eq

    @property
    def calculable(self):
        """Check in the class is fully defined"""
        # Check mix state
        self._multicomponent = False
        if len(self.kwargs["ids"]) > 1:
            self._multicomponent = True

        # Check supported fluid
        COOLPROP_available = True
        for id in self.kwargs["ids"]:
            if id not in __all__ and id not in noIds:
                COOLPROP_available = False
                if not COOLPROP_available:
                    raise(ValueError)

        # Check correct fluid definition
        if self._multicomponent:
            if self.kwargs["ids"] and len(self.kwargs["fraccionMolar"]) == \
                    len(self.kwargs["ids"]):
                self._definition = True
            else:
                self._definition = False
        elif self.kwargs["ids"]:
            self._definition = True
        else:
            self._definition = False

        # Update the kwargs with the special coolprop namespace
        if self.kwargs["x"] != CoolProp.kwargs["x"]:
            self.kwargs["Q"] = self.kwargs["x"]
        if self.kwargs["v"] != CoolProp.kwargs["v"]:
            self.kwargs["Dmass"] = 1/self.kwargs["v"]
        if self.kwargs["rho"] != CoolProp.kwargs["rho"]:
            self.kwargs["Dmass"] = self.kwargs["rho"]
        if self.kwargs["s"] != CoolProp.kwargs["s"]:
            self.kwargs["Smass"] = self.kwargs["s"]
        if self.kwargs["h"] != CoolProp.kwargs["h"]:
            self.kwargs["Hmass"] = self.kwargs["h"]
        if self.kwargs["u"] != CoolProp.kwargs["u"]:
            self.kwargs["Umass"] = self.kwargs["u"]

        # Check thermo definition
        self._thermo = ""
        for def_ in ["P-T", "Q-T", "P-Q", "Dmass-T", "Dmass-P", "Hmass-P",
                     "P-Smass", "Hmass-Smass", "Hmass-T", "T-Umass", "P-Umass",
                     "Dmass-Hmass", "Dmass-Smass", "Dmass-Umass",
                     "Hmass-Umass", "Smass-Umass"]:
            inputs = def_.split("-")
            if self.kwargs[inputs[0]] != CoolProp.kwargs[inputs[0]] and \
                    self.kwargs[inputs[1]] != CoolProp.kwargs[inputs[1]]:
                self._thermo = def_.replace("-", "")
                self._inputs = inputs
                self._par = CP.__getattribute__("%s_INPUTS" % self._thermo)
                break

        return self._definition and self._thermo

    def args(self):
        var1 = self.kwargs[self._inputs[0]]
        var2 = self.kwargs[self._inputs[1]]

        args = [var1, var2]
        return args

    def _name(self):
        lst = []
        for fld in self.kwargs["ids"]:
            if fld in __all__:
                lst.append(__all__[fld])
            elif fld in noIds:
                lst.append(fld)
        name = "&".join(lst)
        return name

    def calculo(self):
        fluido = self._name()
        args = self.args()
        estado = CP.AbstractState("HEOS", fluido)
        self.eq = self._limit(fluido, estado)
        if self._multicomponent:
            estado.set_mole_fractions(self.kwargs["fraccionMolar"])
        estado.update(self._par, *args)

        self.M = unidades.Dimensionless(estado.molar_mass()*1000)

        if self._multicomponent:
            # Disabled CoolProp critical properties for multicomponent,
            # see issue #1087

            # Calculate critical properties with mezcla method
            # Coolprop for mixtures can fail and it's slow
            Cmps = [Componente(int(i)) for i in self.kwargs["ids"]]

            # Calculate critic temperature, API procedure 4B1.1 pag 304
            V = sum([xi*cmp.Vc for xi, cmp in
                     zip(self.kwargs["fraccionMolar"], Cmps)])
            k = [xi*cmp.Vc/V for xi, cmp in
                 zip(self.kwargs["fraccionMolar"], Cmps)]
            Tcm = sum([ki*cmp.Tc for ki, cmp in zip(k, Cmps)])
            self.Tc = unidades.Temperature(Tcm)

            # Calculate pseudocritic temperature
            tpc = sum([x*cmp.Tc for x, cmp in
                       zip(self.kwargs["fraccionMolar"], Cmps)])

            # Calculate pseudocritic pressure
            ppc = sum([x*cmp.Pc for x, cmp in
                       zip(self.kwargs["fraccionMolar"], Cmps)])

            # Calculate critic pressure, API procedure 4B2.1 pag 307
            sumaw = 0
            for xi, cmp in zip(self.kwargs["fraccionMolar"], Cmps):
                sumaw += xi*cmp.f_acent
            pc = ppc+ppc*(5.808+4.93*sumaw)*(self.Tc-tpc)/tpc
            self.Pc = unidades.Pressure(pc)

            # Calculate critic volume, API procedure 4B3.1 pag 314
            sumaxvc23 = sum([xi*cmp.Vc**(2./3) for xi, cmp in
                             zip(self.kwargs["fraccionMolar"], Cmps)])
            k = [xi*cmp.Vc**(2./3)/sumaxvc23 for xi, cmp in
                 zip(self.kwargs["fraccionMolar"], Cmps)]

            # TODO: Calculate C value from component type.
            # For now it suppose all are hidrycarbon (C=0)
            C = 0

            V = [[-1.4684*abs((cmpi.Vc-cmpj.Vc)/(cmpi.Vc+cmpj.Vc))+C
                  for cmpj in Cmps] for cmpi in Cmps]
            v = [[V[i][j]*(cmpi.Vc+cmpj.Vc)/2. for j, cmpj in enumerate(
                Cmps)] for i, cmpi in enumerate(Cmps)]
            suma1 = sum([ki*cmp.Vc for ki, cmp in zip(k, Cmps)])
            suma2 = sum([ki*kj*v[i][j] for j, kj in enumerate(k)
                         for i, ki in enumerate(k)])
            self.rhoc = unidades.Density((suma1+suma2)*self.M)

        else:
            self.Tc = unidades.Temperature(estado.T_critical())
            self.Pc = unidades.Pressure(estado.p_critical())
            self.rhoc = unidades.Density(estado.rhomass_critical())

        self.R = unidades.SpecificHeat(estado.gas_constant()/self.M)
        self.Tt = unidades.Temperature(estado.Ttriple())

        if self._multicomponent:
            estado2 = CP.AbstractState("HEOS", fluido)
            estado2.set_mole_fractions(self.kwargs["fraccionMolar"])
            estado2.update(CP.PQ_INPUTS, 101325, 1)
            self.Tb = unidades.Temperature(estado2.T())
        else:
            Pt = estado.trivial_keyed_output(CP.iP_triple)
            if Pt < 101325:
                estado2 = CP.AbstractState("HEOS", fluido)
                estado2.update(CP.PQ_INPUTS, 101325, 1)
                self.Tb = unidades.Temperature(estado2.T())

        self.f_accent = unidades.Dimensionless(estado.acentric_factor())

        # Dipole moment only available for REFPROP backend
        # self.momentoDipolar(estado.keyed_output(CP.idipole_moment))

        self.phase, x = self.getphase(estado)
        self.x = unidades.Dimensionless(x)
        if self._multicomponent:
            string = fluido.replace("&", " (%0.2f), ")
            string += " (%0.2f)"
            self.name = string % tuple(self.kwargs["fraccionMolar"])
            self.CAS = ""
            self.synonim = ""
            self.formula = ""
        else:
            self.name = fluido
            self.CAS = estado.fluid_param_string("CAS")
            self.synonim = estado.fluid_param_string("aliases")
            self.formula = estado.fluid_param_string("formula")

        self.P = unidades.Pressure(estado.p())
        self.T = unidades.Temperature(estado.T())
        self.Tr = unidades.Dimensionless(self.T/self.Tc)
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.rho = unidades.Density(estado.rhomass())
        self.v = unidades.SpecificVolume(1./self.rho)

        cp0 = self._prop0(estado)
        self._cp0(cp0)

        self.Liquido = ThermoAdvanced()
        self.Gas = ThermoAdvanced()
        if self.x == 0:
            # liquid phase
            self.fill(self.Liquido, estado)
            self.fill(self, estado)
            self.fillNone(self.Gas)
        elif self.x == 1:
            # vapor phase
            self.fill(self.Gas, estado)
            self.fill(self, estado)
            self.fillNone(self.Liquido)
        else:
            # Two phase
            liquido = CP.AbstractState("HEOS", fluido)
            if self._multicomponent:
                xi = estado.mole_fractions_liquid()
                liquido.set_mole_fractions(xi)
            liquido.specify_phase(CP.iphase_liquid)
            liquido.update(CP.QT_INPUTS, 0, self.T)
            self.fill(self.Liquido, liquido)

            vapor = CP.AbstractState("HEOS", fluido)
            if self._multicomponent:
                yi = estado.mole_fractions_vapor()
                vapor.set_mole_fractions(yi)
            vapor.specify_phase(CP.iphase_gas)
            vapor.update(CP.QT_INPUTS, 1, self.T)
            self.fill(self.Gas, vapor)
            self.fill(self, estado)

        # Calculate special properties useful only for one phase
        if self._multicomponent:
            self.sigma = unidades.Tension(None)
        elif x < 1 and self.Tt <= self.T <= self.Tc:
            try:
                self.sigma = unidades.Tension(estado.surface_tension())
            except ValueError:
                self.sigma = unidades.Tension(None)
        else:
            self.sigma = unidades.Tension(None)

        self.virialB = unidades.SpecificVolume(estado.Bvirial())
        self.virialC = unidades.SpecificVolume_square(estado.Cvirial())
        self.invT = unidades.InvTemperature(-1/self.T)

        if 0 < self.x < 1:
            self.Hvap = unidades.Enthalpy(self.Gas.h-self.Liquido.h)
            self.Svap = unidades.SpecificHeat(self.Gas.s-self.Liquido.s)
        else:
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)

    def _prop0(self, estado):
        """Ideal gas properties"""
        tau = self.Tc/self.T
        fio = estado.alpha0()
        fiot = estado.dalpha0_dTau()
        fiott = estado.d2alpha0_dTau2()

        propiedades = {}
        propiedades["v"] = self.R*self.T/self.P.kPa
        propiedades["h"] = self.R*self.T*(1+tau*fiot)*1000
        propiedades["s"] = self.R*(tau*fiot-fio)*1000
        propiedades["cv"] = -self.R*tau**2*fiott*1000
        propiedades["cp"] = self.R*(-tau**2*fiott+1)*1000
        propiedades["w"] = (self.R*self.T*1000/(1+tau**2/abs(fiott)))**0.5
        return propiedades

    def fill(self, fase, estado):
        fase._bool = True
        fase.M = unidades.Dimensionless(estado.molar_mass()*1000)
        fase.rho = unidades.Density(estado.rhomass())
        fase.v = unidades.SpecificVolume(1./fase.rho)
        fase.Z = unidades.Dimensionless(estado.keyed_output(CP.iZ))

        fase.h = unidades.Enthalpy(estado.hmass())
        fase.s = unidades.SpecificHeat(estado.smass())
        fase.u = unidades.Enthalpy(estado.umass())
        fase.a = unidades.Enthalpy(fase.u-self.T*fase.s)
        fase.g = unidades.Enthalpy(fase.h-self.T*fase.s)
        if self._multicomponent:
            fase.fi = []
            fase.f = []
            for i in range(len(self.kwargs["ids"])):
                fase.fi.append(unidades.Dimensionless(
                    estado.fugacity_coefficient(i)))
                fase.f.append(unidades.Pressure(estado.fugacity(i)))
        else:
            fase.fi = [unidades.Dimensionless(1)]
            # print(estado.alphar(), estado.delta(), estado.dalphar_dDelta(),
                    # (1+estado.delta()*estado.dalphar_dDelta()))
            # fase.fi = [unidades.Dimensionless(
                # exp(estado.alphar()+estado.delta()*estado.dalphar_dDelta() -
                    # log(1+estado.delta()*estado.dalphar_dDelta())))]
            fase.f = [unidades.Pressure(self.P*f) for f in fase.fi]

        fase.cv = unidades.SpecificHeat(estado.cvmass())
        fase.cp = unidades.SpecificHeat(estado.cpmass())
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.gamma = fase.cp_cv
        fase.w = unidades.Speed(estado.speed_sound())

        fase.rhoM = unidades.MolarDensity(estado.rhomolar(), "molm3")
        fase.hM = unidades.MolarEnthalpy(estado.hmolar(), "Jmol")
        fase.sM = unidades.MolarSpecificHeat(estado.smolar(), "JmolK")
        fase.uM = unidades.MolarEnthalpy(estado.umolar(), "Jmol")
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(estado.cvmolar(), "JmolK")
        fase.cpM = unidades.MolarSpecificHeat(estado.cpmolar(), "JmolK")

        fase.joule = unidades.TemperaturePressure(
            estado.first_partial_deriv(CP.iT, CP.iP, CP.iHmass))
        fase.Gruneisen = unidades.Dimensionless(
            fase.v/fase.cv*estado.first_partial_deriv(CP.iP, CP.iT, CP.iDmass))
        fase.alfav = unidades.InvTemperature(
            estado.isobaric_expansion_coefficient())
        fase.kappa = unidades.InvPressure(estado.isothermal_compressibility())
        fase.kappas = unidades.InvPressure(
            -1/fase.v*self.derivative("v", "P", "s", fase))
        fase.alfap = unidades.Density(fase.alfav/self.P/fase.kappa)
        fase.betap = unidades.Density(
            -1/self.P*self.derivative("P", "v", "T", fase))
        fase.betas = unidades.TemperaturePressure(
            estado.first_partial_deriv(CP.iT, CP.iP, CP.iSmass))

        fase.kt = unidades.Dimensionless(
            fase.rho/self.P*estado.first_partial_deriv(
                CP.iP, CP.iDmass, CP.iT))
        fase.Ks = unidades.Pressure(
            fase.rho*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iSmass))
        fase.Kt = unidades.Pressure(
            fase.rho*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        fase.ks = unidades.Dimensionless(
            fase.rho/self.P*estado.first_partial_deriv(
                CP.iP, CP.iDmass, CP.iSmass))
        fase.dhdT_rho = unidades.SpecificHeat(
            estado.first_partial_deriv(CP.iHmass, CP.iT, CP.iDmass))
        fase.dhdT_P = unidades.SpecificHeat(
            estado.first_partial_deriv(CP.iHmass, CP.iT, CP.iP))
        fase.dhdP_T = unidades.EnthalpyPressure(
            estado.first_partial_deriv(CP.iHmass, CP.iP, CP.iT))  # deltat
        fase.deltat = fase.dhdP_T
        fase.dhdP_rho = unidades.EnthalpyPressure(
            estado.first_partial_deriv(CP.iHmass, CP.iP, CP.iDmass))
        fase.dhdrho_T = unidades.EnthalpyDensity(
            estado.first_partial_deriv(CP.iHmass, CP.iDmass, CP.iT))
        fase.dhdrho_P = unidades.EnthalpyDensity(
            estado.first_partial_deriv(CP.iHmass, CP.iDmass, CP.iP))
        fase.dpdT_rho = unidades.PressureTemperature(
            estado.first_partial_deriv(CP.iP, CP.iT, CP.iDmass))
        fase.dpdrho_T = unidades.PressureDensity(
            estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        fase.drhodP_T = unidades.DensityPressure(
            estado.first_partial_deriv(CP.iDmass, CP.iP, CP.iT))
        fase.drhodT_P = unidades.DensityTemperature(
            estado.first_partial_deriv(CP.iDmass, CP.iT, CP.iP))

        fase.Z_rho = unidades.SpecificVolume((fase.Z-1)/fase.rho)
        fase.IntP = unidades.Pressure(
            self.T*estado.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)-self.P)
        fase.hInput = unidades.Enthalpy(
            -fase.rho*estado.first_partial_deriv(CP.iHmass, CP.iDmass, CP.iP))

        fase.virialB = unidades.SpecificVolume(estado.Bvirial())
        fase.virialC = unidades.SpecificVolume_square(estado.Cvirial())
        fase.invT = unidades.InvTemperature(-1/self.T)

        try:
            fase.mu = unidades.Viscosity(estado.viscosity())
        except ValueError:
            fase.mu = unidades.Viscosity(None)
            fase.Prandt = unidades.Dimensionless(None)

        try:
            fase.k = unidades.ThermalConductivity(estado.conductivity())
        except ValueError:
            fase.k = unidades.ThermalConductivity(None)

        fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        fase.alfa = unidades.Diffusivity(fase.k/fase.rho/fase.cp)
        fase.fraccion = estado.get_mole_fractions()
        fase.fraccion_masica = estado.get_mass_fractions()
        fase.epsilon = unidades.Dimensionless(None)

    def getphase(self, estado):
        """Return fluid phase with translation support"""
        phase = estado.phase()
        if phase == CP.iphase_supercritical:
            msg = QtWidgets.QApplication.translate("Supercritical fluid")
            x = 1
        elif phase == CP.iphase_supercritical_liquid:
            msg = QtWidgets.QApplication.translate("Supercritical liquid")
            x = 1
        elif phase == CP.iphase_supercritical_gas:
            msg = QtWidgets.QApplication.translate("Supercritical gas")
            x = 1
        elif phase == CP.iphase_gas:
            msg = QtWidgets.QApplication.translate("Vapor")
            x = 1
        elif phase == CP.iphase_liquid:
            msg = QtWidgets.QApplication.translate("Liquid")
            x = 0
        elif phase == CP.iphase_twophase:
            msg = QtWidgets.QApplication.translate("Two phases")
            x = estado.Q()
        elif phase == CP.iphase_critical_point:
            msg = QtWidgets.QApplication.translate("Critical point")
            x = 1

        return msg, x


if __name__ == '__main__':
    fluido = CoolProp(ids=[107], T=300, P=1e5)
    print(fluido.status, fluido.msg)
    print(fluido.rho, fluido.msg)
