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
# Library to multiparameter equation of state calculation using coolprop
# http://coolprop.sourceforge.net/index.html
# optional method to meos tools calculations and to multicomponent streams
###############################################################################


from PyQt5.QtWidgets import QApplication

try:
    import CoolProp as CP
except ImportError as e:
    raise e

from lib import unidades
from lib.thermo import Fluid, Thermo
from lib.compuestos import Componente


noIds = {
    'Deuterium',
    'HeavyWater',
    'ParaDeuterium',
    'ParaHydrogen',
    'OrthoDeuterium',
    'OrthoHydrogen',
    'D4',
    'D5',
    'D6',
    'MD2M',
    'MD3M',
    'MD4M',
    'MDM',
    'MM',
    'HFE143m',
    'Novec649',
    'Krypton',
    'Xenon',
    'MethylLinoleate',
    'MethylLinolenate',
    'MethylOleate',
    'MethylPalmitate',
    'MethylStearate',
    'DimethylCarbonate',
    'R1233zd(E)',
    'R1234yf',
    'R1234ze(Z)',
    'R13I1',
    'R365MFC',
    'R404A',
    'R407C',
    'R410A',
    'R507A',
    'SES36'}


__all__ = {
    24: '1-Butene',
    140: 'Acetone',
    475: 'Air',
    63: 'Ammonia',
    98: 'Argon',
    40: 'Benzene',
    49: 'CarbonDioxide',
    48: 'CarbonMonoxide',
    219: 'CarbonylSulfide',
    25: 'cis-2-Butene',
    38: 'CycloHexane',
    36: 'Cyclopentane',
    258: 'CycloPropane',
    126: 'Dichloroethane',  # or maybe 127, 1,2 dichloro
    162: 'DiethylEther',
    162: 'DimethylEther',
    3: 'Ethane',
    134: 'Ethanol',
    45: 'EthylBenzene',
    22: 'Ethylene',
    129: 'EthyleneOxide',
    208: 'Fluorine',
    212: 'Helium',
    1: 'Hydrogen',
    104: 'HydrogenChloride',
    50: 'HydrogenSulfide',
    5: 'IsoButane',
    27: 'IsoButene',
    52: 'Isohexane',
    7: 'Isopentane',
    43: 'm-Xylene',
    2: 'Methane',
    117: 'Methanol',
    6: 'n-Butane',
    14: 'n-Decane',
    16: 'n-Dodecane',
    11: 'n-Heptane',
    10: 'n-Hexane',
    13: 'n-Nonane',
    12: 'n-Octane',
    8: 'n-Pentane',
    4: 'n-Propane',
    15: 'n-Undecane',
    107: 'Neon',
    9: 'Neopentane',
    46: 'Nitrogen',
    110: 'NitrousOxide',
    42: 'o-Xylene',
    47: 'Oxygen',
    44: 'p-Xylene',
    23: 'Propylene',
    66: 'Propyne',
    217: 'R11',
    232: 'R113',
    231: 'R114',
    229: 'R115',
    236: 'R116',
    216: 'R12',
    1631: 'R123',
    671: 'R1234ze(E)',
    1629: 'R124',
    1231: 'R125',
    215: 'R13',
    1235: 'R134a',
    218: 'R14',
    1633: 'R141b',
    241: 'R142b',
    243: 'R143a',
    245: 'R152A',
    247: 'R161',
    642: 'R21',
    671: 'R218',
    220: 'R22',
    1872: 'R227EA',
    643: 'R23',
    693: 'R236EA',
    1873: 'R236FA',
    1817: 'R245fa',
    645: 'R32',
    115: 'R40',
    225: 'R41',
    692: 'RC318',
    51: 'SulfurDioxide',
    953: 'SulfurHexafluoride',
    41: 'Toluene',
    26: 'trans-2-Butene',
    62: 'Water'
}


class CoolProp(Thermo):
    """Stream class using coolProp external library
    Parameters needed to define it are:

        -ids: index of fluid
        -fraccionMolar: molar fraction

        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -rho: Density, kg/m3
        -H: Enthalpy, J/kg
        -S: Entropy, J/kgK
        -x: Quality, -
    """
    kwargs = {"ids": [],
              "fraccionMolar": [],

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "Q": None,
              "rho": 0.0,
              "Dmass": 0.0,
              "H": None,
              "Hmass": None,
              "S": None,
              "Smass": None}

    __doi__ = [
        {"autor": "Bell, Ian H. and Wronski, Jorrit and Quoilin, Sylvain and Lemort, Vincent",
         "title": "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp",
         "ref": "Ind. Eng. Chem. Res., 2014, 53 (6), pp 2498–2508",
         "doi": "10.1021/ie4033999"}]

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
        else:
            self._definition = True

        # Update the kwargs with the special coolprop namespace
        if self.kwargs["x"] != CoolProp.kwargs["x"]:
            self.kwargs["Q"] = self.kwargs["x"]
        if self.kwargs["rho"] != CoolProp.kwargs["rho"]:
            self.kwargs["Dmass"] = self.kwargs["rho"]
        if self.kwargs["S"] != CoolProp.kwargs["S"]:
            self.kwargs["Smass"] = self.kwargs["S"]
        if self.kwargs["H"] != CoolProp.kwargs["H"]:
            self.kwargs["Hmass"] = self.kwargs["H"]

        # Check thermo definition
        self._thermo = ""
        for def_ in ["P-T", "Q-T", "P-Q", "Dmass-T", "Dmass-P", "Hmass-P",
                     "P-Smass", "Hmass-Smass"]:
            inputs = def_.split("-")
            if self.kwargs[inputs[0]] != CoolProp.kwargs[inputs[0]] and \
                    self.kwargs[inputs[1]] != CoolProp.kwargs[inputs[1]]:
                self._thermo = def_.replace("-", "")
                self._par = CP.__getattribute__("%s_INPUTS" % self._thermo)
                break

        return self._definition and self._thermo

    def args(self):
        var1 = self.kwargs[self._thermo[0]]
        var2 = self.kwargs[self._thermo[1]]

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
        if self._multicomponent:
            estado.set_mole_fractions(self.kwargs["fraccionMolar"])
        estado.update(self._par, *args)

        self.M = unidades.Dimensionless(estado.molar_mass()*1000)

        # Disabled CoolProp, see issue #1087
        # self.Tc = unidades.Temperature(estado.T_critical())
        # self.Pc = unidades.Pressure(estado.p_critical())
        # self.rhoc = unidades.Density(estado.rhomass_critical())

        # Calculate critical properties with mezcla method
        # Coolprop for mixtures can fail and it's slow
        Cmps = [Componente(int(i)) for i in self.kwargs["ids"]]

        # Calculate critic temperature, API procedure 4B1.1 pag 304
        V = sum([xi*cmp.Vc for xi, cmp in zip(self.kwargs["fraccionMolar"], Cmps)])
        k = [xi*cmp.Vc/V for xi, cmp in zip(self.kwargs["fraccionMolar"], Cmps)]
        Tcm = sum([ki*cmp.Tc for ki, cmp in zip(k, Cmps)])
        self.Tc = unidades.Temperature(Tcm)

        # Calculate pseudocritic temperature
        tpc = sum([x*cmp.Tc for x, cmp in zip(self.kwargs["fraccionMolar"], Cmps)])

        # Calculate pseudocritic pressure
        ppc = sum([x*cmp.Pc for x, cmp in zip(self.kwargs["fraccionMolar"], Cmps)])

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

        self.Tt = unidades.Temperature(estado.Ttriple())
        estado2 = CP.AbstractState("HEOS", fluido)
        if self._multicomponent:
            estado2.set_mole_fractions(self.kwargs["fraccionMolar"])
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
        else:
            self.name = fluido
            self.CAS = estado.fluid_param_string("CAS")
            self.synonim = estado.fluid_param_string("aliases")

        self.P = unidades.Pressure(estado.p())
        self.T = unidades.Temperature(estado.T())
        self.Tr = unidades.Dimensionless(self.T/self.Tc)
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.rho = unidades.Density(estado.rhomass())
        self.v = unidades.SpecificVolume(1./self.rho)

        cp0 = self._prop0(estado)
        self._cp0(cp0)

        self.Liquido = Thermo()
        self.Gas = Thermo()
        if self.x == 0:
            # liquid phase
            self.fill(self.Liquido, estado)
            self.fill(self, estado)
        elif self.x == 1:
            # vapor phase
            self.fill(self.Gas, estado)
            self.fill(self, estado)
        else:
            # Two phase
            liquido = CP.AbstractState("HEOS", fluido)
            if self._multicomponent:
                xi = estado.mole_fractions_liquid()
                liquido.set_mole_fractions(xi)
            liquido.specify_phase(CP.iphase_liquid)
            liquido.update(self._par, *args)
            self.fill(self.Liquido, liquido)

            vapor = CP.AbstractState("HEOS", fluido)
            if self._multicomponent:
                yi = estado.mole_fractions_vapor()
                vapor.set_mole_fractions(yi)
            vapor.specify_phase(CP.iphase_gas)
            vapor.update(self._par, *args)
            self.fill(self.Gas, vapor)
            self.fill(self, estado)

            ['Tmax', 'Tmin', 'all_critical_points',
             'apply_simple_mixing_rule', 'build_phase_envelope', 'change_EOS',
             'conductivity_contributions', 'conformal_state', 'criticality_contour_values',
             'd2alpha0_dDelta2', 'd2alpha0_dDelta_dTau', 'd2alpha0_dTau2', 'd2alphar_dDelta2', 'd2alphar_dDelta_dTau',
             'd2alphar_dTau2', 'd3alpha0_dDelta2_dTau', 'd3alpha0_dDelta3', 'd3alpha0_dDelta_dTau2', 'd3alpha0_dTau3', 'd3alphar_dDelta2_dTau',
             'd3alphar_dDelta3', 'd3alphar_dDelta_dTau2', 'd3alphar_dTau3', 'dalpha0_dDelta', 'dalpha0_dTau', 'dalphar_dDelta', 'dalphar_dTau',
             'first_saturation_deriv', 'first_two_phase_deriv', 'first_two_phase_deriv_splined',
             'fluid_param_string', 'get_binary_interaction_double',
             'get_binary_interaction_string', 'get_phase_envelope_data', 'has_melting_line',
             'ideal_curve',
             'melting_line', 'pmax',
             'rhomass', 'rhomass_reducing', 'rhomolar_critical', 'rhomolar_reducing',
             'saturated_liquid_keyed_output', 'saturated_vapor_keyed_output', 'saturation_ancillary', 'second_partial_deriv',
             'second_saturation_deriv', 'second_two_phase_deriv', 'set_binary_interaction_double', 'set_binary_interaction_string',
             'tangent_plane_distance', 'trivial_keyed_output', 'true_critical_point',
             'update', 'update_with_guesses', 'viscosity_contributions']


        # Calculate special properties useful only for one phase
        if self._multicomponent:
            self.Liquido.sigma = unidades.Tension(None)
        elif x < 1 and self.Tt <= self.T <= self.Tc:
            self.Liquido.sigma = unidades.Tension(estado.surface_tension())
        else:
            self.Liquido.sigma = unidades.Tension(None)

        self.virialB = unidades.SpecificVolume(estado.Bvirial())
        self.virialC = unidades.SpecificVolume_square(estado.Cvirial())
        self.invT = unidades.InvTemperature(-1/self.T)

        if self.Tt <= self.T <= self.Tc:
            self.Hvap = unidades.Enthalpy(self.Gas.h-self.Liquido.h)
            self.Svap = unidades.SpecificHeat(self.Gas.s-self.Liquido.s)
        else:
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)

    def _prop0(self, estado):
        """Ideal gas properties"""
        R = estado.gas_constant()/self.M
        tau = self.Tc/self.T
        fio = estado.alpha0()
        fiot = estado.dalpha0_dTau()
        fiott = estado.d2alpha0_dTau2()

        propiedades = Fluid()
        propiedades.v = R*self.T/self.P
        propiedades.h = R*self.T*(1+tau*fiot)
        propiedades.s = R*(tau*fiot-fio)
        propiedades.cv = -R*tau**2*fiott
        propiedades.cp = R*(-tau**2*fiott+1)
        propiedades.alfap = 1/self.T
        propiedades.betap = self.rho
        propiedades.w = (R*self.T*1000/(1+tau**2/abs(fiott)))**0.5
        propiedades.gamma = 1/propiedades.v/self.P/1000 * \
            estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iSmass)
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
            fase.fi = unidades.Dimensionless([1])
            fase.f = unidades.Pressure([self.P])

        fase.cv = unidades.SpecificHeat(estado.cvmass())
        fase.cp = unidades.SpecificHeat(estado.cpmass())
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.w = unidades.Speed(estado.speed_sound())

        fase.rhoM = unidades.MolarDensity(estado.rhomolar())
        fase.hM = unidades.MolarEnthalpy(estado.hmolar())
        fase.sM = unidades.MolarSpecificHeat(estado.smolar())
        fase.uM = unidades.MolarEnthalpy(estado.umolar())
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(estado.cvmolar())
        fase.cpM = unidades.MolarSpecificHeat(estado.cpmolar())

        fase.joule = unidades.TemperaturePressure(
            estado.first_partial_deriv(CP.iT, CP.iP, CP.iHmass))
        fase.Gruneisen = unidades.Dimensionless(
            estado.first_partial_deriv(CP.iP, CP.iT, CP.iDmass))
        fase.alfav = unidades.InvTemperature(
            estado.isobaric_expansion_coefficient())
        fase.kappa = unidades.InvPressure(estado.isothermal_compressibility())
        fase.alfap = unidades.Density(fase.alfav/self.P/fase.kappa)
        fase.betap = unidades.Density(
            fase.rho/self.P*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        fase.betas = unidades.TemperaturePressure(
            estado.first_partial_deriv(CP.iT, CP.iP, CP.iSmass))
        fase.gamma = unidades.Dimensionless(
            fase.rho/self.P*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iSmass))

        fase.kt = unidades.Dimensionless(
            fase.rho/self.P*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        fase.Kt = unidades.Pressure(
            fase.rho*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iSmass))
        fase.Ks = unidades.Pressure(
            fase.rho*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        fase.ks = unidades.InvPressure(1/fase.Ks)
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

        fase.mu = unidades.Viscosity(estado.viscosity())
        fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        fase.k = unidades.ThermalConductivity(estado.conductivity())
        fase.alfa = unidades.Diffusivity(fase.k/1000/fase.rho/fase.cp)
        fase.Prandt = unidades.Dimensionless(estado.Prandtl())
        fase.fraccion = estado.get_mole_fractions()
        fase.fraccion_masica = estado.get_mass_fractions()

    def getphase(self, estado):
        """Return fluid phase with translation support"""
        phase = estado.phase()
        if phase == CP.iphase_supercritical:
            if self.T > self.Tc:
                return QApplication.translate("pychemqt", "Supercritical fluid"), 1.
            else:
                return QApplication.translate("pychemqt", "Compressible liquid"), 1.
        elif phase == CP.iphase_gas:
            # if estado.superheat > 0:
                return QApplication.translate("pychemqt", "Vapor"), 1.
            # else:
                # return QApplication.translate("pychemqt", "Saturated vapor"), 1.
        elif phase == CP.iphase_liquid:
            # if estado.subcooling > 0:
                return QApplication.translate("pychemqt", "Liquid"), 0.
            # else:
                # return QApplication.translate("pychemqt", "Saturated liquid"), 0.
        elif phase == CP.iphase_twophase:
            return QApplication.translate("pychemqt", "Two phases"), estado.Q()

if __name__ == '__main__':
    fluido = CoolProp(ids=[62], fraccionMolar=[1], T=300, P=101325)
    print(fluido.P, fluido.Pc)
