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
# Library to multiparameter equation of state calculation using refprop
# https://github.com/BenThelen/python-refprop
# refprop dll must be installed from NIST package, license requered
# optional method to meos tools calculations and to stream
###############################################################################


# TODO: Don't work when it's used in qt loop, as library work great


from PyQt5.QtWidgets import QApplication
from scipy.constants import R

try:
    # import multiRP
    import refprop
except:
    pass

from lib import unidades
from lib.thermo import ThermoAdvanced
from lib.config import Preferences


__all__ = {212: "helium",
           107: "neon",
           98: "argon",
           1: "hydrogen",
           46: "nitrogen",
           47: "oxygen",
           208: "fluorine",
           62: "water",
           49: "co2",
           48: "co",
           110: "n2o",
           51: "so2",
           219: "cos",
           63: "ammonia",
           50: "h2s",
           2: "methane",
           3: "ethane",
           4: "propane",
           6: "butane",
           5: "isobutan",
           8: "pentane",
           9: "neopentn",
           7: "ipentane",
           10: "hexane",
           52: "ihexane",
           11: "heptane",
           12: "octane",
           13: "nonane",
           14: "decane",
           16: "c12",
           258: "cyclopro",
           38: "cyclohex",
           40: "benzene",
           41: "toluene",
           22: "ethylene",
           23: "propylen",
           24: "1butene",
           27: "ibutene",
           25: "c2butene",
           26: "t2butene",
           66: "propyne",
           117: "methanol",
           134: "ethanol",
           140: "acetone",
           133: "dme",
           951: "nf3",
           971: "krypton",
           994: "xenon",
           953: "sf6",
           645: "cf3i",
           217: "r11",
           216: "r12",
           215: "r13",
           218: "r14",
           642: "r21",
           220: "r22",
           643: "r23",
           645: "r32",
           225: "r41",
           232: "r113",
           231: "r114",
           229: "r115",
           236: "r116",
           1631: "r123",
           1629: "r124",
           1231: "r125",
           1235: "r134a",
           1633: "r141b",
           241: "r142b",
           243: "r143a",
           245: "r152a",
           671: "r218",
           1872: "r227ea",
           1873: "r236fa",
           1817: "r245fa",
           692: "rc318",
           475: "air"}

noIds = ["d2", "parahyd", "d2o", "r365mfc", "r404a", "r410a", "r407c", "r507a"]


class RefProp(ThermoAdvanced):
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
    """
    kwargs = {"ids": [],
              "fraccionMolar": [],

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "rho": 0.0,
              "H": 0.0,
              "S": 0.0,
              "U": 0.0}

    @property
    def calculable(self):
        """Check in the class is fully defined"""
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
        if self.kwargs["x"] != RefProp.kwargs["x"]:
            self.kwargs["Q"] = self.kwargs["x"]
        if self.kwargs["rho"] != RefProp.kwargs["rho"]:
            self.kwargs["D"] = self.kwargs["rho"]
        if self.kwargs["U"] != RefProp.kwargs["U"]:
            self.kwargs["E"] = self.kwargs["U"]

        # # Check thermo definition
        # self._thermo = ""
        # for def_ in ["P-T", "Q-T", "P-Q", "Dmass-T", "Dmass-P", "Hmass-P",
                     # "P-Smass", "Hmass-Smass"]:
            # inputs = def_.split("-")
            # if self.kwargs[inputs[0]] != CoolProp.kwargs[inputs[0]] and \
                    # self.kwargs[inputs[1]] != CoolProp.kwargs[inputs[1]]:
                # self._thermo = def_.replace("-", "")
                # self._par = CP.__getattribute__("%s_INPUTS" % self._thermo)
                # break


        self._thermo = ""
        if self.kwargs["T"] and self.kwargs["P"]:
            self._thermo = "TP"
        elif self.kwargs["T"] and self.kwargs["x"] is not None:
            self._thermo = "TQ"
        elif self.kwargs["P"] and self.kwargs["x"] is not None:
            self._thermo = "PQ"
        elif self.kwargs["T"] and self.kwargs["rho"]:
            self._thermo = "TD"
        elif self.kwargs["P"] and self.kwargs["rho"]:
            self._thermo = "PD"
        elif self.kwargs["P"] and self.kwargs["H"]:
            self._thermo = "PH"
        elif self.kwargs["P"] and self.kwargs["S"]:
            self._thermo = "PS"
        elif self.kwargs["H"] and self.kwargs["S"]:
            self._thermo = "HS"
        elif self.kwargs["T"] and self.kwargs["H"]:
            self._thermo = "TH"
        elif self.kwargs["T"] and self.kwargs["S"]:
            self._thermo = "TS"
        elif self.kwargs["T"] and self.kwargs["U"]:
            self._thermo = "TE"
        elif self.kwargs["P"] and self.kwargs["U"]:
            self._thermo = "PE"
        elif self.kwargs["S"] and self.kwargs["U"]:
            self._thermo = "ES"
        elif self.kwargs["H"] and self.kwargs["rho"]:
            self._thermo = "DH"
        elif self.kwargs["S"] and self.kwargs["rho"]:
            self._thermo = "DS"
        elif self.kwargs["U"] and self.kwargs["rho"]:
            self._thermo = "DE"

        return self._definition and self._thermo

    def args(self):
        var1 = self.kwargs[self._thermo[0]]
        var2 = self.kwargs[self._thermo[1]]

        # unit conversion to refprop accepted units
        # P in kPa, U,H in kJ/kg, S in kJ/kgK
        if self._thermo[0] == "P":
            var1 /= 1000.
        if self._thermo[0] in ("E", "H", "S"):
            var1 /= 1000. * self.M
        if self._thermo[1] == "P":
            var2 /= 1000.
        if self._thermo[1] in ("E", "H", "S"):
            var2 /= 1000. * self.M

        return self._thermo, var1, var2, self.kwargs["fraccionMolar"]

    def _name(self):
        name = []
        for fld in self.kwargs["ids"]:
            if fld in __all__:
                name.append(__all__[fld])
            elif fld in noIds:
                name.append(fld)
        return name

    def calculo(self):
        preos = Preferences.getboolean("refProp", "preos")
        aga = Preferences.getboolean("refProp", "aga")
        gerg = Preferences.getboolean("refProp", "gerg")
        fluido = self._name()

        # refprop.setmod()
        if gerg:
            gerg04(ixflag=1)
        refprop.setup("def", fluido)
        # refprop.setktv()
        if preos:
            refprop.preos(ixflag=2)
        elif aga:
            refprop.setaga()
        # refprop.setref()

        m = refprop.wmol(self.kwargs["fraccionMolar"])["wmix"]
        self.M = unidades.Dimensionless(m)
        crit = refprop.critp(self.kwargs["fraccionMolar"])
        self.Pc = unidades.Pressure(crit["pcrit"], "kPa")
        self.Tc = unidades.Temperature(crit["tcrit"])
        self.rhoc = unidades.Density(crit["Dcrit"]*self.M)

        args = self.args()
        flash = refprop.flsh(*args)
        self.phase, self.x = self.getphase(flash)
        self.T = unidades.Temperature(flash["t"])
        self.P = unidades.Pressure(flash["p"], "kPa")
        self.Tr = unidades.Dimensionless(self.T/self.Tc)
        self.Pr = unidades.Dimensionless(self.P/self.Pc)
        self.rho = unidades.Density(flash["D"]*self.M)
        self.v = unidades.SpecificVolume(1./self.rho)

        self._cp0(flash)

        if flash["nc"] == 1:
            name = refprop.name(flash["nc"])
            self.name = name["hname"]
            self.synonim = name["hn80"]
            self.CAS = name["hcas"]

            info = refprop.info(flash["nc"])
            self.Tt = unidades.Temperature(info["ttrp"])
            self.Tb = unidades.Temperature(info["tnbpt"])
            self.f_accent = unidades.Dimensionless(info["acf"])
            self.momentoDipolar = unidades.DipoleMoment(info["dip"], "Debye")
        else:
            self.name = ""
            self.synonim = ""
            self.CAS = ""

            self.Tt = unidades.Temperature(None)
            self.Tb = unidades.Temperature(None)
            self.f_accent = unidades.Dimensionless(None)
            self.momentoDipolar = unidades.DipoleMoment(None)

        self.Liquido = ThermoAdvanced()
        self.Gas = ThermoAdvanced()
        if self.x == 0:
            # liquid phase
            self.fill(self.Liquido, flash["t"], flash["Dliq"], flash["xliq"])
            self.fill(self, flash["t"], flash["Dliq"], flash["xliq"])
        elif self.x == 1:
            # vapor phase
            self.fill(self.Gas, flash["t"], flash["Dvap"], flash["xvap"])
            self.fill(self, flash["t"], flash["Dvap"], flash["xvap"])
        else:
            # Two phase
            self.fill(self.Liquido, flash["t"], flash["Dliq"], flash["xliq"])
            self.fill(self.Gas, flash["t"], flash["Dvap"], flash["xvap"])

            self.h = unidades.Enthalpy(self.x*self.Gas.h+(1-self.x)*self.Liquido.h)
            self.s = unidades.SpecificHeat(self.x*self.Gas.s+(1-self.x)*self.Liquido.s)
            self.cp = unidades.SpecificHeat(self.x*self.Gas.cp+(1-self.x)*self.Liquido.cp)
            self.cv = unidades.SpecificHeat(self.x*self.Gas.cv+(1-self.x)*self.Liquido.cv)
            self.cp_cv = unidades.Dimensionless(self.cp/self.cv)
            self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv)

        if self.x < 1 and self.T <= self.Tc:
            surten = refprop.surten(flash["t"], flash["Dliq"], flash["Dvap"],
                                    flash["xliq"], flash["xvap"])
            self.Liquido.sigma = unidades.Tension(surten["sigma"])
        else:
            self.Liquido.sigma = unidades.Tension(None)

        if 0 < self.x < 1:
            self.Hvap = unidades.Enthalpy(self.Gas.h-self.Liquido.h)
            self.Svap = unidades.SpecificHeat(self.Gas.s-self.Liquido.s)
        else:
            self.Hvap = unidades.Enthalpy(None)
            self.Svap = unidades.SpecificHeat(None)

#        cp2=refprop.cv2pk(self.nc, flash["t"], flash["D"])
#        self.cv2p=unidades.SpecificHeat(cp2["cv2p"]/self.M)
#        self.csat=unidades.SpecificHeat(cp2["csat"]/self.M)

    def _cp0(self, flash):
        "Set ideal properties to state"""
        p0 = refprop.therm0(flash["t"], flash["D"], flash["x"])["p"]
        self.v0 = unidades.SpecificVolume(R/self.M*self.T/self.P.kPa)
        self.rho0 = unidades.Density(1/self.v0)
        self.P0 = unidades.Pressure(p0, "kPa")
        self.P_Pideal = unidades.Pressure(self.P-self.P0)

        cp0 = refprop.therm0(float(self.T), self.rho0/self.M, flash["x"])
        self.h0 = unidades.Enthalpy(cp0["h"]/self.M, "kJkg")
        self.u0 = unidades.Enthalpy(cp0["e"]/self.M, "kJkg")
        self.s0 = unidades.SpecificHeat(cp0["s"]/self.M, "kJkgK")
        self.a0 = unidades.Enthalpy(cp0["A"]/self.M, "kJkg")
        self.g0 = unidades.Enthalpy(cp0["G"]/self.M, "kJkg")
        self.cp0 = unidades.SpecificHeat(cp0["cp"]/self.M, "kJkgK")
        self.cv0 = unidades.SpecificHeat(cp0["cv"]/self.M, "kJkgK")
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv0)
        self.w0 = unidades.Speed(cp0["w"])
        # self.gamma0 = unidades.Dimensionless(cp0.gamma)

    def fill(self, fase, T, rho, x):
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
#         if self._multicomponent:
            # fase.fi = []
            # fase.f = []
            # for i in range(len(self.kwargs["ids"])):
                # fase.fi.append(unidades.Dimensionless(
                    # estado.fugacity_coefficient(i)))
                # fase.f.append(unidades.Pressure(estado.fugacity(i)))
        # else:
            # fase.fi = unidades.Dimensionless([1])
            # fase.f = unidades.Pressure([self.P])

        fase.cv = unidades.SpecificHeat(thermo["cv"]/fase.M, "JgK")
        fase.cp = unidades.SpecificHeat(thermo["cp"]/fase.M, "JgK")
        fase.cp_cv = unidades.Dimensionless(fase.cp/fase.cv)
        fase.w = unidades.Speed(thermo["w"])

        fase.rhoM = unidades.MolarDensity(fase.rho*self.M)
        fase.hM = unidades.MolarEnthalpy(fase.h*self.M)
        fase.sM = unidades.MolarSpecificHeat(fase.s*self.M)
        fase.uM = unidades.MolarEnthalpy(fase.u*self.M)
        fase.aM = unidades.MolarEnthalpy(fase.a*self.M)
        fase.gM = unidades.MolarEnthalpy(fase.g*self.M)
        fase.cvM = unidades.MolarSpecificHeat(fase.cv*self.M)
        fase.cpM = unidades.MolarSpecificHeat(fase.cp*self.M)

        '''Compute miscellaneous thermodynamic properties
        betas--Adiabatic compressibility
        thrott--Isothermal throttling coefficient
        '''

        fase.joule = unidades.TemperaturePressure(thermo["hjt"], "KkPa")
        # fase.Gruneisen = unidades.Dimensionless(
            # estado.first_partial_deriv(CP.iP, CP.iT, CP.iDmass))
        fase.alfav = unidades.InvTemperature(thermo["beta"])
        fase.kappa = unidades.InvPressure(thermo["xkappa"], "kPa")
        # fase.alfap = unidades.Density(fase.alfav/self.P/fase.kappa)
        # fase.betap = unidades.Density(
            # fase.rho/self.P*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iT))
        # fase.betas = unidades.TemperaturePressure(
            # estado.first_partial_deriv(CP.iT, CP.iP, CP.iSmass))
        # fase.gamma = unidades.Dimensionless(
            # fase.rho/self.P*estado.first_partial_deriv(CP.iP, CP.iDmass, CP.iSmass))

        fase.kt = unidades.Dimensionless(thermo3["xkt"])
        fase.Kt = unidades.Pressure(thermo3["xkkt"], "kPa")
        fase.Ks = unidades.Pressure(thermo3["bs"], "kPa")
        fase.ks = unidades.InvPressure(thermo3["xisenk"], "kPa")

        dh = refprop.dhd1(T, rho, x)
        fase.dhdT_rho = unidades.SpecificHeat(dh["dhdt_D"]/fase.M, "kJkgK")
        fase.dhdT_P = unidades.SpecificHeat(dh["dhdt_p"]/fase.M, "kJkgK")
        fase.dhdP_T = unidades.EnthalpyPressure(dh["dhdp_t"]/fase.M, "kJkgkPa")
        # dhdtp_D : fix in library
        fase.dhdP_rho = unidades.EnthalpyPressure(dh["dhdtp_D"]/fase.M, "kJkgkPa")
        fase.dhdrho_T = unidades.EnthalpyDensity(dh["dhdD_t"]/fase.M**2, "kJkgkgm3")
        fase.dhdrho_P = unidades.EnthalpyDensity(dh["dhdD_p"]/fase.M**2, "kJkgkgm3")

        fase.dpdT_rho = unidades.PressureTemperature(thermo["dpdt"], "kPaK")
        fase.dpdrho_T = unidades.PressureDensity(thermo["dpdD"]/fase.M, "kPakgm3")
        # TODO: Add unit for derivative d^2p/dD^2 [kPa-L^2/mol^2]
        fase.d2pdrho2 = unidades.Dimensionless(thermo["d2pdD2"]/fase.M**2/1000)  # MPa·m⁶/kg²
        fase.drhodP_T = unidades.DensityPressure(thermo["dDdp"]*fase.M, "kgm3kPa")
        fase.drhodT_P = unidades.DensityTemperature(thermo["dDdt"]*fase.M)


        # fase.Z_rho = unidades.SpecificVolume((fase.Z-1)/fase.rho)
        fase.IntP = unidades.Pressure(thermo3["pint"], "kPa")
        fase.hInput = unidades.Enthalpy(thermo3["spht"]/fase.M)
        fase.invT = unidades.InvTemperature(-1/self.T)

        fpv = refprop.fpv(T, rho, self.P.kPa, x)["Fpv"]
        # supercompressibility factor, Fpv.
        fase.fpv = unidades.Dimensionless(fpv)

        chempot = refprop.chempot(T, rho, x)["u"]
        fase.chempot = [unidades.Enthalpy(c/fase.M) for c in chempot]
        fi = refprop.fugcof(T, rho, x)["f"]
        fase.fi = [unidades.Dimensionless(f) for f in fi]
        f = refprop.fgcty(T, rho, x)["f"]
        fase.f = [unidades.Pressure(f_i, "kPa") for f_i in f]

        b = refprop.virb(T, x)["b"]
        fase.virialb = unidades.SpecificVolume(b/self.M)
        c = refprop.virc(T, x)["c"]
        fase.virialc = unidades.SpecificVolume_square(c/self.M**2)
        # viriald don't supported for windows
        d = refprop.vird(T, x)["d"]
        fase.viriald = unidades.Dimensionless(d/self.M**3)
        ba = refprop.virba(T, x)["ba"]
        fase.virialba = unidades.SpecificVolume(ba/self.M)
        ca = refprop.virca(T, x)["ca"]
        fase.virialca = unidades.SpecificVolume_square(ca/self.M**2)
        dcdt = refprop.dcdt(T, x)["dct"]
        fase.dcdt = unidades.Dimensionless(dcdt/self.M**2)
        dcdt2 = refprop.dcdt2(T, x)["dct2"]
        fase.dcdt2 = unidades.Dimensionless(dcdt2/self.M**2)
        dbdt = refprop.dbdt(T, x)["dbt"]
        fase.dbdt = unidades.Dimensionless(dbdt/self.M)

        fase.fraccion = [unidades.Dimensionless(xi) for xi in x]
        xg = refprop.xmass(x)["xkg"]
        fase.fraccion_masica = [unidades.Dimensionless(xi) for xi in xg]

        transport = refprop.trnprp(T, rho, x)
        fase.mu = unidades.Viscosity(transport["eta"], "muPas")
        fase.nu = unidades.Diffusivity(fase.mu/fase.rho)
        fase.k = unidades.ThermalConductivity(transport["tcx"])
        fase.alfa = unidades.Diffusivity(fase.k/1000/fase.rho/fase.cp)
        fase.Prandt = unidades.Dimensionless(fase.mu*fase.cp/fase.k)

        dielec = refprop.dielec(T, rho, x)
        fase.epsilon = unidades.Dimensionless(dielec["de"])

    def getphase(self, fld):
        """Return fluid phase
        Override refprop original function with translation support"""
        # check if fld above critical pressure
        if fld['p'] > self.Pc.kPa:
            # check if fld above critical pressure
            if fld['t'] > self.Tc:
                return QApplication.translate("pychemqt", "Supercritical fluid"), 1.
            else:
                return QApplication.translate("pychemqt", "Compressible liquid"), 1.
        # check if fld above critical pressure
        elif fld['t'] > self.Tc:
            return QApplication.translate("pychemqt", "Gas"), 1.
        # check if ['q'] in fld
        if 'q' not in list(fld.keys()):
            if 'h' in list(fld.keys()):
                fld['q'] = refprop.flsh('ph', fld['p'], fld['h'], fld['x'])['q']
            elif 's' in list(fld.keys()):
                fld['q'] = refprop.flsh('ps', fld['p'], fld['s'], fld['x'])['q']
        # check q
        if fld['q'] > 1:
            return QApplication.translate("pychemqt", "Vapor"), 1.
        elif fld['q'] == 1:
            return QApplication.translate("pychemqt", "Saturated vapor"), 1.
        elif 0 < fld['q'] < 1:
            return QApplication.translate("pychemqt", "Two phases"), fld['q']
        elif fld['q'] == 0:
            return QApplication.translate("pychemqt", "Saturated liquid"), 0.
        elif fld['q'] < 0:
            return QApplication.translate("pychemqt", "Liquid"), 0.


if __name__ == '__main__':
    # import sys
    # from PyQt5 import QtWidgets
    # app = QtWidgets.QApplication(sys.argv)

    # fluido = RefProp(ids=[62], T=300, P=1e6)
    # print("%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g"
    #       % (fluido.rho, fluido.u.kJkg, fluido.h.kJkg, fluido.s.kJkgK,
    #          fluido.cv.kJkgK, fluido.cp.kJkgK, fluido.cp0.kJkgK,
    #          fluido.cp_cv, fluido.w, fluido.Z, fluido.a, fluido.g))
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


#    fluido=RefProp(fluido=[u'hydrogen', u'nitrogen', u'oxygen', u'water'], fraccionMolar=[0.03, 0.95, 0.01, 0.01], T=300, P=1e5)
#    print(fluido.rho, fluido.cp.kJkgK, fluido.cv.kJkgK)

    from lib.corriente import Corriente
    kw = {"MEoS": True, "refprop": True, "ids": [5, 6, 7, 8, 10]}
    entrada = Corriente(T=300, P=1e6, caudalMasico=0.01,
                        fraccionMolar=[.3, 0.25, 0.05, 0.15, 0.25], **kw)

    # fluido = RefProp(ids=[5, 6, 7, 8, 10], fraccionMolar=[0.3, 0.25, 0.05, 0.15, 0.25], T=300, P=1e6)
    # import pprint
    # pprint.pprint(fluido.__dict__)
