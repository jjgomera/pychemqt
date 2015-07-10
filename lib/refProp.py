#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to multiparameter equation of state calculation using refprop
# https://github.com/BenThelen/python-refprop
# refprop dll must be installed from NIST package, license requered
# optional method to meos tools calculations and to stream
###############################################################################

#TODO: Don't work when it's used in qt loop, as library work great

from PyQt5.QtWidgets import QApplication

try:
    import refprop
except:
    pass

from lib import unidades
from .config import Fluid


class RefProp(object):
    """
    Stream class using refProp external library
    Parameters needed to define it are:

        -ref: reference state
        -fluido: index of fluid
        -fraccionMolar: molar fraction

        -T: Temperature, Kelvin
        -P: Pressure, Pa
        -rho: Density, kg/m3
        -H: Enthalpy, J/kg
        -S: Entropy, J/kgK
        -U: Internal energy, J/kg
        -x: Quality, -
    """
    kwargs = {"ref": "def",
              "fluido": None,
              "fraccionMolar": None,

              "T": 0.0,
              "P": 0.0,
              "x": None,
              "rho": 0.0,
              "H": 0.0,
              "S": 0.0,
              "U": 0.0}

    status = 0
    msg = "Unknown variables"

    def __init__(self, **kwargs):
        self.kwargs = RefProp.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        if len(kwargs["fluido"]) == 1:
            kwargs["fraccionMolar"] = [1.]
        self.kwargs.update(kwargs)

        if self.calculable:
            self.status = 1
            self.calculo()

    @property
    def calculable(self):
        if self.kwargs["ref"] and self.kwargs["fluido"] \
           and self.kwargs["fraccionMolar"]:
            self._definition = True
        else:
            self._definition = False

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
        # Correct refprop custom namespace versus pychemqt namespace
        if "Q" in self._thermo:
            self.kwargs["Q"] = self.kwargs["x"]
        if "D" in self._thermo:
            self.kwargs["D"] = self.kwargs["rho"]
        if "E" in self._thermo:
            self.kwargs["E"] = self.kwargs["U"]

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

    def calculo(self):
        refprop.setup(self.kwargs["ref"], self.kwargs["fluido"])

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
        self.rho = unidades.Density(flash["D"]*self.M)
        self.v = unidades.SpecificVolume(1./self.rho)
        name = refprop.name(flash["nc"])
        info = refprop.info(flash["nc"])
        if flash["nc"] == 1:
            self.name = name["hname"]
            self.synonim = name["hn80"]
            self.CAS = name["hcas"]

            self.Tt = unidades.Temperature(info["ttrp"])
            self.Tb = unidades.Temperature(info["tnbpt"])
            self.f_accent = unidades.Dimensionless(info["acf"])
            self.momentoDipolar = unidades.DipoleMoment(info["dip"], "Debye")

        self.Liquido = Fluid()
        self.Vapor = Fluid()
        if self.x < 1.:  # Hay fase liquida
            liquido_thermo = refprop.therm2(flash["t"], flash["Dliq"],
                                            flash["xliq"])
            liquido_mol = refprop.wmol(flash["xliq"])
            try:
                liquido_transport = refprop.trnprp(flash["t"], flash["Dliq"],
                                                   flash["xliq"])
            except refprop.RefpropError as e:
                print(e)
                liquido_transport = None
            liquido_dielec = refprop.dielec(flash["t"], flash["Dliq"],
                                            flash["xliq"])
            liquido_thermo0 = refprop.therm0(flash["t"], flash["Dliq"],
                                             flash["xliq"])
            self.fill(self.Liquido, flash, liquido_thermo, liquido_mol,
                      liquido_transport, liquido_dielec, liquido_thermo0)

        if self.x > 0.:  # Hay fase vapor
            vapor_thermo = refprop.therm2(flash["t"], flash["Dvap"],
                                          flash["xvap"])
            vapor_mol = refprop.wmol(flash["xvap"])
            try:
                vapor_transport = refprop.trnprp(flash["t"], flash["Dvap"],
                                                 flash["xvap"])
            except refprop.RefpropError as e:
                print(e)
                vapor_transport = None
            vapor_dielec = refprop.dielec(flash["t"], flash["Dvap"],
                                          flash["xvap"])
            vapor_thermo0 = refprop.therm0(flash["t"], flash["Dvap"],
                                           flash["xvap"])
            self.fill(self.Vapor, flash, vapor_thermo, vapor_mol,
                      vapor_transport, vapor_dielec, vapor_thermo0)

#        crit = multiRP.mRP[u'process'](target=multiRP.critp, args=(self.kwargs["fraccionMolar"], setup, multiRP.mRP))
#        mol = multiRP.mRP[u'process'](target=multiRP.wmol, args=(self.kwargs["fraccionMolar"], setup, multiRP.mRP))
#        processlist = [crit, mol]
#        multiRP.run_mRP(processlist)
#        self.Pc=unidades.Pressure(multiRP.mRP[u'result'][processlist[0].name]["pcrit"], "kPa")
#        self.Tc=unidades.Temperature(multiRP.mRP[u'result'][processlist[0].name]["tcrit"])
#        self.rhoc=unidades.Density(multiRP.mRP[u'result'][processlist[0].name]["Dcrit"]*self.M)

#        args=self.args()
#        flash = multiRP.mRP[u'process'](target=multiRP.flsh, args=args, kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#        name = multiRP.mRP[u'process'](target=multiRP.name, args=(1, setup, multiRP.mRP))
#        info = multiRP.mRP[u'process'](target=multiRP.info, args=(1, setup, multiRP.mRP))
#        processlist = [flash, name, info]
#        multiRP.run_mRP(processlist)
#        flash=multiRP.mRP[u'result'][processlist[0].name]
#        self.phase, self.x=self.getphase(flash)
#        self.nc=flash["nc"]
#        self.fraccion=self.kwargs["fraccionMolar"]
#        if flash["nc"] ==1:
#            self.name=multiRP.mRP[u'result'][processlist[1].name]["hname"]
#            self.synonim=multiRP.mRP[u'result'][processlist[1].name]["hn80"]
#            self.CAS=multiRP.mRP[u'result'][processlist[1].name]["hcas"]
#
#            self.Tt=unidades.Temperature(multiRP.mRP[u'result'][processlist[2].name]["ttrp"])
#            self.Tb=unidades.Temperature(multiRP.mRP[u'result'][processlist[2].name]["tnbpt"])
#            self.f_accent=unidades.Dimensionless(multiRP.mRP[u'result'][processlist[2].name]["acf"])
#            self.momentoDipolar=unidades.DipoleMoment(multiRP.mRP[u'result'][processlist[2].name]["dip"], "Debye")
#            self.Rgas=unidades.SpecificHeat(multiRP.mRP[u'result'][processlist[2].name]["Rgas"]/self.M)
#        else:
#            self.Rgas=unidades.SpecificHeat(refprop.rmix2(self.kwargs["fraccionMolar"])["Rgas"]/self.M)

#        self.T=unidades.Temperature(flash["t"])
#        self.P=unidades.Pressure(flash["p"], "kPa")
#        self.rho=unidades.Density(flash["D"]*self.M)
#        self.v=unidades.SpecificVolume(1./self.rho)

#        self.Liquido=Fluid()
#        self.Vapor=Fluid()
#        if self.x<1.: #Hay fase liquida
#            liquido_thermo = multiRP.mRP[u'process'](target=multiRP.therm2, args=(flash["t"], flash["Dliq"], flash["xliq"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            liquido_mol = multiRP.mRP[u'process'](target=multiRP.wmol, args=(flash["xliq"], setup, multiRP.mRP))
#            liquido_transport = multiRP.mRP[u'process'](target=multiRP.trnprp, args=(flash["t"], flash["Dliq"], flash["xliq"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            liquido_dielec = multiRP.mRP[u'process'](target=multiRP.dielec, args=(flash["t"], flash["Dliq"], flash["xliq"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            liquido_thermo0 = multiRP.mRP[u'process'](target=multiRP.therm0, args=(flash["t"], flash["Dliq"], flash["xliq"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            processlist = [liquido_thermo, liquido_mol, liquido_transport, liquido_dielec, liquido_thermo0]
#            try:
#                multiRP.run_mRP(processlist)
#                transport=multiRP.mRP[u'result'][processlist[2].name]
#            except refprop.RefpropError as e:
#                print e
#                transport=None
#            self.fill(self.Liquido, flash, multiRP.mRP[u'result'][processlist[0].name], multiRP.mRP[u'result'][processlist[1].name],
#                transport, multiRP.mRP[u'result'][processlist[3].name], multiRP.mRP[u'result'][processlist[4].name])

#        if self.x>0.: #Hay fase vapor
#            vapor_thermo = multiRP.mRP[u'process'](target=multiRP.therm2, args=(flash["t"], flash["Dvap"], flash["xvap"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            vapor_mol = multiRP.mRP[u'process'](target=multiRP.wmol, args=(flash["xvap"], setup, multiRP.mRP))
#            vapor_transport = multiRP.mRP[u'process'](target=multiRP.trnprp, args=(flash["t"], flash["Dvap"], flash["xvap"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            vapor_dielec = multiRP.mRP[u'process'](target=multiRP.dielec, args=(flash["t"], flash["Dvap"], flash["xvap"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            vapor_thermo0 = multiRP.mRP[u'process'](target=multiRP.therm0, args=(flash["t"], flash["Dvap"], flash["xvap"]), kwargs={u'prop': setup, u'mRP': multiRP.mRP})
#            processlist = [vapor_thermo, vapor_mol, vapor_transport, vapor_dielec, vapor_thermo0]
#            try:
#                multiRP.run_mRP(processlist)
#                transport=multiRP.mRP[u'result'][processlist[2].name]
#            except refprop.RefpropError as e:
#                print e
#                transport=None
#
#            self.fill(self.Vapor, flash, multiRP.mRP[u'result'][processlist[0].name], multiRP.mRP[u'result'][processlist[1].name],
#                transport, multiRP.mRP[u'result'][processlist[3].name], multiRP.mRP[u'result'][processlist[4].name])

        self.h = unidades.Enthalpy(self.x*self.Vapor.h+(1-self.x)*self.Liquido.h)
        self.s = unidades.SpecificHeat(self.x*self.Vapor.s+(1-self.x)*self.Liquido.s)
        self.cp = unidades.SpecificHeat(self.x*self.Vapor.cp+(1-self.x)*self.Liquido.cp)
        self.cp0 = unidades.SpecificHeat(self.x*self.Vapor.cp0+(1-self.x)*self.Liquido.cp0)
        self.cv = unidades.SpecificHeat(self.x*self.Vapor.cv+(1-self.x)*self.Liquido.cv)

        self.cp_cv = unidades.Dimensionless(self.cp/self.cv)
        self.cp0_cv = unidades.Dimensionless(self.cp0/self.cv)

        if self.T <= self.Tc:
            surten = refprop.surten(flash["t"], flash["Dliq"], flash["Dvap"],
                                    flash["xliq"], flash["xvap"])
            self.surten = unidades.Tension(surten["sigma"])
        else:
            self.surten = unidades.Tension(None)

#        cp2=refprop.cv2pk(self.nc, flash["t"], flash["D"])
#        self.cv2p=unidades.SpecificHeat(cp2["cv2p"]/self.M)
#        self.csat=unidades.SpecificHeat(cp2["csat"]/self.M)

    def fill(self, fase, flash, thermo, mol, transport, dielec, thermo0):
        fase.update(Fluid(thermo))
        fase.fraccion = flash["xliq"]
        fase.M = unidades.Dimensionless(mol["wmix"])
        fase.rho = unidades.Density(flash["Dliq"]*fase.M)

        fase.u = unidades.Enthalpy(fase["e"]/fase.M, "Jg")
        fase.cv = unidades.SpecificHeat(fase["cv"]/fase.M, "JgK")
        fase.cp = unidades.SpecificHeat(fase["cp"]/fase.M, "JgK")
        fase.h = unidades.Enthalpy(fase["h"]/fase.M, "Jg")
        fase.s = unidades.SpecificHeat(fase["s"]/fase.M, "JgK")
        fase.w = unidades.Speed(fase["w"])
        fase.joule = unidades.TemperaturePressure(fase["hjt"], "KkPa")
        fase.Z = unidades.Dimensionless(fase["Z"])
        fase.A = unidades.Enthalpy(fase["A"]/fase.M, "Jg")
        fase.G = unidades.Enthalpy(fase["G"]/fase.M, "Jg")
        fase.xkappa = unidades.InvPressure(fase["xkappa"], "kPa")
        fase.alfav = unidades.InvTemperature(fase["beta"])
#            fase.dpdD = fase["dpdD"]      #derivative dP/dD [kPa-L/mol]
#            fase.d2pdD2 = fase["d2pdD2"]  #derivative d^2p/dD^2 [kPa-L^2/mol^2]
#            fase.dpdt = unidades.PressureTemperature(fase["dpdt"], "kPaK")
#            fase.dDdt = fase["dDdt"]      #derivative dD/dt [mol/(L-K)]
#            fase.dDdp = fase["dDdp"]      #derivative dD/dp [mol/(L-kPa)]
#
#            fluido2=refprop.therm3(flash["t"], flash["Dliq"], flash["xliq"])
#            fase.xisenk = fluido2["xisenk"]
#            fase.xkt = fluido2["xkt"]
#            fase.betas = fluido2["betas"]
#            fase.bs = fluido2["bs"]
#            fase.xkkt = fluido2["xkkt"]
#            fase.thrott = fluido2["thrott"]
#            fase.pint = fluido2["pint"]
#            fase.spht = fluido2["spht"]

#            fase.fpv = refprop.fpv(flash["t"], flash["Dliq"], flash["p"], flash["xliq"])
#            fase.chempot = refprop.chempot(flash["t"], flash["Dliq"], flash["xliq"])
#            fase.fgcty = refprop.fgcty(flash["t"], flash["Dliq"], flash["xliq"])
#            fase.fugcof = refprop.fugcof(flash["t"], flash["Dliq"], flash["xliq"])

#            fase.virb = refprop.virb(flash["t"], flash["xliq"])["b"]
#            fase.virc = refprop.virc(flash["t"], flash["xliq"])["c"]
#            fase.vird = refprop.vird(flash["t"], flash["xliq"])["d"]
#            fase.virba = refprop.virba(flash["t"], flash["xliq"])["ba"]
#            fase.virca = refprop.virca(flash["t"], flash["xliq"])["ca"]

        if transport:
            fase.mu = unidades.Viscosity(transport["eta"], "muPas")
            fase.k = unidades.ThermalConductivity(transport["tcx"])
            fase.Prandt = unidades.Dimensionless(fase.mu*fase.cp/fase.k)
        else:
            fase.mu = unidades.Viscosity(None)
            fase.k = unidades.ThermalConductivity(None)
            fase.Prandt = unidades.Dimensionless(None)

        fase.dielec = unidades.Dimensionless(dielec["de"])
        fase.cp0 = unidades.SpecificHeat(thermo0["cp"]/fase.M)
        fase.cp0_cv = unidades.Dimensionless(fase.cp0/fase.cv)

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

noId = ["d2", "parahyd", "d2o", "r365mfc", "r404a", "r410a", "r407c", "r507a"]

if __name__ == '__main__':
    import sys
    from PyQt4 import QtGui
#    app = QtGui.QApplication(sys.argv)
    from ctypes import CDLL, RTLD_GLOBAL
    _rp = CDLL("/usr/local/lib/librefprop.so", mode=RTLD_GLOBAL)
#    fluido=RefProp(fluido=[u"water"], T=303.15, P=101325.)
#    print fluido.Liquido.rho, fluido.Liquido.h, fluido.Liquido.s, fluido.Liquido.cv, fluido.Liquido.cp, fluido.Liquido.Z, fluido.Liquido.G

    prop = refprop.setup('def', 'air')
#    print prop
#    prop = refprop.wmol(prop[u'x'])
#    print prop

#    fluido=RefProp(fluido=[u'hydrogen', u'nitrogen', u'oxygen', u'water'], fraccionMolar=[0.03, 0.95, 0.01, 0.01], T=300, P=1e5)
#    print fluido.rho, fluido.cp.kJkgK, fluido.cv.kJkgK
