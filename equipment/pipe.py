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
#   Library for pipe equipment definition
###############################################################################


import os

from tools.qt import tr
from scipy.constants import g, pi

from lib import unidades
from lib.friction import f_friccion
from lib.adimensional import Re
from equipment.parents import equipment
from equipment.heatExchanger import Heat_Exchanger


class Pipe(equipment):
    """Class to model a pipe segment

    Parameters:
        entrada: Corriente instance to define the input stream to equipment
        metodo: Calculate method
            0   -   Single phase flow
            1   -   Water flow (Hazen-Williams equation)
            2   -   Steam flow (Fritzsche equation)
            3   -   Isothermic gas flow
            4   -   Multiphase flow (Método Baker)
            5   -   Multiphase flow (Método de Beggs and Brill)
        thermal: Thermal mode:
            0   -   Adiabatic
            1   -   Fix heat flow
            2   -   Heat exchanger calculation
        Material: Array with material properties
            0   -   Name
            1   -   Class
            2   -   Roughness, mm
            3   -   Catalog
            4   -   Di, mm
            5   -   Width, mm
            6   -   De, mm
            7   -   Weight, kg/m
            8   -   Volume, m³/100m
            9   -   Surface, m²/100m
            10 -    Index of material in catalog
            11 -    Index of dimenion in catalog
        Accesorios: Array with each fitting in pipe:
            Index with fitting type
            Index diameter
            K
            Count
            Type
            Diameter in mm
            Nominal diameter, in inch
            Text

        l: Length of pipe
        h: Head increase between input and output pipe
        K: Fittings total coefficient
        C: Parameter for Hazen-Williams equation
        T_ext: External temperatura of pipe for heat exchanger calculation
        U: Global heat transfer coeficient for pipe wall
        Q: Heat transfered by pipe wall

    Coste:
        Only available for steal pipes, Ref Darby pag 217

    >>> from lib.corriente import Corriente
    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    >>> material=["Cast Iron", "Class A", 0.12192, '6"', 152.4, 11.43, \
                  175.26, 42.07, 1.824, 55.06]
    >>> pipe=Pipe(entrada=agua, metodo=0, l=5, material=material)
    >>> print("%0.4f %6g %6g %6g" % (pipe.Di, pipe.V, pipe.Re, pipe.DeltaP))
    0.1524 0.0626815 9427.99 1.83620
    """
    title = tr("pychemqt", "Pipe")
    help = ""
    kwargs = {
        "entrada": None,
        "metodo": 0,
        "thermal": 0,
        "material": [],
        "accesorios": [],
        "l": 0.0,
        "h": 0.0,
        "K": 0.0,
        "C": 100.,
        "T_ext": 0.0,
        "U": 0.0,
        "Q": 0.0,

        "f_install": 2.8,
        "Base_index": 0.0,
        "Current_index": 0.0}
    kwargsInput = ("entrada", )
    kwargsValue = ("l", "h", "C")
    kwargsList = ("metodo", "thermal")
    calculateValue = ("DeltaP", "DeltaP_f", "DeltaP_ac", "DeltaP_h",
                      "DeltaP_v", "DeltaP_100ft", "V", "f", "Re", "Tout")
    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 5
    salida = [None]

    TEXT_METODO = (
        tr("pychemqt", "Single Phase flow"),
        tr("pychemqt", "Water (Hazen-Williams)"),
        tr("pychemqt", "Steam (Fritzsche)"),
        tr("pychemqt", "Isotermic gas flow"),
        tr("pychemqt", "Two Phase flow (Baker method)"),
        tr("pychemqt",
                               "Two Phase flow (Beggs and Brill method)"))
    TEXT_THERMAL = (
        tr("pychemqt", "Adiabatic"),
        tr("pychemqt", "Heat flux"),
        tr("pychemqt", "Heat transfer"))

    @property
    def isCalculable(self):
        materialCoste = ['Stainless Steel (ANSI)', 'Steel Galvanised (ANSI)',
                         'Steel (ANSI)']
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"] and self.kwargs["material"] and \
                self.kwargs["material"][0] in materialCoste:
            self.statusCoste = True
        else:
            self.statusCoste = False
            self.C_adq = unidades.Currency(None)
            self.C_inst = unidades.Currency(None)

        if not self.kwargs["entrada"]:
            self.msg = tr("pychemqt", "undefined input")
            self.status = 0
            return
        if not self.kwargs["l"]:
            self.msg = tr("pychemqt",
                                              "undefined pipe length")
            self.status = 0
            return
        if not self.kwargs["material"]:
            self.msg = tr("pychemqt", "undefined material")
            self.status = 0
            return

        if self.kwargs["thermal"] == 1 and not self.kwargs["Q"]:
            self.msg = tr("pychemqt",
                                              "undefined heat flux")
            self.status = 0
            return
        elif self.kwargs["thermal"] == 2 and \
                (not self.kwargs["T_ext"] or not self.kwargs["U"]):
            self.msg = tr(
                "pychemqt", "undefined heat transfer conditions")
            self.status = 0
            return

        if self.kwargs["metodo"] == 1 and not self.kwargs["C"]:
            self.msg = tr(
                "pychemqt", "undefined C William Factor")
            self.status = 0
            return

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        self.L = unidades.Length(self.kwargs["l"])

        if self.kwargs["entrada"].x == 0:
            self.rho = self.kwargs["entrada"].Liquido.rho
            self.mu = self.kwargs["entrada"].Liquido.mu
        else:
            self.rho = self.kwargs["entrada"].Gas.rho
            self.mu = self.kwargs["entrada"].Gas.mu

        self.material = self.kwargs["material"][0] + " " + \
            self.kwargs["material"][1]
        self.Dn = self.kwargs["material"][3]
        self.rugosidad = unidades.Length(self.kwargs["material"][2], "mm")
        self.De = unidades.Length(self.kwargs["material"][6], "mm")
        self.w = unidades.Length(self.kwargs["material"][5], "mm")
        self.Di = unidades.Length((self.De-2*self.w))
        self.eD = unidades.Dimensionless(self.rugosidad/self.Di)
        self.seccion = unidades.Area(pi/4*self.Di**2)
        self.A = unidades.Area(pi*self.De*self.L)
        self.V = unidades.Speed(self.kwargs["entrada"].Q/self.seccion)
        self.Re = Re(self.Di, self.V, self.rho, self.mu)
        K = 0
        for accesorio in self.kwargs["accesorios"]:
            K += accesorio[2]*accesorio[3]
        self.K = unidades.Dimensionless(K)
        self.DeltaP_h = unidades.Pressure(g*self.kwargs["h"]*self.rho)
        self.DeltaP_ac = unidades.Pressure(self.K*self.V**2/2*self.rho)

        self.f = f_friccion(self.Re, self.eD)
        self.DeltaP_f = self.__DeltaP_friccion()
        # TODO:
        self.DeltaP_v = unidades.Pressure(0)

        self.DeltaP = unidades.Pressure(self.DeltaP_h + self.DeltaP_ac +
                                        self.DeltaP_f + self.DeltaP_v)
        self.DeltaP_100ft = self.DeltaP*100/self.L.ft
        self.Pout = unidades.Pressure(self.kwargs["entrada"].P-self.DeltaP)

        if self.kwargs["thermal"] == 0:
            self.Tout = self.kwargs["entrada"].T
            self.Heat = unidades.Power(0)
        else:
            ch = Heat_Exchanger()
            ch(entrada=self.kwargs["entrada"], modo=self.kwargs["thermal"],
               Heat=self.kwargs["Q"], deltaP=self.DeltaP, A=self.A,
               U=self.kwargs["U"], Text=self.kwargs["Text"])
            self.Tout = ch.salida[0].T
            self.Heat = ch.Heat

        self.salida = [self.kwargs["entrada"].clone(T=self.Tout, P=self.Pout)]
        self.Pin = self.kwargs["entrada"].P
        self.Pout = self.salida[0].P

    def __DeltaP_friccion(self):
        """Método para el calculo de la perdida de presión"""
        if self.kwargs["metodo"] == 0:
            dp = unidades.DeltaP(self.L*self.V**2/self.Di*self.f*self.rho/2)
        elif self.kwargs["metodo"] == 1:
            p = (self.kwargs["entrada"].Q.galUSmin*self.L.ft**0.54/0.442 /
                 self.Di.inch**2.63/self.kwargs["C"])**(1./0.54)
            dp = unidades.DeltaP(p, "psi")
        elif self.kwargs["metodo"] == 2:
            q = self.kwargs["entrada"].caudalmasico.lbh
            p = 2.1082*self.L.ft*q**1.85/self.rho.lbft3/1e7/self.Di.inch**4.97
            dp = unidades.Pressure(p, "psi")
        elif self.kwargs["metodo"] == 3:
            pass

        elif self.kwargs["metodo"] == 4:
            pass

        elif self.kwargs["metodo"] == 5:
            pass

        return dp

    def coste(self):
        """
        Coste solo disponible para las tuberías de acero
        Ref Darby pag 217

        kwargs:
            schedule: Clase de acero
        """
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]
        codigo = str(self.kwargs["material"][1])
        if codigo in ('Sch.  40', 'Sch.  5S'):
            a = 30.
            p = 1.31
        elif codigo in ('Sch.  80', 'Sch.  10S'):
            a = 38.1
            p = 1.35
        elif codigo in ('Sch.  160', 'Sch.  40S'):
            a = 55.3
            p = 1.39
        else:
            a = 0
            p = 1

        self.C_adq = unidades.Currency(a*self.Di.ft**p*self.L*CI/BI)
        self.C_inst = unidades.Currency(self.C_adq*self.kwargs["f_install"])

    def writeListtoJSON(self, kwarg, key, value):
        """Personalizar en el caso de equipos con listas complejas"""
        kwarg_list = {}
        if key == "material":
            kwarg_list["material"] = value[0]
            kwarg_list["class"] = value[1]
            kwarg_list["rugosidad"] = value[2]
            kwarg_list["nominal"] = value[3]
            kwarg_list["Di"] = value[4]
            kwarg_list["w"] = value[5]
            kwarg_list["De"] = value[6]
            kwarg_list["weight"] = value[7]
            kwarg_list["volume"] = value[8]
            kwarg_list["area"] = value[9]
            kwarg_list["ind_material"] = value[10]
            kwarg_list["ind_size"] = value[11]
        elif key == "accesorios":
            for i, accesorio in enumerate(value):
                ac = {}
                ac["ind_equip"] = value[i][0]
                ac["ind_size"] = value[i][1]
                ac["K"] = value[i][2]
                ac["count"] = value[i][3]
                ac["type"] = value[i][4]
                ac["D_mm"] = value[i][5]
                ac["D_in"] = value[i][6]
                ac["description"] = value[i][7]
                kwarg_list[i] = ac
        kwarg[key] = kwarg_list

    def readListFromJSON(self, data, key):
        """Read list from file, customize in entities with complex list"""
        kwarg = []
        if key == "material":
            kwarg.append(data[key]["material"])
            kwarg.append(data[key]["class"])
            kwarg.append(data[key]["rugosidad"])
            kwarg.append(data[key]["nominal"])
            kwarg.append(data[key]["Di"])
            kwarg.append(data[key]["w"])
            kwarg.append(data[key]["De"])
            kwarg.append(data[key]["weight"])
            kwarg.append(data[key]["volume"])
            kwarg.append(data[key]["area"])
            kwarg.append(data[key]["ind_material"])
            kwarg.append(data[key]["ind_size"])

        elif key == "accesorios":
            for i, accesorio in data[key].items():
                ac = []
                ac.append(accesorio["ind_equip"])
                ac.append(accesorio["ind_size"])
                ac.append(accesorio["K"])
                ac.append(accesorio["count"])
                ac.append(accesorio["type"])
                ac.append(accesorio["D_mm"])
                ac.append(accesorio["D_in"])
                ac.append(accesorio["description"])
                kwarg.append(ac)
        return kwarg

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["L"] = self.L
        state["rho"] = self.rho
        state["mu"] = self.mu
        state["material"] = self.material
        state["Dn"] = self.Dn
        state["rugosidad"] = self.rugosidad
        state["De"] = self.De
        state["w"] = self.w
        state["Di"] = self.Di
        state["eD"] = self.eD
        state["seccion"] = self.seccion
        state["A"] = self.A
        state["V"] = self.V
        state["Re"] = self.Re
        state["K"] = self.K
        state["DeltaP_h"] = self.DeltaP_h
        state["DeltaP_ac"] = self.DeltaP_ac
        state["f"] = self.f
        state["DeltaP_f"] = self.DeltaP_f
        state["DeltaP_v"] = self.DeltaP_v
        state["DeltaP"] = self.DeltaP
        state["DeltaP_100ft"] = self.DeltaP_100ft
        state["Pout"] = self.Pout
        state["Tout"] = self.Tout
        state["Heat"] = self.Heat
        state["Pin"] = self.Pin
        state["Pout"] = self.Pout
        state["statusCoste"] = self.statusCoste
        if self.statusCoste:
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.L = unidades.Length(state["L"])
        self.rho = unidades.Density(state["rho"])
        self.mu = unidades.Viscosity(state["mu"])
        self.material = state["material"]
        self.Dn = state["Dn"]
        self.rugosidad = unidades.Length(state["rugosidad"])
        self.De = unidades.Length(state["De"])
        self.w = unidades.Length(state["w"])
        self.Di = unidades.Length(state["Di"])
        self.eD = unidades.Dimensionless(state["eD"])
        self.seccion = unidades.Area(state["seccion"])
        self.A = unidades.Area(state["A"])
        self.V = unidades.Speed(state["V"])
        self.Re = unidades.Dimensionless(state["Re"])
        self.K = unidades.Dimensionless(state["K"])
        self.DeltaP_h = unidades.DeltaP(state["DeltaP_h"])
        self.DeltaP_ac = unidades.DeltaP(state["DeltaP_ac"])
        self.f = unidades.Dimensionless(state["f"])
        self.DeltaP_f = unidades.DeltaP(state["DeltaP_f"])
        self.DeltaP_v = unidades.DeltaP(state["DeltaP_v"])
        self.DeltaP = unidades.DeltaP(state["DeltaP"])
        self.DeltaP_100ft = unidades.Dimensionless(state["DeltaP_100ft"])
        self.Tout = unidades.Temperature(state["Tout"])
        self.Heat = unidades.Power(state["Heat"])
        self.Pin = unidades.Pressure(state["Pin"])
        self.Pout = unidades.Pressure(state["Pout"])
        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.C_adq = unidades.Currency(state["C_adq"])
            self.C_inst = unidades.Currency(state["C_inst"])
        self.salida = [None]

    def propTxt(self):
        txt = "#---------------"
        txt += tr("pychemqt", "Catalog")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(9))

        if self.kwargs["accesorios"]:
            txt += os.linesep + "#---------------"
            txt += tr("pychemqt", "Fittings")
            txt += "-----------------#"+os.linesep
            txt += self.propertiesToText(9)
            for accesorio in self.kwargs["accesorios"]:
                txt += "%5i %-22s\t %s" % (accesorio[3], accesorio[7],
                                           accesorio[2]) + os.linesep

        txt += os.linesep + "#---------------"
        txt += tr("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(11, 24))

        if self.kwargs["thermal"]:
            txt += self.propertiesToText(range(24, 26))

        if self.statusCoste:
            txt += os.linesep+"#---------------"
            txt += tr(
                "pychemqt", "Preliminary Cost Estimation")
            txt += "-----------------#" + os.linesep
            txt += self.propertiesToText(range(26, 31))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(tr("pychemqt", "Material"), "material", str),
             (tr("pychemqt", "Nominal Diameter"),
              "Dn", str),
             (tr("pychemqt", "Length"), "L",
              unidades.Length),
             (tr("pychemqt", "Roughness"), "rugosidad",
              unidades.Length),
             (tr("pychemqt", "Internal Diamter"), "Di",
              unidades.Length),
             (tr("pychemqt", "External Diamter"), "De",
              unidades.Length),
             (tr("pychemqt", "Thickness"), "w",
              unidades.Length),
             (tr("pychemqt", "Transversal section"),
              "seccion", unidades.Area),
             (tr("pychemqt", "External Area"), "A",
              unidades.Area),
             (tr("pychemqt", "K total"), "K",
              unidades.Dimensionless),
             (tr("pychemqt", "Fittings"), "accesorios",
              None),
             (tr("pychemqt", "Method"),
              ("TEXT_METODO", "metodo"),  str),
             (tr("pychemqt", "Input Pressure"), "Pin",
              unidades.Pressure),
             (tr("pychemqt", "Output Pressure"), "Pout",
              unidades.Pressure),
             (tr("pychemqt", "ΔP Total"), "DeltaP",
              unidades.DeltaP),
             (tr("pychemqt", "ΔP friction"), "DeltaP_f",
              unidades.DeltaP),
             (tr("pychemqt", "ΔP fittings"), "DeltaP_ac",
              unidades.DeltaP),
             (tr("pychemqt", "ΔP elevation"), "DeltaP_h",
              unidades.DeltaP),
             (tr("pychemqt", "ΔP acceleration"),
              "DeltaP_v", unidades.DeltaP),
             (tr("pychemqt", "Thermal Condition"),
              ("TEXT_THERMAL", "thermal"),  str),
             (tr("pychemqt", "Fluid Speed"), "V",
              unidades.Speed),
             (tr("pychemqt", "Reynolds number"), "Re",
              unidades.Dimensionless),
             (tr("pychemqt", "Relative roughness"), "eD",
              unidades.Dimensionless),
             (tr("pychemqt", "Factor Friction"), "f",
              unidades.Dimensionless),
             (tr("pychemqt", "Output Temperature"), "Tout",
              unidades.Temperature),
             (tr("pychemqt", "Heat Transfer"), "Heat",
              unidades.Power),
             (tr("pychemqt", "Base index"),
              "Base_index", float),
             (tr("pychemqt", "Current index"),
              "Current_index", float),
             (tr("pychemqt", "Install factor"),
              "f_install", float),
             (tr("pychemqt", "Purchase Cost"), "C_adq",
              unidades.Currency),
             (tr("pychemqt", "Installed Cost"), "C_inst",
              unidades.Currency)]
        return l

    def propertiesListTitle(self, index):
        lista = []
        for accesorio in self.kwargs["accesorios"]:
            lista.append("%3i %s" % (accesorio[3], accesorio[7]))
        return lista


if __name__ == '__main__':
    import doctest
    doctest.testmod()
