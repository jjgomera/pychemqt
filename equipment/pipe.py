#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
#   Library for pipe equipment definition
###############################################################################

import os

from PyQt5.QtWidgets import QApplication
from scipy.constants import g, pi

from lib import unidades
from lib.utilities import representacion
from lib.corriente import Corriente
from lib.friction import f_friccion
from lib.adimensional import Re
from equipment.parents import equipment
from equipment.heatExchanger import Heat_Exchanger


class Pipe(equipment):
    """Clase que modela un una tubería

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        metodo: Método de cálculo
            0   -   Flujo unifásico
            1   -   Flujo de agua (Eq. de Hazen-Williams)
            2   -   Flujo de vapor (Eq de Fritzsche)
            3   -   Flujo de gas isotérmico
            4   -   Flujo bifásico (Método Baker)
            5   -   Flujo bifásico (Método de Beggs and Brill)
        thermal: Método de funcionamiento adiabático:
            0   -   Adiabático
            1   -   Flujo de calor fíjo
            2   -   Calcular el intercambio de calor conocidas U y Tª externa
        Material: Array con los datos del material
            0   -   Nombre
            1   -   Clase
            2   -   Rugosidad, mm
            3   -   Catalogo
            4   -   Di, mm
            5   -   Espesor, mm
            6   -   De, mm
            7   -   Peso, kg/m
            8   -   Volumen, m³/100m
            9   -   Superficie, m²/100m
            10 -    Indice de la lista de materiales
            11 -    Indice de la lista de de dimensiones
        Accesorios: Array en el que cada entrada representa un equipo con perdida de carga:
            Indice tipo equipo
            Indice diametro
            K
            Numero de equipos
            Tipo
            Diametro en mm
            Diametro en pulgadas
            Texto explicativo

        l: Longitud en metros de la tubería
        h: diferencia de altura entre la entrada y la salida de la tubería, m
        K: Coeficiente del los accesorios
        C:Coeficiente para el método de Hazen-Williams
        T_ext: Temperatura en el exterior de la tubería
        U: Coeficiente global de transimisión de calor entre la tubería y el exterior
        Q: Calor transmitido a través de las paredes de la tubería

    Coste:
        solo disponible para las tuberías de acero, Ref Darby pag 217

    >>> agua=agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    >>> tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Cast Iron", "Class A", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06])
    >>> print tuberia.Di, tuberia.V,  tuberia.Re, tuberia.DeltaP
    0.1524 0.0626814744044 9427.99142792 1.8362005711
    """
    title = QApplication.translate("pychemqt", "Pipe")
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
        QApplication.translate("pychemqt", "Single Phase flow"),
        QApplication.translate("pychemqt", "Water (Hazen-Williams)"),
        QApplication.translate("pychemqt", "Steam (Fritzsche)"),
        QApplication.translate("pychemqt", "Isotermic gas flow"),
        QApplication.translate("pychemqt", "Two Phase flow (Baker method)"),
        QApplication.translate("pychemqt", "Two Phase flow (Beggs and Brill method)"))
    TEXT_THERMAL = (
        QApplication.translate("pychemqt", "Adiabatic"),
        QApplication.translate("pychemqt", "Heat flux"),
        QApplication.translate("pychemqt", "Heat transfer"))

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"] and self.kwargs["material"] and \
                self.kwargs["material"][0] in ['Stainless Steel (ANSI)', 'Steel Galvanised (ANSI)', 'Steel (ANSI)']:
            self.statusCoste = True
        else:
            self.statusCoste = False
            self.C_adq = unidades.Currency(None)
            self.C_inst = unidades.Currency(None)

        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            return
        if not self.kwargs["l"]:
            self.msg = QApplication.translate("pychemqt", "undefined pipe length")
            self.status = 0
            return
        if not self.kwargs["material"]:
            self.msg = QApplication.translate("pychemqt", "undefined material")
            self.status = 0
            return

        if self.kwargs["thermal"] == 1 and not self.kwargs["Q"]:
            self.msg = QApplication.translate("pychemqt", "undefined heat flux")
            self.status = 0
            return
        elif self.kwargs["thermal"] == 2 and (not self.kwargs["T_ext"] or not self.kwargs["U"]):
            self.msg = QApplication.translate("pychemqt", "undefined heat transfer conditions")
            self.status = 0
            return

        if self.kwargs["metodo"] == 1 and not self.kwargs["C"]:
            self.msg = QApplication.translate("pychemqt", "undefined C William Factor")
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

        self.material = self.kwargs["material"][0] + " " + self.kwargs["material"][1]
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

        self.DeltaP = unidades.Pressure(self.DeltaP_f+self.DeltaP_ac+self.DeltaP_h)
        self.DeltaP_100ft = self.DeltaP*100/self.L.ft
        self.Pout = unidades.Pressure(self.kwargs["entrada"].P-self.DeltaP)

        if self.kwargs["thermal"] == 0:
            self.Tout = self.kwargs["entrada"].T
            self.Heat = unidades.Power(0)
        else:
            cambiador = Heat_Exchanger()
            cambiador.calculo(entrada=self.kwargs["entrada"], modo=self.kwargs["thermal"], Heat=self.kwargs["Q"], deltaP=self.DeltaP, A=self.A, U=self.kwargs["U"], Text=self.kwargs["Text"])
            self.Tout = cambiador.salida[0].T
            self.Heat = cambiador.Heat

        self.salida = [self.kwargs["entrada"].clone(T=self.Tout, P=self.Pout)]
        self.Pin = self.kwargs["entrada"].P
        self.Pout = self.salida[0].P

    def __DeltaP_friccion(self):
        """Método para el calculo de la perdida de presión"""
        if self.kwargs["metodo"] == 0:
            delta = unidades.Pressure(self.L*self.V**2/self.Di*self.f*self.rho/2)
        elif self.kwargs["metodo"] == 1:
            delta = unidades.Pressure((self.kwargs["entrada"].Q.galUSmin*self.L.ft**0.54/0.442/self.Di.inch**2.63/self.kwargs["C"])**(1./0.54), "psi")
        elif self.kwargs["metodo"] == 2:
            delta = unidades.Pressure(2.1082*self.L.ft*self.kwargs["entrada"].caudalmasico.lbh**1.85/self.rho.lbft3/1e7/self.Di.inch**4.97, "psi")
        elif self.kwargs["metodo"] == 3:
            pass

        elif self.kwargs["metodo"] == 4:
            pass

        elif self.kwargs["metodo"] == 5:
            pass

        return delta

    def coste(self):
        """
        Coste solo disponible para las tuberías de acero
        Ref Darby pag 217

        kwargs:
            schedule: Clase de acero
        """
        codigo = str(self.kwargs["material"][1])
        if codigo in ('Sch. 40', 'Sch. 5S'):
            a = 30.
            p = 1.31
        elif codigo in ('Sch. 80', 'Sch. 10S'):
            a = 38.1
            p = 1.35
        elif codigo in ('Sch. 160', 'Sch.  40S'):
            a = 55.3
            p = 1.39
        else:
            a = 0
            p = 1

        self.C_adq = unidades.Currency(a*self.Di.ft**p*self.L * self.kwargs["Current_index"] / self.kwargs["Base_index"])
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
        txt="#---------------"+QApplication.translate("pychemqt", "Catalog")+"-----------------#"+os.linesep
        txt+="%-25s\t %s %s" %(QApplication.translate("pychemqt", "Material"), self.kwargs["material"][0], self.kwargs["material"][1])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Nominal Diameter"), self.kwargs["material"][3])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Length"), self.L.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Roughness"), self.rugosidad.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Internal Diamter"), self.Di.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Diamter"), self.De.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thickness"), self.w.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Transversal section"), self.seccion.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Area"), self.A.str)+os.linesep

        if self.kwargs["accesorios"]:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Fittings")+"-----------------#"+os.linesep
            txt+="%-25s\t %s" %("K "+QApplication.translate("pychemqt", "Total"), self.K)+os.linesep
            for accesorio in self.kwargs["accesorios"]:
                txt+="%5i %-22s\t %s" %(accesorio[3], accesorio[7], accesorio[2])+os.linesep

        txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Method"), self.TEXT_METODO[self.kwargs["metodo"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.kwargs["entrada"].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP Total", None), self.DeltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP friction", None), self.DeltaP_f.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP fittings", None), self.DeltaP_ac.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP elevation", None), self.DeltaP_h.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP acceleration", None), self.DeltaP_v.str)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Thermal Condition"), self.TEXT_THERMAL[self.kwargs["thermal"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Fluid Speed"), self.V.str)+os.linesep
        txt+="%-25s\t %s" %("Reynolds", self.Re)+os.linesep
        txt+="%-25s\t %s" %("ε/D", self.eD)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Factor Friction"), self.f)+os.linesep

        if self.kwargs["thermal"]:
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.Tout.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Heat Transfer"), self.Heat.str)+os.linesep

        if self.statusCoste:
            txt += os.linesep
            txt += "#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt += "%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Material"), "material", str),
             (QApplication.translate("pychemqt", "Nominal Diameter"), "Dn", str),
             (QApplication.translate("pychemqt", "Length"), "L", unidades.Length),
             (QApplication.translate("pychemqt", "Roughness"), "rugosidad", unidades.Length),
             (QApplication.translate("pychemqt", "Internal Diamter"), "Di", unidades.Length),
             (QApplication.translate("pychemqt", "External Diamter"), "De", unidades.Length),
             (QApplication.translate("pychemqt", "Thickness"), "w", unidades.Length),
             (QApplication.translate("pychemqt", "Transversal section"), "seccion", unidades.Area),
             (QApplication.translate("pychemqt", "External Area"), "A", unidades.Area),
             (QApplication.translate("pychemqt", "K total"), "K", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Fittings"), "accesorios", None),
             (QApplication.translate("pychemqt", "Method"), ("TEXT_METODO", "metodo"),  str),
             (QApplication.translate("pychemqt", "Input Pressure"), "Pin", unidades.Pressure),
             (QApplication.translate("pychemqt", "Output Pressure"), "Pout", unidades.Pressure),
             (QApplication.translate("pychemqt", "ΔP Total", None), "DeltaP", unidades.DeltaP),
             (QApplication.translate("pychemqt", "ΔP friction", None), "DeltaP_f", unidades.DeltaP),
             (QApplication.translate("pychemqt", "ΔP fittings", None), "DeltaP_ac", unidades.DeltaP),
             (QApplication.translate("pychemqt", "ΔP elevation", None), "DeltaP_h", unidades.DeltaP),
             (QApplication.translate("pychemqt", "ΔP acceleration", None), "DeltaP_v", unidades.DeltaP),
             (QApplication.translate("pychemqt", "Thermal Condition"), ("TEXT_THERMAL", "thermal"),  str),
             (QApplication.translate("pychemqt", "Fluid Speed"), "V", unidades.Speed),
             (QApplication.translate("pychemqt", "Reynolds number"), "Re", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Relative roughness"), "eD", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Factor Friction"), "f", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Output Temperaturet"), "Tout", unidades.Temperature),
             (QApplication.translate("pychemqt", "Heat Transfer"), "Heat", unidades.Power),
             (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", unidades.Currency),
             (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", unidades.Currency)]
        return l

    def propertiesListTitle(self, index):
        lista = []
        for accesorio in self.kwargs["accesorios"]:
            lista.append("%3i %s" % (accesorio[3], accesorio[7]))
        return lista


if __name__ == '__main__':
#    import doctest
#    doctest.testmod()


#    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
#    tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Cast Iron", "Class A", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, 0, 2])
#    print tuberia.Di, tuberia.V,  tuberia.Re, tuberia.DeltaP
    # import config
    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Steel (ANSI)", "Sch. 40", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, -1, 2], notas="Tuberia")
    print((tuberia.DeltaP))
