#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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
# library for heat exchanger equipment calculation
# - Heat_Exchanger
# - Fired_Heater
# - Shell_Tube
# - Hairpin
###############################################################################


from math import sqrt, exp, log
import os

from numpy import arccos, sin, cos, tanh
from scipy.constants import g, pi
from scipy.optimize import fsolve
from tools.qt import translate

from equipment.parents import equipment
from lib import unidades
from lib.adimensional import Re, Pr, Gr, Gz
from lib.friction import f_friccion
from lib.heatTransfer import *  # noqa


class Heat_Exchanger(equipment):
    """Define a simple heat exchanger, only make energy balance

    Parameters:
        entrada: Corriente instance to define input stream
        Heat: global heat value exchanged
        Tout: Output temperature
        DeltaT: Increase temperature in stream
        A: area for heat exchange
        U: Global coefficient of heat transmision
        Text: Ambient external temperature
        DeltaP: Opcional pressure losses of equipment

    >>> from lib.corriente import Corriente
    >>> agua = Corriente(T=300, P=101325., caudalMasico=1., fraccionMolar=[1.])
    >>> Cambiador = Heat_Exchanger(entrada=agua, Tout=350)
    >>> print("%6g" % Cambiador.HeatCalc.MJh)
    753.063
    """
    title = translate("equipment", "Heat Exchanger")
    help = ""
    kwargs = {
        "entrada": None,
        "Heat": 0.0,
        "Tout": 0.0,
        "DeltaT": 0.0,
        "A": 0.0,
        "U": 0.0,
        "Text": 0.0,
        "DeltaP": 0.0}

    kwargsInput = ("entrada", )
    kwargsValue = ("Heat", "Tout", "DeltaT", "A", "U", "Text", "DeltaP")
    calculateValue = ("HeatCalc", "ToutCalc")

    def cleanOldValues(self, **kwargs):
        """Clean kwargs values for old heat exchanger definition"""
        if kwargs.get("Tout", 0):
            self.kwargs["DeltaT"] = 0
            self.kwargs["Heat"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        elif kwargs.get("DeltaT", 0):
            self.kwargs["Tout"] = 0
            self.kwargs["Heat"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        elif kwargs.get("Heat", 0):
            self.kwargs["Tout"] = 0
            self.kwargs["DeltaT"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        self.kwargs.update(kwargs)

    @property
    def isCalculable(self):
        """
        modo: unknown variable to calculate
            1 - Known output temperature, calculate other variables
            2 - known heat exchange
            3 - known heat exchanger characteristic, calculate output stream
        """
        if not self.kwargs["entrada"]:
            self.msg = translate("equipment", "undefined input")
            self.status = 0
            self.modo = 0
            return

        if self.kwargs["A"] and self.kwargs["U"] and self.kwargs["Text"]:
            self.modo = 3
        elif self.kwargs["Heat"]:
            self.modo = 2
        elif self.kwargs["Tout"] or self.kwargs["DeltaT"]:
            self.modo = 1
        else:
            self.msg = translate(
                "equipment", "undefined output temperature specification")
            self.status = 0
            self.modo = 0

        if self.modo:
            self.msg = ""
            self.status = 1
            return True

    def calculo(self):
        entrada = self.kwargs["entrada"]
        self.deltaP = unidades.DeltaP(self.kwargs["DeltaP"])
        self.HeatCalc = unidades.Power(self.kwargs["Heat"])
        if self.kwargs["Tout"]:
            Tout = unidades.Temperature(self.kwargs["Tout"])
        elif self.kwargs["DeltaT"]:
            Tout = unidades.Temperature(entrada.T+self.kwargs["DeltaT"])
        A = unidades.Area(self.kwargs["A"])
        U = unidades.HeatTransfCoef(self.kwargs["U"])
        Text = unidades.Temperature(self.kwargs["Text"])

        if self.modo == 1:
            self.salida = [entrada.clone(T=Tout, P=entrada.P-self.deltaP)]
            self.HeatCalc = unidades.Power(self.salida[0].h-entrada.h)
        else:
            if self.modo == 2:
                self.HeatCalc = unidades.Power(0)
            else:
                self.HeatCalc = unidades.Power(A*U*(Text-entrada.T))

            def f():
                output = entrada.clone(T=T, P=entrada.P-self.deltaP)
                return output.h-entrada.h-self.HeatCalc
            T = fsolve(f, entrada.T)[0]
            if T > max(Text, entrada.T) or T < min(Text, entrada.T):
                T = self.Text
            self.salida = [entrada.clone(T=T, P=entrada.P-self.deltaP)]

        self.Tin = entrada.T
        self.ToutCalc = self.salida[0].T
        self.deltaT = unidades.DeltaT(self.ToutCalc-entrada.T)

    def propTxt(self):
        txt = "#---------------"
        txt += translate("equipment", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(5))
        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(translate("equipment", "Input Temperature"),
              "Tin", unidades.Temperature),
             (translate("equipment", "Output Temperature"),
              "ToutCalc", unidades.Temperature),
             (translate("equipment", "Temperature increase"),
              "deltaT", unidades.DeltaT),
             (translate("equipment", "Pressure increase"), "deltaP", unidades.DeltaP),
             (translate("equipment", "Heat"), "HeatCalc", unidades.Power)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Tin"] = self.Tin
        state["ToutCalc"] = self.ToutCalc
        state["deltaT"] = self.deltaT
        state["deltaP"] = self.deltaP
        state["HeatCalc"] = self.HeatCalc

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Tin = unidades.Temperature(state["Tin"])
        self.ToutCalc = unidades.Temperature(state["ToutCalc"])
        self.deltaT = unidades.DeltaT(state["deltaT"])
        self.deltaP = unidades.DeltaP(state["deltaP"])
        self.HeatCalc = unidades.Power(state["HeatCalc"])
        self.salida = [None]


class Fired_Heater(equipment):
    """Class to model a fire heater

    Parameters
    ----------
    entrada : Corriente
        Input stream to equipment
    Tout : float
        Desired output stream temperature, [K]
    deltaP: float, optional
        Pressure loss in equipment, [Pa]
    Hmax: float, optional
        Design heat transfer of equipment, maximum equipment capacity, [kJ/s]
    eficiencia: float, optional, default 0.75
        Thermic efficiency in combustion process, [-]
    poderCalorifico: float, optional, default 900 Btu/stdft³
        fuel heat power, [-]

    Coste:
        tipo:
            0   -   Box type
            1   -   Cylindrical type
        subtipoBox:
            0   -   Process heater
            1   -   Pyrolysis
            2   -   Reformer without catalyst
        subtipoCylindrical:
            0   -   Cylindrical
            1   -   Dowtherm
        material:
            0   -   Carbon steel
            1   -   CrMo steel
            2   -   Stainless
        P_dis: Design pressure of equipment, if no specified use the input
            stream pressure

    >>> from lib.corriente import Corriente
    >>> agua = Corriente(T=400, P=101325., caudalMasico=1., fraccionMolar=[1.])
    >>> cambiador = Fired_Heater(entrada=agua, Tout=450)
    >>> print("%6g" % cambiador.Heat.MJh)
    357.685
    """
    title = translate("equipment", "Fired Heater")
    help = os.environ["pychemqt"] + "doc/fireHeater.htm"
    kwargs = {
        "entrada": None,
        "Tout": 0.0,
        "deltaP": 0.0,
        "Hmax": 0.0,
        "eficiencia": 0.0,
        "poderCalorifico": 0.0,

        "f_install": 1.3,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "tipo": 0,
        "subtipoBox": 0,
        "subtipoCylindrical": 0,
        "material": 0,
        "P_dis": 0.0}

    kwargsInput = ("entrada", )
    kwargsValue = ("Tout", "deltaP", "Hmax", "eficiencia", "poderCalorifico",
                   "P_dis")
    kwargsList = ("tipo", "subtipoBox", "subtipoCylindrical", "material")
    calculateValue = ("CombustibleRequerido", "Heat")
    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 3
    salida = [None]

    TEXT_TIPO = [translate("equipment", "Box Type"),
                 translate("equipment", "Cylindrical type")]
    TEXT_SUBTIPOBOX = [
        translate("equipment", "Process heater"),
        translate("equipment", "Pyrolysis"),
        translate("equipment", "Reforme without catalysis")]
    TEXT_SUBTIPOCYLINDRICAL = [
        translate("equipment", "Cylindrical"),
        translate("equipment", "Dowtherm")]
    TEXT_MATERIAL = [translate("equipment", "Carbon steel"),
                     translate("equipment", "Cr-Mo steel"),
                     translate("equipment", "Stainless steel")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = translate("equipment", "undefined input")
            self.status = 0
            return

        if not self.kwargs["Tout"]:
            self.msg = translate("equipment", "undefined output temperature condition")
            self.status = 0
            return
        if self.kwargs["Tout"] <= self.kwargs["entrada"].T:
            self.msg = translate("equipment", "bad output temperature condition")
            self.status = 0
            return

        if not self.kwargs["eficiencia"]:
            self.msg = translate("equipment", "using default efficiency")
            self.status = 3
            return True
        if not self.kwargs["poderCalorifico"]:
            self.msg = translate("equipment", "using default fuel calorific value")
            self.status = 3
            return True

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        entrada = self.kwargs["entrada"]
        self.Tout = unidades.Temperature(self.kwargs["Tout"])
        self.deltaP = unidades.DeltaP(self.kwargs["deltaP"])
        self.Hmax = unidades.Power(self.kwargs["Hmax"])
        if self.kwargs["eficiencia"]:
            eficiencia = self.kwargs["eficiencia"]
        else:
            eficiencia = 0.75

        if self.kwargs["poderCalorifico"]:
            poderCalorifico = self.kwargs["poderCalorifico"]
        else:
            poderCalorifico = 900

        salida = entrada.clone(T=self.Tout, P=entrada.P-self.deltaP)
        Ho = entrada.h
        H1 = salida.h
        Heat = unidades.Power(H1-Ho)

        if self.Hmax and Heat > self.Hmax:
            self.Heat = unidades.Power(self.Hmax)
            To = (entrada.T+self.Tout)/2
            T = fsolve(lambda T: entrada.clone(
                T=T, P=entrada.P-self.deltaP).h-Ho-self.Hmax, To)[0]
            self.salida = [entrada.clone(T=T, P=entrada.P-self.deltaP)]
        else:
            self.Heat = Heat
            self.salida = [salida]

        fuel = self.Heat.Btuh/poderCalorifico/eficiencia
        self.CombustibleRequerido = unidades.VolFlow(fuel, "ft3h")
        self.deltaT = unidades.DeltaT(self.salida[0].T-entrada.T)
        self.eficiencia = unidades.Dimensionless(eficiencia)
        self.poderCalorifico = unidades.Dimensionless(poderCalorifico)
        self.Tin = entrada.T
        self.Tout = self.salida[0].T

    def coste(self):
        """
        tipo:
            0   -   Box type
            1   -   Cylindrical type
        subtipoBox:
            0   -   Process heater
            1   -   Pyrolysis
            2   -   Reformer without catalyst
        subtipoCylindrical:
            0   -   Cylindrical
            1   -   Dowtherm
        material:
            0   -   Carbon steel
            1   -   CrMo steel
            2   -   Stainless
        P_dis: Presión de diseño del equipo
        """
        if self.kwargs["P_dis"]:
            P_dis = unidades.Pressure(self.kwargs["P_dis"])
        else:
            P_dis = self.kwargs["entrada"].P
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]

        if self.kwargs["tipo"] == 0:  # Boxtype
            k = [25.5, 33.8, 45.][self.kwargs["material"]]
            Fd = [0, 0.1, 0.35][self.kwargs["subtipoBox"]]

            if P_dis.psi <= 500:
                Fp = 0
            elif P_dis.psi <= 1000:
                Fp = 0.1
            elif P_dis.psi <= 1500:
                Fp = 0.15
            elif P_dis.psi <= 2000:
                Fp = 0.25
            elif P_dis.psi <= 2500:
                Fp = 0.4
            else:
                Fp = 0.6

            C = k*(1+Fd+Fp)*self.Heat.MBtuh**0.86*1000

        else:  # Cylindrical
            k = [27.3, 40.2, 42.][self.kwargs["material"]]
            Fd = [0, 0.33][self.kwargs["subtipoCylindrical"]]

            if P_dis.psi <= 500:
                Fp = 0
            elif P_dis.psi <= 1000:
                Fp = 0.15
            else:
                Fp = 0.2

            C = k*(1+Fd+Fp)*self.Heat.MBtuh**0.82*1000

        self.P_dis = P_dis
        self.C_adq = unidades.Currency(C * CI / BI)
        self.C_inst = unidades.Currency(self.C_adq * self.kwargs["f_install"])

    def propTxt(self):
        txt = "#---------------"
        txt += translate("equipment", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(5))

        stimated = " (%s)" % translate("equipment", "stimated")
        txt += self.propertiesToText(5, kwCheck=True, kwSuffix=stimated,
                                     kwKey="eficiencia", kwValue=0.0)
        txt += self.propertiesToText(6, kwCheck=True, kwSuffix=stimated,
                                     kwKey="poderCalorifico", kwValue=0.0)
        txt += self.propertiesToText(7)

        if self.statusCoste:
            txt += os.linesep+"#---------------"
            txt += translate("equipment", "Preliminary Cost Estimation")
            txt += "-----------------#" + os.linesep
            txt += self.propertiesToText(range(8, 18))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(translate("equipment", "Input Temperature"),
              "Tin", unidades.Temperature),
             (translate("equipment", "Output Temperature"),
              "Tout", unidades.Temperature),
             (translate("equipment", "Temperature increase"),
              "deltaT", unidades.DeltaT),
             (translate("equipment", "Pressure increase"), "deltaP", unidades.DeltaP),
             (translate("equipment", "Maximum heat"), "Hmax", unidades.Power),
             (translate("equipment", "Thermal Efficiency"),
              "eficiencia", unidades.Dimensionless),
             (translate("equipment", "Fuel Heating Value"),
              "poderCalorifico", unidades.Dimensionless),
             (translate("equipment", "Required Fuel"),
              "CombustibleRequerido", unidades.VolFlow),
             (translate("equipment", "Base index"), "Base_index", float),
             (translate("equipment", "Current index"), "Current_index", float),
             (translate("equipment", "Install factor"), "f_install", float),
             (translate("equipment", "FireHeater type"), ("TEXT_TIPO", "tipo"), str),
             (translate("equipment", "Cylindrical type"),
              ("TEXT_SUBTIPOCYLINDRICAL", "subtipoCylindrical"), str),
             (translate("equipment", "Box type"),
              ("TEXT_SUBTIPOBOX", "subtipoBox"), str),
             (translate("equipment", "Material"), ("TEXT_MATERIAL", "material"), str),
             (translate("equipment", "Design Pressure"), "P_dis", unidades.Pressure),
             (translate("equipment", "Purchase Cost"), "C_adq", unidades.Currency),
             (translate("equipment", "Installed Cost"), "C_inst", unidades.Currency)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Tout"] = self.Tout
        state["deltaP"] = self.deltaP
        state["Hmax"] = self.Hmax
        state["Heat"] = self.Heat
        state["CombustibleRequerido"] = self.CombustibleRequerido
        state["deltaT"] = self.deltaT
        state["eficiencia"] = self.eficiencia
        state["poderCalorifico"] = self.poderCalorifico
        state["Tin"] = self.Tin
        state["statusCoste"] = self.statusCoste
        if self.statusCoste:
            state["P_dis"] = self.P_dis
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Tout = unidades.Temperature(state["Tout"])
        self.deltaP = unidades.DeltaP(state["deltaP"])
        self.Hmax = unidades.Power(state["Hmax"])
        self.Heat = unidades.Power(state["Heat"])
        self.CombustibleRequerido = unidades.VolFlow(
            state["CombustibleRequerido"])
        self.deltaT = unidades.DeltaT(state["deltaT"])
        self.eficiencia = unidades.Dimensionless(state["eficiencia"])
        self.poderCalorifico = unidades.Dimensionless(state["poderCalorifico"])
        self.Tin = unidades.Temperature(state["Tin"])
        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.P_dis = unidades.Pressure(state["P_dis"])
            self.C_adq = unidades.Currency(state["C_adq"])
            self.C_inst = unidades.Currency(state["C_inst"])
        self.salida = [None]


class Hairpin(equipment):
    """Class to model double pipe section heat exchanger (hairpin)

    Parameters:
        entradaInterior: Corriente instance to define the fluid stream that
            flow at internal section
        entradaExterior: Corriente instance to define the fluid stream that
            flow at external (anulli) section

        modo: Calculate module
            0 - Design
            1 - Rating
        flujo: Flow type
            0 - Counterflow
            1 - Parallelflow
        orientacion: Pipe orientation
            0 - Horizontal
            1 - Vertical, internal down
            2 - Vertical, internal up
        metodo:
            0 - Mean temperature
            1 - Split the pipe in segments
        tubesideLaminar: Method to calculate the global heat transfer
            coefficient in laminar flow for tubeside
            0 - Eubank-Proctor
            1 - VDI mean Nusselt
            2 - Hausen
            3 - Sieder-Tate
        tubesideTurbulent: Method to calculate the global heat transfer
            coefficient in turbulent flow for tubeside
            0 - Sieder-Tate
            1 - Colburn
            2 - Dittus-Boelter
            3 - ESDU
            4 - Gnielinski
            5 - VDI mean Nusselt

        LTube: Tube length
        DeeTube: External diameter of annulli
        DeTube: External diameter of internal pipe
        DiTube: Internal diameter of internal pipe
        wTube: Pipe width
        rTube: Internal roughness of pipe
        kTube: Thermal conductivity

        tubeFouling: Fouling at tubeside
        annulliFouling: Fouling at annulliside

        tubeFinned: Boolean to use finned tube
        hFin: Finned height
        thicknessBaseFin: Thickness of the bottom of fin
        thicknessTopFin: Thickness of the top of fin
        rootDoFin: External diameter in the bottom of fin
        kFin: Thermal conductiviti of material of fin
        nFin: Fin count per meter of pipe

        tubeTout: Output temperature of fluid in tubeside
        annulliTout: Output temperature of fluid in annulliside

    Coste:
        material:
            0 - Carbon steel/carbon steel
            1 - Carbon steel/304 stainless
            2 - Carbon steel/316 stainless
        P_dis: Design pressure

    >>> from lib.corriente import Corriente
    >>> kw = {"ids": [62], "fraccionMolar": [1.]}
    >>> caliente = Corriente(T=90+273.15, P=361540., caudalMasico=0.36, **kw)
    >>> fria = Corriente(T=20+273.15, P=101325., caudalMasico=500/3600., **kw)
    >>> Cambiador = Hairpin(entradaTubo=caliente, entradaExterior=fria, \
                            modo=1, DiTube=0.0525, DeTube=0.0603, LTube=2.5, \
                            DeeTube=0.0779, kTube=54, rTube=0.0459994e-3, \
                            annulliFouling=0.000352, tubeFouling=0.000176)
    >>> print("%6g %6g" % (Cambiador.ReTube, Cambiador.ReAnnulli))
    27783.3 1277.55
    >>> print("%6g %6g" % (Cambiador.hTube.kWm2K, Cambiador.hAnnulli.kWm2K))
    1555.53 52.8267
    """
    title = translate("equipment", "Hairpin Heat Exchanger")
    help = ""
    kwargs = {
        "entradaTubo": None,
        "entradaExterior": None,

        "modo": 0,
        "flujo": 0,
        "orientacion": 0,
        "tubesideLaminar": 0,
        "tubesideTurbulent": 0,
        "metodo": 0,
        "phase": 0,

        "DeeTube": 0.0,
        "DeTube": 0.0,
        "DiTube": 0.0,
        "wTube": 0.0,
        "rTube": 0.0,
        "kTube": 0.0,
        "LTube": 0.0,
        "nTube": 0.0,

        "tubeFouling": 0.0,
        "annulliFouling": 0.0,

        "tubeFinned": 0,
        "hFin": 0.0,
        "thicknessBaseFin": 0.0,
        "thicknessTopFin": 0.0,
        "rootDoFin": 0.0,
        "kFin": 0.0,
        "nFin": 0,

        "tubeTout": 0.0,
        "tubeXout": -1.0,
        "annulliTout": 0.0,
        "annulliXout": -1.0,

        "f_install": 3.,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "material": 0,
        "P_dis": 0}

    kwargsInput = ("entradaTubo", "entradaExterior")
    kwargsValue = ("DeTube", "DiTube", "wTube", "rTube", "kTube", "LTube",
                   "nTube", "tubeFouling", "annulliFouling", "P_dis",
                   "tubeTout", "annulliTout")
    kwargsList = ("modo", "flujo", "orientacion")
    kwargsCheck = ("tubeFinned", )
    calculateValue = ("Q", "ToutAnnulli", "ToutTube", "U", "A", "L",
                      "deltaPTube", "deltaPAnnulli", "CF")
    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 2

    TEXT_MODO = [
        translate("equipment", "Design"),
        translate("equipment", "Rating")]
    TEXT_FLUJO = [
        translate("equipment", "Counterflow"),
        translate("equipment", "Parallelflow")]
    TEXT_ORIENTACION = [
        translate("equipment", "Horizontal"),
        translate("equipment", "Vertical, (in down)"),
        translate("equipment", "Vertical, (in up)")]
    TEXT_MATERIAL = [
        translate("equipment", "Carbon steel/carbon steel"),
        translate("equipment", "Carbon steel/304 stainless"),
        translate("equipment", "Carbon steel/316 stainless")]
    CODE_FLUJO = ("CF", "PF")

    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entradaTubo"]:
            self.msg = translate("equipment", "undefined internal stream input")
            self.status = 0
            return
        if not self.kwargs["entradaExterior"]:
            self.msg = translate("equipment", "undefined external stream input")
            self.status = 0
            return

        if not self.kwargs["DeeTube"]:
            self.msg = translate("equipment", "undefined pipe external diameter")
            self.status = 0
            return

        self.statusPipe = 0
        if self.kwargs["DeTube"] and self.kwargs["DiTube"]:
            self.statusPipe = 1
        elif self.kwargs["DeTube"] and self.kwargs["wTube"]:
            self.statusPipe = 2
        elif self.kwargs["DiTube"] and self.kwargs["wTube"]:
            self.statusPipe = 3
        else:
            self.msg = translate("equipment", "undefined pipe diameters")
            self.status = 0
            return

        if not self.kwargs["kTube"]:
            self.msg = translate(
                "equipment", "undefined pipe material thermal conductivity")
            self.status = 0
            return

        self.statusFinned = 0
        self.tubefinned = translate("equipment", "Bare Tube")
        if self.kwargs["tubeFinned"]:
            self.tubefinned = translate("equipment", "Finned Tube")
            if self.kwargs["hFin"] and (self.kwargs["thicknessBaseFin"] or
                                        self.kwargs["thicknessTopFin"]):
                self.statusFinned = 1
                self.msg = ""
            else:
                self.msg = translate("equipment", "fin not specified, using bare tube")
                self.status = 3

        if self.kwargs["modo"]:
            if not self.kwargs["LTube"]:
                self.msg = translate("equipment", "undefined pipe length")
                self.status = 0
                return
        else:
            self.statusOut = 0
            o1 = self.kwargs["tubeTout"] or self.kwargs["tubeXout"] != -1
            o2 = self.kwargs["annulliTout"] or self.kwargs["annulliXout"] != -1
            if o1 and o2:
                self.statusOut = 1
            elif o1:
                self.statusOut = 2
            elif o2:
                self.statusOut = 3
            else:
                self.msg = translate("equipment", "undefined output condition")
                self.status = 0
                return
        return True

    def calculo(self):
        # Define tipo de flujo
        if self.kwargs["flujo"]:
            self.flujo = "PF"
        else:
            self.flujo = "CF"

        # Calculate pipe dimension
        if self.statusPipe == 1:
            self.De = unidades.Length(self.kwargs["DeTube"])
            self.Di = unidades.Length(self.kwargs["DiTube"])
            self.w = unidades.Length((self.De-self.Di)/2)
            if self.kwargs["wTube"] and w != self.kwargs["wTube"]:
                self.msg = translate("equipment", "Pipe thickness discard")
                self.status = 3
        elif self.statusPipe == 2:
            self.De = unidades.Length(self.kwargs["DeTube"])
            self.w = unidades.Length(self.kwargs["wTube"])
            self.Di = unidades.Length(self.De-w*2)
        else:
            self.Di = unidades.Length(self.kwargs["DiTube"])
            self.w = unidades.Length(self.kwargs["wTube"])
            self.De = unidades.Length(self.Di+w*2)
        self.Dee = unidades.Length(self.kwargs["DeeTube"])
        self.Dei = unidades.Length(self.De+self.w)
        self.rugosidad = unidades.Length(self.kwargs["rTube"])
        self.k = unidades.ThermalConductivity(self.kwargs["kTube"])
        self.fi = unidades.Fouling(self.kwargs["tubeFouling"])
        self.fo = unidades.Fouling(self.kwargs["annulliFouling"])

        if self.kwargs["modo"]:
            self.rating()
        else:
            self.design()

        eD = unidades.Dimensionless(self.kwargs["rTube"]/self.Di)
        f = f_friccion(self.ReTube, eD)
        dp_tube = self.L*self.VTube**2/self.Di*f*self.rhoTube/2
        self.deltaPTube = unidades.DeltaP(dp_tube)

        f_a = f_friccion(self.ReAnnulli, geometry=6, Di=self.Dei, Do=self.Dee)
        dp_annulli = self.L*self.VAnnulli**2/self.De*f_a*self.rhoAnnulli/2
        self.deltaPAnnulli = unidades.DeltaP(dp_annulli)

        self.salida = [
            self.outTube.clone(P=self.outTube.P-self.deltaPTube),
            self.outAnnulli.clone(P=self.outAnnulli.P-self.deltaPAnnulli)]
        self.ToutTube = self.salida[0].T
        self.ToutAnnulli = self.salida[1].T
        self.XoutTube = self.salida[0].x
        self.XoutAnnulli = self.salida[1].x
        self.TinTube = self.kwargs["entradaTubo"].T
        self.XinTube = self.kwargs["entradaTubo"].x
        self.TinAnnulli = self.kwargs["entradaExterior"].T
        self.XinAnnulli = self.kwargs["entradaExterior"].x

    def rating(self):
        """Rating of a specified pipe"""
        # Input stream
        inTube = self.kwargs["entradaTubo"]
        inAnnulli = self.kwargs["entradaExterior"]

        self.L = unidades.Length(self.kwargs["LTube"])

        # Mean temperature method
        if self.kwargs["metodo"] == 0:
            self.phaseTube = self.ThermalPhase(inTube, inTube)
            self.phaseAnnulli = self.ThermalPhase(inAnnulli, inAnnulli)

            Ci = inTube.Liquido.cp*inTube.caudalmasico
            Co = inAnnulli.Liquido.cp*inAnnulli.caudalmasico
            Cmin = min(Ci, Co)
            Cmax = max(Ci, Co)
            C_ = Cmin/Cmax

            self.A = unidades.Area(self.L*pi*self.De)
            hi = self._hTube(inTube)
            ho = self._hAnnulli(inAnnulli)
            ni, no = self.rendimientoAletas(hi, ho)
            self.Ug(hi, ni, ho, no)

            NTU = self.A*self.U/Cmin
            ep = effectiveness(NTU, C_, self.CODE_FLUJO[self.kwargs["flujo"]])
            self.Q = unidades.Power(ep*Cmin*abs(inTube.T-inAnnulli.T))

            if inTube.T > inAnnulli.T:
                QTube = self.Q
                QAnnulli = -self.Q
            else:
                QTube = -self.Q
                QAnnulli = self.Q

            def f(T):
                return inTube.clone(T=T).h-inTube.h+QTube
            T = fsolve(f, inTube.T)[0]
            self.outTube = inTube.clone(T=T)

            def f(T):
                return inAnnulli.clone(T=T).h-inAnnulli.h+QAnnulli
            T = fsolve(f, inAnnulli.T)[0]
            self.outAnnulli = inAnnulli.clone(T=T)

    def design(self):
        """Design a pipe to meet the specified heat transfer requeriments"""
        # Input stream
        inTube = self.kwargs["entradaTubo"]
        inAnnulli = self.kwargs["entradaExterior"]

        # Metodo temperaturas medias
        if self.kwargs["metodo"] == 0:

            # Calculate output condition and sensible/latent thermal situation,
            # global thermal balance
            if self.statusOut == 1:
                if self.kwargs["tubeTout"]:
                    self.outTube = inTube.clone(T=self.kwargs["tubeTout"])
                else:
                    self.outTube = inTube.clone(x=self.kwargs["tubeXout"])
                if self.kwargs["annulliTout"]:
                    Tout = self.kwargs["annulliTout"]
                    self.outAnnulli = inAnnulli.clone(T=Tout)
                else:
                    Xout = self.kwargs["annulliXout"]
                    self.outAnnulli = inAnnulli.clone(x=Xout)

                Qo = abs(self.outAnnulli.h-inAnnulli.h)
                Qi = abs(self.outTube.h-inTube.h)
                self.Q = unidades.Power((Qo+Qi)/2.)

            elif self.statusOut == 2:
                if self.kwargs["tubeTout"]:
                    self.outTube = inTube.clone(T=self.kwargs["tubeTout"])
                else:
                    self.outTube = inTube.clone(x=self.kwargs["tubeXout"])

                Qi = abs(self.outTube.h-inTube.h)
                self.Q = unidades.Power(Qi)

                def f(T):
                    return inAnnulli.clone(T=T).h-inAnnulli.h-Qi
                T = fsolve(f, inAnnulli.T)[0]
                self.outAnnulli = inAnnulli.clone(T=T)

            elif self.statusOut == 3:
                if self.kwargs["annulliTout"]:
                    Tout = self.kwargs["annulliTout"]
                    self.outAnnulli = inAnnulli.clone(T=Tout)
                else:
                    Xout = self.kwargs["annulliXout"]
                    self.outAnnulli = inAnnulli.clone(x=Xout)

                Qo = abs(self.outAnnulli.h-inAnnulli.h)
                self.Q = unidades.Power(Qo)

                def f(T):
                    return inTube.clone(T=T).h-inTube.h-Qi
                T = fsolve(f, inTube.T)[0]
                self.outTube = inTube.clone(T=T)

            self.phaseTube = self.ThermalPhase(inTube, self.outTube)
            self.phaseAnnulli = self.ThermalPhase(inAnnulli, self.outAnnulli)

            fluidTube = inTube.clone(T=(inTube.T+self.outTube.T)/2.)
            T = (inAnnulli.T+self.outAnnulli.T)/2.
            fluidAnnulli = inAnnulli.clone(T=T)

            hi = self._hTube(fluidTube)
            ho = self._hAnnulli(fluidAnnulli)
            ni, no = self.rendimientoAletas(hi, ho)
            self.Ug(hi, ni, ho, no)

            if self.kwargs["flujo"]:
                DTin = abs(inAnnulli.T-inTube.T)
                DTout = abs(self.kwargs["tubeTout"]-self.kwargs["annulliTout"])
            else:
                DTin = abs(self.kwargs["tubeTout"]-inAnnulli.T)
                DTout = abs(self.kwargs["annulliTout"]-inTube.T)
            if DTin == DTout:
                DTm = DTin
            else:
                DTm = (DTin-DTout)/log(DTin/DTout)

        self.A = unidades.Area(self.Q/self.U/DTm)
        self.L = unidades.Length(self.A/2/pi)

    def Ug(self, hi, ni, ho, no):
        """Calculate global heat transfer coefficient"""
        Ui = self.De/self.Di/hi/ni
        Ufi = self.De*self.fi/self.Di/ni
        k = self.De*log(self.De/self.Di)/2/self.k
        U = 1/(Ui+Ufi+k+self.fo/no+1/ho/no)
        Uc = 1/(Ui+k+1/ho/no)
        self.hTube = unidades.HeatTransfCoef(hi)
        self.hAnnulli = unidades.HeatTransfCoef(ho)
        self.U = unidades.HeatTransfCoef(U)
        self.CF = unidades.Dimensionless(U/Uc)
        self.OS = unidades.Dimensionless(Uc*(self.fi+self.fo))

    def ThermalPhase(self, input, output):
        # Calculate thermal fundamentals
        if input.x == output.x:
            if input.x == 0:
                phase = "Latent-Liquid"
            else:
                phase = "Latent-Vapor"
        elif input.x < output.x:
            phase = "Evaporator"
        else:
            phase = "Condenser"
        return phase

    def rendimientoAletas(self, hi, ho):
        """Calculate thermal efficiency of fins"""
        self.hFin = unidades.Length(self.kwargs["hFin"])
        self.thicknessBaseFin = unidades.Length(self.kwargs["thicknessBaseFin"])
        self.thicknessTopFin = unidades.Length(self.kwargs["thicknessTopFin"])
        self.rootDoFin = unidades.Length(self.kwargs["rootDoFin"])
        self.kFin = unidades.ThermalConductivity(self.kwargs["kFin"])
        self.nFin = unidades.Dimensionless(self.kwargs["nFin"])
        if self.kwargs["tubeFinned"]:
            if self.statusFinned:
                # For now only use circular fin
                do = self.kwargs["rootDoFin"]
                D = do + self.kwargs["hFin"] * 2
                phi = (D/do-1)*(1+0.35*log(D/do))
                w_b = self.kwargs["thicknessBaseFin"]
                w_t = self.kwargs["thicknessTopFin"]
                if w_b and w_t:
                    delta = (w_b + w_t) / 2
                else:
                    delta = w_b + w_t
                if self.kwargs["kFin"]:
                    kf = self.kwargs["kFin"]
                else:
                    kf = self.kwargs["kTube"]

                X = phi*do/2*sqrt(2*ho/kf/delta)
                no = tanh(X)/X
            else:
                no = 1
            # TODO: Define internal fins
            ni = 1
        else:
            ni = 1
            no = 1
        return ni, no

    def _hTube(self, fluidTube):
        """Calculate convection heat trasnfer coefficient in tubeside"""
        if self.phaseTube[:6] == "Latent":
            if fluidTube.x == 0:
                fluido = fluidTube.Liquido
            else:
                fluido = fluidTube.Vapor

            rho_i = fluido.rho
            mu = fluido.mu
            k = fluido.k
            v_i = fluidTube.Q*4/pi/self.Di**2
            re_i = Re(D=self.Di, V=v_i, rho=rho_i, mu=mu)
            self.VTube = unidades.Speed(v_i)
            self.rhoTube = rho_i
            self.ReTube = unidades.Dimensionless(re_i)
            pr = fluido.Prandt

            if re_i < 2300:
                L = self.L
                cp = fluido.cp
                w = fluido.caudalmasico
                gz = Gz(w=w, cp=cp, k=k, L=L)
                beta = fluido.alfav
                gr = Gr(beta=beta, T=fluidTube.T, To=fluidTube.T, L=L, mu=mu)
                if self.kwargs["tubesideLaminar"] == 0:
                    Nu = h_tubeside_laminar_Eubank_Proctor(
                        Pr=pr, Gz=gz, Gr=gr, D=self.Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 1:
                    Nu = h_tubeside_laminar_VDI(Re=re_i, Pr=pr, D=self.Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 2:
                    Nu = h_tubeside_laminar_Hausen(Gz=gz)
                elif self.kwargs["tubesideLaminar"] == 3:
                    Nu = h_tubeside_laminar_Sieder_Tate(Gz=gz, Gr=gr)
            else:
                if self.kwargs["tubesideTurbulent"] == 0:
                    Nu = h_tubeside_turbulent_Sieder_Tate(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 1:
                    Nu = h_tubeside_turbulent_Colburn(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 2:
                    frio = self.kwargs["entradaCarcasa"].T > fluidTube.T
                    Nu = h_tubeside_turbulent_Dittus_Boelter(
                        Re=re_i, Pr=pr, calentamiento=frio)
                elif self.kwargs["tubesideTurbulent"] == 3:
                    Nu = h_tubeside_turbulent_ESDU(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 4:
                    Nu = h_tubeside_turbulent_Gnielinski(
                        Re=re_i, Pr=pr, D=self.Di, L=L)
                elif self.kwargs["tubesideTurbulent"] == 5:
                    line = self.kwargs["distribucionTube"] == 3
                    filas = self.kwargs["NTube"]**0.5
                    Nu = h_tubeside_turbulent_VDI(
                        Re=re_i, Pr=pr, filas_tubos=filas, alineados=line)

        if self.phaseTube == "Condenser":
            if self.kwargs["orientation"] == 0:
                # Condensation in horizontal tubes
                if 0 < fluidTube.x < 1:
                    X_lockhart = ((1-fluidTube.x)/fluidTube.x)**0.9 * \
                        (fluidTube.Vapor.rho/fluidTube.Liquido.rho)**0.5 * \
                        (fluidTube.Liquido.mu/fluidTube.Vapor.mu)**0.1
                    G = fluidTube.caudalmasico*4/pi/self.Di**2
                    j = fluidTube.x*G/(
                        g*self.Di*fluidTube.Vapor.rho *
                        (fluidTube.Liquido.rho-fluidTube.Vapor.rho))**0.5
                    print((j, X_lockhart))

            else:
                # Condensation in vertical tubes
                pass

        return unidades.HeatTransfCoef(Nu*k/self.Di)

    def _hAnnulli(self, fluidAnnulli):
        """Calculate convection heat transfer coefficient in annulliside"""
        a = self.Dee/self.De
        dh = self.Dee-self.De

        rho = fluidAnnulli.Liquido.rho
        mu = fluidAnnulli.Liquido.mu
        k = fluidAnnulli.Liquido.k
        v = fluidAnnulli.Q*4/pi/(self.Dee**2-self.De**2)
        re = Re(D=dh, V=v, rho=rho, mu=mu)
        self.VAnnulli = unidades.Speed(v)
        self.rhoAnnulli = rho
        self.ReAnnulli = unidades.Dimensionless(re)
        pr = fluidAnnulli.Liquido.Prandt

        if re <= 2300:
            Nu = h_anulli_Laminar(re, pr, a)
        elif re >= 1e4:
            Nu = h_anulli_Turbulent(re, pr, a)
        else:
            Nu = h_anulli_Transition(re, pr, a)

        return unidades.HeatTransfCoef(Nu*k/self.Di)

    def coste(self):
        self.material = self.kwargs["material"]
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]

        if self.kwargs["P_dis"]:
            self.P_dis = unidades.Pressure(self.kwargs["P_dis"])
        else:
            self.P_dis = max(self.kwargs["entradaTubo"].P,
                             self.kwargs["entradaExterior"].P)

        Pd = self.P_dis.psi
        Fm = [1., 1.9, 2.2][self.kwargs["material"]]

        if Pd < 4:
            Fp = 1.
        elif Pd < 6:
            Fp = 1.1
        else:
            Fp = 1.25

        C = Fm*Fp*900*self.A.ft2**0.18
        self.C_adq = unidades.Currency(C * CI / BI)
        self.C_inst = unidades.Currency(self.C_adq * self.kwargs["f_install"])

    def propTxt(self):
        txt = "#---------------"
        txt += translate("equipment", "Catalog")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(11))

        if self.kwargs["tubeFinned"]:
            txt += "\t" + self.propertiesToText(range(11, 17))

        txt += os.linesep + "#---------------"
        txt += translate("equipment", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(17, 20)) + os.linesep
        txt += self.propertiesToText(range(20, 29)) + os.linesep  # Tube
        txt += self.propertiesToText(range(29, 38)) + os.linesep  # Annulli
        txt += self.propertiesToText(range(38, 42)) + os.linesep

        if self.statusCoste:
            txt += os.linesep+"#---------------"
            txt += translate("equipment", "Preliminary Cost Estimation")
            txt += "-----------------#" + os.linesep
            txt += self.propertiesToText(range(42, 49))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(translate("equipment", "Length"), "L", unidades.Length),
             (translate("equipment", "Pipe Internal Diameter"), "Di", unidades.Length),
             (translate("equipment", "Pipe External Diameter"), "De", unidades.Length),
             (translate("equipment", "Annulli External Diameter"),
              "Dee", unidades.Length),
             (translate("equipment", "Thickness"), "w", unidades.Length),
             (translate("equipment", "Roughness"), "rugosidad", unidades.Length),
             (translate("equipment", "External Area"), "A", unidades.Area),
             (translate("equipment", "Thermal Conductivity"), "k",
              unidades.ThermalConductivity),
             (translate("equipment", "Internal Fouling"), "fi", unidades.Fouling),
             (translate("equipment", "External Fouling"), "fo", unidades.Fouling),
             (translate("equipment", "Finned Tube"), "tubefinned", str),
             (translate("equipment", "Fin height"), "hFin", unidades.Length),
             (translate("equipment", "Thickness at bottom of fin"),
              "thicknessBaseFin", unidades.Length),
             (translate("equipment", "Thickness at top of fin"),
              "thicknessTopFin", unidades.Length),
             (translate("equipment", "External diameter at bottom of fin"),
              "rootDoFin", unidades.Length),
             (translate("equipment", "Fin thermal conductivity"),
              "kFin", unidades.ThermalConductivity),
             (translate("equipment", "Fin count per meter"),
              "nFin", unidades.Dimensionless),
             (translate("equipment", "Mode"), ("TEXT_MODO", "modo"), str),
             (translate("equipment", "Arrangement Flow"),
              ("TEXT_FLUJO", "flujo"), str),
             (translate("equipment", "Layout"),
              ("TEXT_ORIENTACION", "orientacion"), str),
             (translate("equipment", "Tube Mechanism"), "phaseTube", str),
             (translate("equipment", "Tube Fluid Speed"), "VTube", unidades.Speed),
             (translate("equipment", "Tube Reynolds"), "ReTube",
              unidades.Dimensionless),
             (translate("equipment", "Tube In Temperature"),
              "TinTube", unidades.Temperature),
             (translate("equipment", "Tube In Quality"),
              "XinTube", unidades.Dimensionless),
             (translate("equipment", "Tube Out Temperature"),
              "ToutTube", unidades.Temperature),
             (translate("equipment", "Tube Out Quality"),
              "XoutTube", unidades.Dimensionless),
             (translate("equipment", "ΔP Tube"), "deltaPTube", unidades.DeltaP),
             (translate("equipment", "Tube heat transfer"),
              "hTube", unidades.HeatTransfCoef),
             (translate("equipment", "Annulli Mechanism"), "phaseAnnulli", str),
             (translate("equipment", "Annulli Fluid Speed"),
              "VAnnulli", unidades.Speed),
             (translate("equipment", "Annulli Reynolds"),
              "ReAnnulli", unidades.Dimensionless),
             (translate("equipment", "Annulli In Temperature"),
              "TinAnnulli", unidades.Temperature),
             (translate("equipment", "Annulli In Quality"),
              "XinAnnulli", unidades.Dimensionless),
             (translate("equipment", "Annulli Out Temperature"),
              "ToutAnnulli", unidades.Temperature),
             (translate("equipment", "Annulli Out Quality"),
              "XoutAnnulli", unidades.Dimensionless),
             (translate("equipment", "ΔP Annulli", None),
              "deltaPAnnulli", unidades.DeltaP),
             (translate("equipment", "Annulli heat transfer"),
              "hAnnulli", unidades.HeatTransfCoef),
             (translate("equipment", "U"), "U", unidades.HeatTransfCoef),
             (translate("equipment", "Q"), "Q", unidades.Power),
             (translate("equipment", "Clean Factor"), "CF", unidades.Dimensionless),
             (translate("equipment", "Over Surface"), "OS", unidades.Dimensionless),
             (translate("equipment", "Base index"), "Base_index", float),
             (translate("equipment", "Current index"), "Current_index", float),
             (translate("equipment", "Install factor"), "f_install", float),
             (translate("equipment", "Material"), ("TEXT_MATERIAL", "material"), str),
             (translate("equipment", "Design Pressure"), "P_dis", unidades.Pressure),
             (translate("equipment", "Purchase Cost"), "C_adq", unidades.Currency),
             (translate("equipment", "Installed Cost"), "C_inst", unidades.Currency)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["L"] = self.L
        state["Di"] = self.Di
        state["De"] = self.De
        state["Dee"] = self.Dee
        state["w"] = self.w
        state["rugosidad"] = self.rugosidad
        state["A"] = self.A
        state["k"] = self.k
        state["fi"] = self.fi
        state["fo"] = self.fo
        state["tubefinned"] = self.tubefinned

        state["hFin"] = self.hFin
        state["thicknessBaseFin"] = self.thicknessBaseFin
        state["thicknessTopFin"] = self.thicknessTopFin
        state["rootDoFin"] = self.rootDoFin
        state["kFin"] = self.kFin
        state["nFin"] = self.nFin

        state["phaseTube"] = self.phaseTube
        state["VTube"] = self.VTube
        state["ReTube"] = self.ReTube
        state["TinTube"] = self.TinTube
        state["XinTube"] = self.XinTube
        state["ToutTube"] = self.ToutTube
        state["XoutTube"] = self.XoutTube
        state["deltaPTube"] = self.deltaPTube
        state["hTube"] = self.hTube
        state["phaseAnnulli"] = self.phaseAnnulli
        state["VAnnulli"] = self.VAnnulli
        state["ReAnnulli"] = self.ReAnnulli
        state["TinAnnulli"] = self.TinAnnulli
        state["XinAnnulli"] = self.XinAnnulli
        state["ToutAnnulli"] = self.ToutAnnulli
        state["XoutAnnulli"] = self.XoutAnnulli
        state["deltaPAnnulli"] = self.deltaPAnnulli
        state["hAnnulli"] = self.hAnnulli
        state["U"] = self.U
        state["Q"] = self.Q
        state["CF"] = self.CF
        state["OS"] = self.OS

        state["statusCoste"] = self.statusCoste
        if self.statusCoste:
            state["P_dis"] = self.P_dis
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.L = unidades.Length(state["L"])
        self.Di = unidades.Length(state["Di"])
        self.De = unidades.Length(state["De"])
        self.Dee = unidades.Length(state["Dee"])
        self.w = unidades.Length(state["w"])
        self.rugosidad = unidades.Length(state["rugosidad"])
        self.A = unidades.Area(state["A"])
        self.k = unidades.ThermalConductivity(state["k"])
        self.fi = unidades.Fouling(state["fi"])
        self.fo = unidades.Fouling(state["fo"])
        self.tubefinned = state["tubefinned"]

        self.hFin = unidades.Length(state["hFin"])
        self.thicknessBaseFin = unidades.Length(state["thicknessBaseFin"])
        self.thicknessTopFin = unidades.Length(state["thicknessTopFin"])
        self.rootDoFin = unidades.Length(state["rootDoFin"])
        self.kFin = unidades.ThermalConductivity(state["kFin"])
        self.nFin = unidades.Dimensionless(state["nFin"])

        self.phaseTube = state["phaseTube"]
        self.VTube = unidades.Speed(state["VTube"])
        self.ReTube = unidades.Dimensionless(state["ReTube"])
        self.TinTube = unidades.Temperature(state["TinTube"])
        self.XinTube = unidades.Dimensionless(state["XinTube"])
        self.ToutTube = unidades.Temperature(state["ToutTube"])
        self.XoutTube = unidades.Dimensionless(state["XoutTube"])
        self.deltaPTube = unidades.DeltaP(state["deltaPTube"])
        self.hTube = unidades.HeatTransfCoef(state["hTube"])
        self.phaseAnnulli = state["phaseAnnulli"]
        self.VAnnulli = unidades.Speed(state["VAnnulli"])
        self.ReAnnulli = unidades.Dimensionless(state["ReAnnulli"])
        self.TinAnnulli = unidades.Temperature(state["TinAnnulli"])
        self.XinAnnulli = unidades.Dimensionless(state["XinAnnulli"])
        self.ToutAnnulli = unidades.Temperature(state["ToutAnnulli"])
        self.XoutAnnulli = unidades.Dimensionless(state["XoutAnnulli"])
        self.deltaPAnnulli = unidades.DeltaP(state["deltaPAnnulli"])
        self.hAnnulli = unidades.HeatTransfCoef(state["hAnnulli"])
        self.U = unidades.HeatTransfCoef(state["U"])
        self.Q = unidades.Power(state["Q"])
        self.CF = unidades.Dimensionless(state["CF"])
        self.OS = unidades.Dimensionless(state["OS"])

        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.P_dis = unidades.Pressure(state["P_dis"])
            self.C_adq = unidades.Currency(state["C_adq"])
            self.C_inst = unidades.Currency(state["C_inst"])
        self.salida = [None]


class Air_Cooler(equipment):
    """Clase que modela los enfriadores de aire"""
    title = "Enfriador de aire"
    help = ""

    def coste(self, *args):
        self._indicesCoste(*args)

        C = 24.6*self.area.ft2**0.4*1000

        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)


class Evaporator(equipment):
    """Clase que modela los evaporadores"""
    title = "Evaporador"
    help = ""

    def coste(self, *args, **kwargs):
        """
        tipo:
            0   -   Forced circulation
            1   -   Long tube
            2   -   Falling film
        material:
        forced circulation materials shell/tube
            0   -  Steel/copper
            1   -  Monel/cupronickel
            2   -  Nickel/nickel
        Long tube evaporator materials shell/tube
            0   -  Steel/copper
            1   -  Steel/steel
            2   -  Steel/aluminum
            3   -  Nickel/nickel
        """
        self._indicesCoste(*args)

        self.tipo = kwargs.get("tipo", 0)
        self.material = kwargs.get("material", 0)

        A = self.area.ft2

        if self.tipo == 0:
            C_base = exp(5.9785-0.6056*log(A)+0.08514*log(A)**2)*1000
            if self.material == 0:
                Fm = 1.
            elif self.material == 1:
                Fm = 1.35
            else:
                Fm = 1.8
        elif self.tipo == 1:
            C_base = 0.36*A**0.85*1000
            if self.material == 0:
                Fm = 1.
            elif self.material == 1:
                Fm = 0.6
            elif self.material == 2:
                Fm = 0.7
            else:
                Fm = 3.3
        else:
            C_base = exp(3.2362-0.0126*log(A)+0.0244*log(A)**2)*1000
            Fm = 1

        C = Fm*C_base
        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)


class Refrigeration(equipment):
    """Clase que modela las unidades de refrigeración"""
    title = "Refrigerador"
    help = ""

    def Coste(self, *args, **kwargs):
        self._indicesCoste(*args)

        Tmax = self.Tmax.C
        Q = self.calor.MBtuh

        if Tmax >= 0:
            F = 1.
        elif Tmax >= -10:
            F = 1.55
        elif Tmax >= -20:
            F = 2.1
        elif Tmax >= -30:
            F = 2.65
        elif Tmax >= -3:
            F = 3.2
        else:
            F = 4.

        C = 146*F*Q**0.65*1000
        C_adq = C*self.Current_index/self.Base_index
        C_inst = C_adq*self.kwargs["f_install"]

        self.C_adq = C_adq
        self.C_inst = C_inst


def unsteady():
#========================================================================
#Input constants
#========================================================================
    dtau = 0.001 # Set dimensionless time increments
    dx = 0.05 # Set dimensionless length increments
    Tmax = 0.95 # Set maximum dimensionless temperature
    M = 21 # Counter for length discretization
#========================================================================
#Calculate parameters
#========================================================================
    dx = 1.0/(M-1)
    dx_x = 1.0/(M-1)
    ratio = dtau/(dx**2)
    const = 1.0 - 2.0*ratio
#========================================================================
#Set counters to zero
#========================================================================
    i = 0
    tau = 0.0
#========================================================================
# Set up arrays for solution and print
#========================================================================
    Tnew = zeros(M, dtype = float)
    T = zeros(M, dtype = float)
    T[0] = 1.0
    T[-1] = 1.0
    print(("T initial = ", T))
#========================================================================
# I just pick 400 on trial and error for the total array
#========================================================================
    T_sol = zeros((400,M), dtype = float)
    T_sol[:,0] = 1.0
    T_sol[:,-1] = 1.0
#========================================================================
# While loop to iterate until mid-point temperature reaches Tmax
#========================================================================
    while T[10] < Tmax:
        i = i + 1
        tau = tau + dtau
#========================================================================
# Calculate new tempertures
#========================================================================
        for j in range(1,M-1):
            Tnew[j] = ratio*(T[j-1] + T[j+1]) + const*T[j]
#========================================================================
# Substitute new Temperatures in array for T
#========================================================================
        for k in range(1,M-1):
            T[k] = Tnew[k]
            T_sol[i,k] = T[k]

    print(("Tau and T_final =", tau, T_sol[i,:]))
#========================================================================
# Set up array for spatial values of x to plot
#========================================================================
    x = [i*dx_x for i in range(M-1)]
    x.append(1.0)
#========================================================================
# Plot the solutions
#=======================================================================
    plot(x,T_sol[50,:])
    plot(x,T_sol[100,:])
    plot(x,T_sol[150,:])
    plot(x,T_sol[250,:])
    plot(x,T_sol[i,:])
#legend(['Tau = 0.5','Tau = 0.1','Tau = 0.15','Tau = 0.25',
#'Tau = final time'])
    title('Normalized Slab Temperatures')
    plt.show()
    grid()

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    # from lib.corriente import Corriente
    # aguaTubo=Corriente(T=283.15, P=101325., ids=[62], caudalMasico=10., fraccionMolar=[1.])
    # aguaCarcasa=Corriente(T=370, P=101325., ids=[62], caudalMasico=1000., fraccionMolar=[1.])
    # Ds=unidades.Length(19.25, "inch")
    # Do=unidades.Length(1, "inch")
    # L=unidades.Length(14, "ft")
    # pt=unidades.Length(1.25, "inch")
    # B=unidades.Length(3.85, "inch")
    # dotl=unidades.Length(0.5*(19.25-17.91), "inch")
    # cambiador = Shell_Tube(entradaTubo=aguaTubo, entradaCarcasa=aguaCarcasa, shellsideSensible=1, DShell=Ds, NTube=124, DeTube=Do, wTube=0.006, LTube=L, pitch=pt, distribucionTube=3, baffleSpacing=B, baffleSpacingIn=B, baffleSpacingOut=B, baffleCut=0.2, clearanceTubeBaffle=0.0004, clearanceShellBaffle=0.0025279, clearanceShellBundle=dotl, sealingStrips=0.1*9.24)


#    from scipy import *
#    from pylab import *
#    import matplotlib.pyplot as plt
#    unsteady()
