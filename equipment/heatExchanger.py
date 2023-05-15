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
# library for heat exchanger equipment calculation
# - Heat_Exchanger
# - Fired_Heater
# - Shell_Tube
# - Hairpin
###############################################################################


from math import sqrt, exp, log
import os

from tools.qt import tr
from scipy import pi, arccos, sin, cos, tanh
from scipy.constants import g
from scipy.optimize import fsolve

from lib import unidades
from lib.adimensional import Re, Pr, Gr, Gz
from lib.friction import f_friccion
from lib.heatTransfer import *  # noqa
from equipment.parents import equipment


# Equipment
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
    title = tr("pychemqt", "Heat Exchanger")
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
            self.msg = tr("pychemqt", "undefined input")
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
            self.msg = tr(
                "pychemqt", "undefined output temperature specification")
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
        txt += tr("pychemqt", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(5))
        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(tr("pychemqt", "Input Temperature"),
              "Tin", unidades.Temperature),
             (tr("pychemqt", "Output Temperature"),
              "ToutCalc", unidades.Temperature),
             (tr("pychemqt", "Temperature increase"),
              "deltaT", unidades.DeltaT),
             (tr("pychemqt", "Pressure increase"),
              "deltaP", unidades.DeltaP),
             (tr("pychemqt", "Heat"), "HeatCalc",
              unidades.Power)]
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
    title = tr("pychemqt", "Fired Heater")
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

    TEXT_TIPO = [tr("pychemqt", "Box Type"),
                 tr("pychemqt", "Cylindrical type")]
    TEXT_SUBTIPOBOX = [
        tr("pychemqt", "Process heater"),
        tr("pychemqt", "Pyrolysis"),
        tr("pychemqt", "Reforme without catalysis")]
    TEXT_SUBTIPOCYLINDRICAL = [
        tr("pychemqt", "Cylindrical"),
        tr("pychemqt", "Dowtherm")]
    TEXT_MATERIAL = [tr("pychemqt", "Carbon steel"),
                     tr("pychemqt", "Cr-Mo steel"),
                     tr("pychemqt", "Stainless steel")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = tr("pychemqt", "undefined input")
            self.status = 0
            return

        if not self.kwargs["Tout"]:
            self.msg = tr(
                "pychemqt", "undefined output temperature condition")
            self.status = 0
            return
        if self.kwargs["Tout"] <= self.kwargs["entrada"].T:
            self.msg = tr(
                "pychemqt", "bad output temperature condition")
            self.status = 0
            return

        if not self.kwargs["eficiencia"]:
            self.msg = tr(
                "pychemqt", "using default efficiency")
            self.status = 3
            return True
        if not self.kwargs["poderCalorifico"]:
            self.msg = tr(
                "pychemqt", "using default fuel calorific value")
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
        txt += tr("pychemqt", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(5))

        stimated = " (%s)" % tr("pychemqt", "stimated")
        txt += self.propertiesToText(5, kwCheck=True, kwSuffix=stimated,
                                     kwKey="eficiencia", kwValue=0.0)
        txt += self.propertiesToText(6, kwCheck=True, kwSuffix=stimated,
                                     kwKey="poderCalorifico", kwValue=0.0)
        txt += self.propertiesToText(7)

        if self.statusCoste:
            txt += os.linesep+"#---------------"
            txt += tr(
                "pychemqt", "Preliminary Cost Estimation")
            txt += "-----------------#" + os.linesep
            txt += self.propertiesToText(range(8, 18))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(tr("pychemqt", "Input Temperature"),
              "Tin", unidades.Temperature),
             (tr("pychemqt", "Output Temperature"),
              "Tout", unidades.Temperature),
             (tr("pychemqt", "Temperature increase"),
              "deltaT", unidades.DeltaT),
             (tr("pychemqt", "Pressure increase"),
              "deltaP", unidades.DeltaP),
             (tr("pychemqt", "Maximum heat"),
              "Hmax", unidades.Power),
             (tr("pychemqt", "Thermal Efficiency"),
              "eficiencia", unidades.Dimensionless),
             (tr("pychemqt", "Fuel Heating Value"),
              "poderCalorifico", unidades.Dimensionless),
             (tr("pychemqt", "Required Fuel"),
              "CombustibleRequerido", unidades.VolFlow),
             (tr("pychemqt", "Base index"),
              "Base_index", float),
             (tr("pychemqt", "Current index"),
              "Current_index", float),
             (tr("pychemqt", "Install factor"),
              "f_install", float),
             (tr("pychemqt", "FireHeater type"),
              ("TEXT_TIPO", "tipo"), str),
             (tr("pychemqt", "Cylindrical type"),
              ("TEXT_SUBTIPOCYLINDRICAL", "subtipoCylindrical"), str),
             (tr("pychemqt", "Box type"),
              ("TEXT_SUBTIPOBOX", "subtipoBox"), str),
             (tr("pychemqt", "Material"),
              ("TEXT_MATERIAL", "material"), str),
             (tr("pychemqt", "Design Pressure"),
              "P_dis", unidades.Pressure),
             (tr("pychemqt", "Purchase Cost"),
              "C_adq", unidades.Currency),
             (tr("pychemqt", "Installed Cost"),
              "C_inst", unidades.Currency)]
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
    title = tr("pychemqt", "Hairpin Heat Exchanger")
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
        tr("pychemqt", "Design"),
        tr("pychemqt", "Rating")]
    TEXT_FLUJO = [
        tr("pychemqt", "Counterflow"),
        tr("pychemqt", "Parallelflow")]
    TEXT_ORIENTACION = [
        tr("pychemqt", "Horizontal"),
        tr("pychemqt", "Vertical, (in down)"),
        tr("pychemqt", "Vertical, (in up)")]
    TEXT_MATERIAL = [
        tr("pychemqt", "Carbon steel/carbon steel"),
        tr("pychemqt", "Carbon steel/304 stainless"),
        tr("pychemqt", "Carbon steel/316 stainless")]
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
            self.msg = tr(
                "pychemqt", "undefined internal stream input")
            self.status = 0
            return
        if not self.kwargs["entradaExterior"]:
            self.msg = tr(
                "pychemqt", "undefined external stream input")
            self.status = 0
            return

        if not self.kwargs["DeeTube"]:
            self.msg = tr(
                "pychemqt", "undefined pipe external diameter")
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
            self.msg = tr(
                "pychemqt", "undefined pipe diameters")
            self.status = 0
            return

        if not self.kwargs["kTube"]:
            self.msg = tr(
                "pychemqt", "undefined pipe material thermal conductivity")
            self.status = 0
            return

        self.statusFinned = 0
        self.tubefinned = tr("pychemqt", "Bare Tube")
        if self.kwargs["tubeFinned"]:
            self.tubefinned = tr("pychemqt", "Finned Tube")
            if self.kwargs["hFin"] and (self.kwargs["thicknessBaseFin"] or
                                        self.kwargs["thicknessTopFin"]):
                self.statusFinned = 1
                self.msg = ""
            else:
                self.msg = tr(
                    "pychemqt", "fin not specified, using bare tube")
                self.status = 3

        if self.kwargs["modo"]:
            if not self.kwargs["LTube"]:
                self.msg = tr(
                    "pychemqt", "undefined pipe length")
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
                self.msg = tr(
                    "pychemqt", "undefined output condition")
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
                self.msg = tr(
                    "pychemqt", "Pipe thickness discard")
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

        f_a = f_friccion(self.ReAnnulli, geometry=6)
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
            ep = efectividad(NTU, C_, self.CODE_FLUJO[self.kwargs["flujo"]])
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
        txt += tr("pychemqt", "Catalog")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(11))

        if self.kwargs["tubeFinned"]:
            txt += "\t" + self.propertiesToText(range(11, 17))

        txt += os.linesep + "#---------------"
        txt += tr("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(17, 20)) + os.linesep
        txt += self.propertiesToText(range(20, 29)) + os.linesep  # Tube
        txt += self.propertiesToText(range(29, 38)) + os.linesep  # Annulli
        txt += self.propertiesToText(range(38, 41)) + os.linesep

        if self.statusCoste:
            txt += os.linesep+"#---------------"
            txt += tr(
                "pychemqt", "Preliminary Cost Estimation")
            txt += "-----------------#" + os.linesep
            txt += self.propertiesToText(range(41, 48))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(tr("pychemqt", "Length"), "L",
              unidades.Length),
             (tr("pychemqt", "Pipe Internal Diameter"),
              "Di", unidades.Length),
             (tr("pychemqt", "Pipe External Diameter"),
              "De", unidades.Length),
             (tr("pychemqt", "Annulli External Diameter"),
              "Dee", unidades.Length),
             (tr("pychemqt", "Thickness"), "w",
              unidades.Length),
             (tr("pychemqt", "Roughness"), "rugosidad",
              unidades.Length),
             (tr("pychemqt", "External Area"), "A",
              unidades.Area),
             (tr("pychemqt", "Thermal Conductivity"), "k",
              unidades.ThermalConductivity),
             (tr("pychemqt", "Internal Fouling"), "fi",
              unidades.Fouling),
             (tr("pychemqt", "External Fouling"), "fo",
              unidades.Fouling),
             (tr("pychemqt", "Finned Tube"), "tubefinned",
              str),
             (tr("pychemqt", "Fin height"), "hFin",
              unidades.Length),
             (tr("pychemqt", "Thickness at bottom of fin"),
              "thicknessBaseFin", unidades.Length),
             (tr("pychemqt", "Thickness at top of fin"),
              "thicknessTopFin", unidades.Length),
             (tr("pychemqt",
                                     "External diameter at bottom of fin"),
              "rootDoFin", unidades.Length),
             (tr("pychemqt", "Fin thermal conductivity"),
              "kFin", unidades.ThermalConductivity),
             (tr("pychemqt", "Fin count per meter"),
              "nFin", unidades.Dimensionless),
             (tr("pychemqt", "Mode"),
              ("TEXT_MODO", "modo"), str),
             (tr("pychemqt", "Arrangement Flow"),
              ("TEXT_FLUJO", "flujo"), str),
             (tr("pychemqt", "Layout"),
              ("TEXT_ORIENTACION", "orientacion"), str),
             (tr("pychemqt", "Tube Mechanism"),
              "phaseTube", str),
             (tr("pychemqt", "Tube Fluid Speed"), "VTube",
              unidades.Speed),
             (tr("pychemqt", "Tube Reynolds"), "ReTube",
              unidades.Dimensionless),
             (tr("pychemqt", "Tube In Temperature"),
              "TinTube", unidades.Temperature),
             (tr("pychemqt", "Tube In Quality"),
              "XinTube", unidades.Dimensionless),
             (tr("pychemqt", "Tube Out Temperature"),
              "ToutTube", unidades.Temperature),
             (tr("pychemqt", "Tube Out Quality"),
              "XoutTube", unidades.Dimensionless),
             (tr("pychemqt", "ΔP Tube"),
              "deltaPTube", unidades.DeltaP),
             (tr("pychemqt", "Tube heat transfer"),
              "hTube", unidades.HeatTransfCoef),
             (tr("pychemqt", "Annulli Mechanism"),
              "phaseAnnulli", str),
             (tr("pychemqt", "Annulli Fluid Speed"),
              "VAnnulli", unidades.Speed),
             (tr("pychemqt", "Annulli Reynolds"),
              "ReAnnulli", unidades.Dimensionless),
             (tr("pychemqt", "Annulli In Temperature"),
              "TinAnnulli", unidades.Temperature),
             (tr("pychemqt", "Annulli In Quality"),
              "XinAnnulli", unidades.Dimensionless),
             (tr("pychemqt", "Annulli Out Temperature"),
              "ToutAnnulli", unidades.Temperature),
             (tr("pychemqt", "Annulli Out Quality"),
              "XoutAnnulli", unidades.Dimensionless),
             (tr("pychemqt", "ΔP Annulli", None),
              "deltaPAnnulli", unidades.DeltaP),
             (tr("pychemqt", "Annulli heat transfer"),
              "hAnnulli", unidades.HeatTransfCoef),
             (tr("pychemqt", "U"), "U",
              unidades.HeatTransfCoef),
             (tr("pychemqt", "Clean Factor"), "CF",
              unidades.Dimensionless),
             (tr("pychemqt", "Over Surface"), "OS",
              unidades.Dimensionless),
             (tr("pychemqt", "Base index"), "Base_index",
              float),
             (tr("pychemqt", "Current index"),
              "Current_index", float),
             (tr("pychemqt", "Install factor"),
              "f_install", float),
             (tr("pychemqt", "Material"),
              ("TEXT_MATERIAL", "material"), str),
             (tr("pychemqt", "Design Pressure"), "P_dis",
              unidades.Pressure),
             (tr("pychemqt", "Purchase Cost"), "C_adq",
              unidades.Currency),
             (tr("pychemqt", "Installed Cost"), "C_inst",
              unidades.Currency)]
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
        self.CF = unidades.Dimensionless(state["CF"])
        self.OS = unidades.Dimensionless(state["OS"])

        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.P_dis = unidades.Pressure(state["P_dis"])
            self.C_adq = unidades.Currency(state["C_adq"])
            self.C_inst = unidades.Currency(state["C_inst"])
        self.salida = [None]


class Shell_Tube(equipment):
    """Clase que define un cambiador de calor de carcasa y tubos

    Parámetros:
        entrada: Array con dos Instancia de clase corriente que define las corrientes que fluye por el equipo, en el orden [tubo, carcasa]
        entradaTubo: Instancia de clase corriente que define la corriente que pasa por los tubos
        entradaCarcasa: Instancia de calse corriente que define la corriente que pasa por la carcasa

    Standard:
        class_: Clase del standard TEMA:
            0   -   Clase R
            1   -   Clase B
            2   -   Clase C
        frontHead: Tipo de cabezal inicial
            0   -   Channel & Removable Cover
            1   -   Bonnet
            2   -   Removable Bundle
            3   -   Special High Pressure Closure
            4   -   Channel with Tubesheet & Removable Cover
        shell: Tipo de carcasa
            0   -   One Pass
            1   -   Two Pass
            2   -   Split Flow
            3   -   Double Split Flow
            4   -   Divided Flow
            5   -   Kettle Reboiler
            6   -   Cross Flow
        rearHead: Tipo de cabezal final
            0   -   Fixed Tubesheet (A head)
            1   -   Fixed Tubesheet (B head)
            2   -   Fixed Tubesheet (N head)
            3   -   Outside Packed Flt Head
            4   -   Flt Head with Backing Dev
            5   -   Pull Throught Flt Heat
            6   -   U-Tube Bundle
            7   -   Exit Sealed Flt Tubesheet
        orientation: Orientaction del cambiador
            0   -   Horizontal
            1   -   Vertical

    Métodos:
        tubesideLaminar: Método de cálculo de h en el lado del tubo en regimen laminar
            0   -   Eubank-Proctor
            1   -   VDI mean Nusselt
            2   -   Hausen
            3   -   Sieder-Tate
        tubesideTurbulent: Método de cálculo de h en el lado del tubo en regimen turbulento
            0   -   Sieder-Tate
            1   -   Colburn
            2   -   Dittus-Boelter
            3   -   ESDU
            4   -   Gnielinski
            5   -   VDI mean Nusselt
        shellsideSensible: Método de cálculo de h en el lado de la carcasa
            0   -   Stream analysis
            1   -   Bell-Delaware
            2   -   Kern

    Tubo:
        NTubes: Número de tubos
        NPases: Número de pasos de los tubos por la carcasa
        LTube: Lóngitud de tubos
        DeTube: Diametro externo
        wTube: Espesor de la tubería
        rTube: rugosidad interna de los tubos
        kTube: Conductividad térmica
        distribucionTube: Distribucion de tubos
            0   -   Triangular, 30º
            1   -   Diamante, 45º
            2   -   Rotated Triangular, 60º
            3   -   Square, 90º
        pitch: Espacio entre tuberias
        finned: boolean que indica que la tubería tiene alerones
            0   -   tubería lisa
            1   -   tubería con alerones
        Nfin: numero de aletas por metro de tubería
        heightFin:
        foulingTube: resistencia por depositos en la parte del tubo

    Carcasa:
        parallel: Número de intercambiadores en paralelo
        serie: Número de intercambiadores en serie
        DShell: Diematro de la carcasa
        foulingShell: resistencia por depositos en la parte de la carcasa

    Baffle:
        typeBaffle: Tipo de baffle
            0   -   Single segmental
            1   -   Double segmental
            2   -   Triple segmental
            3   -   No tubes in window
            4   -   Disk & donut
            5   -   Rod
        baffleSpacingIn: Espacio de separación anterior al primer bafle
        baffleSpacing: Espacio de separación entre baffles
        baffleSpacingOut: Espacio de separación posterior al último bafle
        baffleThickness: Espesor de los baffles
        BaffleOrientation: Orientación de los baffles
            0   -   Horizontal
            1   -   Vertical
        baffleCut: Porcentaje de corte de los baffles
        baffleCutBase: Base de cálculo del porcentaje de corte
            0   -   Diámetro
            1   -   Área

    Clearances:
        clearanceTubeBaffle
        clearanceShellBaffle
        clearanceShellBundle
        sealingStrips
    Coste:
        tipo: tipo de cambiador
            0   -   Fired head
            1   -   Kettle reboiler
            2   -   U-tubes
        material:
            0   -   Carbon Steel
            1   -   Stainless steel 316
            2   -   Stainless steel 304
            3   -   Stainless steel 347
            4   -   Nickel 200
            5   -   Monel 400
            6   -   Inconel 600
            7   -   Incoloy 825
            8   -   Titanium
            9   -   Hastelloy
        P_dis: Presión de diseño, si no se especifica se usará la máxima presión del las corrientes del proceso
    """

    title = tr("pychemqt", "Shell and Tube Heat Exchanger")
    help = ""
    kwargs = {
        "entrada": [],
        "entradaTubo": None,
        "entradaCarcasa": None,

        "class_": 0,
        "frontHead": 0,
        "shell": 0,
        "rearHead": 0,
        "orientation": 0,

        "tubesideLaminar": 0,
        "tubesideTurbulent": 0,
        "shellsideSensible": 0,

        "NTube": 0,
        "NPases": 0,
        "LTube": 0.0,
        "DeTube": 0.0,
        "wTube": 0.0,
        "rTube": 0.0,
        "kTube": 0.0,
        "distribucionTube": 0,
        "pitch": 0,
        "finned": 0,
        "Nfin": 0,
        "heightFin": 0.0,
        "foulingTube": 0.0,

        "parallel": 0,
        "serie": 0,
        "DShell": 0.0,
        "foulingShell": 0.0,

        "baffleType": 0,
        "baffleSpacingIn": 0.0,
        "baffleSpacing": 0.0,
        "baffleSpacingOut": 0.0,
        "baffleThickness": 0.0,
        "BaffleOrientation": 0,
        "baffleCut": 0.0,
        "baffleCutBase": 0,

        "nozzleInTubesideDiameter": 0.0,
        "nozzleOutTubesideDiameter": 0.0,
        "nozzleInShellsideDiameter": 0.0,
        "nozzleOutShellsideDiameter": 0.0,

        "clearanceTubeBaffle": 0.0,
        "clearanceShellBaffle": 0.0,
        "clearanceShellBundle": 0.0,
        "sealingStrips": 0.0,

        "modo": 0,

        "f_install": 3.,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "tipoCoste": 0,
        "materialCoste": 0,
        "P_dis": 0.0}

    indiceCostos = 2

    TEXT_METHOD_TUBE_LAMINAR = ["Eubank-Proctor", "VDI mean Nusselt",
                                "Hausen", "Sieder-Tate"]
    TEXT_METHOD_TUBE_TURBULENT = ["Sieder-Tate", "Colburn", "Dittus-Boelter",
                                  "ESDU", "Gnielinski", "VDI mean Nusselt"]
    TEXT_METHOD_SHELL = ["Stream analysis", "Bell-Delaware", "Kern"]
    TEXT_CLASS = ["TEMA R", "TEMA B", "TEMA C"]
    TEXT_FRONTHEAD = [
        "A - "+tr("pychemqt", "Channel & Removable Cover"),
        "B - "+tr("pychemqt", "Bonnet"),
        "C - "+tr("pychemqt", "Removable Bundle"),
        "D - "+tr("pychemqt", "Special High Pressure Closure"),
        "N - "+tr("pychemqt", "Channel with Tubesheet & Removable Cover")]
    TEXT_SHELL = [
        "E - "+tr("pychemqt", "One Pass"),
        "F - "+tr("pychemqt", "Two Pass"),
        "G - "+tr("pychemqt", "Split Flow"),
        "H - "+tr("pychemqt", "Double Split Flow"),
        "J - "+tr("pychemqt", "Divided Flow"),
        "K - "+tr("pychemqt", "Kettle Reboiler"),
        "X - "+tr("pychemqt", "Cross Flow")]
    TEXT_REARHEAD = [
        "L - "+tr("pychemqt", "Fixed Tubesheet (A head)"),
        "M - "+tr("pychemqt", "Fixed Tubesheet (B head)"),
        "N - "+tr("pychemqt", "Fixed Tubesheet (N head)"),
        "P - "+tr("pychemqt", "Outside Packed Flt Head"),
        "S - "+tr("pychemqt", "Flt Head with Backing Dev"),
        "T - "+tr("pychemqt", "Pull Throught Flt Heat"),
        "U - "+tr("pychemqt", "U-Tube Bundle"),
        "W - "+tr("pychemqt", "Exit Sealed Flt Tubesheet")]
    TEXT_ORIENTATION = [
        tr("pychemqt", "Horizontal"),
        tr("pychemqt", "Vertical")]
    TEXT_DISTRIBUTION_TUBE = [
        tr("pychemqt", "Triangular")+", 30º",
        tr("pychemqt", "Diamond")+", 45º",
        tr("pychemqt", "Rotated Triangular")+", 60º",
        tr("pychemqt", "Square")+", 90º"]
    TEXT_BAFFLE_TYPE = [
        tr("pychemqt", "Single segmental"),
        tr("pychemqt", "Double segmental"),
        tr("pychemqt", "Triple segmental"),
        tr("pychemqt", "No tubes in window"),
        tr("pychemqt", "Disk & donut"),
        tr("pychemqt", "Rod")]
    TEXT_COST_TYPE = [
        tr("pychemqt", "Fixed Head"),
        tr("pychemqt", "Kettle Reboiler"),
        tr("pychemqt", "U-Tube")]
    TEXT_COST_MATERIAL = [
        tr("pychemqt", "Carbon Steel"),
        tr("pychemqt", "Stainless Steel 316"),
        tr("pychemqt", "Stainless Steel 304"),
        tr("pychemqt", "Stainless Steel 347"),
        tr("pychemqt", "Nickel 200"),
        tr("pychemqt", "Monel 400"),
        tr("pychemqt", "Inconel 600"),
        tr("pychemqt", "Incoloy 825"),
        tr("pychemqt", "Titanium"),
        tr("pychemqt", "Hastelloy")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entradaTubo"]:
            self.msg = tr("pychemqt", "undefined tubeside input")
            self.status = 0
            return
        if not self.kwargs["entradaCarcasa"]:
            self.msg = tr("pychemqt", "undefined shellside input")
            self.status = 0
            return

        return True

    def calculo(self):
        if self.kwargs["modo"]:
            #Diseño
            pass

        else:  # Evaluación
            N = self.kwargs["NTube"]
            De = unidades.Length(self.kwargs["DeTube"])
            Di = unidades.Length(De-2*self.kwargs["wTube"])
            L = unidades.Length(self.kwargs["LTube"])

            #TubeSide
            rho = self.kwargs["entradaTubo"].Liquido.rho
            mu = self.kwargs["entradaTubo"].Liquido.mu
            cp = self.kwargs["entradaTubo"].Liquido.cp
            k = self.kwargs["entradaTubo"].Liquido.k
            w = self.kwargs["entradaTubo"].Liquido.caudalmasico
            beta = self.kwargs["entradaTubo"].Liquido.alfav
            v = self.kwargs["entradaTubo"].Q/N*4/pi/Di
            re = Re(D=Di, V=v, rho=rho, mu=mu)
            pr = self.kwargs["entradaTubo"].Liquido.Pr
            gz = Gz(w=w, cp=cp, k=k, L=L)
            gr = Gr(beta=beta, T=self.kwargs["entradaTubo"].T, To=self.kwargs["entradaTubo"].T, L=L, mu=mu)

            if re < 2300:
                if self.kwargs["tubesideLaminar"] == 0:
                    Nu = h_tubeside_laminar_Eubank_Proctor(Pr=pr, Gz=gz, Gr=gr, D=Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 1:
                    Nu = h_tubeside_laminar_VDI(Re=re, Pr=pr, D=Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 2:
                    Nu = h_tubeside_laminar_Hausen(Gz=gz)
                elif self.kwargs["tubesideLaminar"] == 3:
                    Nu = h_tubeside_laminar_Sieder_Tate(Gz=gz, Gr=gr)
            else:
                if self.kwargs["tubesideTurbulent"] == 0:
                    Nu = h_tubeside_turbulent_Sieder_Tate(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 1:
                    Nu = h_tubeside_turbulent_Colburn(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 2:
                    frio = self.kwargs["entradaCarcasa"] > self.kwargs["entradaTubo"]
                    Nu = h_tubeside_turbulent_Dittus_Boelter(Re=re, Pr=pr, calentamiento=frio)
                elif self.kwargs["tubesideTurbulent"] == 3:
                    Nu = h_tubeside_turbulent_ESDU(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 4:
                    Nu = h_tubeside_turbulent_Gnielinski(Re=re, Pr=pr, D=Di, L=L)
                elif self.kwargs["tubesideTurbulent"] == 5:
                    line = self.kwargs["distribucionTube"] == 3
                    filas = self.kwargs["NTube"]**0.5
                    Nu = h_tubeside_turbulent_VDI(Re=re, Pr=pr, filas_tubos=filas, alineados=line)

            hi = unidades.HeatTransfCoef(Nu*k/Di)

            # ShellSide
            if self.kwargs["shellsideSensible"] == 0:
                h, DP = self.h_shellside_turbulent_Stream_Analysis()
            elif self.kwargs["shellsideSensible"] == 1:
                self.h_shellside_turbulent_Bell_Delaware()
            else:
                h = self.h_shelside_turbulent_Kern()

            # Fouling
            fi = self.kwargs["foulingShell"]
            fo = self.kwargs["foulingTube"]

        #F: Heat exchanger Design, Operation, Maintenance and Enhancement - Ali A. Rabah
#        U=1/((1/ho+fo)/Ef+fw+fi(Ao/Ai)+Ao/Ai/hi)
#        deltaT=(deltaT2-deltaT1)/log(deltaT2/deltaT1)
#        U=1/(Do/hi/Di+Do*log(Do/Di)/2/k+1/ho)
#        #TODO: añadir resistencias de depositos en la pared
#        """Serth - Process heat transfer_ principles and applications pag 102"""
#        q=U*Ao*deltaT
        self.area = unidades.Area(25)


    def fw(self):
        if self.kwargs["finned"]:
            fw=self.kwargs["wTube"]/self.kwargs["kTube"]*((self.kwargs["DeTube"]+2*self.kwargs["Nfin"]*self.kwargs["heightFin"]*(self.kwargs["DeTube"]+self.kwargs["heightFin"]))/(self.kwargs["DeTube"]-self.kwargs["wTube"]))
        else:
            fw=self.kwargs["DeTube"]/2/self.kwargs["kTube"]*log(self.kwargs["DeTube"]/(self.kwargs["DeTube"]-2*self.kwargs["wTube"]))
        return fw


    @staticmethod
    def h_shellside_turbulent_Stream_Analysis():
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        Serth - Process Heat Transfer - Principles and applications Cap. 7
        """


    def h_shellside_turbulent_Bell_Delaware(self):
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        Serth - Process Heat Transfer - Principles and applications Cap. 6
        """
        fi = 1

        P = self.kwargs["pitch"]/self.kwargs["DeTube"]
        mo = self.kwargs["entradaCarcasa"].caudalmasico
        if self.kwargs["distribucionTube"] == 2:
            Ptef = self.kwargs["pitch"]/2**0.5
        else:
            Ptef = self.kwargs["pitch"]

        if self.kwargs["distribucionTube"] < 3:
            tita_tp = unidades.Angle([30, 45, 60][self.kwargs["distribucionTube"]], "deg")
            Pt_ = self.kwargs["pitch"]*cos(tita_tp)
        else:
            Pt_ = self.kwargs["pitch"]
        Nc = self.kwargs["DShell"]*(1-2*self.kwargs["baffleCut"])/Pt_
        Ncw = 0.8*self.kwargs["baffleCut"]*self.kwargs["DShell"]/Pt_

        Dotl = self.kwargs["DShell"]-2*self.kwargs["clearanceShellBundle"]
        Sm = self.kwargs["baffleSpacing"]*(self.kwargs["DShell"]-Dotl+(Dotl-self.kwargs["DeTube"])/Ptef*(self.kwargs["pitch"]-self.kwargs["DeTube"]))

        G = mo/Sm
        Re = self.kwargs["DeTube"]*G/self.kwargs["entradaCarcasa"].Liquido.mu
        Pr = self.kwargs["entradaCarcasa"].Liquido.Pr
        rho = self.kwargs["entradaCarcasa"].Liquido.rho
        cp = self.kwargs["entradaCarcasa"].Liquido.cp

        Dctl = Dotl-self.kwargs["DeTube"]
        tita_ctl = 2*arccos(self.kwargs["DShell"]*(1-2*self.kwargs["baffleCut"])/Dctl)
        Fc = 1+1./pi*(sin(tita_ctl)-tita_ctl)
        Fw = 1./2/pi*(tita_ctl-sin(tita_ctl))
        Stb = pi/8*((self.kwargs["DeTube"]+2*self.kwargs["clearanceTubeBaffle"])**2-self.kwargs["DeTube"]**2)*self.kwargs["NTube"]*(1+Fc)

        tita_ds = 2*arccos(1-2*self.kwargs["baffleCut"])
        Ssb = self.kwargs["DShell"]*self.kwargs["clearanceShellBaffle"]*(pi-0.5*tita_ds)

        Sb = self.kwargs["baffleSpacing"]*(self.kwargs["DShell"]-Dotl)
        Sw = 1./8*self.kwargs["DShell"]**2*(tita_ds-sin(tita_ds))-1./4*self.kwargs["NTube"]*Fw*pi*self.kwargs["DeTube"]**2
        Dw = 4*Sw/(pi*self.kwargs["DeTube"]*self.kwargs["NTube"]*0.5*(1-Fc)+self.kwargs["DShell"]*tita_ds)

        Jc = 0.55+0.72*Fc

        rs = Ssb/(Ssb+Stb)
        rl = (Ssb+Stb)/Sm
        Jl = 0.44*(1-rs)+(1-0.44*(1-rs))*exp(-2.2*rl)
        Rl = exp(-1.33*(1+rs)*rl**(0.8-0.15*(1+rs)))

        rss = self.kwargs["sealingStrips"]/Nc
        if rss < 0.5:
            if Re < 100:
                Cj = 1.35
                Cr = 4.5
            else:
                Cj = 1.25
                Cr = 3.7
            Jb = exp(-Cj*Sb/Sm*(1-(2*rss)**(1./3)))
            Rb = exp(-Cr*Sb/Sm*(1-(2*rss)**(1./3)))
        else:
            Jb = 1.
            Rb = 1.

        if Re < 100:
            n1 = 1./3
            n2 = 1.
        else:
            n1 = 0.6
            n2 = 0.2
        nb = self.kwargs["LTube"]/self.kwargs["baffleSpacing"]
        if self.kwargs["baffleSpacingIn"] > self.kwargs["baffleSpacing"]:
            nb += 1
        if self.kwargs["baffleSpacingOut"] > self.kwargs["baffleSpacing"]:
            nb += 1
        Js = (nb-1+(self.kwargs["baffleSpacingIn"]/self.kwargs["baffleSpacing"])**(1-n1)+(self.kwargs["baffleSpacingOut"]/self.kwargs["baffleSpacing"])**(1-n1))/(nb-1+(self.kwargs["baffleSpacingIn"]/self.kwargs["baffleSpacing"])+self.kwargs["baffleSpacingOut"]/self.kwargs["baffleSpacing"])
        Rs = 0.5*((self.kwargs["baffleSpacing"]/self.kwargs["baffleSpacingIn"])**(2-n2)+(self.kwargs["baffleSpacing"]/self.kwargs["baffleSpacingOut"])**(2-n2))

        Nct = (nb+1)*(Nc+Ncw)
        if Re <= 20:
            Jr = (10/Nct)**0.18
        elif Re >= 100:
            Jr = 1.
        else:  # Interpolacion entre los valores de arriba
            Jr = 0.853379+0.0014662*Re

        if self.kwargs["distribucionTube"] == 1:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 1.55
                a2 = -0.667
            elif Re < 100:
                a1 = 0.498
                a2 = -0.656
            elif Re < 1000:
                a1 = 0.73
                a2 = -0.500
            elif Re < 10000:
                a1 = 0.37
                a2 = -0.396
            else:
                a1 = 0.37
                a2 = -0.396
                a3 = 1.93
                a4 = 0.5
        elif self.kwargs["distribucionTube"] == 3:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 0.97
                a2 = -0.667
            elif Re < 100:
                a1 = 0.9
                a2 = -0.631
            elif Re < 1000:
                a1 = 0.408
                a2 = -0.46
            elif Re < 10000:
                a1 = 0.107
                a2 = -0.266
            else:
                a1 = 0.37
                a2 = -0.395
                a3 = 1.187
                a4 = 0.37
        else:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 1.4
                a2 = -0.667
            elif Re < 100:
                a1 = 1.36
                a2 = -0.657
            elif Re < 1000:
                a1 = 0.593
                a2 = -0.477
            elif Re < 10000:
                a1 = 0.321
                a2 = -0.388
            else:
                a1 = 0.321
                a2 = -0.388
                a3 = 1.45
                a4 = 0.519

        if self.kwargs["distribucionTube"] == 1:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 32
                b2 = -1.
            elif Re < 100:
                b1 = 26.2
                b2 = -0.913
            elif Re < 1000:
                b1 = 3.5
                b2 = -0.476
            elif Re < 10000:
                b1 = 0.333
                b2 = -0.136
            else:
                b1 = 0.303
                b2 = -0.126
                b3 = 6.59
                b4 = 0.52
        elif self.kwargs["distribucionTube"] == 3:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 35.0
                b2 = -1.
            elif Re < 100:
                b1 = 32.1
                b2 = -0.963
            elif Re < 1000:
                b1 = 6.09
                b2 = -0.602
            elif Re < 10000:
                b1 = 0.0815
                b2 = 0.022
            else:
                b1 = 0.391
                b2 = -0.148
                b3 = 6.3
                b4 = 0.378
        else:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 48.0
                b2 = -1.
            elif Re < 100:
                b1 = 45.1
                b2 = -0.973
            elif Re < 1000:
                b1 = 4.570
                b2 = -0.476
            elif Re < 10000:
                b1 = 0.486
                b2 = -0.152
            else:
                b1 = 0.372
                b2 = -0.123
                b3 = 7.0
                b4 = 0.5

        a = a3/(1+0.14*Re**a4)
        b = b3/(1+0.14*Re**b4)

        j = a1*(1.33/P)**a*Re**a2
        f = b1*(1.33/P)**b*Re**b2
        hid = j*cp*G*fi/Pr**(2./3)
        h = unidades.HeatTransfCoef(hid*Jc*Jl*Jb*Jr*Js)

        DPideal = 2*f*Nc*G**2/g/rho/fi
        DPc = (nb-1)*DPideal*Rl*Rb
        if Re > 100:
            DPwideal = (2+0.6*Ncw)*mo**2/2/g/rho/Sm/Sw
        else:
            DPwideal = 26*nu*mo/g/(Sm*Sw)**0.5*(Ncw/P+self.kwargs["baffleCut"]*self.kwargs["DShell"]/Dw**2)+mo**2/g/rho/Sm/Sw
        DPw = nb*DPwideal*Rl
        DPe = 2*DPideal*(1+Ncw/Nc)*Rb*Rs

        mon = mo/self.kwargs["parallel"]

        DPn = 0
        for D in self.kwargs["nozzleInShellsideDiameter"], self.kwargs["nozzleOutShellsideDiameter"]:
            Ren = D*mon/self.kwargs["entradaCarcasa"].Liquido.mu
            Sn = pi/4*D**2
            if Ren > 100:
                DPn += 2e-13*self.kwargs["serie"]*mon**2/Sn
            else:
                DPn += 4e-13*self.kwargs["serie"]*mon**2/Sn
        DP = unidades.DeltaP(DPc+DPw+DPe+DPn)
        return h, DP

    @staticmethod
    def h_shelside_turbulent_Kern(Re, Pr):
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        10<Re<1e6
        kern pag 137
        """
        return 0.36*Re**0.55*Pr**(1./3)*(mu/muw)**0.14


    def h_tubeside_laminar_condensation_Kern(self):
        return 0.815*(k**3*rho_l*(rho_l-rho_g)*g*l/(pi*mu_l*Do*(T-Tw)))**0.25


    def h_tubeside_laminar_condensation_Nusselt(self):
        return 0.72*eg**0.75*(k**3*rho_l*(rho_l-rho_g)*g*hlg/(mu_l*Do*(T-Tw)))**0.25

    def h_tubeSide_fined_Young(self):
        """Briggs, Katz, and Young, Chem.Eng. Prog., 59(11), 49–59 (1963)"""
        return 0.1378*Re**0.718*Pr**(1./3)*(finSpacing/finHeight)**0.296

    def coste(self):
        if self.kwargs["P_dis"]:
            Pd = unidades.Pressure(self.kwargs["P_dis"])
        else:
            Pd = unidades.Pressure(max(self.kwargs["entradaTubo"].P, self.kwargs["entradaCarcasa"].P))

        if self.kwargs["tipoCoste"] == 0:  # Fired head
            Fd = exp(-1.1156+0.09060*log(self.area.ft2))
        elif self.kwargs["tipoCoste"] == 1:  # Kettle reboiler
            Fd = 1.35
        else:  # U-tubes
            Fd = exp(-0.9816+0.0803*log(self.area.ft2))

        g1 = [0., 0.8603, 0.8193, 0.6116, 1.5092, 1.2989, 1.204, 1.1854, 1.5420, 0.1549][self.kwargs["materialCoste"]]
        g2 = [1., 0.23296, 0.15984, 0.22186, 0.60859, 0.43377, 0.50764, 0.49706, 0.42913, 0.51774][self.kwargs["materialCoste"]]
        Fm = g1+g2*log(self.area.ft2)

        if Pd.psi <= 300:
            Fp = 0.771+0.04981*log(self.area.ft2)
        elif Pd.psi <= 600:
            Fp = 1.0305+0.0714*log(self.area.ft2)
        else:
            Fp = 1.14+0.12088*log(self.area.ft2)

        C_base = exp(8.821-0.30863*log(self.area.ft2)+0.0681*log(self.area.ft2)**2)
        C = Fd*Fm*Fp*C_base

        self.C_adq = unidades.Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst = unidades.Currency(self.C_adq*self.kwargs["f_install"])


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
