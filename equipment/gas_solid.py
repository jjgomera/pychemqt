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
# Library with gas-solid separation equipment
#  - Gravty settling chamber
#  - Cyclon
#  - Baghouse
#  - Electrostatic precipitator
###############################################################################


import os

from PyQt5.QtWidgets import QApplication
from scipy import exp, sqrt, ceil, roots
from scipy.constants import pi, g, e, epsilon_0
from scipy.optimize import fsolve

from lib.unidades import (Length, Pressure, DeltaP, Speed, Time, Area, VolFlow,
                          PotencialElectric, Currency, Dimensionless, MassFlow)
from lib.datasheet import pdf
from lib.corriente import Corriente
from equipment.parents import equipment


class Separador_SolidGas(equipment):
    """Generic class with common functionality of gas-solid equipment"""

    def calcularRendimiento(self, rendimientos):
        entrada = self.kwargs["entrada"]
        rendimiento_global = 0
        for i, fraccion in enumerate(entrada.solido.fracciones):
            rendimiento_global += rendimientos[i]*fraccion
        return Dimensionless(rendimiento_global)

    def CalcularSalidas(self):
        entrada = self.kwargs["entrada"]
        unfiltered, filtered = entrada.solido.Separar(self.rendimiento_parcial)
        Pout = entrada.P-self.deltaP
        self.salida = []
        if self.rendimiento == 0:
            self.salida.append(entrada.clone(P=Pout))
            self.salida.append(entrada.clone(split=0, P=Pout))
        else:
            self.salida.append(entrada.clone(solido=unfiltered, P=Pout))
            self.salida.append(
                Corriente(P=Pout, T=entrada.T, solido=filtered))

        self.Pin = entrada.P
        self.Pout = self.salida[0].P
        self.Min = entrada.solido.caudal
        self.Dmin = entrada.solido.diametro_medio
        self.Mr = self.salida[0].solido.caudal
        self.Dmr = self.salida[0].solido.diametro_medio
        self.Ms = self.salida[1].solido.caudal
        self.Dms = self.salida[1].solido.diametro_medio

    def propTxt(self, i):
        """i: index of common properties in equipment subclas list"""
        entrada = self.kwargs["entrada"]

        txt = os.linesep + "#---------------"
        txt += QApplication.translate("pychemqt", "Separation efficiency")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(i, i+6)) + os.linesep
        if len(entrada.solido.diametros) >= 1:
            txt += "%12s %11s %9s %9s %9s" % ("Dp", "Xi", "ηi", "Xis", "Xig")
            txt += os.linesep
            for d, X, eta, Xs, Xg in zip(
                    entrada.solido.diametros, entrada.solido.fracciones,
                    self.rendimiento_parcial, self.salida[1].solido.fracciones,
                    self.salida[0].solido.fracciones):
                txt += "%15s %9.4f %9.5f %9.5f %9.5f" % (d.str, X, eta, Xs, Xg)
                txt += os.linesep
        return txt

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Pin"] = self.Pin
        state["Pout"] = self.Pout
        state["Min"] = self.Min
        state["Dmin"] = self.Dmin
        state["Mr"] = self.Mr
        state["Dmr"] = self.Dmr
        state["Ms"] = self.Ms
        state["Dms"] = self.Dms
        state["rendimiento_parcial"] = self.rendimiento_parcial
        state["rendimiento"] = self.rendimiento

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Pin = Pressure(state["Pin"])
        self.Pout = Pressure(state["Pout"])
        self.Min = MassFlow(state["Min"])
        self.Dmin = Length(state["Dmin"])
        self.Mr = MassFlow(state["Mr"])
        self.Dmr = Length(state["Dmr"])
        self.Ms = MassFlow(state["Ms"])
        self.Dms = Length(state["Dms"])
        eta = [Dimensionless(nu) for nu in state["rendimiento_parcial"]]
        self.rendimiento_parcial = eta
        self.rendimiento = Dimensionless(state["rendimiento"])


class GravityChamber(Separador_SolidGas):
    """Class to define a gravity chamber for separation of solid from gaseous
    stream

    Parameters:
        entrada: Corriente instance to define the input stream to equipment
        metodo: Integer to define the calculate type
            0 - Rating of a equipment with known geometry and input stream
            1 - Design of a equipment to meet a required efficiency
        modelo: Integer to define the method to calculate
            0 - Plug Flow model without vertical mixing
            1 - Perfect mix
        W, H, L: dimensión of equipment (width, height, length)
        rendimientoAdmisible: Required efficiency of equipment, mandatory if
            chosen method of calculate is design
        velocidadAdmisible: Maximum speed of gas in chamber, default 1 m/s
        deltaP: Pressure loss

    >>> from lib.corriente import Corriente
    >>> from lib.solids import Solid
    >>> dm = [17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, \
        54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, \
        0.1, 0.05, 0.03, 0.02]
    >>> solido = Solid(caudalSolido=[0.1], distribucion_diametro=dm, \
        distribucion_fraccion=fracciones, solids=[638])
    >>> kw = {"ids": [475], "fraccionMolar": [1.], "MEoS": True}
    >>> entrada = Corriente(T=300, P=1e5, caudalMasico=1, solido=solido, **kw)
    >>> camara = GravityChamber(entrada=entrada, modelo=0, H=1, W=1, L=1)
    >>> print("%0.4f %0.4f" % (camara.Vgas, camara.rendimiento))
    0.4249 0.5588
    """

    title = QApplication.translate("pychemqt", "Gravity settling chamber")
    kwargs = {"entrada": None,
              "metodo": 0,
              "modelo": 0,
              "W": 0.0,
              "H": 0.0,
              "L": 0.0,
              "rendimientoAdmisible": 0.0,
              "velocidadAdmisible": 0.0,
              "deltaP": 0.0}

    kwargsInput = ("entrada", )
    kwargsValue = ("W", "H", "L", "rendimientoAdmisible", "velocidadAdmisible",
                   "deltaP")
    kwargsList = ("metodo", "modelo")
    calculateValue = ("Q", "LCalc", "WCalc", "HCalc", "Vgas",  "rendimiento")

    TEXT_TIPO = [QApplication.translate("pychemqt", "Rating"),
                 QApplication.translate("pychemqt", "Design")]
    TEXT_MODEL = [
        QApplication.translate("pychemqt", "Plug flow without vertical mix"),
        QApplication.translate("pychemqt", "Vertical mix")]

    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            return

        if self.kwargs["metodo"]:
            if self.kwargs["rendimientoAdmisible"]:
                self.msg = ""
                self.status = 1
                return True
            else:
                self.msg = QApplication.translate(
                    "pychemqt", "undefined efficiency")
                self.status = 0
        else:
            if not self.kwargs["W"]:
                self.msg = QApplication.translate(
                    "pychemqt", "undefined width")
                self.status = 0
            elif self.kwargs["W"] and self.kwargs["H"] and self.kwargs["L"]:
                self.msg = ""
                self.status = 1
                return True
            elif not self.kwargs["H"]:
                self.msg = QApplication.translate(
                    "pychemqt", "height undefined, using default")
                self.status = 3
                return True
            else:
                self.msg = ""
                self.status = 1
                return True

    def calculo(self):
        entrada = self.kwargs["entrada"]
        W = Length(self.kwargs["W"])
        L = Length(self.kwargs["L"])
        self.deltaP = DeltaP(self.kwargs["deltaP"])
        if self.kwargs["H"]:
            self.H = Length(self.kwargs["H"])
        else:
            self.H = Length(1)
        if self.kwargs["velocidadAdmisible"]:
            velocidadAdmisible = self.kwargs["velocidadAdmisible"]
        else:
            velocidadAdmisible = 1.

        if self.kwargs["metodo"] == 0:  # Calculo
            self.Vgas = Speed(entrada.Q/self.H/W)
            self.rendimiento_parcial = self.calcularRendimientos_parciales(L)
            rendimiento = self.calcularRendimiento(self.rendimiento_parcial)
            self.rendimiento = Dimensionless(rendimiento)
        else:  # Diseño
            self.Vgas = Speed(velocidadAdmisible)
            W = Length(entrada.Q/velocidadAdmisible/self.H)

            def f(longitud):
                eta = self.calcularRendimientos_parciales(longitud)
                eta_adm = self.kwargs["rendimientoAdmisible"]
                return self.calcularRendimiento(eta)-eta_adm

            if f(0.1) > self.kwargs["rendimientoAdmisible"]:
                longitud = 1
            else:
                longitud_ = 0.1
                longitud = fsolve(f, longitud_)[0]
            L = Length(float(longitud))
            eta_i = self.calcularRendimientos_parciales(L)
            self.rendimiento_parcial = eta_i
            self.rendimiento = Dimensionless(self.calcularRendimiento(eta_i))

        self.LCalc = L
        self.WCalc = W
        self.HCalc = self.H
        self.Q = entrada.Q

        self.CalcularSalidas()

    def calcularRendimientos_parciales(self, L):
        """Calculate the efficiency of separation process for each diameter of
        solid fraction"""
        entrada = self.kwargs["entrada"]
        rhoS = entrada.solido.rho
        rhoG = entrada.Gas.rho
        muG = entrada.Gas.mu
        rendimientos = []
        for d in entrada.solido.diametros:
            Ar = d**3*rhoG*(rhoS-rhoG)*g/muG**2
            Vt = muG/d*rhoG*((14.42+1.827*Ar**0.5)**0.5-3.798)**2
            if self.kwargs["modelo"] == 0:
                r = Vt*L/self.Vgas/self.H
            else:
                r = 1-exp(-Vt*L/self.Vgas/self.H)
            if r > 1:
                rendimientos.append(1)
            else:
                rendimientos.append(r)
        return rendimientos

    def propTxt(self):
        txt = os.linesep + "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(11))

        txt += Separador_SolidGas.propTxt(self, 12)
        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Mode"),
              ("TEXT_TIPO", "metodo"), str),
             (QApplication.translate("pychemqt", "Model"),
              ("TEXT_MODEL", "modelo"), str),
             (QApplication.translate("pychemqt", "Input Pressure"), "Pin",
              Pressure),
             (QApplication.translate("pychemqt", "Output Pressure"), "Pout",
              Pressure),
             (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP",
              DeltaP),
             (QApplication.translate("pychemqt", "Height"), "HCalc", Length),
             (QApplication.translate("pychemqt", "Width"), "WCalc", Length),
             (QApplication.translate("pychemqt", "Length"), "LCalc", Length),
             (QApplication.translate("pychemqt", "Gas Speed"), "Vgas", Speed),
             (QApplication.translate("pychemqt", "Gas Volumetric Flow"), "Q",
              VolFlow),
             (QApplication.translate("pychemqt", "Efficiency"), "rendimiento",
              Dimensionless),
             (QApplication.translate("pychemqt", "Partial Efficiency"),
              "rendimiento_parcial", Dimensionless),
             (QApplication.translate("pychemqt", "Input Solid Mass Flow"),
              "Min", MassFlow),
             (QApplication.translate("pychemqt", "Input Solid Mean Diameter"),
              "Dmin", Length),
             (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"),
              "Mr", MassFlow),
             (QApplication.translate(
                 "pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
             (QApplication.translate("pychemqt", "Solid Output Mass Flow"),
              "Ms", MassFlow),
             (QApplication.translate("pychemqt", "Solid Output Mean Diameter"),
              "Dms", Length)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Vgas"] = self.Vgas
        state["LCalc"] = self.LCalc
        state["WCalc"] = self.WCalc
        state["HCalc"] = self.HCalc
        state["Q"] = self.Q
        state["deltaP"] = self.deltaP
        Separador_SolidGas.writeStatetoJSON(self, state)

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Vgas = Speed(state["Vgas"])
        self.LCalc = Length(state["LCalc"])
        self.WCalc = Length(state["WCalc"])
        self.HCalc = Length(state["HCalc"])
        self.Q = VolFlow(state["Q"])
        self.deltaP = DeltaP(state["deltaP"])
        Separador_SolidGas.readStatefromJSON(self, state)

        self.salida = [None]


class Ciclon(Separador_SolidGas):
    """Class to model a cyclone equipment

    Parameters:
        Entrada: Corriente instance to define the input stream to equipment
        tipo_calculo:
            0 - Rating
            1 - Design
        modelo_rendimiento: Model to simulate the equipment:
            0 - Rosin-Rammler-Intelmann
            1 - Leith-Licht
        modelo_DeltaP: Method to calculate the pressure loss:
            0 - Simple, only cinetic loss
            1 - Casal-Martinez-Benet
            2 - Leith-Licht
            3 - Sheferd, Lapple y Ter Linden
        modelo_ciclón: Use a standard model to dimension the equipment
            0 - Stairmand (High η)
            1 - Swift (High η)
            2 - Lapple (Low η)
            3 - Swift (Low η)
            4 - Peterson/Whitby (Low η)
            5 - Lorenz I
            6 - Lorenz II
            7 - Lorenz III
            8 - Custom
        diametro: Cyclone diameter
        num_ciclones: Cyclone count working in parallel
        dimensiones: If use a custom cyclone model, it can be defined here with
            array with the custom dimension in the order:
                Inlet Height, Hc
                Inlet Width, Bc
                Solid Output Diameter, Jc
                Cylinder Cyclone Section Length, Lc
                Conical Cyclone Section Length, Zc
                Clean Gas Output Diameter, De
                Clean Gas Inlet Orifice Length, Sc

        rendimientoAdmisible: Required efficiency orcyclone (design)
        DeltaPAdmisible: Maximum pressure loss permisible of cyclone
        velocidadAdmisible: Input gas speed to equipment

    Coste
        tipo_costo:
            0 - Heavy duty
            1 - Standart duty
            2 - Multiciclone

    >>> from lib.corriente import Corriente
    >>> from lib.solids import Solid
    >>> dm = [17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, \
        54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, \
        0.1, 0.05, 0.03, 0.02]
    >>> sol = Solid(caudalSolido=[0.001], distribucion_diametro=dm, \
        distribucion_fraccion=fracciones, solids=[638])
    >>> kw = {"ids": [475], "fraccionMolar": [1.], "MEoS": True}
    >>> entrada = Corriente(T=300, P=1e5, caudalMasico=0.3, solido=sol, **kw)
    >>> ciclon = Ciclon(entrada=entrada, tipo_calculo=1, \
        rendimientoAdmisible=0.95, velocidadAdmisible=5)
    >>> print("%0.2f %0.2f" % (ciclon.C_instTotal, ciclon.C_adqTotal))
    2302.76 1644.83
    """
    title = QApplication.translate("pychemqt", "Cyclone")
    help = os.environ["pychemqt"] + "doc/Ciclones.htm"
    kwargs = {"entrada": None,
              "tipo_calculo": 0,
              "modelo_rendimiento": 0,
              "modelo_DeltaP": 0,
              "modelo_ciclon": 0,
              "Dc": 0.0,
              "num_ciclones": 0,
              "dimensiones": [],
              "rendimientoAdmisible": 0.0,
              "DeltaPAdmisible": 0.0,
              "velocidadAdmisible": 0.0,

              "f_install": 1.4,
              "Base_index": 0.0,
              "Current_index": 0.0,
              "tipo_costo": 0}
    kwargsInput = ("entrada", )
    kwargsValue = ("Dc", "num_ciclones", "rendimientoAdmisible",
                   "velocidadAdmisible", "DeltaPAdmisible")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP",
                  "modelo_ciclon", "tipo_costo")
    calculateValue = ("deltaP", "V", "rendimientoCalc", "NCalc", "Dcc", "Hc",
                      "Bc", "Jc", "Lc", "Zc", "De", "Sc")
    calculateCostos = ("C_adq", "C_inst", "num_ciclonesCoste", "Q")
    indiceCostos = 2

    TEXT_TIPO = [
        QApplication.translate("pychemqt", "Rating"),
        QApplication.translate("pychemqt", "Design")]
    TEXT_MODEL = ["Rossin, Rammler & Intelmann", "Leith & Licht"]
    TEXT_MODEL_DELTAP = [
        QApplication.translate("pychemqt", "Standart"),
        "Casal & Martinez-Benet", "Leith & Licht",
        "Sheferd, Lapple & Ter Linden"]
    TEXT_MODEL_CICLON = [
        "Stairmand ("+QApplication.translate("pychemqt", "High η")+")",
        "Swift ("+QApplication.translate("pychemqt", "High η")+")",
        "Lapple ("+QApplication.translate("pychemqt", "Low η")+")",
        "Swift ("+QApplication.translate("pychemqt", "Low η")+")",
        "Peterson/Whitby ("+QApplication.translate("pychemqt", "Low η")+")",
        "Lorenz I", "Lorenz II", "Lorenz III",
        QApplication.translate("pychemqt", "Custom")]
    TEXT_COST = [
        QApplication.translate("pychemqt", "Heavy duty"),
        QApplication.translate("pychemqt", "Standard dury"),
        QApplication.translate("pychemqt", "Multicyclone")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            return

        if self.kwargs["tipo_calculo"]:
            if self.kwargs["modelo_ciclon"] == 8:
                if self.kwargs["Dc"] and self.kwargs["Hc"] and \
                        self.kwargs["Bc"] and self.kwargs["De"]:
                    self.msg = ""
                    self.status = 1
                    return True
                else:
                    self.msg = QApplication.translate(
                        "pychemqt", "undefined cyclone dimension")
                    self.status = 0
            else:
                if (self.kwargs["DeltaPAdmisible"] or
                    self.kwargs["velocidadAdmisible"]) and \
                        self.kwargs["rendimientoAdmisible"]:
                    self.msg = ""
                    self.status = 1
                    return True
                elif self.kwargs["rendimientoAdmisible"]:
                    self.msg = QApplication.translate(
                        "pychemqt", "undefined efficiency")
                    self.status = 0
                else:
                    self.msg = QApplication.translate(
                        "pychemqt", "undefined loss pressure specification")
                    self.status = 0

        else:
            if self.kwargs["Dc"] and self.kwargs["num_ciclones"]:
                self.msg = ""
                self.status = 1
                return True
            elif self.kwargs["Dc"]:
                self.msg = QApplication.translate(
                    "pychemqt", "undefined cyclone number")
                self.status = 0
            else:
                self.msg = QApplication.translate(
                    "pychemqt", "undefined cyclone diameter")
                self.status = 0

    def calculo(self):
        entrada = self.kwargs["entrada"]
        self.Dc = Length(self.kwargs["Dc"])
        self.num_ciclones = int(self.kwargs["num_ciclones"])
        self.dimensiones = self.kwargs["dimensiones"]
        self.rendimientoAdmisible = self.kwargs["rendimientoAdmisible"]
        self.DeltaPAdmisible = Pressure(self.kwargs["DeltaPAdmisible"])

        rhoS = entrada.solido.rho
        rhoG = entrada.Gas.rho
        muG = entrada.Gas.mu

        if self.kwargs["tipo_calculo"] == 0:  # Calculo
            if self.kwargs["modelo_ciclon"] != 8:
                # Hc, Bc, Jc, Lc, Zc, De, Sc,kf, G
                self.dimensiones = self.dimensionado(Dc=self.Dc)
            else:
                dimensiones = self.dimensionado(dimensiones=self.dimensiones)
                self.dimensiones = dimensiones
            self.Hc = Length(self.dimensiones[0])
            self.Bc = Length(self.dimensiones[1])
            self.Jc = Length(self.dimensiones[2])
            self.Lc = Length(self.dimensiones[3])
            self.Zc = Length(self.dimensiones[4])
            self.De = Length(self.dimensiones[5])
            self.Sc = Length(self.dimensiones[6])
            self.kf = Length(self.dimensiones[7])
            self.G = Length(self.dimensiones[8])

            self.V = Speed(entrada.Q/(self.Bc*self.Hc*self.num_ciclones))
            if self.kwargs["modelo_rendimiento"] == 0:
                self.N = 11.3*(self.Hc*self.Bc/self.Sc**2)**2+3.33
                dc = sqrt(9*self.Bc*muG/(2*pi*self.N*self.V*(rhoS-rhoG)))
                self.dc = Length(dc, "ParticleDiameter")
            else:
                self.N = self.V*(0.1079-0.00077*self.V+1.924e-6*self.V**2)

            eta_i = self.calcularRendimientos_parciales()
            self.rendimiento_parcial = eta_i
            self.rendimiento = self.calcularRendimiento(eta_i)

        else:
            if self.DeltaPAdmisible and self.kwargs["velocidadAdmisible"]:
                V_fpresion = self.velocidad_f_presion()
                if self.kwargs["velocidadAdmisible"] > V_fpresion:
                    self.kwargs["velocidadAdmisible"] = V_fpresion
            elif self.DeltaPAdmisible:
                self.kwargs["velocidadAdmisible"] = self.velocidad_f_presion()

            def f(diametro):
                self.dimensiones = self.dimensionado(Dc=diametro)
                self.Dc = Length(diametro)
                self.Hc = Length(self.dimensiones[0])
                self.Bc = Length(self.dimensiones[1])
                self.Jc = Length(self.dimensiones[2])
                self.Lc = Length(self.dimensiones[3])
                self.Zc = Length(self.dimensiones[4])
                self.De = Length(self.dimensiones[5])
                self.Sc = Length(self.dimensiones[6])
                self.kf = Length(self.dimensiones[7])
                self.G = Length(self.dimensiones[8])
                n = entrada.Q/self.Bc/self.Hc/self.kwargs["velocidadAdmisible"]
                self.num_ciclones = int(ceil(n))
                self.V = Speed(entrada.Q/(self.Bc*self.Hc*ceil(n)))
                if self.kwargs["modelo_rendimiento"] == 0:
                    N = 11.3*(self.Hc*self.Bc/self.Sc**2)**2+3.33
                    self.N = Dimensionless(N)
                    dc = sqrt(9*self.Bc*muG/(2*pi*N*self.V*(rhoS-rhoG)))
                    self.dc = Length(dc, magnitud="ParticleDiameter")
                else:
                    N = self.V*(0.1079-0.00077*self.V+1.924e-6*self.V**2)
                    self.N = Dimensionless(N)
                rendimiento_parcial = self.calcularRendimientos_parciales()
                rendimiento = self.calcularRendimiento(rendimiento_parcial)
                return rendimiento-self.kwargs["rendimientoAdmisible"]

            Dc0 = 1
            diametro = fsolve(f, Dc0)
            self.Dc = Length(diametro[0])
            self.dimensiones = self.dimensionado(self.Dc)
            self.Hc = Length(self.dimensiones[0])
            self.Bc = Length(self.dimensiones[1])
            self.Jc = Length(self.dimensiones[2])
            self.Lc = Length(self.dimensiones[3])
            self.Zc = Length(self.dimensiones[4])
            self.De = Length(self.dimensiones[5])
            self.Sc = Length(self.dimensiones[6])
            self.kf = Length(self.dimensiones[7])
            self.G = Length(self.dimensiones[8])
            eta_i = self.calcularRendimientos_parciales()
            self.rendimiento_parcial = eta_i
            self.rendimiento = self.calcularRendimiento(eta_i)

        self.deltaP = self.PerdidaPresion()
        self.num_ciclonesCoste = self.num_ciclones
        self.NCalc = self.num_ciclones
        self.Q = entrada.Q
        self.rendimientoCalc = self.rendimiento
        self.Dcc = self.Dc
        self.CalcularSalidas()

    def dimensionado(self, Dc=0, dimensiones=[]):
        coef = [
            # Stairmand (Alta η)
            [0.5, 0.2, 0.375, 1.5, 2.5, 0.5, 0.5, 6.4, 551.3],
            # Swift (Alta η)
            [0.44, 0.21, 0.4, 1.4, 2.5, 0.4, 0.5, 9.2, 699.2],
            # Lapple (Baja η)
            [0.5, 0.25, 0.25, 2, 2, 0.5, 0.625, 8, 402.9],
            # Swift (Baja η)
            [0.5, 0.25, 0.4, 1.75, 2, 0.5, 0.6, 7.6, 381.8],
            # Peterson/Whitby (Baja η)
            [0.583, 0.208, 0.5, 1.333, 1.837, 0.5, 0.583, 7.760896, 0],
            # Lorenz I
            [0.533, 0.133, 0.333, 0.693, 1.887, 0.333, 0.733, 11.187981, 0],
            # Lorenz II
            [0.533, 0.133, 0.333, 0.693, 1.887, 0.233, 0.733, 22.85221684, 0],
            # Lorenz III
            [0.4, 0.1, 0.333, 0.693, 1.887, 0.233, 0.733, 11.78876, 0]]

        if self.kwargs["modelo_ciclon"] == 8:
            coef = dimensiones
            dimensiones.append(16*coef[0]*coef[1]/coef[5]**2)
            dimensiones.append(0)
        else:
            coef = coef[self.kwargs["modelo_ciclon"]]

        Hc = Dc*coef[0]
        Bc = Dc*coef[1]
        Jc = Dc*coef[2]
        Lc = Dc*coef[3]
        Zc = Dc*coef[4]
        De = Dc*coef[5]
        Sc = Dc*coef[6]
        kf = coef[7]
        G = coef[8]

        return Hc, Bc, Jc, Lc, Zc, De, Sc, kf, G

    def calcularRendimientos_parciales(self):
        entrada = self.kwargs["entrada"]
        rhoS = entrada.solido.rho
        muG = entrada.Gas.mu

        rendimiento_parcial = []
        if self.kwargs["modelo_rendimiento"]:
            # modelo Leith-Licht
            Vo = entrada.Q/self.num_ciclones
            n = 1-(1-(12*self.Dc.ft)**0.14/2.5)*((entrada.T.F+460)/530)**0.3
            if self.G == 0:
                Vs = entrada.solido.caudal/entrada.solido.rho/self.num_ciclones
                self.G = 4*self.Dc*(2*Vs+Vo)/self.num_ciclones / \
                    self.Hc**2/self.Bc**2
            for d in entrada.solido.diametros:
                t = rhoS*(d)**2/18/muG
                rendimiento_parcial.append(
                    1-exp(-2*(self.G*t*Vo/self.Dc**3*(n+1))**(0.5/(n+1))))
        else:
            # model Rosin-Rammler-Intelmann
            for d in entrada.solido.diametros:
                rendimiento_parcial.append((d/self.dc)**2/(1+(d/self.dc)**2))
        return rendimiento_parcial

    def PerdidaPresion(self):
        entrada = self.kwargs["entrada"]
        rhoS = entrada.solido.rho
        rhoG = entrada.Gas.rho
        QS = entrada.solido.caudal
        Q = entrada.caudalmasico

        if self.kwargs["modelo_DeltaP"] == 0:
            # Cinetic loss
            DeltaP = Pressure(self.kf/2.*rhoG*self.V**2)
        elif self.kwargs["modelo_DeltaP"] == 1:
            # Casal-Martinez-Benet
            DeltaP = Pressure(0.003*self.N*rhoG*self.V**2, "inH2O")
        elif self.kwargs["modelo_DeltaP"] == 2:
            # Leith-Licht
            DeltaP = Pressure(0.003*rhoG*16*(entrada.Q/self.num_ciclones)**2 /
                              self.De**2/self.Bc/self.Hc, "mmHg")
        else:
            # Sheferd, Lapple y Ter Linden
            Ae = self.Bc*self.Hc
            As = pi/4*self.De**2
            rhom = rhoG+QS/rhoG/(Q+QS/rhoG)*(rhoS-rhoG)
            DeltaP = Pressure(1.078*(Ae/As)**1.21*rhom*self.V**2, "mmH2O")
        return DeltaP

    def velocidad_f_presion(self):
        entrada = self.kwargs["entrada"]
        rhoS = entrada.solido.rho
        rhoG = entrada.Gas.rho
        dP = self.DeltaPAdmisible
        QS = entrada.solido.caudal
        Q = entrada.caudalmasico

        if self.kwargs["modelo_DeltaP"] == 0:
            velocidad = sqrt(dP*2/self.kf/rhoG)
        elif self.kwargs["modelo_DeltaP"] == 1:
            velocidad = sqrt(dP.inH2O/self.N/0.003/rhoG)
        elif self.kwargs["modelo_DeltaP"] == 2:
            velocidad = sqrt(dP*self.De**2*self.Bc*self.Hc/0.003/rhoG/16)
        else:
            Ae = self.Bc*self.Hc
            As = pi/4*self.De**2
            rhom = rhoG+QS/rhoG/(Q+QS/rhoG)*(rhoS-rhoG)
            velocidad = sqrt(dP.mmH2O/1.078/(Ae/As)**1.21/rhom)
        return velocidad

    def coste(self):
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]
        # TODO: Q is at standard condition
        Q = self.kwargs["entrada"].Q.kft3min
        if self.num_ciclones != 1:
            Q = Q/self.num_ciclones

        if self.kwargs["tipo_costo"] == 0:
            C = 1.39*Q**0.98*1000  # Heavy duty
        elif self.kwargs["tipo_costo"] == 1:
            C = 0.65*Q**0.91*1000  # Standart duty
        else:
            C = 1.56*Q**0.68*1000  # multiciclone

        self.C_adq = Currency(C*CI/BI)
        self.C_inst = Currency(self.C_adq * self.kwargs["f_install"])
        self.C_adqTotal = Currency(self.C_adq*self.num_ciclones)
        self.C_instTotal = Currency(self.C_inst*self.num_ciclones)

    def Pdf(self):
        archivo = pdf("CICLÓN")
        archivo.ciclon("Ciclon limpieza polvo", self)
        archivo.dibujar()
        # if we sant to open with a external editor
        # os.system("xpdf datasheet.pdf")

    def propTxt(self):
        txt = os.linesep + "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(12))

        txt += os.linesep + "#---------------"
        txt += QApplication.translate("pychemqt", "Cyclone Dimensions")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(12, 20))

        txt += Separador_SolidGas.propTxt(self, 20)

        if self.statusCoste:
            txt += os.linesep + "#---------------"
            txt += QApplication.translate(
                "pychemqt", "Preliminary Cost Estimation")
            txt += "-----------------#"+os.linesep
            txt += self.propertiesToText(range(26, 32))

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Mode"),
              ("TEXT_TIPO", "tipo_calculo"), str),
             (QApplication.translate("pychemqt", "Model"),
              ("TEXT_MODEL", "modelo_rendimiento"), str),
             (QApplication.translate("pychemqt", "Pressure Loss Model"),
              ("TEXT_MODEL_DELTAP", "modelo_DeltaP"), str),
             (QApplication.translate("pychemqt", "Ciclon Model"),
              ("TEXT_MODEL_CICLON", "modelo_ciclon"), str),
             (QApplication.translate("pychemqt", "Input Pressure"), "Pin",
              Pressure),
             (QApplication.translate("pychemqt", "Output Pressure"), "Pout",
              Pressure),
             (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP",
              DeltaP),
             (QApplication.translate("pychemqt", "Critic Particle Diameter"),
              "dc", Length),
             (QApplication.translate("pychemqt", "Gas Internal Cycles"), "N",
              Dimensionless),
             (QApplication.translate("pychemqt", "Gas Speed"), "V", Speed),
             (QApplication.translate("pychemqt", "Efficiency"), "rendimiento",
              Dimensionless),
             (QApplication.translate("pychemqt", "Cyclone number"),
              "num_ciclones", int),
             (QApplication.translate("pychemqt", "Ciclon Diameter"), "Dc",
              Length),
             (QApplication.translate("pychemqt", "Inlet Height"), "Hc",
              Length),
             (QApplication.translate("pychemqt", "Inlet Width"), "Bc", Length),
             (QApplication.translate("pychemqt", "Solid Output Diameter"),
              "Jc", Length),
             (QApplication.translate(
                 "pychemqt", "Cylinder Cyclone Section Length"), "Lc", Length),
             (QApplication.translate(
                 "pychemqt", "Conical Cyclone Section Length"), "Zc", Length),
             (QApplication.translate("pychemqt", "Clean Gas Output Diameter"),
              "De", Length),
             (QApplication.translate(
                 "pychemqt", "Clean Gas Inlet Orifice Length"), "Sc", Length),
             (QApplication.translate("pychemqt", "Input Solid Mass Flow"),
              "Min", MassFlow),
             (QApplication.translate("pychemqt", "Input Solid Mean Diameter"),
              "Dmin", Length),
             (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"),
              "Mr", MassFlow),
             (QApplication.translate(
                 "pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
             (QApplication.translate("pychemqt", "Solid Output Mass Flow"),
              "Ms", MassFlow),
             (QApplication.translate("pychemqt", "Solid Output Mean Diameter"),
              "Dms", Length),
             (QApplication.translate("pychemqt", "Base index"),
              "Base_index", float),
             (QApplication.translate("pychemqt", "Current index"),
              "Current_index", float),
             (QApplication.translate("pychemqt", "Install factor"),
              "f_install", float),
             (QApplication.translate("pychemqt", "Cost Mode"),
              ("TEXT_COST", "tipo_costo"), str),
             (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq",
              Currency),
             (QApplication.translate("pychemqt", "Installed Cost"), "C_inst",
              Currency)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Pin"] = self.Pin
        state["Pout"] = self.Pout
        state["deltaP"] = self.deltaP
        state["dc"] = self.dc
        state["N"] = self.N
        state["V"] = self.V
        state["num_ciclones"] = self.num_ciclones
        state["Dc"] = self.Dc
        state["Hc"] = self.Hc
        state["Bc"] = self.Bc
        state["Jc"] = self.Jc
        state["Lc"] = self.Lc
        state["Zc"] = self.Zc
        state["De"] = self.De
        state["Sc"] = self.Sc
        Separador_SolidGas.writeStatetoJSON(self, state)
        state["statusCoste"] = self.statusCoste
        if self.statusCoste:
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Pin = Pressure(state["Pin"])
        self.Pout = Pressure(state["Pout"])
        self.deltaP = DeltaP(state["deltaP"])
        self.dc = Length(state["dc"])
        self.N = Dimensionless(state["N"])
        self.V = Speed(state["V"])
        self.num_ciclones = state["num_ciclones"]
        self.Dc = Length(state["Dc"])
        self.Hc = Length(state["Hc"])
        self.Bc = Length(state["Bc"])
        self.Jc = Length(state["Jc"])
        self.Lc = Length(state["Lc"])
        self.Zc = Length(state["Zc"])
        self.De = Length(state["De"])
        self.Sc = Length(state["Sc"])
        Separador_SolidGas.readStatefromJSON(self, state)
        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.C_adq = Currency(state["C_adq"])
            self.C_inst = Currency(state["C_inst"])
        self.salida = [None]


class Baghouse(Separador_SolidGas):
    """Clase que define los filtros de mangas

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la cámara
        metodo: Integer que indica el tipo de cálculo a realizar:
            0   -   Calcular caida de presión
            1   -   Calcular el tiempo de filtración
            2   -   Calcular el número de filtros necesarios
        num_filtros: número de filtros
        tiempo: tiempo de filtración que lleva operando
        deltaP: perdida de presión a traves de la membrana
        resistenciaFiltro: coeficiente de perdida de presión del filtros, (in water)/(cP)(ft/min)
        resistenciaTorta: coeficiente de perdida de presión de la torta,  (in water)/(cP)(gr/ft2)(ft/min)
        limpieza: indica el número de filtros que están en limpieza
        membranasFiltro: número de membranas que tiene cada filtro
        diametroMembrana: diametro  de cada membrana
        areaMembrana: area de filtración que tiene cada membrana
        rendimientos: array con los rendimientos de fabrica de la membrana

    >>> diametros=[17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> solido=Solid(caudalSolido=[12426.28], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    >>> filtro=Baghouse(entrada=corriente, metodo=1, num_filtros=4, tiempo=3600, deltaP=0.1)
    >>> print filtro.floorArea, filtro.Vgas.ftmin
    7.24643712 0.480966666862
    """
    title=QApplication.translate("pychemqt", "Baghouse")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "num_filtros": 0,
                    "tiempo": 0.0,
                    "deltaP": 0.0,
                    "resistenciaFiltro": 0.0,
                    "resistenciaTorta": 0.0,
                    "limpieza": 0,
                    "membranasFiltro": 0,
                    "diametroMembrana": 0.0,
                    "areaMembrana": 0.0,
                    "rendimientos": []}
    kwargsInput=("entrada", )
    kwargsValue=("num_filtros", "tiempo", "deltaP", "resistenciaFiltro", "resistenciaTorta", "limpieza", "membranasFiltro", "diametroMembrana", "areaMembrana")
    kwargsList=("metodo", )
    calculateValue=("floorArea", "rendimiento", "Vgas", "num_filtrosCalc", "tiempoCalc", "deltaPCalc")

    TEXT_TIPO=[QApplication.translate("pychemqt", "Calculate Pressure drop"),
                        QApplication.translate("pychemqt", "Calculate time of filtration"),
                        QApplication.translate("pychemqt", "Calculate number of cells")]


    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return

        if self.kwargs["metodo"]==0 and (not self.kwargs["num_filtros"] or not self.kwargs["tiempo"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return
        elif self.kwargs["metodo"]==1 and (not self.kwargs["num_filtros"] or not self.kwargs["deltaP"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return
        elif self.kwargs["metodo"]==2 and (not self.kwargs["tiempo"] or not self.kwargs["deltaP"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return

        if self.kwargs["metodo"]==2 and self.kwargs["limpieza"]>=self.kwargs["num_filtros"]:
            self.msg=QApplication.translate("pychemqt", "All filters cleaned")
            self.status=0
            return

        if not self.kwargs["rendimientos"]:
            self.msg=QApplication.translate("pychemqt", "using default efficiency")
            self.status=3
            return True
        if len(self.kwargs["rendimientos"])!=self.kwargs["entrada"].solido.diametros:
            self.msg=QApplication.translate("pychemqt", "using default efficiency")
            self.status=3
            return True

        self.msg=""
        self.status=1
        return True


    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.metodo=self.kwargs["metodo"]
        self.num_filtros=self.kwargs["num_filtros"]
        self.tiempo=Time(self.kwargs["tiempo"])
        self.deltaP=Pressure(self.kwargs["deltaP"])

        if self.kwargs["resistenciaFiltro"]:
            self.resistenciaFiltro=Dimensionless(self.kwargs["resistenciaFiltro"])
        else:
            self.resistenciaFiltro=Dimensionless(0.84)
        if self.kwargs["resistenciaTorta"]:
            self.resistenciaTorta=Dimensionless(self.kwargs["resistenciaTorta"])
        else:
            self.resistenciaTorta=Dimensionless(0.1)
        if self.kwargs["limpieza"]:
            self.limpieza=self.kwargs["limpieza"]
        else:
            self.limpieza=1
        if self.kwargs["membranasFiltro"]:
            self.membranasFiltro=self.kwargs["membranasFiltro"]
        else:
            self.membranasFiltro=78

        if self.kwargs["diametroMembrana"]:
            self.diametroMembrana=Length(self.kwargs["diametroMembrana"])
        else:
            self.diametroMembrana=Length(0.5, "ft")

        if self.kwargs["areaMembrana"]:
            self.areaMembrana=Area(self.kwargs["areaMembrana"])
        else:
            self.areaMembrana=Area(16, "ft2")

        if any(self.kwargs["rendimientos"]):
            self.rendimiento_parcial=self.kwargs["rendimientos"]
        else:
            self.rendimiento_parcial=self.defaultRendimiento()
        self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)

        if self.kwargs["metodo"]==0:
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            ms=self.entrada.solido.caudal.gmin*self.tiempo.min/self.area.ft2
            self.Vgas=Speed(self.entrada.Q/self.area)
            self.deltaP=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin+self.resistenciaTorta*self.entrada.Gas.mu.cP*ms*self.Vgas.ftmin, "inH2O")
        elif self.kwargs["metodo"]==1:
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            self.Vgas=Speed(self.entrada.Q/self.area)
            deltaPMinimo=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin, "inH2O")
            if deltaPMinimo>self.deltaP:
                self.deltaP=deltaPMinimo
                self.tiempo=Time(0)
            else:
                ms=(self.deltaP.inH2O-self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin)/(self.resistenciaTorta*self.entrada.Gas.mu.cP*self.Vgas.ftmin)
                self.tiempo=Time(ms*self.area.ft2/self.entrada.solido.caudal.gmin, "min")
        else:
            N=roots([self.deltaP.inH2O, -self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.entrada.Q.ft3min/self.areaMembrana.ft2/self.membranasFiltro,-self.resistenciaTorta*self.entrada.Gas.mu.cP*self.entrada.Q.ft3min*self.entrada.solido.caudal.gmin*self.tiempo.min/self.areaMembrana.ft2**2/self.membranasFiltro**2])
            if N[0]>0:
                self.num_filtros=ceil(N[0])+self.limpieza
            else:
                self.num_filtros=ceil(N[1])+self.limpieza
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            ms=self.entrada.solido.caudal.gmin*self.tiempo.min/self.area.ft2
            self.Vgas=Speed(self.entrada.Q/self.area)
            self.deltaP=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin+self.resistenciaTorta*self.entrada.Gas.mu.cP*ms*self.Vgas.ftmin, "inH2O")

        self.floorArea=Area(self.num_filtros*self.membranasFiltro*self.diametroMembrana**2)
        self.num_filtrosCalc=self.num_filtros
        self.tiempoCalc=self.tiempo
        self.deltaPCalc=self.deltaP

        self.CalcularSalidas()


    def defaultRendimiento(self):
        """Modelo de rendimiento sylvan, se usara si no se indica el rendimiento"""
        rendimiento=[(d/0.3177e-6)/(1+(d/0.3177e-6)) for d in self.entrada.solido.diametros]
        return rendimiento


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["metodo"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t %i" %(QApplication.translate("pychemqt", "Filter Number"), self.num_filtros)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Operation Time"), self.tiempo.str)

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Cloth resistence"), self.resistenciaFiltro.str)
        if not self.kwargs["resistenciaFiltro"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Cake resistence"), self.resistenciaTorta.str)
        if not self.kwargs["resistenciaTorta"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t %i" %(QApplication.translate("pychemqt", "Cells cleaned"), self.limpieza)
        if not self.kwargs["limpieza"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t %i" %(QApplication.translate("pychemqt", "Bags per cell"), self.membranasFiltro)
        if not self.kwargs["membranasFiltro"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Bag diameter"), self.diametroMembrana.str)
        if not self.kwargs["diametroMembrana"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Area per bag"), self.areaMembrana.str)
        if not self.kwargs["areaMembrana"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Speed"), self.Vgas.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Surface"), self.floorArea.str)+os.linesep

        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep
        if self.entrada.solido.diametros:
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
            if len(self.entrada.solido.diametros)>=1:
                txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Partial Efficiency")
                if not any(self.kwargs["rendimientos"]):
                    txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")

                txt+=os.linesep+"%10s %8s %10s %10s %10s" %("Dp("+Length.text("ParticleDiameter")+")", "Xi", "ŋi", "Xis", "Xig")+os.linesep
                for i in range(len(self.rendimiento_parcial)):
                    txt+="%10.4f %10.4f %10.4f %10.4f %10.4f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "metodo"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Filter Number"), "num_filtros", int),
                (QApplication.translate("pychemqt", "Operation Time"), "tiempo", Time),
                (QApplication.translate("pychemqt", "Cloth resistence"), "resistenciaFiltro", Dimensionless),
                (QApplication.translate("pychemqt", "Cake resistence"), "resistenciaTorta", Dimensionless),
                (QApplication.translate("pychemqt", "Cells cleaned"), "limpieza", int),
                (QApplication.translate("pychemqt", "Bags per cell"), "membranasFiltro", int),
                (QApplication.translate("pychemqt", "Bag diameter"), "diametroMembrana", Length),
                (QApplication.translate("pychemqt", "Area per bag"), "areaMembrana", int),
                (QApplication.translate("pychemqt", "Speed"), "Vgas", Speed),
                (QApplication.translate("pychemqt", "Surface"), "floorArea", Area),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length)]
        return list


class ElectricPrecipitator(Separador_SolidGas):
    """ Clase que define el precipitador electrostático

    Parámetros:
        entrada: Instancia de clase corriente que define el gas que fluye por la cámara
        metodo: indice que indica el tipo de cálculo a realizar
            0   -   Evaluación, conocidas las dimensiones, calcular el rendimiento del equipo
            1   -   Diseño, imponiendo el rendimiento a obtener calcular el area y potencia del carga
        potencialCarga: potencial de carga del equipo en V/m
        potencialDescarga: potencial de descarga del equipo en V/m
        area: area de recolección del equipo, m2
        epsilon: constante dieléctrica relativa del material
        rendimientoAdmisible: rendimiento exigido al equipo (modo diseño)
        deltaP: parámetro opcional que indica la perdida de presión del equipo

    >>> diametros=[17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> solido=Solid(caudalSolido=[12426.28], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    >>> precipitador=ElectricPrecipitator(entrada=corriente, metodo=1, rendimientoAdmisible=0.9)
    >>> print precipitador.area, precipitador.rendimiento
    2389.41624164 0.9
    """
    title=QApplication.translate("pychemqt", "Electrostatic precipitator")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "potencialCarga": 0.0,
                    "potencialDescarga": 0.0,
                    "area": 0.0,
                    "epsilon": 0.0,
                    "rendimientoAdmisible": 0.0,
                    "deltaP": 0.0}
    kwargsInput=("entrada", )
    kwargsValue=("potencialCarga", "potencialDescarga", "area", "rendimientoAdmisible", "epsilon", "deltaP")
    kwargsList=("metodo", )
    calculateValue=("areaCalculada", "rendimiento")

    TEXT_TIPO=[QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
                        QApplication.translate("pychemqt", "Design: Calculate dimensions to fit a requerid efficiency")]


    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return

        if self.kwargs["metodo"]==0 and not self.kwargs["area"]:
            self.msg=QApplication.translate("pychemqt", "undefined area")
            self.status=0
            return
        elif self.kwargs["metodo"]==1 and not self.kwargs["rendimientoAdmisible"]:
            self.msg=QApplication.translate("pychemqt", "undefined efficiency")
            self.status=0
            return


        if not self.kwargs["potencialCarga"]:
            self.msg=QApplication.translate("pychemqt", "using default charging field")
            self.status=3
            return True
        if not self.kwargs["potencialDescarga"]:
            self.msg=QApplication.translate("pychemqt", "using default collecting field")
            self.status=3
            return True
        if not self.kwargs["epsilon"]:
            self.msg=QApplication.translate("pychemqt", "using default dielectric constant")
            self.status=3
            return True

        self.msg=""
        self.status=1
        return True


    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.areaCalculada=Area(self.kwargs["area"])
        self.deltaP=Pressure(self.kwargs["deltaP"])

        if self.kwargs["potencialCarga"]:
            self.potencialCarga=PotencialElectric(self.kwargs["potencialCarga"])
        else:
            self.potencialCarga=PotencialElectric(24000)
        if self.kwargs["potencialDescarga"]:
            self.potencialDescarga=PotencialElectric(self.kwargs["potencialDescarga"])
        else:
            self.potencialDescarga=PotencialElectric(24000)
        if self.kwargs["epsilon"]:
            self.epsilon=Dimensionless(self.kwargs["epsilon"])
        else:
            self.epsilon=Dimensionless(4.)

        if self.kwargs["metodo"]==0:
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.areaCalculada)
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)
        else:
            def f(area):
                rendimientos=self.calcularRendimientos_parciales(area)
                return self.calcularRendimiento(rendimientos)-self.kwargs["rendimientoAdmisible"]

            self.areaCalculada=Area(fsolve(f, 100)[0])
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.areaCalculada)
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)

        self.CalcularSalidas()


    def calcularRendimientos_parciales(self, A):
        """Método que calcula el rendimiento de separación, desglosado según la distribución de sólidos"""
        rendimiento_parcial=[]
        for dp in self.entrada.solido.diametros:
            if dp <= 1e-6:
                q=dp*e*1e8
            else:
                q=pi*epsilon_0*self.potencialCarga*dp**2*(1+2*(self.epsilon-1.)/(self.epsilon+2.))
            U=q*self.potencialDescarga/(3*pi*dp*self.entrada.Gas.mu)
            rendimiento_parcial.append(1-exp(-U*A/self.entrada.Q))
        return rendimiento_parcial


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["metodo"]].split(":")[0])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Charging field"), self.potencialCarga.str)
        if self.kwargs["potencialCarga"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Collecting field"), self.potencialDescarga.str)
        if self.kwargs["potencialDescarga"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Dielectric constant"), self.epsilon.str)
        if self.kwargs["epsilon"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Area"), self.areaCalculada.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
        if len(self.entrada.solido.diametros)>=1:
            txt+="%10s %10s %10s %10s %10s" %("Dp(µm)", "Xi", "ηi", "Xis", "Xig")+os.linesep
            for i in range(len(self.rendimiento_parcial)):
                txt+="%10.1f %10.2f %10.3f %10.3f %10.3f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "metodo"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Charging field"), "potencialCarga", PotencialElectric),
                (QApplication.translate("pychemqt", "Collecting field"), "potencialDescarga", PotencialElectric),
                (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", Dimensionless),
                (QApplication.translate("pychemqt", "Area"), "areaCalculada", Area),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length)]
        return list


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    from lib.solids import Solid
#    dm = [17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4,
#          60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
#    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#                  0.05, 0.03, 0.02]
#    solido = Solid(caudalSolido=[0.1], distribucion_diametro=dm,
#                   distribucion_fraccion=fracciones, solids=[638])
#    kw = {"ids": [475], "fraccionMolar": [1.], "MEoS": True}
#    corriente = Corriente(T=300, P=1e5, caudalMasico=1, solido=solido, **kw)
#

#
#    filtro=Baghouse(entrada=corriente, metodo=1, num_filtros=4, tiempo=3600, deltaP=0.1)
#    print((filtro.floorArea, filtro.Vgas.ftmin))


#    precipitador=ElectricPrecipitator(entrada=corriente, metodo=1, rendimientoAdmisible=0.9)
#    print precipitador.msg, precipitador.status
#    print precipitador.areaCalculada, precipitador.rendimiento
