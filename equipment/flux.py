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
#   Library for flux equipment definition
#   * Divider: Simple divider equipment
#   * Mixer: Simple mixer equipment
#   * Valve: Simple valve equipment
###############################################################################


import os

from PyQt5.QtWidgets import QApplication
from scipy.optimize import fsolve

from lib.corriente import Corriente
from lib import unidades
from equipment.parents import equipment


class Divider(equipment):
    """Class to define a simple divider equipment, only splliting input stream

    Parameters:
        entrada: Corriente instance to define de input stream
        salidas: Number of output streams
        criterio: split ratio
            0   -   fractions, ratio of input stream in each output stream
            1   -   molar flow. Overwrite the input stream molar flow.
        fracciones: array with the split values, it depend of criterio,
                fractions of input stream or molar flow in output streams
        deltaP: Pressure loss in equipment, optional

    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, ids=[62],
                       fraccionMasica=[1])
    >>> divisor=Divider(entrada=agua, criterio=0, salidas=3,
                        split=[0.7, 0.2, 0.1])
    >>> print(divisor.entrada.caudalmasico.kgh)
    3600.0
    >>> for salida in divisor.salida: print("%0.1f" % salida.caudalmasico.kgh)
    2520.0
    720.0
    360.0
    """

    title = QApplication.translate("pychemqt", "Divider")
    kwargs = {"entrada": None,
              "criterio": 0,
              "salidas": 0,
              "split": [],
              "deltaP": 0.0}
    kwargsInput = ("entrada", )
    kwargsList = ("criterio", )
    kwargsValue = ("deltaP", )

    TEXT_CRITERIO = [
        QApplication.translate("pychemqt", "Flux ratio"),
        QApplication.translate("pychemqt", "Flowrate (overwrite input flow)")]

    @property
    def isCalculable(self):
        """Check equipment input parameter

        Mandatory parameter:
            entrada, split
        Incompatibilities:
            salidas must be equal to len of split
        """
        # Auto select split for trivial divider with one output stream
        if self.kwargs["salidas"] == 1:
            self.kwargs["split"] = [1.]

        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
        elif not self.kwargs["split"]:
            self.msg = QApplication.translate("pychemqt",
                                              "undefined split fraction")
            self.status = 0
        elif self.kwargs["salidas"] != len(self.kwargs["split"]):
            self.msg = QApplication.translate("pychemqt",
                                              "incomplete split fraction")
            self.status = 0
        else:
            self.status = 1
            self.msg = ""
            return True

    def calculo(self):
        """Calculate procedure, only a mass balance"""
        self.entrada = self.kwargs["entrada"]
        self.criterio = self.kwargs["criterio"]
        self.split = [unidades.Dimensionless(i) for i in self.kwargs["split"]]
        self.deltaP = unidades.DeltaP(self.kwargs["deltaP"])
        # Normalize fractions
        if self.criterio == 0:
            if sum(self.split) != 1:
                fracciones_normalizadas = []
                total = sum(self.split)
                for i in self.split:
                    fracciones_normalizadas.append(i/total)
                self.split = fracciones_normalizadas

        # Calculate output stream
        self.salida = []
        if self.criterio == 0:
            for i in self.split:
                self.salida.append(self.entrada.clone(
                    P=self.entrada.P-self.deltaP, split=i))
        else:
            self.entrada = self.entrada.clone(caudalmasico=sum(self.split))
            for i in self.split:
                self.salida.append(self.entrada.clone(
                    P=self.entrada.P-self.deltaP, split=i))

        # Calculate other properties
        self.inputMolarFlow = self.entrada.caudalmolar
        self.inputMassFlow = self.entrada.caudalmasico
        self.inputVolFlow = self.entrada.Q
        self.inputT = self.entrada.T
        self.inputP = self.entrada.P
        self.output = unidades.Dimensionless(self.kwargs["salidas"])

    def propTxt(self):
        """Text format for report"""
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(0)
        for i, salida in enumerate(self.salida):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Output stream")+str(i),
                salida.caudalmolar.str)+os.linesep
        txt += os.linesep
        txt += self.propertiesToText(1)
        for i, salida in enumerate(self.salida):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Output stream")+str(i),
                salida.caudalmasico.str)+os.linesep
        txt += os.linesep
        txt += self.propertiesToText(2)
        for i, salida in enumerate(self.salida):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Output stream")+str(i),
                salida.Q.str)+os.linesep
        txt += os.linesep
        txt += self.propertiesToText(7)
        return txt

    @classmethod
    def propertiesEquipment(cls):
        """Properties availables to show in report"""
        l = [(QApplication.translate("pychemqt", "Feed Molar Flow"),
              "inputMolarFlow", unidades.MolarFlow),
             (QApplication.translate("pychemqt", "Feed Mass Flow"),
              "inputMassFlow", unidades.MassFlow),
             (QApplication.translate("pychemqt", "Feed Volumetric Flow"),
              "inputVolFlow", unidades.VolFlow),
             (QApplication.translate("pychemqt", "Feed Temperature"),
              "inputT", unidades.Temperature),
             (QApplication.translate("pychemqt", "Feed Pressure"), "inputP",
              unidades.Pressure),
             (QApplication.translate("pychemqt", "Number of Product Streams"),
              "output", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Split Fractions"), "split",
              unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP",
              unidades.DeltaP)]
        return l

    def propertiesListTitle(self, index):
        """Titles of properties of type list"""
        lista = []
        for i in range(self.kwargs["salidas"]):
            lista.append("Output %i" % (i+1))
        return lista

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["criterio"] = self.criterio
        state["split"] = self.split
        state["deltaP"] = self.deltaP
        state["inputMolarFlow"] = self.inputMolarFlow
        state["inputMassFlow"] = self.inputMassFlow
        state["inputVolFlow"] = self.inputVolFlow
        state["inputT"] = self.inputT
        state["inputP"] = self.inputP
        state["output"] = self.output

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.criterio = state["criterio"]
        self.split = (unidades.Dimensionless(x) for x in state["split"])
        self.deltaP = unidades.DeltaP(state["deltaP"])
        self.inputMolarFlow = unidades.MolarFlow(state["inputMolarFlow"])
        self.inputMassFlow = unidades.MassFlow(state["inputMassFlow"])
        self.inputVolFlow = unidades.VolFlow(state["inputVolFlow"])
        self.inputT = unidades.Temperature(state["inputT"])
        self.inputP = unidades.Pressure(state["inputP"])
        self.output = unidades.Dimensionless(state["output"])
        self.salida = [None]*self.kwargs["salidas"]

    def ajustState(self, stream):
        print(self.kwargs)


class Mixer(equipment):
    """Class to define a simple stream mixer, only apply a mass balance and
    specified the output pressure

    Parameters
        n_entradas: Input number
        entrada: Corriente instance to define a input stream to equipment,
            if use a list of them, it fix all the input stream
        id_entrada: id of defined input
        criterio: mode of define output pressure
            0 - Minimum pressure of input streams
            1 - Mean pressure of input streams
            2 - Specified pressure, must be input in Pout parameter
        Pout: Output pressure

    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, ids=[62],
                       fraccionMasica=[1.])
    >>> agua2=Corriente(T=350, P=101325, caudalMasico=5, ids=[62],
                        fraccionMasica=[1.])
    >>> mezclador=Mixer(entrada=[agua, agua2], criterio=0)
    >>> print(mezclador.salida[0].T)
    341.6818851215901
    """
    title = QApplication.translate("pychemqt", "Mixer")
    help = ""
    kwargs = {"entrada": [],
              "id_entrada": 0,
              "criterio": 0,
              "Pout": 0.0}
    kwargs_forbidden = ["entrada", "id_entrada"]
    kwargsValue = ("Pout", )
    kwargsList = ("criterio", )

    TEXT_METODO = [
        QApplication.translate("pychemqt", "Inputs minimum pressure"),
        QApplication.translate("pychemqt", "Inputs mean pressure"),
        QApplication.translate("pychemqt", "Custom")]

    @property
    def isCalculable(self):
        """Check equipment input parameter

        Mandatory parameter:
            entrada
            Pout if criterio is setted to 2
        Incompatibilities:
            All corriente instance in entrada must be defined fine
            salidas must be equal to len of split
        Warnings:
            Some input stream have a warning status definitions
        """
        if self.kwargs["criterio"] == 2 and not self.kwargs["Pout"]:
            self.msg = QApplication.translate("pychemqt",
                                              "pressure output not defined")
            self.status = 0
        elif not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
        elif sum([s.status for s in self.kwargs["entrada"]]) == 0:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
        elif len(self.kwargs["entrada"]) != sum(
                [s.status for s in self.kwargs["entrada"]]):
            self.msg = QApplication.translate(
                "pychemqt", "some input stream isn't defined")
            self.status = 3
            return True
        else:
            self.msg = ""
            self.status = 1
            return True

    def cleanOldValues(self, **kwargs):
        """Clean incompatible kwargs parameters
            New defined entrada delete old entrada
            Entrada as list not need the id_entrada
        """
        if "entrada" in kwargs:
            if isinstance(kwargs["entrada"], list):
                kwargs["id_entrada"] = None
            else:
                corriente = kwargs["entrada"]
                kwargs["entrada"] = self.kwargs["entrada"][:]
                while len(kwargs["entrada"]) < kwargs["id_entrada"]+1:
                    kwargs["entrada"].append(Corriente())
                kwargs["entrada"][kwargs["id_entrada"]] = corriente
                kwargs["id_entrada"] = None
        self.kwargs.update(kwargs)

    def calculo(self):
        """Calculate procedure, only a mass balance"""
        self.entrada = self.kwargs["entrada"]
        self.criterio = self.kwargs["criterio"]

        # Define output pressure
        if self.criterio == 2:
            Pout = self.kwargs["Pout"]
        else:
            lst = []
            for entrada in self.entrada:
                if entrada.status:
                    lst.append(entrada.P)
            if self.criterio == 0:
                Pout = min(lst)
            else:
                Pout = sum(lst, 0.0) / len(lst)
        self.Pout = unidades.Pressure(Pout)

        # Heat balance for calculate output temperature
        h_in = 0
        To = 0
        massUnitFlow = [0]*len(self.entrada[0].fraccion)
        for entrada in self.entrada:
            if entrada.status:
                h_in += entrada.h
                To += entrada.T*entrada.caudalmasico
                for i, caudal in enumerate(entrada.caudalunitariomasico):
                    massUnitFlow[i] += caudal
        To /= sum(massUnitFlow)

        def f(T):
            output = Corriente(T=T, P=self.Pout,
                               caudalUnitarioMasico=massUnitFlow)
            return output.h-h_in
        T = fsolve(f, To)[0]

        # TODO: Add solid mixer capability
        if self.entrada[0].solido:
            pass

        salida = Corriente(T=T, P=self.Pout, caudalUnitarioMasico=massUnitFlow)
        self.salida = [salida]

        # Calculate other properties
        self.outT = salida.T
        self.outX = salida.x
        self.outMolarFlow = salida.caudalmolar
        self.outMassFlow = salida.caudalmasico
        self.outVolFlow = salida.Q

    def propTxt(self):
        """Text format for report"""
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(3))

        txt += os.linesep
        txt += self.propertiesToText(3)
        for i, entrada in enumerate(self.kwargs["entrada"]):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Input stream")+str(i),
                entrada.caudalmolar.str)+os.linesep
        txt += os.linesep
        txt += self.propertiesToText(4)
        for i, entrada in enumerate(self.kwargs["entrada"]):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Input stream")+str(i),
                entrada.caudalmasico.str)+os.linesep
        txt += os.linesep
        txt += self.propertiesToText(5)
        for i, entrada in enumerate(self.kwargs["entrada"]):
            txt += "%-25s\t%s" % (
                QApplication.translate("pychemqt", "Input stream")+str(i),
                entrada.Q.str)+os.linesep

        txt += os.linesep+"#"
        txt += QApplication.translate("pychemqt", "Output Molar Composition")
        txt += os.linesep
        for comp, x in zip(self.salida[0].componente, self.salida[0].fraccion):
            txt += "%-25s\t %0.4f" % (comp.nombre, x)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        """Properties availables to show in report"""
        l = [(QApplication.translate("pychemqt", "Output Temperature"),
              "outT", unidades.Temperature),
             (QApplication.translate("pychemqt", "Output Pressure"),
              "Pout", unidades.Pressure),
             (QApplication.translate("pychemqt", "Output vapor fraction"),
              "outX", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Output Molar Flow"),
              "outMolarFlow", unidades.MolarFlow),
             (QApplication.translate("pychemqt", "Output Mass Flow"),
              "outMassFlow", unidades.MassFlow),
             (QApplication.translate("pychemqt", "Output Volumetric Flow"),
              "outVolFlow", unidades.VolFlow)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["criterio"] = self.criterio
        state["Pout"] = self.Pout
        state["outT"] = self.outT
        state["outX"] = self.outX
        state["outMolarFlow"] = self.outMolarFlow
        state["outMassFlow"] = self.outMassFlow
        state["outVolFlow"] = self.outVolFlow

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.criterio = state["criterio"]
        self.Pout = unidades.Pressure(state["Pout"])
        self.outT = unidades.Temperature(state["outT"])
        self.outX = unidades.Dimensionless(state["outX"])
        self.outMolarFlow = unidades.MolarFlow(state["outMolarFlow"])
        self.outMassFlow = unidades.MassFlow(state["outMassFlow"])
        self.outVolFlow = unidades.VolFlow(state["outVolFlow"])
        self.salida = [None]


class Valve(equipment):
    """Class to define a simple valve, only output pressure calculation

    Parameters:
        entrada: Instance of Corriente to define the valve input stream
        off: state of valve
            0 - open
            1 - open partially
            2 - closed
        Pout: Output pressure
        DeltaP: Pressure loss in equipment
        Dew: Dew point temperature
        Bubble: Bubble point temperature

    >>> agua=Corriente(T=300, P=202650, caudalMasico=1, ids=[62],
                       fraccionMasica=[1])
    >>> valvula=Valve(entrada=agua, Pout=101325, off=1)
    >>> print(agua.P.atm, valvula.salida[0].P.atm)
    2.0 1.0
    >>> valvula(DeltaP=2650)
    >>> print(valvula.salida[0].P.atm)
    1.9738465334320257
    >>> valvula(off=0)
    >>> print(valvula.salida[0].P.atm)
    2.0
    """

    title = QApplication.translate("pychemqt", "Valve")
    help = ""
    kwargs = {"entrada": None,
              "off": 0,
              "Pout": 0.0,
              "DeltaP": 0.0,
              "Dew": 0.0,
              "Bubble": 0.0}
    kwargsInput = ("entrada", )
    kwargsValue = ("Pout", "DeltaP", "Dew", "Bubble")
    kwargsList = ("off", )

    TEXT_WORKING = [QApplication.translate("pychemqt", "Totally open"),
                    QApplication.translate("pychemqt", "Partially open"),
                    QApplication.translate("pychemqt", "Close")]

    @property
    def isCalculable(self):
        """Check equipment input parameter

        Mandatory parameter:
            entrada if valve is working, open or partially open
            Pout (DeltaP, Dew, Bubble) if valve is partially open
        """
        if self.kwargs["off"] == 1:
            if not self.kwargs["entrada"]:
                self.msg = QApplication.translate("pychemqt",
                                                  "undefined input")
                self.status = 0
            elif self.kwargs["Pout"] or self.kwargs["DeltaP"] or \
                    self.kwargs["Dew"] or self.kwargs["Bubble"]:
                self.status = 1
                self.msg = ""
                return True
            else:
                self.msg = QApplication.translate("pychemqt",
                                                  "undefined exit condition")
                self.status = 0
        elif self.kwargs["off"] == 2:
            self.msg = ""
            self.status = 1
            return True
        else:
            if not self.kwargs["entrada"]:
                self.msg = QApplication.translate("pychemqt",
                                                  "undefined input")
                self.status = 0
            else:
                self.msg = ""
                self.status = 1
                return True

    def cleanOldValues(self, **kwargs):
        """Clean incompatible kwargs parameters
            Each output pressure definition disabled any old definition
        """
        if "Pout" in kwargs:
            self.kwargs["DeltaP"] = 0
            self.kwargs["Dew"] = 0
            self.kwargs["Bubble"] = 0
        elif "DeltaP" in kwargs:
            self.kwargs["Pout"] = 0
            self.kwargs["Dew"] = 0
            self.kwargs["Bubble"] = 0
        elif "Dew" in kwargs:
            self.kwargs["Pout"] = 0
            self.kwargs["DeltaP"] = 0
            self.kwargs["Bubble"] = 0
        elif "Bubble" in kwargs:
            self.kwargs["Pout"] = 0
            self.kwargs["DeltaP"] = 0
            self.kwargs["Dew"] = 0
        self.kwargs.update(kwargs)

    def calculo(self):
        """Calculate procedure, only apply pressure loss method chosen"""
        self.entrada = self.kwargs["entrada"]
        Pout = self.kwargs["Pout"]
        DeltaP = self.kwargs["DeltaP"]
        Dew = self.kwargs["Dew"]
        Bubble = self.kwargs["Bubble"]

        if self.kwargs["off"] == 1:
            if Pout:
                self.Pout = unidades.Pressure(Pout)
            elif DeltaP:
                self.Pout = unidades.Pressure(self.entrada.P-DeltaP)
            elif Dew:
                corriente = self.entrada.clone(T=Dew)
                self.Pout = corriente.eos._Dew_P()
            elif Bubble:
                corriente = self.entrada.clone(T=Bubble)
                self.Pout = corriente.eos._Bubble_P()
            self.salida = [self.entrada.clone(P=self.Pout)]

        elif self.kwargs["off"] == 2:
            self.entrada = Corriente()
            self.salida = [self.entrada]

        else:
            self.salida = [self.entrada]
            self.Pout = unidades.Pressure(self.salida[0].P)

        # Calculate other properties
        self.outT = self.salida[0].T
        self.outX = self.salida[0].x

    def propTxt(self):
        """Text format for report"""
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#"+os.linesep
        txt += self.propertiesToText(range(4))
        return txt

    @classmethod
    def propertiesEquipment(cls):
        """Properties availables to show in report"""
        l = [(QApplication.translate("pychemqt", "Output Temperature"),
              "outT", unidades.Temperature),
             (QApplication.translate("pychemqt", "Output Pressure"),
              "Pout", unidades.Pressure),
             (QApplication.translate("pychemqt", "Output vapor fraction"),
              "outX", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Working Condition"),
              ("TEXT_WORKING", "off"), str)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Pout"] = self.Pout
        state["outT"] = self.outT
        state["outX"] = self.outX

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Pout = unidades.Pressure(state["Pout"])
        self.outT = unidades.Temperature(state["outT"])
        self.outX = unidades.Dimensionless(state["outX"])
        self.salida = [None]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
