#!/usr/bin/python
# -*- coding: utf-8 -*-


###############################################################################
# Module to define chemical reaction functionality
###############################################################################


from math import exp, log
import sqlite3

from numpy import polyval
from scipy.optimize import fsolve
from PyQt4.QtGui import QApplication

from lib import unidades
from lib.sql import databank_name


class Reaction(object):
    """Chemical reaction object"""

    status = 0
    msg = QApplication.translate("pychemqt", "undefined")
    error = 0

    kwargs = {"componentes": [],
              "coeficientes": [],
              "tipo": 0,
              "fase": 0,
              "key": 0,
              "base": 0,
              "customHr": False,
              "Hr": 0.0,
              "formula": False,
              "estequiometria": [],
              "keq": None}
    kwargsValue = ("Hr",)
    kwargsList = ("tipo", "fase", "key", "base")
    kwargsCheck = ("customHr", "formula")
    calculateValue = ("DeltaP", "DeltaP_f", "DeltaP_ac", "DeltaP_h",
                      "DeltaP_v", "DeltaP_100ft", "V", "f", "Re", "Tout")

    TEXT_TYPE = [QApplication.translate("pychemqt", "Estequiometric"),
                 QApplication.translate("pychemqt", "Equilibrium"),
                 QApplication.translate("pychemqt", "Kinetic"),
                 QApplication.translate("pychemqt", "Catalitic")]
    TEXT_PHASE = [QApplication.translate("pychemqt", "Global"),
                  QApplication.translate("pychemqt", "Liquid"),
                  QApplication.translate("pychemqt", "Gas")]
    TEXT_BASE = [QApplication.translate("pychemqt", "Mole"),
                 QApplication.translate("pychemqt", "Mass"),
                 QApplication.translate("pychemqt", "Partial pressure")]

    def __init__(self, **kwargs):
        """constructor, kwargs keys can be:
            componentes: array with index of reaction components
            coeficientes: array with stequiometric coefficient for each component
            tipo: Kind of reaction
                0   -   Stequiometric, without equilibrium or kinetic calculations
                1   -   Equilibrium, without kinetic calculation
                2   -   Equilibrium by minimization of Gibbs free energy
                3   -   Kinetic
                4   -   Catalytic
            fase: Phase where reaction work
                0   -   Global
                1   -   Liquid
                2   -   Gas
            key: Index of key component
            base
                0   -   Mol
                1   -   Mass
                2   -   Partial pressure
            Hr: Heat of reaction, calculate from heat of formation if no input
            formula: boolean to show compound names in formules
            estequiometria: array with conversion reaction values
            keq: equilibrium constant
                -it is float if it don't depend with temperature
                -it is array if it depends with temperature
        """
        self.kwargs = Reaction.kwargs.copy()
        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        oldkwargs = self.kwargs.copy()
        self.kwargs.update(kwargs)

        if oldkwargs != self.kwargs and self.isCalculable:
            self.calculo()

    @property
    def isCalculable(self):
        if not self.kwargs["componentes"]:
            self.msg = QApplication.translate("pychemqt", "undefined components")
            self.status = 0
            return
        if not self.kwargs["coeficientes"]:
            self.msg = QApplication.translate("pychemqt", "undefined stequiometric")
            self.status = 0
            return

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        self.componentes = self.kwargs["componentes"]
        self.coeficientes = self.kwargs["coeficientes"]

        self.tipo = self.kwargs["tipo"]
        self.base = self.kwargs["base"]
        self.fase = self.kwargs["fase"]
        self.estequiometria = self.kwargs["estequiometria"]
        self.calor = self.kwargs["Hr"]
        self.formulas = self.kwargs["formula"]
        self.keq = self.kwargs["keq"]

        databank = sqlite3.connect(databank_name).cursor()
        databank.execute("select nombre, peso_molecular, formula, \
                         calor_formacion_gas from compuestos where id  IN \
                         %s" % str(tuple(self.componentes)))
        nombre = []
        peso_molecular = []
        formula = []
        calor_reaccion = 0
        check_estequiometria = 0
        for i, compuesto in enumerate(databank):
            nombre.append(compuesto[0])
            peso_molecular.append(compuesto[1])
            formula.append(compuesto[2])
            calor_reaccion += compuesto[3]*self.coeficientes[i]
            check_estequiometria += self.coeficientes[i]*compuesto[1]
        self.nombre = nombre
        self.peso_molecular = peso_molecular
        self.formula = formula
        if self.calor:
            self.Hr = self.kwargs.get("Hr", 0)
        else:
            self.Hr = unidades.MolarEnthalpy(calor_reaccion/abs(
                self.coeficientes[self.base]), "Jkmol")
        self.error = round(check_estequiometria, 1)
        self.state = self.error == 0
        self.text = self._txt(self.formulas)

    def conversion(self, corriente, T):
        """Calculate reaction conversion
        corriente: Corriente instance for reaction
        T: Temperature of reaction"""
        if self.tipo == 0:
            # Material balance without equilibrium or kinetics considerations
            alfa = polyval(self.estequiometria, T)
            if alfa < 0:
                alfa = 0
            elif alfa > 1:
                alfa = 1

        elif self.tipo == 1:
            # Chemical equilibrium without kinetics
            if isinstance(self.keq, list):
                A, B, C, D, E, F, G, H = self.keq
                keq = exp(A+B/T+C*log(T)+D*T+E*T**2+F*T**3+G*T**4+H*T**5)
            else:
                keq = self.keq

            def f(alfa):
                conc_out = [
                    (corriente.caudalunitariomolar[i]+alfa*self.coeficientes[i])
                    / corriente.Q.m3h for i in range(len(self.componentes))]
                productorio = 1
                for i in range(len(self.componentes)):
                    productorio *= conc_out[i]**self.coeficientes[i]
                return keq-productorio

            alfa = fsolve(f, 0.5)
            print alfa, f(alfa)

        avance = alfa*self.coeficientes[self.base]*corriente.caudalunitariomolar[self.base]
        Q_out = [corriente.caudalunitariomolar[i]+avance*self.coeficientes[i] /
                 self.coeficientes[self.base] for i in range(len(self.componentes))]
        minimo = min(Q_out)
        if minimo < 0:
            # The key component is not correct, redo the result
            indice = Q_out.index(minimo)
            avance = self.coeficientes[indice]*corriente.caudalunitariomolar[indice]
            Q_out = [corriente.caudalunitariomolar[i]+avance*self.coeficientes[i] /
                     self.coeficientes[indice] for i in range(len(self.componentes))]
            h = unidades.Power(self.Hr*self.coeficientes[self.base] /
                               self.coeficientes[indice]*avance, "Jh")
        else:
            h = unidades.Power(self.Hr*avance, "Jh")

        caudal = sum(Q_out)
        fraccion = [caudal_i/caudal for caudal_i in Q_out]
        return fraccion, h


#    def cinetica(self, tipo, Ko, Ei):
#        """Método que define la velocidad de reacción"""
#
#
    def _txt(self, nombre=False):
        """Function to get text representation for reaction"""
        if nombre:
            txt = self.nombre
        else:
            txt = self.formula

        reactivos = []
        productos = []
        for i in range(len(self.componentes)):
            if self.coeficientes[i] == int(self.coeficientes[i]):
                self.coeficientes[i] = int(self.coeficientes[i])
            if self.coeficientes[i] < -1:
                reactivos.append(str(-self.coeficientes[i])+txt[i])
            elif self.coeficientes[i] == -1:
                reactivos.append(txt[i])
            elif -1 < self.coeficientes[i] < 0:
                reactivos.append(str(-self.coeficientes[i])+txt[i])
            elif 0 < self.coeficientes[i] < 1:
                productos.append(str(self.coeficientes[i])+txt[i])
            elif self.coeficientes[i] == 1:
                productos.append(txt[i])
            elif self.coeficientes[i] > 1:
                productos.append(str(self.coeficientes[i])+txt[i])
        return " + ".join(reactivos)+" ---> "+" + ".join(productos)

    def __repr__(self):
        if self.status:
            eq = self._txt()
            return eq + "   " + "Hr= %0.4e Jkmol" % self.Hr
        else:
            return str(self.msg)


if __name__ == "__main__":
#    from lib.corriente import Corriente, Mezcla
#    mezcla=Corriente(300, 1, 1000, Mezcla([1, 46, 47, 62], [0.03, 0.01, 0.96, 0]))
#    reaccion=Reaction([1, 46, 47, 62], [-2, 0, -1, 2], base=2, estequiometria=[0, 0, 0.5])
#    reaccion.conversion(mezcla)
#    print reaccion

    reaccion = Reaction(componentes=[1, 47, 62], coeficientes=[-1, -0.5, 1],
                        base=0, estequiometria=[0, 0, 0.5])
    print reaccion
