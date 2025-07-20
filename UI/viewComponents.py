#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Module to view/edit component in database
#   -View_Component: Dialog to view/edit compounds of database
#   -DIPPR_widget: DIPPR equation composite widget
#   -Parametric_widget: Parametric equation composite widget
###############################################################################


from functools import partial
import os
import tempfile

from numpy import array, linspace
from scipy import optimize
from tools.qt import QtCore, QtGui, QtWidgets

from lib.plot import PlotDialog
from lib.compuestos import (Componente, MuL_Parametric, Pv_Antoine, Pv_Wagner,
                            Tension_Parametric, Henry, DIPPR)
from lib import unidades, sql
from lib.config import IMAGE_PATH, Preferences
from lib.EoS.Cubic.MSRK import dat as MSRKdat
from lib.EoS.Cubic.PRMelhem import dat as Melhemdat
from lib.EoS.Cubic.SRKAPI import S1, S2
from lib.EoS.Cubic.PRMathiasCopeman import dat as Mathiasdat
from lib.EoS.Cubic.PRSV import dat as PRSVdat
from lib.EoS.Cubic.PRZV import dat as PRZVdat
from lib.EoS.Cubic.PRAlmeida import dat as Almeidadat
from lib.EoS.Cubic.PT import dat as PTdat
from lib.EoS.Cubic.TB import dat as TBdat
from lib.EoS.Cubic.TBS import dat as TBSdat
from lib.openbabel import imageFromSmile

from UI.inputTable import InputTableDialog, eqDIPPR
from UI.delegate import SpinEditor
from UI.widgets import Entrada_con_unidades, Tabla, okToContinue, QLabelMath


class DIPPR_widget(QtWidgets.QGroupBox):
    """Composite widget to edit/view a DIPPR coefficients"""
    valueChanged = QtCore.pyqtSignal()
    dict_proptex = {
        "rhoS": r"\rho_s",
        "rhoL": r"\rho_l",
        "Hv": r"H_v",
        "Pv": r"P_v",
        "cpS": r"Cp_s",
        "cpL": r"Cp_l",
        "cpG": r"Cp_g",
        "muL": r"\mu_l",
        "muG": r"\mu_g",
        "kL": r"\lambda_l",
        "kG": r"\lambda_v",
        "sigma": r"\sigma"}
    latex = {
        1: r"A + BT + CT^2 + DT^3 + ET^4",
        2: r"e^{A + BT + C\ln(T) + DT^E}",
        3: r"\frac{AT^B}{1+\frac{C}{T}+\frac{D}{T^2}}",
        4: r"A+B e^{-\frac{C}{T^D}}",
        5: r"A + BT + CT^3 + DT^8 + ET^9",
        6: r"\frac{A}{B^{1+(1-T/C)^D}}",
        7: r"A(1-T_r)^{B+CT_r+dT_r^2+ET_r^3}",
        8: r"A+B\frac{C/T}{\sinh(C/T)}+D\frac{E/T}{\cosh(E/T)}",
        9: r"A^2/Tr+B-2ACTr-ADTr^2-C^2Tr^3/3-CDTr^4/2-D^2Tr^5/5"}

    def __init__(self, title, unit, id=0, prop="", array=None, parent=None):
        """Constructor of widget

        Parameters
        ----------
        title : string
            Title so show of qgroupbox
        unit : string
            Aditional string with unit representation
        id : integer
            Index of compound in database to show
        prop : string
            Latex representation of property to show in matplotlib
        array : list
            List of DIPPR equation representation in format
            [eq, A, B, C, D, E, Tmin, Tmax]
        """
        super().__init__(title, parent)
        self.prop = prop

        self.proptex = self.dict_proptex[prop]

        self.unit = unit
        self.parent = parent
        self.t = []
        self.data = []

        lyt = QtWidgets.QGridLayout(self)
        self.btnFit = QtWidgets.QToolButton()
        self.btnFit.setToolTip(
            self.tr("Fit parameters from experimental data"))
        self.btnFit.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "fit.png"))))
        self.btnFit.setIconSize(QtCore.QSize(32, 32))
        self.btnFit.setFixedSize(QtCore.QSize(32, 32))
        self.btnFit.clicked.connect(self.fit)
        lyt.addWidget(self.btnFit, 1, 1, 2, 1)
        self.btnPlot = QtWidgets.QToolButton()
        self.btnPlot.setToolTip(self.tr("Plot equation vs temperature"))
        self.btnPlot.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "plot.png"))))
        self.btnPlot.setIconSize(QtCore.QSize(32, 32))
        self.btnPlot.setFixedSize(QtCore.QSize(32, 32))
        self.btnPlot.setEnabled(False)
        self.btnPlot.clicked.connect(self.plot)
        lyt.addWidget(self.btnPlot, 1, 2, 2, 1)
        lyt.addWidget(QtWidgets.QLabel("T<sub>min</sub>"), 4, 1)
        self.tmin = Entrada_con_unidades(unidades.Temperature, width=70)
        self.tmin.valueChanged.connect(self.valueChanged.emit)
        lyt.addWidget(self.tmin, 4, 2)
        lyt.addWidget(QtWidgets.QLabel("T<sub>max</sub>"), 5, 1)
        self.tmax = Entrada_con_unidades(unidades.Temperature, width=70)
        self.tmax.valueChanged.connect(self.valueChanged.emit)
        lyt.addWidget(self.tmax, 5, 2)

        txt = ["A", "B", "C", "D", "E"]
        self.coeff = []
        for i in range(5):
            lyt.addWidget(QtWidgets.QLabel(txt[i]), 1+i, 4)
            self.coeff.append(Entrada_con_unidades(float))
            self.coeff[-1].valueChanged.connect(self.valueChanged.emit)
            lyt.addWidget(self.coeff[-1], 1+i, 5)
        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 3)

        self.eq = eqDIPPR(1)
        lyt.addWidget(self.eq, 0, 0, 1, 6)
        self.eqformula = QLabelMath("")
        lyt.addWidget(self.eqformula, 0, 0, 1, 6)

        self.changeId(id)
        if array:
            self.fill(array)

    def changeId(self, id):
        """Common action to do when the component id change, change from user
        to system compound, editable to readOnly compound"""
        self.id = id
        if id == 0 or id > 10000:
            self.eqformula.setVisible(False)
            self.setReadOnly(False)
            self.btnFit.setEnabled(True)
        else:
            self.setReadOnly(True)
            self.eq.setVisible(True)
            self.btnFit.setEnabled(False)

    def fill(self, array):
        """Populate the widgets with the DIPPR coefficient in array in format
        [eq, A, B, C, D, E, Tmin, Tmax]"""
        if array[0] != 0:
            for valor, entrada in zip(array[1:6], self.coeff):
                entrada.setValue(valor)
            self.tmin.setValue(array[6])
            self.tmax.setValue(array[7])
            self.eq.setValue(array[0])
            latex = self.latex[array[0]]
            tex = "$%s = %s$" % (self.proptex, latex)
            self.eqformula.setTex(tex)
            self.eq.setVisible(False)
            self.eqformula.setVisible(True)
            self.btnPlot.setEnabled(True)
            self.equation = array[0]
        self.btnFit.setEnabled(self.parent.index >= 10001)

    def setReadOnly(self, bool):
        """Set widget readOnly state"""
        for entrada in self.coeff:
            entrada.setReadOnly(bool)
        self.tmin.setReadOnly(bool)
        self.tmax.setReadOnly(bool)
        self.btnFit.setDisabled(bool)

    @property
    def value(self):
        """Return the DIPPR equation parameters in the format:
        [eq, A, B, C, D, E, Tmin, Tmax]"""
        valor = [self.eq.value()]
        for entrada in self.coeff:
            elemento = entrada.value
            if elemento:
                valor.append(elemento)
            else:
                valor.append(0)
        if valor[1:] == [0]*5:
            return []
        valor.append(self.tmin.value)
        valor.append(self.tmax.value)
        return valor

    def clear(self):
        """Clear widget contents"""
        self.tmin.clear()
        self.tmax.clear()
        self.eq.clear()
        self.eqformula.clear()
        for entrada in self.coeff:
            entrada.clear()
        self.t = []
        self.data = []
        self.btnPlot.setEnabled(False)

    def formula_DIPPR(self, eq, args):
        """Calculate the formula of DIPPR equation in a latex format

        Parameters
        ----------
        eq : integer
            Index of DIPPR equation
        args : list
            Coefficient of DIPPR equation
        """
        if eq == 1:
            string = "$%s = " % self.proptex
            if args[0]:
                string += "%0.3g" % args[0]
            if args[1]:
                string += "%+0.3gT" % args[1]
            if args[2]:
                string += "%+0.3gT^2" % args[2]
            if args[3]:
                string += "%+0.3gT^3" % args[3]
            if args[4]:
                string += "%+0.3gT^4" % args[4]
            string += "$"

        elif eq == 2:
            string = "$%s = e^{" % self.proptex
            if args[0]:
                string += "%0.3g" % args[0]
            if args[1]:
                string += "%+0.3fgT" % args[1]
            if args[2]:
                string += r"%+0.3g \ln(T)" % args[2]
            if args[3]:
                string += "%+0.3gT^{%0.3g}" % (args[3], args[4])
            string += "}$"

        elif eq == 3:
            string = "$%s = %0.3g T^{\\frac{%0.3g}{1%+0.3g T" % (
                    self.proptex, args[0], args[1], args[2])
            if args[3]:
                string += "%+0.3g^2" % args[3]
            string += "}}$"

        elif eq == 4:
            string = "$%s = %0.3g%+0.3g*exp(" % (
                self.proptex, args[0], args[1])
            if args[2] < 0:
                string += "-"
            string += "%0.3g/T^{%0.3g})$" % (args[2], args[3])

        elif eq == 5:
            string = "$%s = " % self.proptex
            if args[0]:
                string += "%0.3g" % args[0]
            if args[1]:
                string += "+\\frac{%0.3g}{T}" % args[1]
            if args[2]:
                string += "+\\frac{%0.3g}{T^3}" % args[2]
            if args[3]:
                string += "+\\frac{%0.3g}{T^8}" % args[3]
            if args[4]:
                string += "+\\frac{%0.3g}{T^9}" % args[4]
            string += "$"

        elif eq == 6:
            string = "$%s = \\frac{%0.3g}{%0.3g^{1+\\left(1-T/%0.3g" % (
                    self.proptex, args[0], args[1], args[2])
            string += "\\right)^{%0.3g}}}$" % args[3]

        elif eq == 7:
            string = "$%s = %0.3g (1-T_r)^{" % (self.proptex, args[0])
            if args[1]:
                string += "%0.3g" % args[1]
            if args[2]:
                string += "%+0.3gT_r" % args[2]
            if args[3]:
                string += "%+0.3gT_r^2" % args[3]
            if args[4]:
                string += "%+0.3gT_r^3" % args[4]
            string += "}$"

        elif eq == 8:
            string = "$%s = " % self.proptex
            if args[0]:
                string += "%0.3g" % args[0]
            if args[1]:
                string += "%+0.3g\\left(\\frac{%0.3g/T}{sinh(%0.3g/T)}" % (
                    args[1], args[2], args[2])
                string += "\\right)^2"
            if args[3]:
                string += "%+0.3g\\left(\\frac{%0.3g/T}{cosh(%0.3g/T)}" % (
                        args[3], args[4], args[4])
                string += "\\right)^2"
            string += "$"

        elif eq == 9:
            string = "$%s = " % self.proptex
            if args[0]:
                string += "\\frac{%0.3g}{T_r}" % args[0]
            if args[1]:
                string += "%+0.3g" % args[1]
            if args[0] and args[2]:
                string += "%+0.3gT_r" % -2*args[0]*args[2]
            if args[0] and args[3]:
                string += "%+0.3gT_r^2" % -args[0]*args[3]
            if args[2]:
                string += "%+0.3gT_r^3" % -args[2]**2/3
            if args[2] and args[3]:
                string += "%+0.3gT_r^4" % -args[2]*args[3]/2
            if args[3]:
                string += "%+0.3gT_r^5" % -args[3]**2/5
            string += "$"

        return string

    def plot(self, r=None):
        """Plot the current DIPPR correlation equation
        r correlation coefficient, optional parameter for reuse the code in the
        fit data procedure
        """
        array = self.value
        t = list(linspace(array[-2], array[-1], 100))
        kw = {}
        kw["Tc"] = self.parent.Tc.value
        kw["M"] = self.parent.M.value
        var = [DIPPR(self.prop, ti, array[:-2], **kw) for ti in t]
        dialog = PlotDialog()
        dialog.addData(t, var)
        if self.data:
            dialog.addData(self.t, self.data, "ro")
        dialog.plot.ax.grid(True)
        dialog.plot.ax.set_title(self.title(), size="14")

        # Annotate the equation
        formula = self.formula_DIPPR(array[0], array[1:-2])
        dialog.plot.ax.annotate(formula, (0.05, 0.9), xycoords='axes fraction',
                                size="10", va="center")

        # Add correlation coefficient
        if r:
            dialog.plot.ax.annotate(
                "$r^2=%0.5f$" % r, (0.05, 0.8), xycoords='axes fraction',
                size="10", va="center")

        dialog.plot.ax.set_xlabel("T, K", ha='right', size="12")
        ylabel = "$"+self.proptex+r",\;"+self.unit.text()+"$"
        dialog.plot.ax.set_ylabel(ylabel, ha='right', size="12")
        dialog.exec()

    def fit(self):
        """Fit experimental data to a DIPPR equation"""
        uniT = unidades.Temperature.text()
        hHeader = ["T, %s" % uniT, "%s, %s" % (self.title(), self.unit.text())]
        dlg = InputTableDialog(
                title=self.title(), DIPPR=True, horizontalHeader=hHeader,
                hasTc=True, Tc=self.parent.Tc.value, t=self.t,
                property=self.data, eq=self.eq.value())

        if dlg.exec():
            t = array(dlg.widget.column(0, unidades.Temperature))
            p = array(dlg.widget.column(1, self.unit))
            eq = dlg.widget.eqDIPPR.value()

            def errf(parametros, eq, t, f):
                var = array([eq]+list(parametros))
                return f-array([DIPPR(self.prop, ti, var) for ti in t])

            # Do the least square fitting
            p0 = [1, 1, 1, 1, 1]
            coeff, cov, info, mesg, ier = optimize.leastsq(
                    errf, p0, args=(eq, t, p), full_output=True)

            ss_err = (info['fvec']**2).sum()
            ss_tot = ((p-p.mean())**2).sum()
            r = 1-(ss_err/ss_tot)

            if ier in range(1, 5):
                self.fill([eq]+list(coeff)+[min(t), max(t)])
                self.t = list(t)
                self.data = list(p)
                self.plot(r)
                self.valueChanged.emit()
            else:
                title = self.tr("Warning")
                msg = self.tr("Fit unsuccessfully")
                QtWidgets.QMessageBox.warning(self, title, msg)


class Parametric_widget(QtWidgets.QGroupBox):
    """Composite widget to edit/view a Parametric equation"""
    valueChanged = QtCore.pyqtSignal()

    dict_latex = {
        "henry": r"$\ln H = \frac{A}{T}+B\ln T+CT+D$",
        "tension": r"$\sigma=A(1-T_r)^B$",
        "antoine": r"$P_v=e^{A-\frac{B}{C+T}}$",
        "wagner": r"$P_v=\frac{P_c}{Tr}e^{a\tau+b\tau^{1.5}+c\tau^3+d\tau^6}$",
        "viscosity": r"$\mu=e^{A\left(\frac{1}{T}-\frac{1}{B}\right)}$"}
    dict_prop = {
        "henry": "K_H",
        "tension": r"\sigma",
        "antoine": "P_v",
        "wagner": "P_v",
        "viscosity": r"\mu_g"}
    dict_count = {
        "henry": 4,
        "tension": 2,
        "antoine": 3,
        "wagner": 4,
        "viscosity": 2}

    def __init__(self, title, unit, id=0, prop="", array=None, parent=None):
        """Constructor of widget

        Parameters
        ----------
        title : string
            Title so show of qgroupbox
        unit : string
            Aditional string with unit representation
        id : integer
            Index of compound in database to show
        prop : string
            code name of property to show
        array : list
            List of DIPPR equation representation in format
            [eq, A, B, C, D, E, Tmin, Tmax]
        """
        super().__init__(title, parent)
        self.prop = prop
        self.unit = unit
        self.id = id
        self.parent = parent
        self.t = []
        self.data = []

        self.latex = self.dict_latex[prop]
        self.prop_latex = self.dict_prop[prop]
        self.count = self.dict_count[prop]

        layout = QtWidgets.QGridLayout(self)
        self.formula = QLabelMath(self.latex)
        layout.addWidget(self.formula, 0, 1, 1, 5)
        self.btnFit = QtWidgets.QToolButton()
        self.btnFit.setToolTip(self.tr("Fit parameters from experimental data"))
        self.btnFit.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "fit.png"))))
        self.btnFit.setIconSize(QtCore.QSize(32, 32))
        self.btnFit.setFixedSize(QtCore.QSize(32, 32))
        self.btnFit.clicked.connect(self.fit)
        layout.addWidget(self.btnFit, 1, 1, 2, 1)
        self.btnPlot = QtWidgets.QToolButton()
        self.btnPlot.setToolTip(self.tr("Plot equation vs temperature"))
        self.btnPlot.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "plot.png"))))
        self.btnPlot.setIconSize(QtCore.QSize(32, 32))
        self.btnPlot.setFixedSize(QtCore.QSize(32, 32))
        self.btnPlot.clicked.connect(self.plot)
        self.btnPlot.setEnabled(False)
        layout.addWidget(self.btnPlot, 1, 2, 2, 1)

        txt = ["A", "B", "C", "D"]
        self.coeff = []
        for i in range(self.count):
            layout.addWidget(QtWidgets.QLabel(txt[i]), 1+i, 4)
            self.coeff.append(Entrada_con_unidades(float))
            self.coeff[-1].valueChanged.connect(self.valueChanged.emit)
            layout.addWidget(self.coeff[-1], 1+i, 5)

        if prop == "antoine":
            self.advancedCoeff = []
            self.advancedLabel = []
            self.advancedLabel.append(QtWidgets.QLabel("to"))
            layout.addWidget(self.advancedLabel[-1], 4, 1)
            self.advancedCoeff.append(Entrada_con_unidades(float))
            layout.addWidget(self.advancedCoeff[-1], 4, 2)
            self.advancedLabel.append(QtWidgets.QLabel("n"))
            layout.addWidget(self.advancedLabel[-1], 5, 1)
            self.advancedCoeff.append(Entrada_con_unidades(float))
            layout.addWidget(self.advancedCoeff[-1], 5, 2)
            self.advancedLabel.append(QtWidgets.QLabel("E"))
            layout.addWidget(self.advancedLabel[-1], 4, 4)
            self.advancedCoeff.append(Entrada_con_unidades(float))
            layout.addWidget(self.advancedCoeff[-1], 4, 5)
            self.advancedLabel.append(QtWidgets.QLabel("F"))
            layout.addWidget(self.advancedLabel[-1], 5, 4)
            self.advancedCoeff.append(Entrada_con_unidades(float))
            layout.addWidget(self.advancedCoeff[-1], 5, 5)

        layout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)

        self.changeId(id)
        if array:
            self.fill(array)

    def changeId(self, id):
        """Common action to do when the component id change, change from user
        to system compound, editable to readOnly compound"""
        self.id = id
        if id == 0 or id > 10000:
            self.setReadOnly(False)
            self.btnFit.setEnabled(True)
        else:
            self.setReadOnly(True)
            self.btnFit.setEnabled(False)

        if self.prop == "antoine":
            for entrada in self.advancedCoeff:
                entrada.setVisible(False)
            for label in self.advancedLabel:
                label.setVisible(False)

    def fill(self, array):
        """Populate the widgets with the parametric coefficient in array"""
        # Use local variable enabled to enable the btnPlot button if the
        # parameters are available
        enabled = False
        self.btnFit.setEnabled(self.parent.index >= 10001)

        self.array = array
        for input, value in zip(self.coeff, array):
            if value:
                enabled = True
                input.setValue(value)
        self.btnPlot.setEnabled(enabled)

        if self.prop == "antoine" and len(array) > 3:
            for entrada, value in zip(self.advancedCoeff, array[3:]):
                entrada.setVisible(True)
                entrada.setValue(value)
            for label in self.advancedLabel:
                label.setVisible(True)

    def setReadOnly(self, bool):
        """Set widget readOnly state"""
        for entrada in self.coeff:
            entrada.setReadOnly(bool)
        self.btnFit.setDisabled(bool)

        if self.prop == "antoine":
            for entrada in self.advancedCoeff:
                entrada.setReadOnly(bool)

    @property
    def value(self):
        valor = []
        for entrada in self.coeff:
            elemento = entrada.value
            if elemento:
                valor.append(elemento)
            else:
                valor.append(0)
        if valor != [0]*self.count:
            return valor
        else:
            return []

    @property
    def advancedValue(self):
        valor = []
        for entrada in self.advancedCoeff:
            elemento = entrada.value
            if elemento:
                valor.append(elemento)
            else:
                valor.append(0)
        if sum(valor):
            return valor
        else:
            return []

    def clear(self):
        """Clear widget contents"""
        for entrada in self.coeff:
            entrada.clear()
        self.t = []
        self.data = []
        self.btnPlot.setEnabled(False)

        if self.prop == "antoine":
            for entrada in self.advancedCoeff:
                entrada.clear()

    def formula_Parametric(self, eq, args=None):
        """Calculate the formula of a parametric equation in a latex format

        Parameters
        ----------
        eq : string
            Name of property to calculate
        args : list
            Coefficient of parametric equation
        """
        args = tuple(args)
        if eq == "viscosity":
            if not args:
                args = self.parent.cmp.viscosidad_parametrica
            string = "$\\log\\mu="
            string += "%0.2g\\left(\\frac{1}{T}-\\frac{1}{%0.2g}" % args
            string += "\\right)$"

        elif eq == "antoine":
            if not args:
                args = self.parent.cmp.antoine
            string = "$P_{v}=e^{%0.2g-\\frac{%0.2g}{T%+0.2g}" % args[:3]
            if len(args) > 3:
                string += "+0.43429\\left(T%+0.2g\\right)^{%0.2g}" % args[3:5]
                string += "%+0.2g\\left(T%+0.2g\\right)^{8}" % (
                    args[5], args[3])
                string += "%+0.2g\\left(T%+0.2g\\right)^{12}" % (
                    args[6], args[3])
            string += "}$"

        elif eq == "tension":
            if not args:
                args = self.parent.cmp.tension_superficial_parametrica
            string = "$\\sigma=%0.2g\\left(1-T_{r}\\right)^{%0.2g}$" % args

        elif eq == "henry":
            if not args:
                args = self.parent.cmp.henry
            string = "$\\lnH=\\frac"
            string += "{%0.2g}{T}%+0.2g\\cdot\\lnT%+0.2g\\cdot T%+0.2g$" % args

        elif eq == "wagner":
            string = "$P_{v}=\\frac{P_{c}}{Tr}e^{"
            string += "%0.2g\\tau%+0.2g\\tau^{1.5}%+0.2g\\tau^3%+0.2g" % args
            string += "\\tau^6}$"

        return string

    def function(self, prop):
        function = {
            "viscosity": MuL_Parametric,
            "antoine": Pv_Antoine,
            "wagner": Pv_Wagner,
            "tension": Tension_Parametric,
            "henry": Henry}
        return function[prop]

    def plot(self, r=None):
        """Plot the current parametric correlation equation
        r correlation coefficient, optional parameter for reuse the code in the
        fit data procedure
        """
        if self.parent.cmp and self.parent.cmp.Tf:
            tmin = self.parent.cmp.Tf
        elif self.parent.cmp and self.parent.cmp.presion_vapor != [0]*8:
            tmin = self.parent.cmp.presion_vapor[-2]
        else:
            tmin = 300

        if self.parent.cmp and self.parent and self.parent.cmp.Tb:
            tmax = self.parent.cmp.Tb
        elif self.parent.cmp and self.parent.cmp.presion_vapor != [0]*8:
            tmax = self.parent.cmp.presion_vapor[-1]
        else:
            tmax = 500
        t = linspace(tmin, tmax, 100)

        coeff = self.value
        kw = {}
        if self.prop == "antoine":
            coeff += self.advancedValue
            kw["Tc"] = self.parent.Tc.value
        elif self.prop == "wagner":
            kw["Tc"] = self.parent.Tc.value
            kw["Pc"] = self.parent.Pc.value
        elif self.prop == "tension":
            kw["Tc"] = self.parent.Tc.value

        args = [coeff]
        eq = self.function(self.prop)
        var = [eq(ti, *args, **kw) for ti in t]
        dialog = PlotDialog()
        dialog.addData(t, var)
        if self.t and self.data:
            dialog.addData(self.t, self.data, "ro")
        dialog.plot.ax.grid(True)
        dialog.plot.ax.set_title(self.title(), size="14")

        # Annotate the equation
        formula = self.formula_Parametric(self.prop, coeff)
        dialog.plot.ax.annotate(formula, (0.05, 0.9), xycoords='axes fraction',
                                size="10", va="center")

        # Add correlation coefficient
        if r:
            dialog.plot.ax.annotate(
                "$r^2=%0.5f$" % r, (0.05, 0.8), xycoords='axes fraction',
                size="10", va="center")

        dialog.plot.ax.set_xlabel("T, K", ha='right', size="12")
        ylabel = "$"+self.prop+r",\;"+self.unit.text()+"$"
        dialog.plot.ax.set_ylabel(ylabel, ha='right', size="12")
        dialog.exec()

    def fit(self):
        """Fit experimental data to a parametric equation"""
        hHeader = ["T", "%s" % self.title()]
        dlg = InputTableDialog(
            title=self.title(), horizontalHeader=hHeader,
            hasTc=self.prop == "sigma", Tc=self.parent.Tc.value,
            unit=[unidades.Temperature, self.unit])
        if dlg.exec():
            t = array(dlg.widget.column(0))
            p = array(dlg.widget.column(1))

            eq = self.function(self.prop)
            args = []
            kw = {}
            if self.prop == "antoine":
                kw["Tc"] = self.parent.Tc.value
            elif self.prop == "wagner":
                # Check necessary critical properties
                if self.parent.Tc.value and self.parent.Pc.value:
                    kw["Tc"] = self.parent.Tc.value
                    kw["Pc"] = self.parent.Pc.value
                else:
                    title = self.tr("Warning")
                    msg = self.tr("Tc and Pc necessary for wagner correlation")
                    QtWidgets.QMessageBox.warning(self, title, msg)
                    return
            elif self.prop == "tension":
                if self.parent.Tc.value and self.parent.Pc.value:
                    kw["Tc"] = self.parent.Tc.value
                else:
                    title = self.tr("Warning")
                    msg = self.tr("Tc necessary for tension correlation")
                    QtWidgets.QMessageBox.warning(self, title, msg)
                    return

            def errf(coeff, t, f):
                return f-array([eq(ti, coeff, *args, **kw) for ti in t])

            # Do the least square fitting
            p0 = [1.]*self.count
            coeff, cov, info, mesg, ier = optimize.leastsq(
                    errf, p0, args=(t, p), full_output=True)

            ss_err = (info['fvec']**2).sum()
            ss_tot = ((p-p.mean())**2).sum()
            r = 1-(ss_err/ss_tot)

            if ier in range(1, 5):
                self.fill(list(coeff))
                self.t = list(t)
                self.data = list(p)
                self.plot(r)
                self.valueChanged.emit()
            else:
                title = self.tr("Warning")
                msg = self.tr("Fit unsuccessfully")
                QtWidgets.QMessageBox.warning(self, title, msg)


class View_Component(QtWidgets.QDialog):
    """Dialog to view the properties of compounds in pychemqt database in user
    friendly format. This dialog can be used too to edit a user component, or
    define a new component"""
    cmp = Componente()

    def __init__(self, index=0, parent=None):
        """
        index of compound in database
          - 0 is the option for define a new component, all the widget can be
            editable
          - 1 < 10000: Compound in pychemqt databse, read only element. The user
            can copy the compound and edit it
          - > 10000:  Compound added by user"""

        super().__init__(parent)
        lyt = QtWidgets.QVBoxLayout(self)

        lytTitle = QtWidgets.QHBoxLayout()
        lytTitle.setSpacing(0)
        if index == 0:
            # Add a QLineEdit to let the user to set the name of compound
            self.setWindowTitle(self.tr("Define New Component"))
            labelName = QtWidgets.QLabel(self.tr("Name")+": ", self)
            lytTitle.addWidget(labelName)
            self.name = QtWidgets.QLineEdit(self)
            self.name.editingFinished.connect(self.setDirty)
            lytTitle.addWidget(self.name)

        else:
            # Add a "titlebar" with navigation options
            if index > 10000:
                self.setWindowTitle(self.tr("Custom Component Properties"))
            else:
                self.setWindowTitle(self.tr("Component Properties"))
            self.btnFirst = QtWidgets.QToolButton(self)
            self.btnFirst.setToolTip(self.tr("Go to first element"))
            self.btnFirst.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "button", "arrow-left-double.png"))))
            self.btnFirst.clicked.connect(partial(self.change, 1))
            lytTitle.addWidget(self.btnFirst)
            self.btnPrev = QtWidgets.QToolButton(self)
            self.btnPrev.setToolTip(self.tr("Go to previous element"))
            self.btnPrev.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "button", "arrow-left.png"))))
            self.btnPrev.clicked.connect(partial(self.change, "-1"))
            lytTitle.addWidget(self.btnPrev)
            self.name = QtWidgets.QLabel(self)
            lytTitle.addWidget(self.name)
            self.btnNext = QtWidgets.QToolButton(self)
            self.btnNext.setToolTip(self.tr("Go to next element"))
            self.btnNext.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "button", "arrow-right.png"))))
            self.btnNext.clicked.connect(partial(self.change, "+1"))
            lytTitle.addWidget(self.btnNext)
            self.btnLast = QtWidgets.QToolButton(self)
            self.btnLast.setToolTip(self.tr("Go to last element"))
            self.btnLast.setIcon(QtGui.QIcon(QtGui.QPixmap(
                os.path.join(IMAGE_PATH, "button", "arrow-right-double.png"))))
            self.btnLast.clicked.connect(partial(self.change, "last"))
            lytTitle.addWidget(self.btnLast)
        lyt.addItem(lytTitle)

        tabWidget = QtWidgets.QTabWidget()
        lyt.addWidget(tabWidget)

        self.btnBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Discard |
            QtWidgets.QDialogButtonBox.StandardButton.Save |
            QtWidgets.QDialogButtonBox.StandardButton.Apply |
            QtWidgets.QDialogButtonBox.StandardButton.Close)
        self.btnBox.clicked.connect(self.buttonClicked)
        lyt.addWidget(self.btnBox)

        # General tab
        tab1 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, self.tr("&General"))
        lytGeneral = QtWidgets.QGridLayout(tab1)
        lytGeneral.addWidget(QtWidgets.QLabel(self.tr("Formula")), 1, 1)
        self.formula1 = QtWidgets.QLineEdit()
        self.formula1.editingFinished.connect(self.setDirty)
        lytGeneral.addWidget(self.formula1, 1, 2)
        lytGeneral.addWidget(QtWidgets.QLabel(self.tr("CAS Number")), 2, 1)
        self.CAS = QtWidgets.QLineEdit()
        self.CAS.editingFinished.connect(self.setDirty)
        lytGeneral.addWidget(self.CAS, 2, 2)
        lytGeneral.addWidget(QtWidgets.QLabel(
            self.tr("Alternative Name")), 3, 1)
        self.alternateName = QtWidgets.QLineEdit()
        self.alternateName.editingFinished.connect(self.setDirty)
        lytGeneral.addWidget(self.alternateName, 3, 2)

        labelSmile = QtWidgets.QLabel(self.tr("Smile Code"))
        lytGeneral.addWidget(labelSmile, 4, 1)
        self.smile = QtWidgets.QLineEdit()
        self.smile.editingFinished.connect(self.setDirty)
        lytGeneral.addWidget(self.smile, 4, 2)

        labelFormula2 = QtWidgets.QLabel(self.tr("Expanded Formula"))
        lytGeneral.addWidget(labelFormula2, 5, 1)
        self.formula2 = QtWidgets.QLabel()
        lytGeneral.addWidget(self.formula2, 6, 1, 5, 3,
                             QtCore.Qt.AlignmentFlag.AlignCenter)
        if os.environ["openbabel"] != "True":
            labelFormula2.setVisible(False)
            self.formula2.setVisible(False)

        lytGeneral.addWidget(QtWidgets.QLabel("M"), 1, 4)
        self.M = Entrada_con_unidades(float, textounidad="g/mol")
        self.M.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.M, 1, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("Tc"), 2, 4)
        self.Tc = Entrada_con_unidades(unidades.Temperature)
        self.Tc.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.Tc, 2, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("Pc"), 3, 4)
        self.Pc = Entrada_con_unidades(unidades.Pressure)
        self.Pc.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.Pc, 3, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("Vc"), 4, 4)
        self.Vc = Entrada_con_unidades(unidades.SpecificVolume)
        self.Vc.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.Vc, 4, 5)
        labelw = QtWidgets.QLabel("ω")
        labelw.setToolTip(self.tr("Acentric factor"))
        lytGeneral.addWidget(labelw, 5, 4)
        self.w = Entrada_con_unidades(float)
        self.w.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.w, 5, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("SG 60ºF"), 6, 4)
        self.SG = Entrada_con_unidades(float)
        self.SG.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.SG, 6, 5)
        labelTm = QtWidgets.QLabel("T<sub>m</sub>")
        labelTm.setToolTip(self.tr("Melting Point")),
        lytGeneral.addWidget(labelTm, 7, 4)
        self.Tm = Entrada_con_unidades(unidades.Temperature)
        self.Tm.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.Tm, 7, 5)
        labelTb = QtWidgets.QLabel("T<sub>b</sub>")
        labelTb.setToolTip(self.tr("Boiling Point")),
        lytGeneral.addWidget(labelTb, 8, 4)
        self.Tb = Entrada_con_unidades(unidades.Temperature)
        self.Tb.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.Tb, 8, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("ΔH<sub>f</sub>"), 9, 4)
        self.calorFormacionGas = Entrada_con_unidades(unidades.Enthalpy)
        self.calorFormacionGas.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.calorFormacionGas, 9, 5)
        lytGeneral.addWidget(QtWidgets.QLabel("ΔG<sub>f</sub>"), 10, 4)
        self.energiaGibbsGas = Entrada_con_unidades(unidades.Enthalpy)
        self.energiaGibbsGas.valueChanged.connect(self.setDirty)
        lytGeneral.addWidget(self.energiaGibbsGas, 10, 5)
        lytGeneral.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Expanding), 11, 3)
        lytGeneral.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 12, 1, 1, 6)

        # Cp tab
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(tab2, self.tr("&Cp"))
        lytCp = QtWidgets.QGridLayout(tab2)
        lytCp.setSpacing(25)

        grpCpIdeal = QtWidgets.QGroupBox(
                self.tr("Cp ideal gas") + ", cal/mol·K")
        lytCp.addWidget(grpCpIdeal, 1, 1, 1, 1)
        lytCpIdeal = QtWidgets.QGridLayout(grpCpIdeal)
        mathTex = "$C_p^o=A+BT+CT^2+DT^3+ET^4+FT^5$"
        self.eqCpGasIdeal = QLabelMath(mathTex)
        lytCpIdeal.addWidget(self.eqCpGasIdeal, 0, 1, 1, 5)
        lytCpIdeal.addWidget(QtWidgets.QLabel("A"), 1, 1)
        self.cpa = Entrada_con_unidades(float)
        self.cpa.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpa, 1, 2)
        lytCpIdeal.addWidget(QtWidgets.QLabel("B"), 2, 1)
        self.cpb = Entrada_con_unidades(float)
        self.cpb.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpb, 2, 2)
        lytCpIdeal.addWidget(QtWidgets.QLabel("C"), 3, 1)
        self.cpc = Entrada_con_unidades(float)
        self.cpc.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpc, 3, 2)
        lytCpIdeal.addWidget(QtWidgets.QLabel("D"), 1, 4)
        self.cpd = Entrada_con_unidades(float)
        self.cpd.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpd, 1, 5)
        lytCpIdeal.addWidget(QtWidgets.QLabel("E"), 2, 4)
        self.cpe = Entrada_con_unidades(float)
        self.cpe.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpe, 2, 5)
        lytCpIdeal.addWidget(QtWidgets.QLabel("F"), 3, 4)
        self.cpf = Entrada_con_unidades(float)
        self.cpf.valueChanged.connect(self.setDirty)
        lytCpIdeal.addWidget(self.cpf, 3, 5)
        lytCpIdeal.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 4, 1, 1, 4)

        self.cpGasDIPPR = DIPPR_widget(
            self.tr("Cp ideal gas DIPPR"), unidades.MolarSpecificHeat, index,
            "cpG", parent=self)
        self.cpGasDIPPR.valueChanged.connect(self.setDirty)
        lytCp.addWidget(self.cpGasDIPPR, 1, 2)
        self.cpLiquidDIPPR = DIPPR_widget(
            self.tr("Cp liquid DIPPR"), unidades.MolarSpecificHeat, index,
            "cpL", parent=self)
        self.cpLiquidDIPPR.valueChanged.connect(self.setDirty)
        lytCp.addWidget(self.cpLiquidDIPPR, 2, 1)
        self.cpSolidDIPPR = DIPPR_widget(
            self.tr("Cp solid DIPPR"), unidades.MolarSpecificHeat, index,
            "cpS", parent=self)
        self.cpSolidDIPPR.valueChanged.connect(self.setDirty)
        lytCp.addWidget(self.cpSolidDIPPR, 2, 2)
        lytCp.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 1)

        # Density tab
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(tab3, self.tr("&Density"))
        lytRho = QtWidgets.QGridLayout(tab3)
        lytRho.setSpacing(25)

        self.RhoSolid = DIPPR_widget(
            self.tr("Solid Density DIPPR"), unidades.MolarDensity, index,
            "rhoS", parent=self)
        self.RhoSolid.valueChanged.connect(self.setDirty)
        lytRho.addWidget(self.RhoSolid, 1, 1)
        self.RhoLiquid = DIPPR_widget(
            self.tr("Liquid Density DIPPR"), unidades.MolarDensity, index,
            "rhoL", parent=self)
        self.RhoLiquid.valueChanged.connect(self.setDirty)
        lytRho.addWidget(self.RhoLiquid, 1, 2)
        lytRho.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Viscosity tab
        tab4 = QtWidgets.QWidget()
        tabWidget.addTab(tab4, self.tr("&Viscosity"))
        lytMu = QtWidgets.QGridLayout(tab4)
        lytMu.setSpacing(25)

        self.muLiquid = DIPPR_widget(
            self.tr("Liquid Viscosity DIPPR"), unidades.Viscosity, index,
            "muL", parent=self)
        self.muLiquid.valueChanged.connect(self.setDirty)
        lytMu.addWidget(self.muLiquid, 1, 1)
        self.muGas = DIPPR_widget(
            self.tr("Gas Viscosity DIPPR"), unidades.Viscosity, index,
            "muG", parent=self)
        self.muGas.valueChanged.connect(self.setDirty)
        lytMu.addWidget(self.muGas, 1, 2)
        self.muParametric = Parametric_widget(
            self.tr("Viscosity"), unidades.Viscosity, index, "viscosity",
            parent=self)
        self.muParametric.valueChanged.connect(self.setDirty)
        lytMu.addWidget(self.muParametric, 2, 1)
        lytMu.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 1)

        # Pressure vapor and vaporization heat tab
        tab5 = QtWidgets.QWidget()
        tabWidget.addTab(tab5, self.tr("P&v && Hv"))
        lytVapor = QtWidgets.QGridLayout(tab5)
        lytVapor.setSpacing(25)

        self.Hv = DIPPR_widget(
            self.tr("Heat of vaporization DIPPR"), unidades.MolarEnthalpy,
            index, "Hv", parent=self)
        self.Hv.valueChanged.connect(self.setDirty)
        lytVapor.addWidget(self.Hv, 1, 1)
        self.PvDIPPR = DIPPR_widget(
            self.tr("Vapor Pressure DIPPR"), unidades.Pressure, index,
            "Pv", parent=self)
        self.PvDIPPR.valueChanged.connect(self.setDirty)
        lytVapor.addWidget(self.PvDIPPR, 1, 2)
        self.PvAntoine = Parametric_widget(
            self.tr("Antoine Vapor Pressure"), unidades.Pressure, index,
            "antoine", parent=self)
        self.PvAntoine.valueChanged.connect(self.setDirty)
        lytVapor.addWidget(self.PvAntoine, 2, 1)
        self.PvWagner = Parametric_widget(
            self.tr("Wagner Vapor Pressure"), unidades.Pressure, index,
            "wagner", parent=self)
        self.PvWagner.valueChanged.connect(self.setDirty)
        lytVapor.addWidget(self.PvWagner, 2, 2)
        lytVapor.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 1)

        # Thermal conductiity and surface tension tab
        tab6 = QtWidgets.QWidget()
        tabWidget.addTab(tab6, self.tr("&Tension && Conductivity"))
        lytThermal = QtWidgets.QGridLayout(tab6)
        lytThermal.setSpacing(25)

        self.ThermalLiquid = DIPPR_widget(
            self.tr("Liquid Thermal Conductivity DIPPR"),
            unidades.ThermalConductivity, index, "kL", parent=self)
        self.ThermalLiquid.valueChanged.connect(self.setDirty)
        lytThermal.addWidget(self.ThermalLiquid, 1, 1)
        self.ThermalGas = DIPPR_widget(
            self.tr("Gas Thermal Conductivity DIPPR"),
            unidades.ThermalConductivity, index, "kG", parent=self)
        self.ThermalGas.valueChanged.connect(self.setDirty)
        lytThermal.addWidget(self.ThermalGas, 1, 2)
        self.SigmaDIPPR = DIPPR_widget(
            self.tr("Surface Tension DIPPR"), unidades.Tension, index,
            "sigma", parent=self)
        self.SigmaDIPPR.valueChanged.connect(self.setDirty)
        lytThermal.addWidget(self.SigmaDIPPR, 2, 1)
        self.SigmaParametric = Parametric_widget(
            self.tr("Surface Tension"), unidades.Tension, index, "tension",
            parent=self)
        self.SigmaParametric.valueChanged.connect(self.setDirty)
        lytThermal.addWidget(self.SigmaParametric, 2, 2)
        lytThermal.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 1)

        # EOS tab
        # Add here other EOS parameter when add to database
        tab7 = QtWidgets.QWidget()
        tabWidget.addTab(tab7, self.tr("&EoS"))
        lytEOS = QtWidgets.QGridLayout(tab7)
        lytEOS.setSpacing(25)

        self.Henry = Parametric_widget(
            self.tr("Henry Costant"), unidades.Pressure, index, "henry",
            parent=self)
        self.Henry.valueChanged.connect(self.setDirty)
        lytEOS.addWidget(self.Henry, 1, 1, 2, 1)
        lytEOS.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 1)

        grp = QtWidgets.QGroupBox("Cubic")
        lyt = QtWidgets.QVBoxLayout(grp)
        self.cubicElection = QtWidgets.QComboBox()
        lyt.addWidget(self.cubicElection)
        self.cubicStacked = QtWidgets.QStackedWidget()
        lyt.addWidget(self.cubicStacked)
        lytEOS.addWidget(grp, 1, 2, 2, 1)

        grpMSRK = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpMSRK)
        self.cubicElection.addItem("MSRK")
        lyt_MSRK = QtWidgets.QGridLayout(grpMSRK)
        math = r"$\alpha=1+m\left(1-Tr\right)+n\left(\frac{1}{T_R}-1\right)$"
        formulaMSRK = QLabelMath(math)
        lyt_MSRK.addWidget(formulaMSRK, 1, 1, 1, 3)
        lyt_MSRK.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_MSRK.addWidget(QtWidgets.QLabel("m"), 3, 1)
        self.MSRKa = Entrada_con_unidades(float)
        self.MSRKa.valueChanged.connect(self.setDirty)
        lyt_MSRK.addWidget(self.MSRKa, 3, 2)
        lyt_MSRK.addWidget(QtWidgets.QLabel("n"), 4, 1)
        self.MSRKb = Entrada_con_unidades(float)
        self.MSRKb.valueChanged.connect(self.setDirty)
        lyt_MSRK.addWidget(self.MSRKb, 4, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(partial(self.help, "lib.EoS.Cubic.MSRK.html"))
        lyt_MSRK.addWidget(button, 4, 6)
        lyt_MSRK.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 5)

        grpMelhem = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpMelhem)
        self.cubicElection.addItem("PR_Melhem")
        lyt_Melhem = QtWidgets.QGridLayout(grpMelhem)
        math = r"$\alpha = \exp\left(m\left(1-T_r\right) +"
        math += r"n \left(1-\sqrt{T_r}\right)^2\right)$"
        formulaMelhem = QLabelMath(math)
        lyt_Melhem.addWidget(formulaMelhem, 1, 1, 1, 2)
        lyt_Melhem.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_Melhem.addWidget(QtWidgets.QLabel("m"), 3, 1)
        self.Melhemm = Entrada_con_unidades(float)
        self.Melhemm.valueChanged.connect(self.setDirty)
        lyt_Melhem.addWidget(self.Melhemm, 3, 2)
        lyt_Melhem.addWidget(QtWidgets.QLabel("n"), 4, 1)
        self.Melhemn = Entrada_con_unidades(float)
        self.Melhemn.valueChanged.connect(self.setDirty)
        lyt_Melhem.addWidget(self.Melhemn, 4, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PRMelhem.html"))
        lyt_Melhem.addWidget(button, 4, 6)
        lyt_Melhem.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 5)

        grpSRKAPI = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpSRKAPI)
        self.cubicElection.addItem("SRK-API")
        lyt_SRKAPI = QtWidgets.QGridLayout(grpSRKAPI)
        math = r"$\alpha^{0.5} = 1 + S_1\left(1-\sqrt{Tr}\right) + "
        math += r"S_2\frac{\left(1-\sqrt{Tr}\right)}{\sqrt{T_r}}$"
        formulaSRKAPI = QLabelMath(math)
        lyt_SRKAPI.addWidget(formulaSRKAPI, 1, 1, 1, 2)
        lyt_SRKAPI.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_SRKAPI.addWidget(QLabelMath("$S_1$"), 3, 1)
        self.SRKAPIS1 = Entrada_con_unidades(float)
        self.SRKAPIS1.valueChanged.connect(self.setDirty)
        lyt_SRKAPI.addWidget(self.SRKAPIS1, 3, 2)
        lyt_SRKAPI.addWidget(QLabelMath("$S_2$"), 4, 1)
        self.SRKAPIS2 = Entrada_con_unidades(float)
        self.SRKAPIS2.valueChanged.connect(self.setDirty)
        lyt_SRKAPI.addWidget(self.SRKAPIS2, 4, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.SRKAPI.html"))
        lyt_SRKAPI.addWidget(button, 4, 6)
        lyt_SRKAPI.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 5)

        grpMathias = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpMathias)
        self.cubicElection.addItem("PR-Mathias-Copeman")
        lyt_Mathias = QtWidgets.QGridLayout(grpMathias)
        lyt_Mathias.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_Mathias.addWidget(QLabelMath("$C_1$"), 3, 1)
        self.MathiasCopemanC1 = Entrada_con_unidades(float)
        self.MathiasCopemanC1.valueChanged.connect(self.setDirty)
        lyt_Mathias.addWidget(self.MathiasCopemanC1, 3, 2)
        lyt_Mathias.addWidget(QLabelMath("$C_2$"), 4, 1)
        self.MathiasCopemanC2 = Entrada_con_unidades(float)
        self.MathiasCopemanC2.valueChanged.connect(self.setDirty)
        lyt_Mathias.addWidget(self.MathiasCopemanC2, 4, 2)
        lyt_Mathias.addWidget(QLabelMath("$C_3$"), 5, 1)
        self.MathiasCopemanC3 = Entrada_con_unidades(float)
        self.MathiasCopemanC3.valueChanged.connect(self.setDirty)
        lyt_Mathias.addWidget(self.MathiasCopemanC3, 5, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PRMathiasCopeman.html"))
        lyt_Mathias.addWidget(button, 5, 6)
        lyt_Mathias.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 5)

        grpSV = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpSV)
        self.cubicElection.addItem("PR-Stryjek-Vera")
        lyt_SV = QtWidgets.QGridLayout(grpSV)
        lyt_SV.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_SV.addWidget(QLabelMath("$K_1$"), 3, 1)
        self.SVk0 = Entrada_con_unidades(float)
        self.SVk0.valueChanged.connect(self.setDirty)
        lyt_SV.addWidget(self.SVk0, 3, 2)
        lyt_SV.addWidget(QLabelMath("$K_2$"), 4, 1)
        self.SVk1 = Entrada_con_unidades(float)
        self.SVk1.valueChanged.connect(self.setDirty)
        lyt_SV.addWidget(self.SVk1, 4, 2)
        lyt_SV.addWidget(QLabelMath("$K_3$"), 5, 1)
        self.SVk2 = Entrada_con_unidades(float)
        self.SVk2.valueChanged.connect(self.setDirty)
        lyt_SV.addWidget(self.SVk2, 5, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PRSV2.html"))
        lyt_SV.addWidget(button, 5, 6)
        lyt_SV.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 5)

        grpZV = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpZV)
        self.cubicElection.addItem("PR-Zaboloy-Vera")
        lyt_ZV = QtWidgets.QGridLayout(grpZV)
        lyt_ZV.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_ZV.addWidget(QLabelMath("$K_1$"), 3, 1)
        self.ZVk1 = Entrada_con_unidades(float)
        self.ZVk1.valueChanged.connect(self.setDirty)
        lyt_ZV.addWidget(self.ZVk1, 3, 2)
        lyt_ZV.addWidget(QLabelMath("$K_2$"), 4, 1)
        self.ZVk2 = Entrada_con_unidades(float)
        self.ZVk2.valueChanged.connect(self.setDirty)
        lyt_ZV.addWidget(self.ZVk2, 4, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PRZV.html"))
        lyt_ZV.addWidget(button, 4, 6)
        lyt_ZV.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 5)

        grpAlmeida = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpAlmeida)
        self.cubicElection.addItem("PR-Almeida")
        lyt_Almeida = QtWidgets.QGridLayout(grpAlmeida)
        math = r"$\alpha = \exp\left(m\left(1-T_r\right)\left|1-T_r\right|^"
        math += r"{\Gamma-1} + n\left(\frac{1}{Tr}-\right)\right)$"
        formulaAlmeida = QLabelMath(math)
        lyt_Almeida.addWidget(formulaAlmeida, 1, 1, 1, 2)
        lyt_Almeida.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_Almeida.addWidget(QtWidgets.QLabel("m"), 3, 1)
        self.Almeidam = Entrada_con_unidades(float)
        self.Almeidam.valueChanged.connect(self.setDirty)
        lyt_Almeida.addWidget(self.Almeidam, 3, 2)
        lyt_Almeida.addWidget(QtWidgets.QLabel("n"), 4, 1)
        self.Almeidan = Entrada_con_unidades(float)
        self.Almeidan.valueChanged.connect(self.setDirty)
        lyt_Almeida.addWidget(self.Almeidan, 4, 2)
        lyt_Almeida.addWidget(QLabelMath(r"$\Gamma$"), 5, 1)
        self.AlmeidaG = Entrada_con_unidades(float)
        self.AlmeidaG.valueChanged.connect(self.setDirty)
        lyt_Almeida.addWidget(self.AlmeidaG, 5, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PRAlmeida.html"))
        lyt_Almeida.addWidget(button, 5, 6)
        lyt_Almeida.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 5)

        grpPT = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpPT)
        self.cubicElection.addItem("PR-Patel-Teja")
        lyt_PT = QtWidgets.QGridLayout(grpPT)
        lyt_PT.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 2, 3)
        lyt_PT.addWidget(QLabelMath(r"$\zeta_c$"), 3, 1)
        self.PTZc = Entrada_con_unidades(float)
        self.PTZc.valueChanged.connect(self.setDirty)
        lyt_PT.addWidget(self.PTZc, 3, 2)
        lyt_PT.addWidget(QtWidgets.QLabel("F"), 4, 1)
        self.PTF = Entrada_con_unidades(float)
        self.PTF.valueChanged.connect(self.setDirty)
        lyt_PT.addWidget(self.PTF, 4, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.PT.html"))
        lyt_PT.addWidget(button, 4, 6)
        lyt_PT.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 5)

        grpTB = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpTB)
        self.cubicElection.addItem("Trebble-Bishnoi")
        lyt_TB = QtWidgets.QGridLayout(grpTB)
        lyt_TB.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 3)
        lyt_TB.addWidget(QLabelMath(r"$q_1$"), 4, 1)
        self.TBq1 = Entrada_con_unidades(float)
        self.TBq1.valueChanged.connect(self.setDirty)
        lyt_TB.addWidget(self.TBq1, 4, 2)
        lyt_TB.addWidget(QLabelMath(r"$q_2$"), 5, 1)
        self.TBq2 = Entrada_con_unidades(float)
        self.TBq2.valueChanged.connect(self.setDirty)
        lyt_TB.addWidget(self.TBq2, 5, 2)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.TB.html"))
        lyt_TB.addWidget(button, 5, 6)
        lyt_TB.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 6, 5)

        grpTBS = QtWidgets.QWidget()
        self.cubicStacked.addWidget(grpTBS)
        self.cubicElection.addItem("Trebble-Bishnoi-Salim")
        lyt_TBS = QtWidgets.QGridLayout(grpTBS)
        lyt_TBS.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 3, 3)
        lyt_TBS.addWidget(QLabelMath(r"$\zeta_c$"), 4, 1)
        self.TBSZc = Entrada_con_unidades(float)
        self.TBSZc.valueChanged.connect(self.setDirty)
        lyt_TBS.addWidget(self.TBSZc, 4, 2)
        lyt_TBS.addWidget(QtWidgets.QLabel("m"), 5, 1)
        self.TBSm = Entrada_con_unidades(float)
        self.TBSm.valueChanged.connect(self.setDirty)
        lyt_TBS.addWidget(self.TBSm, 5, 2)
        lyt_TBS.addWidget(QtWidgets.QLabel("p"), 4, 4)
        self.TBSp = Entrada_con_unidades(float)
        self.TBSp.valueChanged.connect(self.setDirty)
        lyt_TBS.addWidget(self.TBSp, 4, 5)
        lyt_TBS.addWidget(QtWidgets.QLabel("d"), 5, 4)
        self.TBSd = Entrada_con_unidades(float)
        self.TBSd.valueChanged.connect(self.setDirty)
        lyt_TBS.addWidget(self.TBSd, 5, 5)
        button = QtWidgets.QToolButton(self)
        button.setIcon(QtGui.QIcon(QtGui.QPixmap(
            os.path.join(IMAGE_PATH, "button", "help.png"))))
        button.clicked.connect(
            partial(self.help, "lib.EoS.Cubic.TBS.html"))
        lyt_TBS.addWidget(button, 5, 7)
        lyt_TBS.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 7, 6)

        lytEOS.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 5, 5)
        self.cubicElection.currentIndexChanged.connect(
            self.cubicStacked.setCurrentIndex)

        # Others tab
        tab8 = QtWidgets.QWidget()
        tabWidget.addTab(tab8, self.tr("&Others"))
        lytOthers = QtWidgets.QGridLayout(tab8)

        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Solubility Parameter")), 1, 1)
        self.SolubilityPar = Entrada_con_unidades(unidades.SolubilityParameter)
        self.SolubilityPar.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.SolubilityPar, 1, 2)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("Dipole Moment")), 2, 1)
        self.Dipole = Entrada_con_unidades(unidades.DipoleMoment)
        self.Dipole.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.Dipole, 2, 2)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Molecular Diameter")), 3, 1)
        self.MolecularDiameter = Entrada_con_unidades(float, textounidad="Å")
        self.MolecularDiameter.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.MolecularDiameter, 3, 2)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Net calorific value")), 4, 1)
        self.NetHeat = Entrada_con_unidades(unidades.MolarEnthalpy)
        self.NetHeat.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.NetHeat, 4, 2)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Gross calorific value")), 5, 1)
        self.GrossHeat = Entrada_con_unidades(unidades.MolarEnthalpy)
        self.GrossHeat.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.GrossHeat, 5, 2)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Liquid Volume Costant")), 6, 1)
        self.VolLiqConstant = Entrada_con_unidades(unidades.MolarVolume)
        self.VolLiqConstant.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.VolLiqConstant, 6, 2)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("API Gravity")), 7, 1)
        self.API = Entrada_con_unidades(float)
        self.API.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.API, 7, 2)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Modified acentric factor")), 8, 1)
        self.wMod = Entrada_con_unidades(float)
        self.wMod.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.wMod, 8, 2)
        lytOthers.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("UNIQUAC area")), 1, 4)
        self.UNIQUACArea = Entrada_con_unidades(float)
        self.UNIQUACArea.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.UNIQUACArea, 1, 5)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("UNIQUAC volume")), 2, 4)
        self.UNIQUACVolume = Entrada_con_unidades(float)
        self.UNIQUACVolume.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.UNIQUACVolume, 2, 5)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("Wilson volume")), 3, 4)
        self.wilson = Entrada_con_unidades(unidades.MolarVolume)
        self.wilson.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.wilson, 3, 5)
        lytOthers.addWidget(QtWidgets.QLabel(
            self.tr("Stiehl polar factor")), 4, 4)
        self.stiel = Entrada_con_unidades(float)
        self.stiel.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.stiel, 4, 5)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("Rackett constant")), 5, 4)
        self.rackett = Entrada_con_unidades(float)
        self.rackett.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.rackett, 5, 5)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("Polar parameter")), 6, 4)
        self.PolarParameter = Entrada_con_unidades(float)
        self.PolarParameter.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.PolarParameter, 6, 5)
        lytOthers.addWidget(QtWidgets.QLabel("Eps/k"), 7, 4)
        self.EpsK = Entrada_con_unidades(float)
        self.EpsK.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.EpsK, 7, 5)
        lytOthers.addWidget(QtWidgets.QLabel(self.tr("Watson factor")), 8, 4)
        self.watson = Entrada_con_unidades(float)
        self.watson.valueChanged.connect(self.setDirty)
        lytOthers.addWidget(self.watson, 8, 5)

        lytOthers.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 9, 1)

        label_UNIFAC = QtWidgets.QLabel(self.tr("UNIFAC Groups"))
        label_UNIFAC.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        lytOthers.addWidget(label_UNIFAC, 10, 1)
        HHeader = [self.tr("Group"), self.tr("Contribution")]
        self.UNIFAC = Tabla(2, filas=1, horizontalHeader=HHeader,
                            verticalHeader=False, readOnly=index)
        self.UNIFAC.setItemDelegateForColumn(0, SpinEditor(self))
        self.UNIFAC.setItemDelegateForColumn(1, SpinEditor(self))
        self.UNIFAC.resizeColumnsToContents()
        self.UNIFAC.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                                  QtWidgets.QSizePolicy.Policy.Fixed)
        lytOthers.addWidget(self.UNIFAC, 10, 2)

        lytOthers.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 11, 3)

        self.change(index)
        self.dirty = False

    def change(self, index):
        if index == "-1":
            if self.index == 10001:
                index = sql.N_comp
            else:
                index = self.index-1
        elif index == "+1":
            if self.index == sql.N_comp:
                index = 10001
            else:
                index = self.index+1
        elif index == "last":
            if sql.N_comp_Custom:
                index = 10000+sql.N_comp_Custom
            else:
                index = sql.N_comp

        btnDiscard = self.btnBox.button(QtWidgets.QDialogButtonBox.StandardButton.Discard)
        btnSave = self.btnBox.button(QtWidgets.QDialogButtonBox.StandardButton.Save)
        btnApply = self.btnBox.button(QtWidgets.QDialogButtonBox.StandardButton.Apply)
        btnClose = self.btnBox.button(QtWidgets.QDialogButtonBox.StandardButton.Close)
        if index == 0:
            btnDiscard.setVisible(True)
            btnSave.setVisible(True)
            btnApply.setVisible(False)
            btnClose.setVisible(False)
        elif index < 10000:
            btnDiscard.setVisible(False)
            btnSave.setVisible(False)
            btnApply.setVisible(False)
            btnClose.setVisible(True)
        else:
            btnDiscard.setVisible(True)
            btnSave.setVisible(False)
            btnApply.setVisible(True)
            btnClose.setVisible(False)

        if index == 0 or index > 10000:
            self.setReadOnly(False)
        else:
            self.setReadOnly(True)

        if index > 0:
            self.clear()
            self.btnFirst.setDisabled(index == 1)
            self.btnPrev.setDisabled(index == 1)
            if sql.N_comp_Custom:
                last = 10000+sql.N_comp_Custom
            else:
                last = sql.N_comp
            self.btnNext.setDisabled(index == last)
            self.btnLast.setDisabled(index == last)
            self.fill(index)
        self.index = index

    def setDirty(self):
        self.dirty = True

    def clear(self):
        self.name.clear()
        self.alternateName.clear()
        self.CAS.clear()
        self.formula1.clear()
        self.smile.clear()
        self.UNIFAC.clear()
        self.M.clear()
        self.Tc.clear()
        self.Pc.clear()
        self.Vc.clear()
        self.Tm.clear()
        self.Tb.clear()
        self.w.clear()
        self.SG.clear()
        self.calorFormacionGas.clear()
        self.energiaGibbsGas.clear()
        self.cpa.clear()
        self.cpb.clear()
        self.cpc.clear()
        self.cpd.clear()
        self.cpe.clear()
        self.cpf.clear()
        self.cpGasDIPPR.clear()
        self.cpLiquidDIPPR.clear()
        self.cpSolidDIPPR.clear()
        self.RhoSolid.clear()
        self.RhoLiquid.clear()
        self.muLiquid.clear()
        self.muGas.clear()
        self.Hv.clear()
        self.PvDIPPR.clear()
        self.ThermalLiquid.clear()
        self.ThermalGas.clear()
        self.SigmaDIPPR.clear()
        self.muParametric.clear()
        self.PvAntoine.clear()
        self.PvWagner.clear()
        self.SigmaParametric.clear()
        self.Henry.clear()
        self.MSRKa.clear()
        self.MSRKb.clear()
        self.Melhemm.clear()
        self.Melhemn.clear()
        self.SRKAPIS1.clear()
        self.SRKAPIS2.clear()
        self.MathiasCopemanC1.clear()
        self.MathiasCopemanC2.clear()
        self.MathiasCopemanC3.clear()
        self.SVk0.clear()
        self.SVk1.clear()
        self.SVk2.clear()
        self.ZVk1.clear()
        self.ZVk2.clear()
        self.Almeidam.clear()
        self.Almeidan.clear()
        self.AlmeidaG.clear()
        self.PTZc.clear()
        self.PTF.clear()
        self.TBq1.clear()
        self.TBq2.clear()
        self.TBSZc.clear()
        self.TBSm.clear()
        self.TBSp.clear()
        self.TBSd.clear()
        self.SolubilityPar.clear()
        self.Dipole.clear()
        self.MolecularDiameter.clear()
        self.NetHeat.clear()
        self.GrossHeat.clear()
        self.VolLiqConstant.clear()
        self.API.clear()
        self.wMod.clear()
        self.UNIQUACArea.clear()
        self.UNIQUACVolume.clear()
        self.wilson.clear()
        self.stiel.clear()
        self.rackett.clear()
        self.EpsK.clear()
        self.watson.clear()

    def fill(self, index):
        self.index = index
        self.cmp = Componente(index)
        self.name.setText("%i - %s" % (self.cmp.id, self.cmp.name))
        self.alternateName.setText(self.cmp.Synonyms)
        self.CAS.setText(self.cmp.CASNumber)
        self.formula1.setText(self.cmp.formula)

        if self.cmp.smile and os.environ["openbabel"] == "True":
            self.smile.setText(self.cmp.smile)

            imageFile = tempfile.NamedTemporaryFile("w", suffix=".png")
            imageFromSmile(self.cmp.smile, imageFile.name, Preferences)

            self.formula2.setPixmap(QtGui.QPixmap(imageFile.name))

        if self.cmp.UNIFAC != []:
            self.UNIFAC.setRowCount(len(self.cmp.UNIFAC))
            for i in range(len(self.cmp.UNIFAC)):
                self.UNIFAC.setRowHeight(i, 25)
                item = QtWidgets.QTableWidgetItem(str(self.cmp.UNIFAC[i][0]))
                item2 = QtWidgets.QTableWidgetItem(str(self.cmp.UNIFAC[i][1]))
                self.UNIFAC.setItem(i, 0, item)
                self.UNIFAC.setItem(i, 1, item2)
                self.UNIFAC.item(i, 1).setTextAlignment(
                    QtCore.Qt.AlignmentFlag.AlignRight
                    | QtCore.Qt.AlignmentFlag.AlignVCenter)

        if self.cmp.M:
            self.M.setValue(self.cmp.M)
        if self.cmp.Tc:
            self.Tc.setValue(self.cmp.Tc)
        if self.cmp.Pc:
            self.Pc.setValue(self.cmp.Pc)
        if self.cmp.Vc:
            self.Vc.setValue(self.cmp.Vc)
        if self.cmp.Tf:
            self.Tm.setValue(self.cmp.Tf)
        if self.cmp.Tb:
            self.Tb.setValue(self.cmp.Tb)
        if self.cmp.f_acent:
            self.w.setValue(self.cmp.f_acent)
        if self.cmp.SG:
            self.SG.setValue(self.cmp.SG)
        if self.cmp.Hf:
            self.calorFormacionGas.setValue(self.cmp.Hf)
        if self.cmp.Gf:
            self.energiaGibbsGas.setValue(self.cmp.Gf)
        self.cpa.setValue(self.cmp.cp[0])
        self.cpb.setValue(self.cmp.cp[1])
        self.cpc.setValue(self.cmp.cp[2])
        self.cpd.setValue(self.cmp.cp[3])
        self.cpe.setValue(self.cmp.cp[4])
        self.cpf.setValue(self.cmp.cp[5])

        self.cpGasDIPPR.fill(self.cmp._dipprCpG)
        self.cpLiquidDIPPR.fill(self.cmp._dipprCpL)
        self.cpSolidDIPPR.fill(self.cmp._dipprCpS)
        self.RhoSolid.fill(self.cmp._dipprRhoS)
        self.RhoLiquid.fill(self.cmp._dipprRhoL)
        self.muLiquid.fill(self.cmp._dipprMuL)
        self.muGas.fill(self.cmp._dipprMuG)
        self.Hv.fill(self.cmp._dipprHv)
        self.PvDIPPR.fill(self.cmp._dipprPv)
        self.ThermalLiquid.fill(self.cmp._dipprKL)
        self.ThermalGas.fill(self.cmp._dipprKG)
        self.SigmaDIPPR.fill(self.cmp._dipprSigma)

        self.muParametric.fill(self.cmp._parametricMu)
        self.PvAntoine.fill(self.cmp.antoine)
        self.PvWagner.fill(self.cmp.wagner)
        self.SigmaParametric.fill(self.cmp._parametricSigma)
        self.Henry.fill(self.cmp.henry)

        if index in MSRKdat:
            self.MSRKa.setValue(MSRKdat[index][0])
            self.MSRKb.setValue(MSRKdat[index][1])

        if index in Melhemdat:
            self.Melhemm.setValue(Melhemdat[index][0])
            self.Melhemn.setValue(Melhemdat[index][1])

        if index in S1:
            self.SRKAPIS1.setValue(S1[index])
        if index in S2:
            self.SRKAPIS2.setValue(S2[index])

        if self.cmp.CASNumber in Mathiasdat:
            self.MathiasCopemanC1.setValue(Mathiasdat[self.cmp.CASNumber][0])
            self.MathiasCopemanC2.setValue(Mathiasdat[self.cmp.CASNumber][1])
            self.MathiasCopemanC3.setValue(Mathiasdat[self.cmp.CASNumber][2])

        if index in PRSVdat:
            self.SVk0.setValue(PRSVdat[index][0])
            self.SVk1.setValue(PRSVdat[index][1])
            self.SVk2.setValue(PRSVdat[index][2])

        if index in PRZVdat:
            self.ZVk1.setValue(PRZVdat[index][0])
            self.ZVk2.setValue(PRZVdat[index][1])

        if index in Almeidadat:
            self.Almeidam.setValue(Almeidadat[index][0])
            self.Almeidan.setValue(Almeidadat[index][1])
            self.AlmeidaG.setValue(Almeidadat[index][2])

        if index in PTdat:
            self.PTZc.setValue(PTdat[index][0])
            self.PTF.setValue(PTdat[index][1])

        if index in TBdat:
            self.TBq1.setValue(TBdat[index][0])
            self.TBq2.setValue(TBdat[index][1])

        if index in TBSdat:
            self.TBSZc.setValue(TBSdat[index][0])
            self.TBSm.setValue(TBSdat[index][1])
            self.TBSp.setValue(TBSdat[index][2])
            self.TBSd.setValue(TBSdat[index][3])

        if self.cmp.SolubilityParameter:
            self.SolubilityPar.setValue(self.cmp.SolubilityParameter)
        if self.cmp.dipole:
            self.Dipole.setValue(self.cmp.dipole)
        if self.cmp.Dm:
            self.MolecularDiameter.setValue(self.cmp.Dm)
        if self.cmp.NetHeating:
            self.NetHeat.setValue(self.cmp.NetHeating)
        if self.cmp.GrossHeating:
            self.GrossHeat.setValue(self.cmp.GrossHeating)
        if self.cmp.Vliq:
            self.VolLiqConstant.setValue(self.cmp.Vliq)
        if self.cmp.API:
            self.API.setValue(self.cmp.API)
        if self.cmp.f_acent_mod:
            self.wMod.setValue(self.cmp.f_acent_mod)
        if self.cmp.UNIQUAC_area:
            self.UNIQUACArea.setValue(self.cmp.UNIQUAC_area)
        if self.cmp.UNIQUAC_volumen:
            self.UNIQUACVolume.setValue(self.cmp.UNIQUAC_volumen)
        if self.cmp.wilson:
            self.wilson.setValue(self.cmp.wilson)
        if self.cmp.stiel:
            self.stiel.setValue(self.cmp.stiel)
        if self.cmp.rackett:
            self.rackett.setValue(self.cmp.rackett)
        if self.cmp.PolarParameter:
            self.PolarParameter.setValue(self.cmp.PolarParameter)
        if self.cmp.ek:
            self.EpsK.setValue(self.cmp.ek)
        if self.cmp.Kw:
            self.watson.setValue(self.cmp.Kw)

    def getComponent(self):
        new = []
        new.append(str(self.formula1.text()))
        if self.index == 0:
            new.append(str(self.name.text()))
        else:
            new.append(str(self.name.text()).split(" - ")[1])

        new.append(self.M.value)
        new.append(self.Tc.value)
        new.append(self.Pc.value)
        new.append(self.Vc.value)
        new.append(self.API.value)

        cpideal = []
        cpideal.append(self.cpa.value)
        cpideal.append(self.cpb.value)
        cpideal.append(self.cpc.value)
        cpideal.append(self.cpd.value)
        cpideal.append(self.cpe.value)
        cpideal.append(self.cpf.value)
        new.append(cpideal)

        new.append(self.PvAntoine.value)
        new.append(self.Henry.value)
        new.append(self.muParametric.value)
        new.append(self.SigmaParametric.value)

        new.append(self.RhoSolid.value)
        new.append(self.RhoLiquid.value)
        new.append(self.PvDIPPR.value)
        new.append(self.Hv.value)
        new.append(self.cpSolidDIPPR.value)
        new.append(self.cpLiquidDIPPR.value)
        new.append(self.cpGasDIPPR.value)
        new.append(self.muLiquid.value)
        new.append(self.muGas.value)
        new.append(self.ThermalLiquid.value)
        new.append(self.ThermalGas.value)
        new.append(self.SigmaDIPPR.value)

        new.append(self.Dipole.value)
        new.append(self.VolLiqConstant.value)
        new.append(self.rackett.value)
        new.append(self.SG.value)
        new.append(self.w.value)
        new.append(self.SolubilityPar.value)
        new.append(self.watson.value)

        # TODO: Add support for save equation parameter in a file
        # msrk = []
        # msrk.append(self.MSRKa.value)
        # msrk.append(self.MSRKb.value)
        # new.append(msrk)

        new.append(self.stiel.value)
        new.append(self.Tb.value)
        new.append(self.Tm.value)
        new.append(str(self.CAS.text()))
        new.append("")

        UNIFAC = self.UNIFAC.getData()
        new.append(UNIFAC)

        new.append(self.MolecularDiameter.value)
        new.append(self.EpsK.value)
        new.append(self.UNIQUACArea.value)
        new.append(self.UNIQUACVolume.value)
        new.append(self.wMod.value)
        new.append(self.calorFormacionGas.value)
        new.append(self.energiaGibbsGas.value)
        new.append(self.wilson.value)
        new.append(self.NetHeat.value)
        new.append(self.GrossHeat.value)

        new.append(str(self.alternateName.text()))
        new.append(0)
        new.append(0)
        new.append(0)

        new.append(self.PolarParameter.value)
        new.append(str(self.smile.text()))

        # Added Antoine and Wagner vapor pressure
        new.append(self.PvAntoine.advancedValue)
        new.append(self.PvWagner.value)

        return new

    def buttonClicked(self, boton):
        role = self.btnBox.buttonRole(boton)
        if role == QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole:
            componente = self.getComponent()
            sql.inserElementsFromArray(sql.databank_Custom_name, [componente])
            self.accept()
        elif role == QtWidgets.QDialogButtonBox.ButtonRole.ApplyRole:
            componente = self.getComponent()
            sql.updateElement(componente, self.index)
            self.accept()
        elif role == QtWidgets.QDialogButtonBox.ButtonRole.DestructiveRole:
            componente = self.getComponent()
            if self.index == 0:
                func = sql.inserElementsFromArray
                arg = (sql.databank_Custom_name, [componente])
            elif self.index > 10000:
                func = sql.updateElement
                arg = (componente, self.index)
            if okToContinue(self, self.dirty, func, arg):
                self.accept()
            else:
                self.reject()
        else:
            self.reject()

    def setReadOnly(self, bool):
        self.formula1.setReadOnly(bool)
        self.CAS.setReadOnly(bool)
        self.alternateName.setReadOnly(bool)
        self.formula1.setReadOnly(bool)
        self.smile.setReadOnly(bool)
        self.M.setReadOnly(bool)
        self.Tc.setReadOnly(bool)
        self.Pc.setReadOnly(bool)
        self.Vc.setReadOnly(bool)
        self.w.setReadOnly(bool)
        self.SG.setReadOnly(bool)
        self.Tb.setReadOnly(bool)
        self.Tm.setReadOnly(bool)
        self.calorFormacionGas.setReadOnly(bool)
        self.energiaGibbsGas.setReadOnly(bool)
        self.cpa.setReadOnly(bool)
        self.cpb.setReadOnly(bool)
        self.cpc.setReadOnly(bool)
        self.cpd.setReadOnly(bool)
        self.cpe.setReadOnly(bool)
        self.cpf.setReadOnly(bool)
        self.MSRKa.setReadOnly(bool)
        self.MSRKb.setReadOnly(bool)
        self.Melhemm.setReadOnly(bool)
        self.Melhemn.setReadOnly(bool)
        self.SRKAPIS1.setReadOnly(bool)
        self.SRKAPIS2.setReadOnly(bool)
        self.MathiasCopemanC1.setReadOnly(bool)
        self.MathiasCopemanC2.setReadOnly(bool)
        self.MathiasCopemanC3.setReadOnly(bool)
        self.SVk0.setReadOnly(bool)
        self.SVk1.setReadOnly(bool)
        self.SVk2.setReadOnly(bool)
        self.ZVk1.setReadOnly(bool)
        self.ZVk2.setReadOnly(bool)
        self.Almeidam.setReadOnly(bool)
        self.Almeidan.setReadOnly(bool)
        self.AlmeidaG.setReadOnly(bool)
        self.PTZc.setReadOnly(bool)
        self.PTF.setReadOnly(bool)
        self.TBq1.setReadOnly(bool)
        self.TBq2.setReadOnly(bool)
        self.TBSZc.setReadOnly(bool)
        self.TBSm.setReadOnly(bool)
        self.TBSp.setReadOnly(bool)
        self.TBSd.setReadOnly(bool)
        self.SolubilityPar.setReadOnly(bool)
        self.Dipole.setReadOnly(bool)
        self.MolecularDiameter.setReadOnly(bool)
        self.NetHeat.setReadOnly(bool)
        self.GrossHeat.setReadOnly(bool)
        self.VolLiqConstant.setReadOnly(bool)
        self.API.setReadOnly(bool)
        self.wMod.setReadOnly(bool)
        self.UNIQUACArea.setReadOnly(bool)
        self.UNIQUACVolume.setReadOnly(bool)
        self.wilson.setReadOnly(bool)
        self.stiel.setReadOnly(bool)
        self.rackett.setReadOnly(bool)
        self.PolarParameter.setReadOnly(bool)
        self.EpsK.setReadOnly(bool)
        self.watson.setReadOnly(bool)

    def help(self, page):
        # First search for a local documentation build
        # By default in docs/_build/html/
        path = os.path.join(os.environ["pychemqt"], "docs") + os.sep
        if os.path.isdir(path):
            indexpath = os.path.join(path, "_build", "html", page)
            url = QtCore.QUrl(indexpath)
        else:
            url = QtCore.QUrl("http://pychemqt.readthedocs.io/")
        QtGui.QDesktopServices.openUrl(url)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = View_Component(623)
    Dialog.show()
    sys.exit(app.exec())
