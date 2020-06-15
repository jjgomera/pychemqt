#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Dialogs for fluid selection:
#   - Ui_ChooseFluid: Dialog to choose fluid for calculations
#   - DialogFilterFluid: Dialog for filter compounds family to show
#   - Dialog_InfoFluid: Dialog to show parameter of element with meos
#       - Widget_MEoS_Data: Widget to show meos data
#   - transportDialog: Dialog for transport and ancillary equations
#       - Widget_Viscosity_Data: Widget to show viscosity data
#       - Widget_Conductivity_Data: Widget to show thermal conductivity data
###############################################################################


import inspect
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from lib import meos, mEoS, unidades
from tools.codeEditor import SimplePythonEditor
from UI.widgets import Entrada_con_unidades, Tabla, QLabelMath


class Ui_ChooseFluid(QtWidgets.QDialog):
    """Dialog to choose fluid for meos plugins calculations"""
    all = True
    group = None

    def __init__(self, config=None, parent=None):
        """config: instance with project config to set initial values"""
        super(Ui_ChooseFluid, self).__init__(parent)
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Choose fluid"))
        layout = QtWidgets.QGridLayout(self)

        self.lista = QtWidgets.QListWidget()
        self.fill(mEoS.__all__)
        self.lista.itemDoubleClicked.connect(self.accept)
        layout.addWidget(self.lista, 1, 1, 5, 1)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            QtCore.Qt.Vertical)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.helpRequested.connect(self.info)
        layout.addWidget(self.buttonBox, 1, 2)

        self.widget = QtWidgets.QWidget(self)
        self.widget.setVisible(False)
        layout.addWidget(self.widget, 6, 1, 1, 2)
        gridLayout = QtWidgets.QGridLayout(self.widget)
        self.radioMEoS = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate("pychemqt", "Use MEoS equation"))
        self.radioMEoS.setChecked(True)
        gridLayout.addWidget(self.radioMEoS, 1, 1, 1, 2)
        gridLayout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Equation")+": "), 2, 1)
        self.eq = QtWidgets.QComboBox()
        gridLayout.addWidget(self.eq, 2, 2)
        self.generalized = QtWidgets.QRadioButton(
            QtWidgets.QApplication.translate(
                "pychemqt", "Use generalizated expression"))
        gridLayout.addWidget(self.generalized, 3, 1, 1, 2)
        self.radioPR = QtWidgets.QRadioButton(QtWidgets.QApplication.translate(
            "pychemqt", "Use Peng-Robinson cubic equation"))
        gridLayout.addWidget(self.radioPR, 4, 1, 1, 2)

        gridLayout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            5, 1)
        gridLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Viscosity")), 6, 1)
        self.visco = QtWidgets.QComboBox()
        gridLayout.addWidget(self.visco, 6, 2)
        gridLayout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Thermal")), 7, 1)
        self.thermal = QtWidgets.QComboBox()
        gridLayout.addWidget(self.thermal, 7, 2)
        gridLayout.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Maximum), 8, 2)

        botonFilter = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "filter.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Filter"))
        botonFilter.clicked.connect(self.filter)
        layout.addWidget(botonFilter, 3, 2, 1, 1)
        botonInfo = QtWidgets.QPushButton(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "helpAbout.png"))),
            QtWidgets.QApplication.translate("pychemqt", "Info"))
        botonInfo.clicked.connect(self.info)
        layout.addWidget(botonInfo, 4, 2, 1, 1)
        self.botonMore = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "More..."))
        self.botonMore.setCheckable(True)
        self.botonMore.clicked.connect(self.widget.setVisible)
        layout.addWidget(self.botonMore, 5, 2, 1, 1)

        self.lista.currentRowChanged.connect(self.update)
        self.radioMEoS.toggled.connect(self.eq.setEnabled)

        if config and config.has_option("MEoS", "fluid"):
            self.lista.setCurrentRow(config.getint("MEoS", "fluid"))
            self.eq.setCurrentIndex(config.getint("MEoS", "eq"))
            self.radioPR.setChecked(config.getboolean("MEoS", "PR"))
            self.generalized.setChecked(
                config.getboolean("MEoS", "Generalized"))
            self.visco.setCurrentIndex(config.getint("MEoS", "visco"))
            self.thermal.setCurrentIndex(config.getint("MEoS", "thermal"))

    def id(self):
        """Return correct id of selected fluid in mEoS.__all__ list"""
        id = self.lista.currentRow()

        # Correct id for hidden classes
        if not self.all:
            hiden = 0
            visible = 0
            for grp, boolean in zip(DialogFilterFluid.classOrder, self.group):
                module = mEoS.__getattribute__(grp)
                if boolean:
                    visible += len(module)
                else:
                    hiden += len(module)

                if visible >= id:
                    break
            # Add the element hidden above the selected one
            id += hiden
        return id

    def fill(self, compounds):
        """Fill list fluid
        compounds: List of MEoS subclasses to show"""
        self.lista.clear()
        for fluido in compounds:
            txt = fluido.name
            if fluido.synonym:
                txt += " ("+fluido.synonym+")"
            self.lista.addItem(txt)

    def filter(self):
        """Show dialog with group compound filter"""
        dlg = DialogFilterFluid(self.all, self.group)
        if dlg.exec_():
            if dlg.showAll.isChecked():
                cmps = mEoS.__all__
                self.all = True
            else:
                self.all = False
                self.group = []
                cmps = []
                for i, key in enumerate(dlg.classOrder):
                    if dlg.groups[i].isChecked():
                        cmps += mEoS.__getattribute__(key)
                        self.group.append(True)
                    else:
                        self.group.append(False)
            self.fill(cmps)

    def info(self):
        """Show info dialog for fluid"""
        dialog = Dialog_InfoFluid(mEoS.__all__[self.lista.currentRow()])
        dialog.exec_()

    def update(self, indice):
        """Update data when selected fluid change"""
        fluido = mEoS.__all__[indice]
        self.eq.clear()
        for eq in fluido.eq:
            self.eq.addItem(eq["__name__"])

        self.visco.clear()
        if fluido._viscosity is not None:
            self.visco.setEnabled(True)
            for eq in fluido._viscosity:
                self.visco.addItem(eq["__name__"])
        else:
            self.visco.addItem(
                QtWidgets.QApplication.translate("pychemqt", "Undefined"))
            self.visco.setEnabled(False)

        self.thermal.clear()
        if fluido._thermal is not None:
            self.thermal.setEnabled(True)
            for eq in fluido._thermal:
                self.thermal.addItem(eq["__name__"])
        else:
            self.thermal.addItem(
                QtWidgets.QApplication.translate("pychemqt", "Undefined"))
            self.thermal.setEnabled(False)


class DialogFilterFluid(QtWidgets.QDialog):
    """Dialog for filter compounds family to show"""
    text = {
        "Nobles": QtWidgets.QApplication.translate("pychemqt", "Noble gases"),
        "Gases": QtWidgets.QApplication.translate("pychemqt", "Gases"),
        "Alkanes": QtWidgets.QApplication.translate("pychemqt", "Alkanes"),
        "Naphthenes": QtWidgets.QApplication.translate(
            "pychemqt", "Naphthenes"),
        "Alkenes": QtWidgets.QApplication.translate("pychemqt", "Alkenes"),
        "Heteroatom": QtWidgets.QApplication.translate(
            "pychemqt", "Heteroatom"),
        "CFCs": QtWidgets.QApplication.translate("pychemqt", "CFCs"),
        "Siloxanes": QtWidgets.QApplication.translate("pychemqt", "Siloxanes"),
        "PseudoCompounds": QtWidgets.QApplication.translate(
            "pychemqt", "Pseudo Compounds")}
    classOrder = ["Nobles", "Gases", "Alkanes", "Naphthenes", "Alkenes",
                  "Heteroatom", "CFCs", "Siloxanes", "PseudoCompounds"]

    def __init__(self, all=True, group=None, parent=None):
        super(DialogFilterFluid, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Filter fluids families to show"))
        layout = QtWidgets.QGridLayout(self)
        self.showAll = QtWidgets.QCheckBox(QtWidgets.QApplication.translate(
            "pychemqt", "Show All"))
        layout.addWidget(self.showAll, 1, 1)

        widget = QtWidgets.QWidget()
        layout.addWidget(widget, 2, 1)
        lyt = QtWidgets.QVBoxLayout(widget)
        self.groups = []
        for name in self.classOrder:
            checkBox = QtWidgets.QCheckBox(self.text[name])
            lyt.addWidget(checkBox)
            self.groups.append(checkBox)

        self.showAll.toggled.connect(widget.setDisabled)

        self.showAll.setChecked(all)
        if group is not None:
            for boolean, checkBox in zip(group, self.groups):
                checkBox.setChecked(boolean)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 3, 1)


class Dialog_InfoFluid(QtWidgets.QDialog):
    """Dialog to show parameter of element with meos"""
    def __init__(self, element, parent=None):
        """element: class of element to show info"""
        super(Dialog_InfoFluid, self).__init__(parent)
        layout = QtWidgets.QGridLayout(self)
        self.element = element

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Name")+":"), 1, 1)
        self.name = QtWidgets.QLabel()
        layout.addWidget(self.name, 1, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "R name")+":"), 2, 1)
        self.r_name = QtWidgets.QLabel()
        layout.addWidget(self.r_name, 2, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Formula")+":"), 3, 1)
        self.formula = QtWidgets.QLabel()
        layout.addWidget(self.formula, 3, 2)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "CAS number")+":"), 4, 1)
        self.CAS = QtWidgets.QLabel()
        layout.addWidget(self.CAS, 4, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            30, 30, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 1, 3, 3, 1)

        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "M")+":"), 1, 4)
        self.M = Entrada_con_unidades(
            float, textounidad="g/mol", readOnly=True)
        layout.addWidget(self.M, 1, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Tc")+":"), 2, 4)
        self.Tc = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tc, 2, 5)
        layout.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Pc")+":"), 3, 4)
        self.Pc = Entrada_con_unidades(unidades.Pressure, readOnly=True)
        layout.addWidget(self.Pc, 3, 5)
        layout.addWidget(QtWidgets.QLabel("ρc"+":"), 4, 4)
        self.rhoc = Entrada_con_unidades(
            unidades.Density, "DenGas", readOnly=True)
        layout.addWidget(self.rhoc, 4, 5)
        layout.addItem(QtWidgets.QSpacerItem(
            30, 30, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 1, 6, 3, 1)

        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "T triple")+":"), 1, 7)
        self.Tt = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tt, 1, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "T boiling")+":"), 2, 7)
        self.Tb = Entrada_con_unidades(unidades.Temperature, readOnly=True)
        layout.addWidget(self.Tb, 2, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Dipole moment")+":"), 3, 7)
        self.momento = Entrada_con_unidades(
            unidades.DipoleMoment, readOnly=True)
        layout.addWidget(self.momento, 3, 8)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "F acentric")+":"), 4, 7)
        self.f_acent = Entrada_con_unidades(float, readOnly=True)
        layout.addWidget(self.f_acent, 4, 8)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            5, 1)
        layout.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Equation")+": "), 6, 1)
        self.eq = QtWidgets.QComboBox()
        layout.addWidget(self.eq, 6, 2, 1, 7)
        self.stacked = QtWidgets.QStackedWidget()
        layout.addWidget(self.stacked, 7, 1, 1, 8)
        self.eq.currentIndexChanged.connect(self.stacked.setCurrentIndex)

        self.moreButton = QtWidgets.QPushButton(
            QtWidgets.QApplication.translate("pychemqt", "Others"))
        self.moreButton.clicked.connect(self.more)
        layout.addWidget(self.moreButton, 9, 1)
        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        btBox.clicked.connect(self.reject)
        layout.addWidget(btBox, 9, 2, 1, 7)

        self.fill(element)

    def fill(self, element):
        """Fill values"""
        self.name.setText(element.name)
        self.r_name.setText(element.synonym)
        self.formula.setText(element.formula)
        self.CAS.setText(element.CASNumber)
        self.M.setValue(element.M)
        self.Tc.setValue(element.Tc)
        self.Pc.setValue(element.Pc)
        self.rhoc.setValue(element.rhoc)
        self.Tb.setValue(element.Tb)
        self.Tt.setValue(element.Tt)
        self.momento.setValue(element.momentoDipolar)
        self.f_acent.setValue(element.f_acent)

        for eq in element.eq:
            widget = Widget_MEoS_Data(eq)
            self.stacked.addWidget(widget)
            self.eq.addItem(eq["__name__"])

    def more(self):
        """Show parameter for transport and ancillary equations"""
        dialog = transportDialog(self.element, parent=self)
        dialog.show()


class Widget_MEoS_Data(QtWidgets.QWidget):
    """Widget to show meos data"""
    def __init__(self, eq, parent=None):
        """eq: dict with equation parameter"""
        super(Widget_MEoS_Data, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)
        txt = " ".join((eq["__doi__"]["autor"], eq["__doi__"]["title"],
                        eq["__doi__"]["ref"]))
        ref = QtWidgets.QLabel(txt)
        ref.setWordWrap(True)
        gridLayout.addWidget(ref, 1, 1)

        tabWidget = QtWidgets.QTabWidget()
        gridLayout.addWidget(tabWidget, 3, 1)

        # Cp tab
        if "ao_log" in eq["cp"]:
            # Cp0 form
            tab1 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab1, QtWidgets.QApplication.translate("pychemqt", "Phi0"))
            gridLayout_Ideal = QtWidgets.QGridLayout(tab1)
            mathTex = r"$\alpha^o=\ln\delta + c_o\ln\tau + \sum c_i\tau^{n_i} "
            mathTex += r"+ \sum_j m_j \ln (1-e^{-\theta_j\tau}) + "
            mathTex += r"\sum_k l_k\ln|\sinh(\psi_k\tau)| - "
            mathTex += r"\sum l_k\ln|\cosh(\psi_k\tau)|$"
            label = QLabelMath(mathTex)
            gridLayout_Ideal.addWidget(label, 1, 1, 1, 3)
            self.Tabla_Cp_poly = Tabla(
                2, horizontalHeader=["n", "d"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_poly, 2, 1)
            self.Tabla_Cp_exp = Tabla(
                2, horizontalHeader=["m", "θ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_exp, 2, 2)
            self.Tabla_Cp_hyp = Tabla(
                2, horizontalHeader=["l", "ψ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_hyp, 2, 3)

        else:
            # Phi0 form
            tab1 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab1, QtWidgets.QApplication.translate("pychemqt", "Cp"))
            gridLayout_Ideal = QtWidgets.QGridLayout(tab1)
            mathTex = r"$\frac{C_p^o}{R}=\sum n_i\tau^{d_i}+"
            mathTex += r"\sum m_j(\theta_j\tau)^2\frac{e^{\theta_j\tau}}"
            mathTex += r"{(e^{\theta_j\tau}-1)^2}"
            mathTex += r"+\sum l_k\left(\frac{\phi_k\tau}"
            mathTex += r"{\sinh(\phi_k\tau)}\right)^2"
            mathTex += r"+\sum l_k\left(\frac{\phi_k\tau}"
            mathTex += r"{\cosh(\phi_k\tau)}\right)^2$"
            label = QLabelMath(mathTex)
            gridLayout_Ideal.addWidget(label, 1, 1, 1, 3)
            self.Tabla_Cp_poly = Tabla(
                2, horizontalHeader=["n", "d"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_poly, 2, 1)
            self.Tabla_Cp_exp = Tabla(
                2, horizontalHeader=["m", "θ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_exp, 2, 2)
            self.Tabla_Cp_hyp = Tabla(
                2, horizontalHeader=["l", "ψ"], stretch=False, readOnly=True)
            gridLayout_Ideal.addWidget(self.Tabla_Cp_hyp, 2, 3)

        if eq["__type__"] == "Helmholtz":
            mathTex = r"$\alpha = \alpha^o+\alpha_{Pol}^r+\alpha_{Exp}^r+"
            mathTex += r"\alpha_{GBS}^r+\alpha_{NA}^r+\alpha_{HE}^r$"
            label = QLabelMath(mathTex)
            gridLayout.addWidget(label, 2, 1)

            # Polinomial tab
            tab2 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab2,
                QtWidgets.QApplication.translate("pychemqt", "Polinomial"))
            gridLayout_pol = QtWidgets.QGridLayout(tab2)
            mathTex = r"$\alpha_{Pol}^r=\sum_i n_i\tau^{t_i}\delta^{d_i}$"
            label = QLabelMath(mathTex)
            gridLayout_pol.addWidget(label, 1, 1)
            self.Tabla_lineal = Tabla(
                3, horizontalHeader=["n", "t", "d"], stretch=False,
                readOnly=True)
            gridLayout_pol.addWidget(self.Tabla_lineal, 2, 1)

            # Exponencial tab
            tab3 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab3,
                QtWidgets.QApplication.translate("pychemqt", "Exponential"))
            gridLayout_Exp = QtWidgets.QGridLayout(tab3)
            mathTex = r"$\alpha_{Exp}^r=\sum_i n_i\tau^{t_i}\delta^{d_i}"
            mathTex += r"e^{-\gamma_i\delta^{c_i}}$"
            label = QLabelMath(mathTex)
            gridLayout_Exp.addWidget(label, 1, 1)
            self.Tabla_exponential = Tabla(
                5, horizontalHeader=["n", "t", "d", "γ", "c"],
                stretch=False, readOnly=True)
            gridLayout_Exp.addWidget(self.Tabla_exponential, 2, 1)

            # Gaussian tab
            tab4 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab4, QtWidgets.QApplication.translate("pychemqt", "Gaussian"))
            gridLayout_gauss = QtWidgets.QGridLayout(tab4)
            mathTex = r"$\alpha_{GBS}^r=\sum_i n_i\tau^{t_i}\delta^{d_i}"
            mathTex += r"e^{-\alpha_i\left(\delta-\epsilon_i\right)^2"
            mathTex += r"-\beta\left(\tau-\gamma_i\right)^2}$"
            label = QLabelMath(mathTex)
            gridLayout_gauss.addWidget(label, 1, 1)
            self.Tabla_gauss = Tabla(
                7, horizontalHeader=["n", "t", "d", "η", "ε", "β", "γ"],
                stretch=False, readOnly=True)
            gridLayout_gauss.addWidget(self.Tabla_gauss, 2, 1)

            # Non analytic tab
            tab5 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab5,
                QtWidgets.QApplication.translate("pychemqt", "Non analytic"))
            gridLayout_NA = QtWidgets.QGridLayout(tab5)
            mathTex = r"$\alpha_{NA}^r=\sum_i n_i\delta\Delta^{b_i}"
            mathTex += r"e^{-C_i\left(\delta-1\right)^2-D_i"
            mathTex += r"\left(\tau-1\right)^2}$"
            label = QLabelMath(mathTex)
            gridLayout_NA.addWidget(label, 1, 1)
            mathTex = r"$\Delta = \left(1-\tau+A_i\left(\left(\delta-1\right)"
            mathTex += r"^2\right)^{1/2\beta_i}\right)^2+B_i\left(\left(\delta"
            mathTex += r"-1\right)^2\right)^{a_i}$"
            label2 = QLabelMath(mathTex)
            gridLayout_NA.addWidget(label2, 2, 1)
            self.Tabla_noanalytic = Tabla(
                8, horizontalHeader=["n", "a", "b", "A", "B", "C", "D", "β"],
                stretch=False, readOnly=True)
            gridLayout_NA.addWidget(self.Tabla_noanalytic, 3, 1)

            # Hard Sphere tab
            tab6 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab6,
                QtWidgets.QApplication.translate("pychemqt", "Hard Sphere"))
            gridLayout_HE = QtWidgets.QGridLayout(tab6)
            mathTex = r"$\alpha_{HE}^r=(\varphi^2-1)\ln(1-\xi)+\frac"
            mathTex += r"{(\varphi^2+3\varphi)\xi-3\varphi\xi^2}{(1-\xi)^2}$"
            label = QLabelMath(mathTex)
            gridLayout_HE.addWidget(label, 1, 1, 1, 2)
            gridLayout_HE.addWidget(QtWidgets.QLabel("φ:"), 2, 1)
            self.fi = Entrada_con_unidades(float, readOnly=True)
            gridLayout_HE.addWidget(self.fi, 2, 2)
            gridLayout_HE.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 3, 1, 1, 2)

        elif eq["__type__"] == "MBWR":
            # Pestaña MBWR
            tab2 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab2, QtWidgets.QApplication.translate("pychemqt", "MBWR"))
            gridLayout_MBWR = QtWidgets.QGridLayout(tab2)
            mathTex = r"$P=\rho RT+\sum_{n=2}^{9}\alpha_n\rho^n + "
            mathTex += r"e^{-\delta^2} \sum_{10}^{15} \alpha_n"
            mathTex += r"\rho^{2n-17}$"
            label = QLabelMath(mathTex)
            gridLayout_MBWR.addWidget(label, 1, 1)
            self.Tabla_MBWR = Tabla(
                1, horizontalHeader=["b"], stretch=False, readOnly=True)
            gridLayout_MBWR.addWidget(self.Tabla_MBWR, 2, 1)

        self.fill(eq)

    def fill(self, eq):
        format = {"format": 1, "total": 5}

        if "ao_log" in eq["cp"]:
            # Phi_o term
            self.Tabla_Cp_poly.setColumn(
                0, eq["cp"]["ao_pow"], **format)
            self.Tabla_Cp_poly.setColumn(1, eq["cp"]["pow"], **format)
            self.Tabla_Cp_poly.resizeColumnsToContents()
            if "ao_exp" in eq["cp"]:
                self.Tabla_Cp_exp.setColumn(0, eq["cp"]["ao_exp"], **format)
                self.Tabla_Cp_exp.setColumn(1, eq["cp"]["titao"], **format)
                self.Tabla_Cp_exp.resizeColumnsToContents()
            if "hyp" in eq["cp"]:
                self.Tabla_Cp_hyp.setColumn(0, eq["cp"]["ao_hyp"], **format)
                self.Tabla_Cp_hyp.setColumn(1, eq["cp"]["hyp"], **format)
                self.Tabla_Cp_hyp.resizeColumnsToContents()
        else:
            # Cp term
            an = eq["cp"].get("an", [])
            t = eq["cp"].get("pow", [])
            ao = eq["cp"].get("ao", 0)
            if ao:
                an.insert(0, ao)
                t.insert(0, 0)

            if an:
                self.Tabla_Cp_poly.setColumn(0, an, **format)
                self.Tabla_Cp_poly.setColumn(1, t, **format)
                self.Tabla_Cp_poly.resizeColumnsToContents()

            if "ao_exp" in eq["cp"]:
                self.Tabla_Cp_exp.setColumn(0, eq["cp"]["ao_exp"], **format)
                self.Tabla_Cp_exp.setColumn(1, eq["cp"]["exp"], **format)
                self.Tabla_Cp_exp.resizeColumnsToContents()

            if "hyp" in eq["cp"]:
                self.Tabla_Cp_hyp.setColumn(0, eq["cp"]["ao_hyp"], **format)
                self.Tabla_Cp_hyp.setColumn(1, eq["cp"]["hyp"], **format)
                self.Tabla_Cp_hyp.resizeColumnsToContents()

        if eq["__type__"] == "Helmholtz":
            if eq.get("nr1", []):
                self.Tabla_lineal.setColumn(0, eq["nr1"], **format)
                self.Tabla_lineal.setColumn(1, eq["t1"], **format)
                self.Tabla_lineal.setColumn(2, eq["d1"], **format)
            if eq.get("nr2", []):
                self.Tabla_exponential.setColumn(0, eq["nr2"], **format)
                self.Tabla_exponential.setColumn(1, eq["t2"], **format)
                self.Tabla_exponential.setColumn(2, eq["d2"], **format)
                self.Tabla_exponential.setColumn(3, eq["gamma2"], **format)
                self.Tabla_exponential.setColumn(4, eq["c2"], **format)
            if eq.get("nr3", []):
                self.Tabla_gauss.setColumn(0, eq["nr3"], **format)
                self.Tabla_gauss.setColumn(1, eq["t3"], **format)
                self.Tabla_gauss.setColumn(2, eq["d3"], **format)
                self.Tabla_gauss.setColumn(3, eq["alfa3"], **format)
                self.Tabla_gauss.setColumn(4, eq["beta3"], **format)
                self.Tabla_gauss.setColumn(5, eq["gamma3"], **format)
                self.Tabla_gauss.setColumn(6, eq["epsilon3"], **format)
            if eq.get("nr4", []):
                self.Tabla_noanalytic.setColumn(0, eq["nr4"], **format)
                self.Tabla_noanalytic.setColumn(1, eq["a4"], **format)
                self.Tabla_noanalytic.setColumn(2, eq["b4"], **format)
                self.Tabla_noanalytic.setColumn(3, eq["A"], **format)
                self.Tabla_noanalytic.setColumn(4, eq["B"], **format)
                self.Tabla_noanalytic.setColumn(5, eq["C"], **format)
                self.Tabla_noanalytic.setColumn(6, eq["D"], **format)
                self.Tabla_noanalytic.setColumn(7, eq["beta4"], **format)
            self.Tabla_lineal.resizeColumnsToContents()
            self.Tabla_exponential.resizeColumnsToContents()
            self.Tabla_gauss.resizeColumnsToContents()
            self.Tabla_noanalytic.resizeColumnsToContents()

        elif eq["__type__"] == "MBWR":
            self.Tabla_MBWR.setColumn(0, eq["b"][1:], **format)
            self.Tabla_MBWR.resizeColumnsToContents()


class transportDialog(QtWidgets.QDialog):
    """Dialog to show parameters for transport and ancillary equations"""
    def __init__(self, element, parent=None):
        super(transportDialog, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)
        self.element = element

        tabWidget = QtWidgets.QTabWidget()
        gridLayout.addWidget(tabWidget, 1, 1)

        # Tab viscosity
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab3, QtWidgets.QApplication.translate("pychemqt", "Viscosity"))
        gridLayout_viscosity = QtWidgets.QGridLayout(tab3)

        gridLayout_viscosity.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate(
                "pychemqt", "Equation")+": "), 1, 1)
        self.eqVisco = QtWidgets.QComboBox()
        gridLayout_viscosity.addWidget(self.eqVisco, 1, 2)
        gridLayout_viscosity.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Fixed), 1, 3)
        self.stackedVisco = QtWidgets.QStackedWidget()
        gridLayout_viscosity.addWidget(self.stackedVisco, 2, 1, 1, 3)
        self.eqVisco.currentIndexChanged.connect(
            self.stackedVisco.setCurrentIndex)

        if element._viscosity is not None:
            for eq in element._viscosity:
                widget = Widget_Viscosity_Data(element, eq)
                self.stackedVisco.addWidget(widget)
                self.eqVisco.addItem(eq["__name__"])
        else:
            self.eqVisco.addItem(QtWidgets.QApplication.translate(
                "pychemqt", "Not Implemented"))

        # Tab thermal conductivity
        tab4 = QtWidgets.QWidget()
        tabWidget.addTab(tab4,
                         QtWidgets.QApplication.translate(
                             "pychemqt", "Thermal Conductivity"))
        gridLayout_conductivity = QtWidgets.QGridLayout(tab4)

        gridLayout_conductivity.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Equation")+": "),
            1, 1)
        self.eqThermo = QtWidgets.QComboBox()
        gridLayout_conductivity.addWidget(self.eqThermo, 1, 2)
        gridLayout_conductivity.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Fixed), 1, 3)
        self.stackedThermo = QtWidgets.QStackedWidget()
        gridLayout_conductivity.addWidget(self.stackedThermo, 2, 1, 1, 3)
        self.eqThermo.currentIndexChanged.connect(
            self.stackedThermo.setCurrentIndex)

        if element._thermal is not None:
            for eq in element._thermal:
                widget = Widget_Conductivity_Data(element, eq)
                self.stackedThermo.addWidget(widget)
                self.eqThermo.addItem(eq["__name__"])
        else:
            self.eqThermo.addItem(QtWidgets.QApplication.translate(
                "pychemqt", "Not Implemented"))

        # Tab dielectric constant
        tab1 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab1, QtWidgets.QApplication.translate("pychemqt", "Dielectric"))
        gridLayout_dielectric = QtWidgets.QGridLayout(tab1)

        if element._Dielectric != meos.MEoS._Dielectric:
            label = QtWidgets.QLabel(element._Dielectric.__doc__)
            label.setWordWrap(True)
            gridLayout_dielectric.addWidget(label, 1, 1)
            self.codigo_Dielectric = SimplePythonEditor()
            self.codigo_Dielectric.setText(
                inspect.getsource(element._Dielectric))
            gridLayout_dielectric.addWidget(self.codigo_Dielectric, 2, 1)
        elif element._dielectric:
            label = QtWidgets.QLabel(element._Dielectric.__doc__)
            label.setWordWrap(True)
            gridLayout_dielectric.addWidget(label, 1, 1)

            self.Table_Dielectric = Tabla(
                1, verticalHeader=True, filas=5, stretch=False, readOnly=True)
            gridLayout_dielectric.addWidget(self.Table_Dielectric, 2, 1)
            i = 0
            for key, valor in element._dielectric.items():
                self.Table_Dielectric.setVerticalHeaderItem(
                    i, QtWidgets.QTableWidgetItem(key))
                self.Table_Dielectric.setItem(
                    0, i, QtWidgets.QTableWidgetItem(str(valor)))
                i += 1
            self.Table_Dielectric.resizeColumnsToContents()
        else:
            gridLayout_dielectric.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_dielectric.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab surface tension
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab2,
            QtWidgets.QApplication.translate("pychemqt", "Surface Tension"))
        gridLayout_surface = QtWidgets.QGridLayout(tab2)

        if element._Surface != meos.MEoS._Surface:
            label = QtWidgets.QLabel(element._Surface.__doc__)
            label.setWordWrap(True)
            gridLayout_surface.addWidget(label, 1, 1)
            self.codigo_Surface = SimplePythonEditor()
            self.codigo_Surface.setText(inspect.getsource(element._Surface))
            gridLayout_surface.addWidget(self.codigo_Surface, 2, 1)
        elif element._surface:
            label = QtWidgets.QLabel(element._Surface.__doc__)
            label.setWordWrap(True)
            gridLayout_surface.addWidget(label, 1, 1)

            self.Table_Surface = Tabla(
                2, horizontalHeader=["σ", "n"], verticalHeader=True,
                stretch=False, readOnly=True)
            self.Table_Surface.setColumn(0, element._surface["sigma"])
            self.Table_Surface.setColumn(1, element._surface["exp"])
            gridLayout_surface.addWidget(self.Table_Surface, 2, 1)
            self.Table_Surface.resizeColumnsToContents()
        else:
            gridLayout_surface.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_surface.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab liquid density
        tab5 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab5,
            QtWidgets.QApplication.translate("pychemqt", "Liquid Density"))
        gridLayout_liquid_density = QtWidgets.QGridLayout(tab5)

        if element._Liquid_Density != meos.MEoS._Liquid_Density:
            label = QtWidgets.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_liquid_density.addWidget(label, 1, 1)
            self.codigo_Liquid_Density = SimplePythonEditor()
            self.codigo_Liquid_Density.setText(
                inspect.getsource(element._Liquid_Density))
            gridLayout_liquid_density.addWidget(
                self.codigo_Liquid_Density, 2, 1)
        elif element._liquid_Density:
            label = QtWidgets.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_liquid_density.addWidget(label, 1, 1)

            self.Table_Liquid_Density = Tabla(
                2, horizontalHeader=["n", "t"],
                verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Liquid_Density.setColumn(
                0, element._liquid_Density["n"])
            self.Table_Liquid_Density.setColumn(
                1, element._liquid_Density["t"])
            gridLayout_liquid_density.addWidget(
                self.Table_Liquid_Density, 2, 1)
            self.Table_Liquid_Density.resizeColumnsToContents()
        else:
            gridLayout_liquid_density.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_liquid_density.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab vapor density
        tab6 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab6,
            QtWidgets.QApplication.translate("pychemqt", "Vapor Density"))
        gridLayout_vapor_density = QtWidgets.QGridLayout(tab6)

        if element._Vapor_Density != meos.MEoS._Vapor_Density:
            label = QtWidgets.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_density.addWidget(label, 1, 1)
            self.codigo_Vapor_Density = SimplePythonEditor()
            self.codigo_Vapor_Density.setText(
                inspect.getsource(element._Vapor_Density))
            gridLayout_vapor_density.addWidget(self.codigo_Vapor_Density, 2, 1)
        elif element._vapor_Density:
            label = QtWidgets.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_density.addWidget(label, 1, 1)

            self.Table_Vapor_Density = Tabla(
                2, horizontalHeader=["n", "t"], verticalHeader=True,
                stretch=False, readOnly=True)
            self.Table_Vapor_Density.setColumn(0, element._vapor_Density["n"])
            self.Table_Vapor_Density.setColumn(
                1, element._vapor_Density["t"])
            gridLayout_vapor_density.addWidget(self.Table_Vapor_Density, 2, 1)
            self.Table_Vapor_Density.resizeColumnsToContents()
        else:
            gridLayout_vapor_density.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_vapor_density.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab vapor presure
        tab7 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab7,
            QtWidgets.QApplication.translate("pychemqt", "Vapor Pressure"))
        gridLayout_vapor_pressure = QtWidgets.QGridLayout(tab7)

        if element._Vapor_Pressure != meos.MEoS._Vapor_Pressure:
            label = QtWidgets.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_pressure.addWidget(label, 1, 1)
            self.codigo_Vapor_Pressure = SimplePythonEditor()
            self.codigo_Vapor_Pressure.setText(
                inspect.getsource(element._Vapor_Pressure))
            gridLayout_vapor_pressure.addWidget(
                self.codigo_Vapor_Pressure, 2, 1)
        elif element._vapor_Pressure:
            label = QtWidgets.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_vapor_pressure.addWidget(label, 1, 1)

            self.Table_Vapor_Pressure = Tabla(
                2, horizontalHeader=["n", "t"], verticalHeader=True,
                stretch=False, readOnly=True)
            self.Table_Vapor_Pressure.setColumn(
                0, element._vapor_Pressure["n"])
            self.Table_Vapor_Pressure.setColumn(
                1, element._vapor_Pressure["t"])
            gridLayout_vapor_pressure.addWidget(
                self.Table_Vapor_Pressure, 2, 1)
            self.Table_Vapor_Pressure.resizeColumnsToContents()
        else:
            gridLayout_vapor_pressure.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_vapor_pressure.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab melting presure
        tab8 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab8,
            QtWidgets.QApplication.translate("pychemqt", "Melting Pressure"))
        gridLayout_melting_pressure = QtWidgets.QGridLayout(tab8)

        if element._Melting_Pressure != meos.MEoS._Melting_Pressure:
            label = QtWidgets.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_melting_pressure.addWidget(label, 1, 1)
            self.codigo_Melting_Pressure = SimplePythonEditor()
            self.codigo_Melting_Pressure.setText(
                inspect.getsource(element._Melting_Pressure))
            gridLayout_melting_pressure.addWidget(
                self.codigo_Melting_Pressure, 2, 1)
        elif element._melting:
            label = QtWidgets.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout_melting_pressure.addWidget(label, 1, 1)

            self.Table_Melting_Pressure = Tabla(
                6, horizontalHeader=["a1", "n1", "a2", "n2", "a3", "n3"],
                verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Melting_Pressure.setColumn(0, element._melting["a1"])
            self.Table_Melting_Pressure.setColumn(1, element._melting["exp1"])
            self.Table_Melting_Pressure.setColumn(2, element._melting["a2"])
            self.Table_Melting_Pressure.setColumn(3, element._melting["exp2"])
            self.Table_Melting_Pressure.setColumn(4, element._melting["a3"])
            self.Table_Melting_Pressure.setColumn(5, element._melting["exp3"])
            gridLayout_melting_pressure.addWidget(
                self.Table_Melting_Pressure, 2, 1)
            self.Table_Melting_Pressure.resizeColumnsToContents()
        else:
            gridLayout_melting_pressure.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout_melting_pressure.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab sublimation presure
        tab9 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab9,
            QtWidgets.QApplication.translate(
                "pychemqt", "Sublimation Pressure"))
        gridLayout__sublimation_pressure = QtWidgets.QGridLayout(tab9)

        if element._Sublimation_Pressure != meos.MEoS._Sublimation_Pressure:
            label = QtWidgets.QLabel(element._Sublimation_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout__sublimation_pressure.addWidget(label, 1, 1)
            self.codigo_Sublimation_Pressure = SimplePythonEditor()
            self.codigo_Sublimation_Pressure.setText(
                inspect.getsource(element._Sublimation_Pressure))
            gridLayout__sublimation_pressure.addWidget(
                self.codigo_Sublimation_Pressure, 2, 1)
        elif element._sublimation:
            label = QtWidgets.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            gridLayout__sublimation_pressure.addWidget(label, 1, 1)

            self.Table_Sublimation_Pressure = Tabla(
                6, horizontalHeader=["a1", "n1", "a2", "n2", "a3", "n3"],
                verticalHeader=True, stretch=False, readOnly=True)
            self.Table_Sublimation_Pressure.setColumn(
                0, element._sublimation["a1"])
            self.Table_Sublimation_Pressure.setColumn(
                1, element._sublimation["exp1"])
            self.Table_Sublimation_Pressure.setColumn(
                2, element._sublimation["a2"])
            self.Table_Sublimation_Pressure.setColumn(
                3, element._sublimation["exp2"])
            self.Table_Sublimation_Pressure.setColumn(
                4, element._sublimation["a3"])
            self.Table_Sublimation_Pressure.setColumn(
                5, element._sublimation["exp3"])
            gridLayout__sublimation_pressure.addWidget(
                self.Table_Sublimation_Pressure, 2, 1)
            self.Table_Sublimation_Pressure.resizeColumnsToContents()
        else:
            gridLayout__sublimation_pressure.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            gridLayout__sublimation_pressure.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        # Tab Peng-Robinson
        tab10 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab10,
            QtWidgets.QApplication.translate("pychemqt", "Peng-Robinson"))
        gridLayout_PengRobinson = QtWidgets.QGridLayout(tab10)

        if element._PR:
            label = QtWidgets.QLabel(element._PengRobinson.__doc__)
            label.setWordWrap(True)
            gridLayout_PengRobinson.addWidget(label, 1, 1, 1, 3)
            gridLayout_PengRobinson.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Fixed,
                QtWidgets.QSizePolicy.Fixed), 2, 1, 1, 3)
            gridLayout_PengRobinson.addWidget(QtWidgets.QLabel("C"), 3, 1)
            self.PR = Entrada_con_unidades(
                float, decimales=6, value=element._PR, readOnly=True)
            gridLayout_PengRobinson.addWidget(self.PR, 3, 2)
            gridLayout_PengRobinson.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 4, 1, 1, 3)
        else:
            gridLayout_PengRobinson.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "No Peneloux correction")), 1, 1)
            gridLayout_PengRobinson.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 2, 1)

        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Close)
        self.buttonBox.clicked.connect(self.reject)
        gridLayout.addWidget(self.buttonBox, 2, 1)


class Widget_Viscosity_Data(QtWidgets.QWidget):
    """Widget to show viscosity data"""
    def __init__(self, element, eq, parent=None):
        """
        element: element class for code extract
        eq: dict with viscosity parameter"""
        super(Widget_Viscosity_Data, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)
        if eq["eq"] == 0:
            txt = element.__getattribute__(element, eq["method"]).__doc__
        else:
            txt = " ".join((eq["__doi__"]["autor"], eq["__doi__"]["title"],
                            eq["__doi__"]["ref"]))
        ref = QtWidgets.QLabel(txt)
        ref.setWordWrap(True)
        gridLayout.addWidget(ref, 1, 1, 1, 3)

        if eq["eq"] == 0:
            # Hardcoded method, show code
            self.codigo_Viscosity = SimplePythonEditor()
            code = ""
            for method in eq.get("__code__", ()):
                code += inspect.getsource(method)
                code += os.linesep
            code += inspect.getsource(
                element.__getattribute__(element, eq["method"]))
            self.codigo_Viscosity.setText(code)
            gridLayout.addWidget(self.codigo_Viscosity, 2, 1, 1, 3)
        elif eq["eq"] == 1:
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            self.ek = Entrada_con_unidades(
                float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(self.ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            self.sigma = Entrada_con_unidades(
                float, value=eq["sigma"], readOnly=True)
            gridLayout.addWidget(self.sigma, 5, 2)
            tab = QtWidgets.QTabWidget()
            gridLayout.addWidget(tab, 6, 1, 1, 3)

            # Integral collision
            self.Tabla_Collision = Tabla(
                1, horizontalHeader=["b"], stretch=False, readOnly=True)
            if "collision" in eq:
                self.Tabla_Collision.setColumn(0, eq["collision"])
                self.Tabla_Collision.resizeColumnsToContents()
            else:
                self.Tabla_Collision.setDisabled(True)
            tab.addTab(
                self.Tabla_Collision,
                QtWidgets.QApplication.translate("pychemqt", "Collision"))

            # Virial
            self.Tabla_Virial = Tabla(
                2, horizontalHeader=["n", "t"], stretch=False, readOnly=True)
            if "n_virial" in eq:
                self.Tabla_Virial.setColumn(0, eq["n_virial"])
                self.Tabla_Virial.setColumn(1, eq["t_virial"])
                self.Tabla_Virial.resizeColumnsToContents()
            else:
                self.Tabla_Virial.setDisabled(True)
            tab.addTab(self.Tabla_Virial,
                       QtWidgets.QApplication.translate("pychemqt", "Virial"))

            # Close-packed
            self.Tabla_Packed = Tabla(
                2, horizontalHeader=["n", "t"], stretch=False, readOnly=True)
            if "n_packed" in eq:
                self.Tabla_Packed.setColumn(0, eq["n_packed"])
                self.Tabla_Packed.setColumn(1, eq["t_packed"])
                self.Tabla_Packed.resizeColumnsToContents()
            else:
                self.Tabla_Packed.setDisabled(True)
            tab.addTab(self.Tabla_Packed, QtWidgets.QApplication.translate(
                "pychemqt", "Close-packed density"))

            # polynomial term
            self.Tabla_Visco1 = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "n_poly" in eq:
                self.Tabla_Visco1.setColumn(0, eq["n_poly"])
                self.Tabla_Visco1.setColumn(1, eq["t_poly"])
                self.Tabla_Visco1.setColumn(2, eq["d_poly"])
                self.Tabla_Visco1.setColumn(3, eq["g_poly"])
                self.Tabla_Visco1.setColumn(4, eq["c_poly"])
                self.Tabla_Visco1.resizeColumnsToContents()
            else:
                self.Tabla_Visco1.setDisabled(True)
            tab.addTab(
                self.Tabla_Visco1,
                QtWidgets.QApplication.translate("pychemqt", "Polinomial"))

            # numerator of rational poly
            self.Tabla_numerator = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "n_num" in eq:
                self.Tabla_numerator.setColumn(0, eq["n_num"])
                self.Tabla_numerator.setColumn(1, eq["t_num"])
                self.Tabla_numerator.setColumn(2, eq["d_num"])
                self.Tabla_numerator.setColumn(3, eq["c_num"])
                self.Tabla_numerator.setColumn(4, eq["g_num"])
                self.Tabla_numerator.resizeColumnsToContents()
            else:
                self.Tabla_numerator.setDisabled(True)
            tab.addTab(
                self.Tabla_numerator,
                QtWidgets.QApplication.translate("pychemqt", "Numerator"))

            # denominator of rational poly
            self.Tabla_denominator = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "n_den" in eq:
                self.Tabla_denominator.setColumn(0, eq["n_den"])
                self.Tabla_denominator.setColumn(1, eq["t_den"])
                self.Tabla_denominator.setColumn(2, eq["d_den"])
                self.Tabla_denominator.setColumn(3, eq["c_den"])
                self.Tabla_denominator.setColumn(4, eq["g_den"])
                self.Tabla_denominator.resizeColumnsToContents()
            else:
                self.Tabla_denominator.setDisabled(True)
            tab.addTab(
                self.Tabla_denominator,
                QtWidgets.QApplication.translate("pychemqt", "Denominator"))
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 10, 3)

        elif eq["eq"] == 2:
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            self.ek = Entrada_con_unidades(
                float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(self.ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            self.sigma = Entrada_con_unidades(
                float, value=eq["sigma"], readOnly=True)
            gridLayout.addWidget(self.sigma, 5, 2)
            self.Tabla_Visco2 = Tabla(
                3, horizontalHeader=["b", "F", "E"], stretch=False,
                readOnly=True)
            if "collision" in eq:
                self.Tabla_Visco2.setColumn(0, eq["collision"])
            self.Tabla_Visco2.setColumn(1, eq["F"])
            self.Tabla_Visco2.setColumn(2, eq["E"])
            self.Tabla_Visco2.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco2, 6, 1, 1, 3)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 10, 3)

        elif eq["eq"] == 3:
            self.Tabla_Visco3 = Tabla(
                8, stretch=False, readOnly=True,
                horizontalHeader=["n-poly", "t-poly", "n-num", "t-num",
                                  "d-num", "n-den", "t-den", "d-den"])
            if "n_poly" in eq:
                self.Tabla_Visco3.setColumn(0, eq["n_poly"])
                self.Tabla_Visco3.setColumn(1, eq["t_poly"])
            if "n_num" in eq:
                self.Tabla_Visco3.setColumn(2, eq["n_num"])
                self.Tabla_Visco3.setColumn(3, eq["t_num"])
                self.Tabla_Visco3.setColumn(4, eq["d_num"])
            if "n_den" in eq:
                self.Tabla_Visco3.setColumn(5, eq["n_den"])
                self.Tabla_Visco3.setColumn(6, eq["t_den"])
                self.Tabla_Visco3.setColumn(7, eq["d_den"])
            self.Tabla_Visco3.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco3, 4, 1, 1, 3)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 10, 3)

        elif eq["eq"] == 4:
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            self.ek = Entrada_con_unidades(
                float, value=eq.get("ek", None), readOnly=True)
            gridLayout.addWidget(self.ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            self.sigma = Entrada_con_unidades(
                float, value=eq.get("sigma", None), readOnly=True)
            gridLayout.addWidget(self.sigma, 5, 2)
            self.Tabla_Visco4 = Tabla(
                7, stretch=False, readOnly=True,
                horizontalHeader=["a", "b", "c", "A", "B", "C", "D"])
            format = {"format": 1, "decimales": 10}
            self.Tabla_Visco4.setColumn(0, eq["a"], **format)
            self.Tabla_Visco4.setColumn(1, eq["b"], **format)
            self.Tabla_Visco4.setColumn(2, eq["c"], **format)
            self.Tabla_Visco4.setColumn(3, eq["A"], **format)
            self.Tabla_Visco4.setColumn(4, eq["B"], **format)
            self.Tabla_Visco4.setColumn(5, eq["C"], **format)
            # self.Tabla_Visco4.setColumn(6, eq["D"], **format)
            self.Tabla_Visco4.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Visco4, 6, 1, 1, 3)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 10, 3)

        elif eq["eq"] == 5:
            gridLayout.addWidget(QtWidgets.QLabel("w"), 4, 1)
            self.w = Entrada_con_unidades(float, value=eq["w"], readOnly=True)
            gridLayout.addWidget(self.w, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("mur"), 5, 1)
            self.mur = Entrada_con_unidades(
                float, value=eq["mur"], readOnly=True)
            gridLayout.addWidget(self.mur, 5, 2)
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 6, 1)
            self.k = Entrada_con_unidades(float, value=eq["k"], readOnly=True)
            gridLayout.addWidget(self.k, 6, 2)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 10, 3)


class Widget_Conductivity_Data(QtWidgets.QWidget):
    """Widget to show thermal conductivity data"""
    def __init__(self, element, eq, parent=None):
        """
        element: element class for code extract
        eq: dict with thermal conductivity parameter"""
        super(Widget_Conductivity_Data, self).__init__(parent)
        gridLayout = QtWidgets.QGridLayout(self)
        if eq["eq"] == 0:
            txt = element.__getattribute__(element, eq["method"]).__doc__
        else:
            txt = " ".join((eq["__doi__"]["autor"], eq["__doi__"]["title"],
                            eq["__doi__"]["ref"]))
        ref = QtWidgets.QLabel(txt)
        ref.setWordWrap(True)
        gridLayout.addWidget(ref, 1, 1, 1, 3)

        if eq["eq"] == 0:
            # Hardcoded method, show code
            self.code = SimplePythonEditor()
            code = ""
            for method in eq.get("__code__", ()):
                code += inspect.getsource(method)
                code += os.linesep
            code += inspect.getsource(
                element.__getattribute__(element, eq["method"]))
            self.code.setText(code)
            gridLayout.addWidget(self.code, 2, 1, 1, 3)

        elif eq["eq"] == 1:
            self.Tabla_Therm1 = Tabla(
                11, stretch=False, readOnly=True,
                horizontalHeader=["no", "co", "noden", "toden", "nb", "tb",
                                  "db", "cb", "nbden", "tbden", "dbden"])
            if "no" in eq:
                self.Tabla_Therm1.setColumn(0, eq["no"])
                self.Tabla_Therm1.setColumn(1, eq["co"])
            if "noden" in eq:
                self.Tabla_Therm1.setColumn(2, eq["noden"])
                self.Tabla_Therm1.setColumn(3, eq["toden"])
            if "nb" in eq:
                self.Tabla_Therm1.setColumn(4, eq["nb"])
                self.Tabla_Therm1.setColumn(5, eq["tb"])
                self.Tabla_Therm1.setColumn(6, eq["db"])
                self.Tabla_Therm1.setColumn(7, eq["cb"])
            if "nbden" in eq:
                self.Tabla_Therm1.setColumn(8, eq["nbden"])
                self.Tabla_Therm1.setColumn(9, eq["tbden"])
                self.Tabla_Therm1.setColumn(10, eq["dbden"])
            self.Tabla_Therm1.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm1, 3, 1, 1, 3)

        elif eq["eq"] == 2:
            self.Tabla_Therm2 = Tabla(
                2, horizontalHeader=["E", "G"], stretch=False, readOnly=True)
            self.Tabla_Therm2.setColumn(0, eq["E"])
            self.Tabla_Therm2.setColumn(1, eq["G"])
            self.Tabla_Therm2.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm2, 3, 1, 1, 3)

        elif eq["eq"] == 3:
            self.Tabla_Therm3 = Tabla(
                3, horizontalHeader=["b", "F", "E"],
                stretch=False, readOnly=True)
            self.Tabla_Therm3.setColumn(0, eq["b"])
            self.Tabla_Therm3.setColumn(1, eq["F"])
            self.Tabla_Therm3.setColumn(2, eq["E"])
            self.Tabla_Therm3.resizeColumnsToContents()
            gridLayout.addWidget(self.Tabla_Therm3, 3, 1, 1, 3)

            parameter = QtWidgets.QWidget()
            gridLayout.addWidget(parameter, 4, 1, 1, 3)
            lyt = QtWidgets.QGridLayout(parameter)
            lyt.addWidget(QtWidgets.QLabel("ε/k"), 1, 1)
            self.ek = Entrada_con_unidades(
                float, value=eq["ek"], readOnly=True)
            lyt.addWidget(self.ek, 1, 2)
            lyt.addWidget(QtWidgets.QLabel("σ"), 2, 1)
            self.sigma = Entrada_con_unidades(
                float, value=eq["sigma"], readOnly=True)
            lyt.addWidget(self.sigma, 2, 2)
            lyt.addWidget(QtWidgets.QLabel("Nchapman"), 3, 1)
            self.Nchapman = Entrada_con_unidades(
                float, value=eq["Nchapman"], readOnly=True)
            lyt.addWidget(self.Nchapman, 3, 2)
            lyt.addWidget(QtWidgets.QLabel("tchapman"), 4, 1)
            self.tchapman = Entrada_con_unidades(
                float, value=eq["tchapman"], readOnly=True)
            lyt.addWidget(self.tchapman, 4, 2)
            lyt.addWidget(QtWidgets.QLabel("rhoc"), 1, 4)
            self.rhoc = Entrada_con_unidades(
                float, value=eq["rhoc"], readOnly=True)
            lyt.addWidget(self.rhoc, 1, 5)
            lyt.addWidget(QtWidgets.QLabel("ff"), 2, 4)
            self.ff = Entrada_con_unidades(
                float, value=eq["ff"], readOnly=True)
            lyt.addWidget(self.ff, 2, 5)
            lyt.addWidget(QtWidgets.QLabel("rm"), 3, 4)
            self.rm = Entrada_con_unidades(
                float, value=eq["rm"], readOnly=True)
            lyt.addWidget(self.rm, 3, 5)

        if "critical" in eq and eq["critical"]:
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Fixed,
                QtWidgets.QSizePolicy.Fixed), 5, 3)
            gridLayout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Critical enhancement")), 6, 1, 1, 2)
            if eq["critical"] == 3:
                gridLayout.addWidget(QtWidgets.QLabel("gnu"), 7, 1)
                self.gnu = Entrada_con_unidades(
                    float, value=eq["gnu"], readOnly=True)
                gridLayout.addWidget(self.gnu, 7, 2)
                gridLayout.addWidget(QtWidgets.QLabel("γ"), 8, 1)
                self.gamma = Entrada_con_unidades(
                    float, value=eq["gamma"], readOnly=True)
                gridLayout.addWidget(self.gamma, 8, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Ro"), 9, 1)
                self.R0 = Entrada_con_unidades(
                    float, value=eq["R0"], readOnly=True)
                gridLayout.addWidget(self.R0, 9, 2)
                gridLayout.addWidget(QtWidgets.QLabel("ξo"), 10, 1)
                self.Xio = Entrada_con_unidades(
                    float, value=eq["Xio"], readOnly=True)
                gridLayout.addWidget(self.Xio, 10, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Γo"), 11, 1)
                self.gam0 = Entrada_con_unidades(
                    float, value=eq["gam0"], readOnly=True)
                gridLayout.addWidget(self.gam0, 11, 2)
                gridLayout.addWidget(QtWidgets.QLabel("qd"), 12, 1)
                self.qd = Entrada_con_unidades(
                    float, value=eq["qd"], readOnly=True)
                gridLayout.addWidget(self.qd, 12, 2)
            elif eq["critical"] == 4:
                gridLayout.addWidget(QtWidgets.QLabel("γ"), 7, 1)
                self.gamma = Entrada_con_unidades(
                    float, value=eq["gamma"], readOnly=True)
                gridLayout.addWidget(self.gamma, 7, 2)
                gridLayout.addWidget(QtWidgets.QLabel("v"), 8, 1)
                self.v = Entrada_con_unidades(
                    float, value=eq["expo"], readOnly=True)
                gridLayout.addWidget(self.v, 8, 2)
                gridLayout.addWidget(QtWidgets.QLabel("α"), 9, 1)
                self.alfa = Entrada_con_unidades(
                    float, value=eq["alfa"], readOnly=True)
                gridLayout.addWidget(self.alfa, 9, 2)
                gridLayout.addWidget(QtWidgets.QLabel("β"), 10, 1)
                self.beta = Entrada_con_unidades(
                    float, value=eq["beta"], readOnly=True)
                gridLayout.addWidget(self.beta, 10, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Γo"), 11, 1)
                self.Xio = Entrada_con_unidades(
                    float, value=eq["Xio"], readOnly=True)
                gridLayout.addWidget(self.Xio, 11, 2)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding), 15, 3)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    SteamTables = Ui_ChooseFluid()
    # SteamTables = Dialog_InfoFluid(mEoS.He)
    # SteamTables = transportDialog(mEoS.__all__[2])

    SteamTables.show()
    sys.exit(app.exec_())
