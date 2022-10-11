#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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

from qt import QtCore, QtGui, QtWidgets

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
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel,
            QtCore.Qt.Orientation.Vertical)
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
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
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
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Maximum), 8, 2)

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
            # Add the hidden element above the selected one
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
        if dlg.exec():
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
        dialog.exec()

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
            QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
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
            30, 30, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 1, 3, 3, 1)

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
            30, 30, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 1, 6, 3, 1)

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
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed),
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
        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
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
        txt = "; ".join((eq["__doi__"]["autor"], eq["__doi__"]["title"],
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
            mathTex += r"\alpha_{GBS}^r+\alpha_{NA}^r$"
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

        elif eq["__type__"] == "MBWR":
            # Pestaña MBWR
            mathTex = r"$P=\rho RT+\sum_{n=2}^{9}\alpha_n\rho^n + "
            mathTex += r"e^{-\delta^2} \sum_{10}^{15} \alpha_n"
            mathTex += r"\rho^{2n-17}$"
            label = QLabelMath(mathTex)
            gridLayout.addWidget(label, 2, 1)

            tab2 = QtWidgets.QWidget()
            tabWidget.addTab(
                tab2, QtWidgets.QApplication.translate("pychemqt", "MBWR"))
            gridLayout_MBWR = QtWidgets.QGridLayout(tab2)
            self.Tabla_MBWR = Tabla(
                1, horizontalHeader=["b"], stretch=False, readOnly=True)
            gridLayout_MBWR.addWidget(self.Tabla_MBWR, 2, 1)

        self.fill(eq)

    def fill(self, eq):
        format = {"fmt": 1, "total": 5}

        if "ao_log" in eq["cp"]:
            # Phi_o term
            self.Tabla_Cp_poly.setColumn(
                0, eq["cp"]["ao_pow"], **format)
            self.Tabla_Cp_poly.setColumn(1, eq["cp"]["pow"], **format)
            self.Tabla_Cp_poly.resizeColumnsToContents()
            if "ao_exp" in eq["cp"]:
                self.Tabla_Cp_exp.setColumn(0, eq["cp"]["ao_exp"], **format)
                self.Tabla_Cp_exp.setColumn(1, eq["cp"]["titao"], **format)
            else:
                self.Tabla_Cp_exp.setEnabled(False)
            self.Tabla_Cp_exp.resizeColumnsToContents()

            if "hyp" in eq["cp"]:
                self.Tabla_Cp_hyp.setColumn(0, eq["cp"]["ao_hyp"], **format)
                self.Tabla_Cp_hyp.setColumn(1, eq["cp"]["hyp"], **format)
            else:
                self.Tabla_Cp_hyp.setEnabled(False)
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
            else:
                self.Tabla_Cp_poly.setEnabled(False)
            self.Tabla_Cp_poly.resizeColumnsToContents()

            if "ao_exp" in eq["cp"]:
                self.Tabla_Cp_exp.setColumn(0, eq["cp"]["ao_exp"], **format)
                self.Tabla_Cp_exp.setColumn(1, eq["cp"]["exp"], **format)
            else:
                self.Tabla_Cp_exp.setEnabled(False)
            self.Tabla_Cp_exp.resizeColumnsToContents()

            if "hyp" in eq["cp"]:
                self.Tabla_Cp_hyp.setColumn(0, eq["cp"]["ao_hyp"], **format)
                self.Tabla_Cp_hyp.setColumn(1, eq["cp"]["hyp"], **format)
            else:
                self.Tabla_Cp_hyp.setEnabled(False)
            self.Tabla_Cp_hyp.resizeColumnsToContents()

        if eq["__type__"] == "Helmholtz":
            if eq.get("nr1", []):
                self.Tabla_lineal.setColumn(0, eq["nr1"], **format)
                self.Tabla_lineal.setColumn(1, eq["t1"], **format)
                self.Tabla_lineal.setColumn(2, eq["d1"], **format)
            else:
                self.Tabla_lineal.setEnabled(False)

            if eq.get("nr2", []):
                self.Tabla_exponential.setColumn(0, eq["nr2"], **format)
                self.Tabla_exponential.setColumn(1, eq["t2"], **format)
                self.Tabla_exponential.setColumn(2, eq["d2"], **format)
                self.Tabla_exponential.setColumn(3, eq["gamma2"], **format)
                self.Tabla_exponential.setColumn(4, eq["c2"], **format)
            else:
                self.Tabla_exponential.setEnabled(False)

            if eq.get("nr3", []):
                self.Tabla_gauss.setColumn(0, eq["nr3"], **format)
                self.Tabla_gauss.setColumn(1, eq["t3"], **format)
                self.Tabla_gauss.setColumn(2, eq["d3"], **format)
                self.Tabla_gauss.setColumn(3, eq["alfa3"], **format)
                self.Tabla_gauss.setColumn(4, eq["beta3"], **format)
                self.Tabla_gauss.setColumn(5, eq["gamma3"], **format)
                self.Tabla_gauss.setColumn(6, eq["epsilon3"], **format)
            else:
                self.Tabla_gauss.setEnabled(False)

            if eq.get("nr4", []):
                self.Tabla_noanalytic.setColumn(0, eq["nr4"], **format)
                self.Tabla_noanalytic.setColumn(1, eq["a4"], **format)
                self.Tabla_noanalytic.setColumn(2, eq["b4"], **format)
                self.Tabla_noanalytic.setColumn(3, eq["A"], **format)
                self.Tabla_noanalytic.setColumn(4, eq["B"], **format)
                self.Tabla_noanalytic.setColumn(5, eq["C"], **format)
                self.Tabla_noanalytic.setColumn(6, eq["D"], **format)
                self.Tabla_noanalytic.setColumn(7, eq["beta4"], **format)
            else:
                self.Tabla_noanalytic.setEnabled(False)

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

        tabWidget = QtWidgets.QTabWidget()
        gridLayout.addWidget(tabWidget, 1, 1)

        # Tab viscosity
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab3, QtWidgets.QApplication.translate("pychemqt", "Viscosity"))
        lyt_viscosity = QtWidgets.QGridLayout(tab3)

        lyt_viscosity.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate(
                "pychemqt", "Equation")+": "), 1, 1)
        eqVisco = QtWidgets.QComboBox()
        lyt_viscosity.addWidget(eqVisco, 1, 2)
        lyt_viscosity.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)
        stackedVisco = QtWidgets.QStackedWidget()
        lyt_viscosity.addWidget(stackedVisco, 2, 1, 1, 3)
        eqVisco.currentIndexChanged.connect(
            stackedVisco.setCurrentIndex)

        if element._viscosity is not None:
            for eq in element._viscosity:
                widget = Widget_Viscosity_Data(element, eq)
                stackedVisco.addWidget(widget)
                eqVisco.addItem(eq["__name__"])
        else:
            eqVisco.addItem(QtWidgets.QApplication.translate(
                "pychemqt", "Not Implemented"))

        # Tab thermal conductivity
        tab4 = QtWidgets.QWidget()
        tabWidget.addTab(tab4,
                         QtWidgets.QApplication.translate(
                             "pychemqt", "Thermal Conductivity"))
        lyt_conductivity = QtWidgets.QGridLayout(tab4)

        lyt_conductivity.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Equation")+": "),
            1, 1)
        eqThermo = QtWidgets.QComboBox()
        lyt_conductivity.addWidget(eqThermo, 1, 2)
        lyt_conductivity.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 3)
        stackedThermo = QtWidgets.QStackedWidget()
        lyt_conductivity.addWidget(stackedThermo, 2, 1, 1, 3)
        eqThermo.currentIndexChanged.connect(
            stackedThermo.setCurrentIndex)

        if element._thermal is not None:
            for eq in element._thermal:
                widget = Widget_Conductivity_Data(element, eq)
                stackedThermo.addWidget(widget)
                eqThermo.addItem(eq["__name__"])
        else:
            eqThermo.addItem(QtWidgets.QApplication.translate(
                "pychemqt", "Not Implemented"))

        # Tab dielectric constant
        tab1 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab1, QtWidgets.QApplication.translate("pychemqt", "Dielectric"))
        lyt_dielectric = QtWidgets.QGridLayout(tab1)

        if element._Dielectric != meos.MEoS._Dielectric:
            label = QtWidgets.QLabel(element._Dielectric.__doc__)
            label.setWordWrap(True)
            lyt_dielectric.addWidget(label, 1, 1)
            codigo_Dielectric = SimplePythonEditor()
            codigo_Dielectric.setText(
                inspect.getsource(element._Dielectric))
            lyt_dielectric.addWidget(codigo_Dielectric, 2, 1)
        elif element._dielectric:
            if element._dielectric["eq"] == 1:
                txt1 = r"$P = \frac{\epsilon-1}{\epsilon+2}$"
            else:
                txt1 = r"$P = \frac{(\epsilon-1)(2\epsilon+1}{9\epsilon}$"
            label1 = QLabelMath(txt1)
            lyt_dielectric.addWidget(label1, 1, 1)
            txt2 = r"$\frac{P}{\rho} = A_{\epsilon} + \frac{A_{\mu}}{T} + "
            txt2 += r"B_{\epsilon}\rho + C\rho^D$"
            label2 = QLabelMath(txt2)
            lyt_dielectric.addWidget(label2, 2, 1)

            Table_Dielectric = Tabla(
                1, verticalHeader=True, filas=5, stretch=False, readOnly=True)
            lyt_dielectric.addWidget(Table_Dielectric, 3, 1)
            i = 0
            for key, valor in element._dielectric.items():
                if key == "eq":
                    continue
                Table_Dielectric.setVerticalHeaderItem(
                    i, QtWidgets.QTableWidgetItem(key))
                Table_Dielectric.setItem(
                    0, i, QtWidgets.QTableWidgetItem(str(valor)))
                i += 1
            Table_Dielectric.resizeColumnsToContents()
        else:
            lyt_dielectric.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_dielectric.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab surface tension
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab2,
            QtWidgets.QApplication.translate("pychemqt", "Surface Tension"))
        lyt_surface = QtWidgets.QGridLayout(tab2)

        if element._Surface != meos.MEoS._Surface:
            label = QtWidgets.QLabel(element._Surface.__doc__)
            label.setWordWrap(True)
            lyt_surface.addWidget(label, 1, 1)
            code_Surface = SimplePythonEditor()
            code_Surface.setText(inspect.getsource(element._Surface))
            lyt_surface.addWidget(code_Surface, 2, 1)
        elif element._surface:
            txt = r"$\sigma=\sum_i \sigma_i\left(1-\frac{T}{T_c}\right)^{n_i}$"
            label = QLabelMath(txt)
            lyt_surface.addWidget(label, 1, 1)

            Table_Surface = Tabla(
                2, horizontalHeader=["σ", "n"], verticalHeader=True,
                stretch=False, readOnly=True)
            Table_Surface.setColumn(0, element._surface["sigma"])
            Table_Surface.setColumn(1, element._surface["exp"])
            lyt_surface.addWidget(Table_Surface, 2, 1)
            Table_Surface.resizeColumnsToContents()
        else:
            lyt_surface.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_surface.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab liquid density
        tab5 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab5,
            QtWidgets.QApplication.translate("pychemqt", "Liquid Density"))
        lyt_liquid_density = QtWidgets.QGridLayout(tab5)

        if element._Liquid_Density != meos.MEoS._Liquid_Density:
            label = QtWidgets.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            lyt_liquid_density.addWidget(label, 1, 1)
            code_Liquid_Density = SimplePythonEditor()
            code_Liquid_Density.setText(
                inspect.getsource(element._Liquid_Density))
            lyt_liquid_density.addWidget(
                code_Liquid_Density, 2, 1)
        elif element._liquid_Density:
            label = QtWidgets.QLabel(element._Liquid_Density.__doc__)
            label.setWordWrap(True)
            lyt_liquid_density.addWidget(label, 1, 1)

            if element._liquid_Density["eq"] in [2, 3]:
                mathTex = r"$\ln \frac{\rho_l}{\rho_c} = "
            else:
                mathTex = r"$\frac{\rho_l}{\rho_c} = "

            if element._liquid_Density["eq"] == 1:
                mathTex += r"1 + "
            elif element._liquid_Density["eq"] == 3:
                mathTex += r"\frac{T_c}{T} "

            mathTex += r"\sum_i n_i\theta^{t_i}$"
            label2 = QLabelMath(mathTex)
            lyt_liquid_density.addWidget(label2, 2, 1)

            Table_Liquid_Density = Tabla(
                2, horizontalHeader=["n", "t"],
                verticalHeader=True, stretch=False, readOnly=True)
            Table_Liquid_Density.setColumn(0, element._liquid_Density["n"])
            Table_Liquid_Density.setColumn(1, element._liquid_Density["t"])
            lyt_liquid_density.addWidget(Table_Liquid_Density, 3, 1)
            Table_Liquid_Density.resizeColumnsToContents()
        else:
            lyt_liquid_density.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_liquid_density.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab vapor density
        tab6 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab6,
            QtWidgets.QApplication.translate("pychemqt", "Vapor Density"))
        lyt_vapor_density = QtWidgets.QGridLayout(tab6)

        if element._Vapor_Density != meos.MEoS._Vapor_Density:
            label = QtWidgets.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            lyt_vapor_density.addWidget(label, 1, 1)
            code_Vapor_Density = SimplePythonEditor()
            code_Vapor_Density.setText(
                inspect.getsource(element._Vapor_Density))
            lyt_vapor_density.addWidget(code_Vapor_Density, 2, 1)
        elif element._vapor_Density:
            label = QtWidgets.QLabel(element._Vapor_Density.__doc__)
            label.setWordWrap(True)
            lyt_vapor_density.addWidget(label, 1, 1)

            if element._vapor_Density["eq"] in [2, 3]:
                mathTex = r"$\ln \frac{\rho_v}{\rho_c} = "
            else:
                mathTex = r"$\frac{\rho_v}{\rho_c} = "

            if element._vapor_Density["eq"] == 1:
                mathTex += r"1 + "
            elif element._vapor_Density["eq"] == 3:
                mathTex += r"\frac{T_c}{T} "

            mathTex += r"\sum_i n_i\theta^{t_i}$"
            label2 = QLabelMath(mathTex)
            lyt_vapor_density.addWidget(label2, 2, 1)

            Table_Vapor_Density = Tabla(
                2, horizontalHeader=["n", "t"], verticalHeader=True,
                stretch=False, readOnly=True)
            Table_Vapor_Density.setColumn(0, element._vapor_Density["n"])
            Table_Vapor_Density.setColumn(1, element._vapor_Density["t"])
            lyt_vapor_density.addWidget(Table_Vapor_Density, 3, 1)
            Table_Vapor_Density.resizeColumnsToContents()
        else:
            lyt_vapor_density.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_vapor_density.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab vapor presure
        tab7 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab7,
            QtWidgets.QApplication.translate("pychemqt", "Vapor Pressure"))
        lyt_vapor_pressure = QtWidgets.QGridLayout(tab7)

        if element._Vapor_Pressure != meos.MEoS._Vapor_Pressure:
            label = QtWidgets.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            lyt_vapor_pressure.addWidget(label, 1, 1)
            code_Vapor_Pressure = SimplePythonEditor()
            code_Vapor_Pressure.setText(
                inspect.getsource(element._Vapor_Pressure))
            lyt_vapor_pressure.addWidget(
                code_Vapor_Pressure, 2, 1)
        elif element._vapor_Pressure:
            label = QtWidgets.QLabel(element._Vapor_Pressure.__doc__)
            label.setWordWrap(True)
            lyt_vapor_pressure.addWidget(label, 1, 1)

            if element._vapor_Pressure["eq"] in [2, 3]:
                mathTex = r"$\ln \frac{P_v}{P_c} = "
            else:
                mathTex = r"$\frac{P_v}{P_c} = "

            if element._vapor_Pressure["eq"] == 1:
                mathTex += r"1 + "
            elif element._vapor_Pressure["eq"] == 3:
                mathTex += r"\frac{T_c}{T} "

            mathTex += r"\sum_i n_i\theta^{t_i}$"
            label2 = QLabelMath(mathTex)
            lyt_vapor_pressure.addWidget(label2, 2, 1)

            Table_Vapor_Pressure = Tabla(
                2, horizontalHeader=["n", "t"], verticalHeader=True,
                stretch=False, readOnly=True)
            Table_Vapor_Pressure.setColumn(0, element._vapor_Pressure["n"])
            Table_Vapor_Pressure.setColumn(1, element._vapor_Pressure["t"])
            lyt_vapor_pressure.addWidget(Table_Vapor_Pressure, 3, 1)
            Table_Vapor_Pressure.resizeColumnsToContents()
        else:
            lyt_vapor_pressure.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_vapor_pressure.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab melting presure
        tab8 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab8,
            QtWidgets.QApplication.translate("pychemqt", "Melting Pressure"))
        lyt_melting = QtWidgets.QGridLayout(tab8)

        if inspect.getsource(element._Melting_Pressure) != \
                inspect.getsource(meos.MEoS._Melting_Pressure):
            label = QtWidgets.QLabel(element._Melting_Pressure.__doc__)
            label.setWordWrap(True)
            lyt_melting.addWidget(label, 1, 1)
            code_Melting = SimplePythonEditor()
            code_Melting.setText(
                inspect.getsource(element._Melting_Pressure))
            lyt_melting.addWidget(
                code_Melting, 2, 1)
        elif element._melting:
            if "__doi__" in element._melting:
                eq = element._melting["__doi__"]
                ref = "; ".join((eq["autor"], eq["title"], eq["ref"]))
                ref += os.linesep
                label2 = QtWidgets.QLabel(ref)
                label2.setWordWrap(True)
                lyt_melting.addWidget(label2, 1, 1)

            if element._melting["eq"] == 1:
                mathTex = r"$\frac{P_m}{P_r} = "
            elif element._melting["eq"] == 2:
                mathTex = r"$P_m - P_r = "
            elif element._melting["eq"] == 3:
                mathTex = r"$\ln \frac{P_m}{P_r} = "

            mathTex += r"a_o + \sum a_1\theta^{n_1} + "
            mathTex += r"\sum a_2\left(\theta^{n_2}-1\right) + "
            mathTex += r"\sum a_3\log \theta^{n_3} + "
            mathTex += r"\sum a_4\left(\theta-1\right)^{n_4}$"
            label = QLabelMath(mathTex)
            lyt_melting.addWidget(label, 2, 1)

            header = ["a1", "n1", "a2", "n2", "a3", "n3", "a4", "n4"]
            Table_Melting = Tabla(
                8, horizontalHeader=header, verticalHeader=True,
                stretch=False, readOnly=True)
            if "a1" in element._melting:
                Table_Melting.setColumn(0, element._melting["a1"])
                Table_Melting.setColumn(1, element._melting["exp1"])
            if "a2" in element._melting:
                Table_Melting.setColumn(2, element._melting["a2"])
                Table_Melting.setColumn(3, element._melting["exp2"])
            if "a3" in element._melting:
                Table_Melting.setColumn(4, element._melting["a3"])
                Table_Melting.setColumn(5, element._melting["exp3"])
            if "a4" in element._melting:
                Table_Melting.setColumn(6, element._melting["a4"])
                Table_Melting.setColumn(7, element._melting["exp4"])
            lyt_melting.addWidget(
                Table_Melting, 3, 1)
            Table_Melting.resizeColumnsToContents()
        else:
            lyt_melting.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_melting.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        # Tab sublimation presure
        tab9 = QtWidgets.QWidget()
        tabWidget.addTab(
            tab9,
            QtWidgets.QApplication.translate(
                "pychemqt", "Sublimation Pressure"))
        lyt_sublimation = QtWidgets.QGridLayout(tab9)

        if inspect.getsource(element._Sublimation_Pressure) != \
                inspect.getsource(meos.MEoS._Sublimation_Pressure):
            label = QtWidgets.QLabel(element._Sublimation_Pressure.__doc__)
            label.setWordWrap(True)
            lyt_sublimation.addWidget(label, 1, 1)
            code_Sublimation = SimplePythonEditor()
            code_Sublimation.setText(
                inspect.getsource(element._Sublimation_Pressure))
            lyt_sublimation.addWidget(
                code_Sublimation, 2, 1)
        elif element._sublimation:
            if "__doi__" in element._sublimation:
                eq = element._sublimation["__doi__"]
                ref = "; ".join((eq["autor"], eq["title"], eq["ref"]))
                ref += os.linesep
                label2 = QtWidgets.QLabel(ref)
                label2.setWordWrap(True)
                lyt_sublimation.addWidget(label2, 1, 1)

            if element._sublimation["eq"] == 1:
                mathTex = r"$\frac{P_s}{P_r} = "
            elif element._sublimation["eq"] == 2:
                mathTex = r"$P_s - P_r = "
            elif element._sublimation["eq"] == 3:
                mathTex = r"$\ln \frac{P_s}{P_r} = "

            mathTex += r"a_o + \sum a_1 \theta^{t_1} + "
            mathTex += r"\sum a_2\left(1-\theta\right)^{t_2} + "
            mathTex += r"\sum a_i\log \theta^{t_i}$"
            label = QLabelMath(mathTex)
            lyt_sublimation.addWidget(label, 2, 1)

            Table_Sublimation = Tabla(
                6, horizontalHeader=["a1", "n1", "a2", "n2", "a3", "n3"],
                verticalHeader=True, stretch=False, readOnly=True)
            if "a1" in element._sublimation:
                Table_Sublimation.setColumn(0, element._sublimation["a1"])
                Table_Sublimation.setColumn(1, element._sublimation["exp1"])
            if "a2" in element._sublimation:
                Table_Sublimation.setColumn(2, element._sublimation["a2"])
                Table_Sublimation.setColumn(3, element._sublimation["exp2"])
            if "a3" in element._sublimation:
                Table_Sublimation.setColumn(4, element._sublimation["a3"])
                Table_Sublimation.setColumn(5, element._sublimation["exp3"])
            lyt_sublimation.addWidget(Table_Sublimation, 3, 1)
            Table_Sublimation.resizeColumnsToContents()
        else:
            lyt_sublimation.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Not Implemented")), 1, 1)
            lyt_sublimation.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 2, 1)

        btBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        btBox.clicked.connect(self.reject)
        gridLayout.addWidget(btBox, 2, 1)


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
            code_Viscosity = SimplePythonEditor()
            code = ""
            for method in eq.get("__code__", ()):
                code += inspect.getsource(method)
                code += os.linesep
            code += inspect.getsource(
                element.__getattribute__(element, eq["method"]))
            code_Viscosity.setText(code)
            gridLayout.addWidget(code_Viscosity, 2, 1, 1, 3)

        elif eq["eq"] == 1:
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            ek = Entrada_con_unidades(
                float, value=eq.get("ek", 0), readOnly=True)
            gridLayout.addWidget(ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            sigma = Entrada_con_unidades(
                float, value=eq.get("sigma", 0), readOnly=True)
            gridLayout.addWidget(sigma, 5, 2)
            tab = QtWidgets.QTabWidget()
            gridLayout.addWidget(tab, 6, 1, 1, 3)

            # Integral collision
            Tabla_Collision = Tabla(
                1, horizontalHeader=["b"], stretch=False, readOnly=True)
            if "collision" in eq:
                Tabla_Collision.setColumn(0, eq["collision"])
            else:
                Tabla_Collision.setDisabled(True)
            Tabla_Collision.resizeColumnsToContents()
            tab.addTab(
                Tabla_Collision,
                QtWidgets.QApplication.translate("pychemqt", "Collision"))

            # Virial
            Tabla_Virial = Tabla(
                2, horizontalHeader=["n", "t"], stretch=False, readOnly=True)
            if "n_virial" in eq:
                Tabla_Virial.setColumn(0, eq["n_virial"])
                Tabla_Virial.setColumn(1, eq["t_virial"])
            else:
                Tabla_Virial.setDisabled(True)
            Tabla_Virial.resizeColumnsToContents()
            tab.addTab(Tabla_Virial,
                       QtWidgets.QApplication.translate("pychemqt", "Virial"))

            # Close-packed
            Tabla_Packed = Tabla(4, horizontalHeader=["f", "g1", "g", "t"],
                                 stretch=False, readOnly=True)
            if "CPf" in eq:
                Tabla_Packed.setColumn(0, [eq["CPf"]])
                Tabla_Packed.setColumn(1, [eq["CPg1"]])
                Tabla_Packed.setColumn(2, eq["CPgi"])
                Tabla_Packed.setColumn(3, eq["CPgi"])
            else:
                Tabla_Packed.setDisabled(True)
            Tabla_Packed.resizeColumnsToContents()
            tab.addTab(Tabla_Packed, QtWidgets.QApplication.translate(
                "pychemqt", "Close-packed density"))

            # polynomial term
            Tabla_Visco1 = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "nr" in eq:
                Tabla_Visco1.setColumn(0, eq["nr"])
                Tabla_Visco1.setColumn(1, eq["tr"])
                Tabla_Visco1.setColumn(2, eq["dr"])
                if "gr" in eq:
                    Tabla_Visco1.setColumn(3, eq["gr"])
                    Tabla_Visco1.setColumn(4, eq["cr"])
            else:
                Tabla_Visco1.setDisabled(True)
            Tabla_Visco1.resizeColumnsToContents()
            tab.addTab(
                Tabla_Visco1,
                QtWidgets.QApplication.translate("pychemqt", "Polinomial"))

            # numerator of rational poly
            Tabla_numerator = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "nr_num" in eq:
                Tabla_numerator.setColumn(0, eq["nr_num"])
                Tabla_numerator.setColumn(1, eq["tr_num"])
                Tabla_numerator.setColumn(2, eq["dr_num"])
                if "gr_num" in eq:
                    Tabla_numerator.setColumn(3, eq["gr_num"])
                    Tabla_numerator.setColumn(4, eq["cr_num"])
            else:
                Tabla_numerator.setDisabled(True)
            Tabla_numerator.resizeColumnsToContents()
            tab.addTab(
                Tabla_numerator,
                QtWidgets.QApplication.translate("pychemqt", "Numerator"))

            # denominator of rational poly
            Tabla_denominator = Tabla(
                5, horizontalHeader=["n", "t", "d", "g", "c"], stretch=False,
                readOnly=True)
            if "nr_den" in eq:
                Tabla_denominator.setColumn(0, eq["nr_den"])
                Tabla_denominator.setColumn(1, eq["tr_den"])
                Tabla_denominator.setColumn(2, eq["dr_den"])
                Tabla_denominator.setColumn(
                    3, eq.get("cr_den", [0]*len(eq["nr_den"])))
                Tabla_denominator.setColumn(
                    4, eq.get("gr_den", [0]*len(eq["nr_den"])))
            else:
                Tabla_denominator.setDisabled(True)
            Tabla_denominator.resizeColumnsToContents()
            tab.addTab(
                Tabla_denominator,
                QtWidgets.QApplication.translate("pychemqt", "Denominator"))
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)

        elif eq["eq"] == 2:
            Tabla_Visco2 = Tabla(
                3, horizontalHeader=["b", "F", "E"], stretch=False,
                readOnly=True)
            if "collision" in eq:
                Tabla_Visco2.setColumn(0, eq["collision"])
            Tabla_Visco2.setColumn(1, eq["F"])
            Tabla_Visco2.setColumn(2, eq["E"])
            Tabla_Visco2.resizeColumnsToContents()
            gridLayout.addWidget(Tabla_Visco2, 6, 1, 1, 3)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)

        elif eq["eq"] == 3:
            # ecs formulation
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            ek = Entrada_con_unidades(
                float, value=eq.get("ek", None), readOnly=True)
            gridLayout.addWidget(ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            sigma = Entrada_con_unidades(
                float, value=eq.get("sigma", None), readOnly=True)
            gridLayout.addWidget(sigma, 5, 2)

            gridLayout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate("pychemqt", "ECS formulation")
                + ": " + eq["ref"]), 5, 1, 1, 2)

        elif eq["eq"] == 4:
            # Quiñones-Cisneros formulation
            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            ek = Entrada_con_unidades(
                float, value=eq.get("ek", None), readOnly=True)
            gridLayout.addWidget(ek, 4, 2)
            gridLayout.addWidget(QtWidgets.QLabel("σ"), 5, 1)
            sigma = Entrada_con_unidades(
                float, value=eq.get("sigma", None), readOnly=True)
            gridLayout.addWidget(sigma, 5, 2)
            Tabla_Visco4 = Tabla(
                8, stretch=False, readOnly=True,
                horizontalHeader=["a", "b", "c", "A", "B", "C", "D", "E"])
            format = {"fmt": 1, "decimales": 10}
            Tabla_Visco4.setColumn(0, eq["a"], **format)
            Tabla_Visco4.setColumn(1, eq["b"], **format)
            Tabla_Visco4.setColumn(2, eq["c"], **format)
            Tabla_Visco4.setColumn(3, eq["A"], **format)
            Tabla_Visco4.setColumn(4, eq["B"], **format)
            Tabla_Visco4.setColumn(5, eq["C"], **format)
            if "D" in eq:
                Tabla_Visco4.setColumn(6, eq["D"], **format)
                Tabla_Visco4.setColumn(7, eq["E"], **format)
            Tabla_Visco4.resizeColumnsToContents()
            gridLayout.addWidget(Tabla_Visco4, 6, 1, 1, 3)
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 10, 3)


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
            code = SimplePythonEditor()
            txt = ""
            for method in eq.get("__code__", ()):
                txt += inspect.getsource(method)
                txt += os.linesep
            txt += inspect.getsource(
                element.__getattribute__(element, eq["method"]))
            code.setText(txt)
            gridLayout.addWidget(code, 2, 1, 1, 3)

        elif eq["eq"] == 1:
            gridLayout.addWidget(QtWidgets.QLabel("No_visco"), 3, 1)
            no = Entrada_con_unidades(
                float, value=eq.get("no_visco", 0), readOnly=True)
            gridLayout.addWidget(no, 3, 2)
            Tabla_Therm1 = Tabla(
                13, stretch=False, readOnly=True,
                horizontalHeader=["no_viscoCP", "to_viscoCP", "no", "to",
                                  "no_num", "to_num", "no_den", "to_den",
                                  "nr", "tr", "dr", "gr", "cr"])
            if "no_viscoCP" in eq:
                Tabla_Therm1.setColumn(0, eq["no_viscoCP"])
                Tabla_Therm1.setColumn(1, eq["to_viscoCP"])
            if "no" in eq:
                Tabla_Therm1.setColumn(2, eq["no"])
                Tabla_Therm1.setColumn(3, eq["to"])
            if "no_num" in eq:
                Tabla_Therm1.setColumn(4, eq["no_num"])
                Tabla_Therm1.setColumn(5, eq["to_num"])
                Tabla_Therm1.setColumn(6, eq["no_num"])
                Tabla_Therm1.setColumn(7, eq["to_num"])
            if "nr" in eq:
                Tabla_Therm1.setColumn(8, eq["nr"])
                Tabla_Therm1.setColumn(9, eq["tr"])
                Tabla_Therm1.setColumn(10, eq["dr"])
            if "gr" in eq:
                Tabla_Therm1.setColumn(11, eq["gr"])
                Tabla_Therm1.setColumn(12, eq["cr"])
            Tabla_Therm1.resizeColumnsToContents()
            gridLayout.addWidget(Tabla_Therm1, 4, 1, 1, 3)

            if "special" in eq:
                # Hardcoded method
                gridLayout.addWidget(QtWidgets.QLabel(
                    QtWidgets.QApplication.translate(
                    "pychemqt", "Special hardcoded term")), 5, 1, 1, 3)
                code_SpecialK = SimplePythonEditor()
                code = inspect.getsource(
                    element.__getattribute__(element, eq["special"]))
                code_SpecialK.setText(code)
                gridLayout.addWidget(code_SpecialK, 6, 1, 1, 3)

        elif eq["eq"] == 2:
            Tabla_Therm2 = Tabla(
                3, horizontalHeader=["G", "F", "E"], stretch=False,
                readOnly=True)
            Tabla_Therm2.setColumn(0, eq["no"])
            Tabla_Therm2.setColumn(1, eq["F"])
            Tabla_Therm2.setColumn(2, eq["E"])
            Tabla_Therm2.resizeColumnsToContents()
            gridLayout.addWidget(Tabla_Therm2, 3, 1, 1, 3)

        elif eq["eq"] == 3:
            Tabla_Therm3 = Tabla(
                2, horizontalHeader=["G", "E"], stretch=False, readOnly=True)
            Tabla_Therm3.setColumn(0, eq["G"])
            Tabla_Therm3.setColumn(1, eq["E"])
            Tabla_Therm3.resizeColumnsToContents()
            gridLayout.addWidget(Tabla_Therm3, 3, 1, 1, 3)

            gridLayout.addWidget(QtWidgets.QLabel("ε/k"), 4, 1)
            ek = Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
            gridLayout.addWidget(ek, 4, 2)

        if "critical" in eq and eq["critical"]:
            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
                QtWidgets.QSizePolicy.Policy.Fixed), 5, 3)
            gridLayout.addWidget(QtWidgets.QLabel(
                QtWidgets.QApplication.translate(
                    "pychemqt", "Critical enhancement")), 7, 1, 1, 3)

            if eq["critical"] == 1:
                gridLayout.addWidget(QtWidgets.QLabel("ek"), 8, 1)
                ek = Entrada_con_unidades(float, value=eq["ek"], readOnly=True)
                gridLayout.addWidget(ek, 8, 2)
                gridLayout.addWidget(QtWidgets.QLabel("f"), 9, 1)
                f = Entrada_con_unidades(float, value=eq["f"], readOnly=True)
                gridLayout.addWidget(f, 9, 2)
                gridLayout.addWidget(QtWidgets.QLabel("gm"), 10, 1)
                gm = Entrada_con_unidades(float, value=eq["gm"], readOnly=True)
                gridLayout.addWidget(gm, 10, 2)

            elif eq["critical"] == 2:
                v = eq["X"]
                gridLayout.addWidget(QtWidgets.QLabel("X₁"), 8, 1)
                X1 = Entrada_con_unidades(float, value=v[0], readOnly=True)
                gridLayout.addWidget(X1, 8, 2)
                gridLayout.addWidget(QtWidgets.QLabel("X₂"), 9, 1)
                X2 = Entrada_con_unidades(float, value=v[1], readOnly=True)
                gridLayout.addWidget(X2, 9, 2)
                gridLayout.addWidget(QtWidgets.QLabel("X₃"), 10, 1)
                X3 = Entrada_con_unidades(float, value=v[2], readOnly=True)
                gridLayout.addWidget(X3, 10, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Xㄙ"), 11, 1)
                X4 = Entrada_con_unidades(float, value=v[3], readOnly=True)
                gridLayout.addWidget(X4, 11, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Z"), 12, 1)
                Z = Entrada_con_unidades(float, value=eq["Z"], readOnly=True)
                gridLayout.addWidget(Z, 12, 2)

            elif eq["critical"] == 3:
                gridLayout.addWidget(QtWidgets.QLabel("gnu"), 8, 1)
                gnu = Entrada_con_unidades(
                    float, value=eq["gnu"], readOnly=True)
                gridLayout.addWidget(gnu, 8, 2)
                gridLayout.addWidget(QtWidgets.QLabel("γ"), 9, 1)
                gamma = Entrada_con_unidades(
                    float, value=eq["gamma"], readOnly=True)
                gridLayout.addWidget(gamma, 9, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Ro"), 10, 1)
                R0 = Entrada_con_unidades(
                    float, value=eq["R0"], readOnly=True)
                gridLayout.addWidget(R0, 10, 2)
                gridLayout.addWidget(QtWidgets.QLabel("ξo"), 11, 1)
                Xio = Entrada_con_unidades(
                    float, value=eq["Xio"], readOnly=True)
                gridLayout.addWidget(Xio, 11, 2)
                gridLayout.addWidget(QtWidgets.QLabel("Γo"), 12, 1)
                gam0 = Entrada_con_unidades(
                    float, value=eq["gam0"], readOnly=True)
                gridLayout.addWidget(gam0, 12, 2)
                gridLayout.addWidget(QtWidgets.QLabel("qd"), 13, 1)
                qd = Entrada_con_unidades(
                    float, value=eq["qd"], readOnly=True)
                gridLayout.addWidget(qd, 13, 2)

            elif eq["critical"] == 4:
                Tabla_Thermk4 = Tabla(
                    5, horizontalHeader=["nc", "αc", "tc", "βc", "dc"],
                    stretch=False, readOnly=True)
                Tabla_Thermk4.setColumn(0, eq["nc"])
                Tabla_Thermk4.setColumn(1, eq["alfac"])
                Tabla_Thermk4.setColumn(2, eq["tc"])
                Tabla_Thermk4.setColumn(3, eq["betac"])
                Tabla_Thermk4.setColumn(4, eq["dc"])
                Tabla_Thermk4.resizeColumnsToContents()
                gridLayout.addWidget(Tabla_Thermk4, 8, 1, 1, 3)

            elif isinstance(eq["critical"], str):
                # Hardcoded method
                code_ThermK = SimplePythonEditor()
                code = inspect.getsource(
                    element.__getattribute__(element, eq["critical"]))
                code_ThermK.setText(code)
                gridLayout.addWidget(code_ThermK, 8, 1, 1, 3)

            gridLayout.addItem(QtWidgets.QSpacerItem(
                10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Expanding), 15, 3)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    # SteamTables = Ui_ChooseFluid()
    # SteamTables = Dialog_InfoFluid(mEoS.He)

    # for i, cmp in enumerate(mEoS.__all__):
        # print(i, cmp.__module__)
    SteamTables = transportDialog(mEoS.__all__[16])

    SteamTables.show()
    sys.exit(app.exec())
