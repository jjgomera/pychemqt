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
# Basic periodic table with properties dialog
#   - qtelemental: Periodic table
#   - boton: Element button in periodic table
#   - ElementDialog: Dialog to show properties of components
###############################################################################


from configparser import ConfigParser
import logging
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from UI.widgets import Entrada_con_unidades, Tabla
from lib import unidades
from lib.elemental import Elemental, _configValues, cleanFloat
from lib.config import conf_dir


font7 = QtGui.QFont()
font7.setPointSize(7.6)
font11 = QtGui.QFont()
font11.setPointSize(11)
font20 = QtGui.QFont()
font20.setPointSize(20)
font_title = QtGui.QFont()
font_title.setWeight(75)
font_title.setPointSize(11)
font_title.setBold(True)
palette = QtGui.QPalette()
palette.setColor(QtGui.QPalette.Window, QtGui.QColor("#c8c8c8"))
alignment = QtCore.Qt.AlignRight | \
    QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter


class qtelemental(QtWidgets.QDialog):
    """Periodic table graph"""
    def __init__(self, parent=None):
        super(qtelemental, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/PeriodicTableIcon.png")))
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Periodic Table"))

        self.Preferences = ConfigParser()
        self.Preferences.read(conf_dir+"pychemqtrc")

        layout = QtWidgets.QGridLayout(self)
        layout.setSpacing(2)

        self.populate()

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            8, 0, 1, 20)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 12, 0, 1, 20)
        asterisco = QtWidgets.QLabel("*")
        asterisco.setFont(font20)
        asterisco.setAlignment(alignment)
        layout.addWidget(asterisco, 6, 3)
        asterisco2 = QtWidgets.QLabel("**")
        asterisco2.setFont(font20)
        asterisco2.setAlignment(alignment)
        layout.addWidget(asterisco2, 7, 3)
        asterisco_ = QtWidgets.QLabel("*")
        asterisco_.setFont(font20)
        asterisco_.setAlignment(alignment)
        layout.addWidget(asterisco_, 10, 2)
        asterisco2_ = QtWidgets.QLabel("**")
        asterisco2_.setFont(font20)
        asterisco2_.setAlignment(alignment)
        layout.addWidget(asterisco2_, 11, 2)

        butonConfig = QtWidgets.QToolButton()
        butonConfig.setIcon(QtGui.QIcon(
            os.environ["pychemqt"] +
            os.path.join("images", "button", "configure.png")))
        butonConfig.clicked.connect(self.configure)
        layout.addWidget(butonConfig, 11, 1)

        self.Info = QtWidgets.QFrame()
        layout.addWidget(self.Info, 0, 5, 3, 3)
        layoutInfo = QtWidgets.QGridLayout(self.Info)
        layoutInfo.setSpacing(1)
        layoutInfo.setContentsMargins(2, 0, 2, 0)
        self.Info.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Info.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Info.setAutoFillBackground(True)
        self.Info.setPalette(palette)
        self.numero_atomico = QtWidgets.QLabel()
        self.numero_atomico.setToolTip(
            QtWidgets.QApplication.translate("pychemqt", "Atomic number"))
        layoutInfo.addWidget(self.numero_atomico, 1, 1)
        self.simbolo = QtWidgets.QLabel()
        self.simbolo.setAlignment(alignment)
        self.simbolo.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Symbol"))
        self.simbolo.setFont(font11)
        layoutInfo.addWidget(self.simbolo, 1, 3)
        self.nombre = QtWidgets.QLabel()
        self.nombre.setAlignment(QtCore.Qt.AlignCenter)
        self.nombre.setFont(font_title)
        layoutInfo.addWidget(self.nombre, 2, 1, 1, 3)
        font8 = QtGui.QFont()
        font8.setPointSize(8)
        self.peso_atomico = QtWidgets.QLabel()
        self.peso_atomico.setFont(font8)
        self.peso_atomico.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic mass, g/mol"))
        layoutInfo.addWidget(self.peso_atomico, 3, 1)
        self.densidad = QtWidgets.QLabel()
        self.densidad.setFont(font8)
        self.densidad.setAlignment(alignment)
        self.densidad.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt",
            "Density:\nBrown: Solid, kg/l\nBlue: Liquid, kg/l\n"
            "Green: Gas, g/l"))
        layoutInfo.addWidget(self.densidad, 3, 3)
        self.Tf = QtWidgets.QLabel()
        self.Tf.setFont(font8)
        self.Tf.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Melting Point, K"))
        layoutInfo.addWidget(self.Tf, 4, 1)
        self.Heat_f = QtWidgets.QLabel()
        self.Heat_f.setFont(font8)
        self.Heat_f.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of fusion, kJmol"))
        self.Heat_f.setAlignment(alignment)
        layoutInfo.addWidget(self.Heat_f, 4, 3)
        self.Tb = QtWidgets.QLabel()
        self.Tb.setFont(font8)
        self.Tb.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Boiling Point, K"))
        layoutInfo.addWidget(self.Tb, 5, 1)
        self.Heat_b = QtWidgets.QLabel()
        self.Heat_b.setFont(font8)
        self.Heat_b.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of vaporization, kJmol"))
        self.Heat_b.setAlignment(alignment)
        layoutInfo.addWidget(self.Heat_b, 5, 3)

        self.configuracion = QtWidgets.QLabel()
        self.configuracion.setFont(font7)
        self.configuracion.setAlignment(QtCore.Qt.AlignCenter)
        self.configuracion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Electronic configuration"))
        layoutInfo.addWidget(self.configuracion, 6, 1, 1, 3)

        self.Info2 = QtWidgets.QFrame()
        layout.addWidget(self.Info2, 0, 8, 3, 3)
        layoutInfo2 = QtWidgets.QGridLayout(self.Info2)
        layoutInfo2.setSpacing(1)
        layoutInfo2.setContentsMargins(2, 0, 2, 0)
        self.Info2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Info2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.Info2.setAutoFillBackground(True)
        self.Info2.setPalette(palette)
        self.atomic_volume = QtWidgets.QLabel()
        self.atomic_volume.setFont(font8)
        self.atomic_volume.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic volume")+", cm³/mol")
        layoutInfo2.addWidget(self.atomic_volume, 1, 1)
        self.atomic_radius = QtWidgets.QLabel()
        self.atomic_radius.setFont(font8)
        self.atomic_radius.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic radius") + ", pm")
        layoutInfo2.addWidget(self.atomic_radius, 2, 1)
        self.covalent_radius = QtWidgets.QLabel()
        self.covalent_radius.setFont(font8)
        self.covalent_radius.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Covalent radius") + ", pm")
        layoutInfo2.addWidget(self.covalent_radius, 3, 1)
        self.vanderWaals_radius = QtWidgets.QLabel()
        self.vanderWaals_radius.setFont(font8)
        self.vanderWaals_radius.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Van der Waals radius")+", pm")
        layoutInfo2.addWidget(self.vanderWaals_radius, 4, 1)
        self.ionic_radii = QtWidgets.QLabel()
        self.ionic_radii.setFont(font7)
        self.ionic_radii.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Ionic radii")+", pm")
        layoutInfo2.addWidget(self.ionic_radii, 5, 1, 1, 3)
        self.electronegativity = QtWidgets.QLabel()
        self.electronegativity.setFont(font8)
        self.electronegativity.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Electronegativity, Pauling scale"))
        self.electronegativity.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.electronegativity, 1, 3)
        self.Cp = QtWidgets.QLabel()
        self.Cp.setFont(font8)
        self.Cp.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Specific heat capacitiy") + ", kJ/kgK")
        self.Cp.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Cp, 2, 3)
        self.k = QtWidgets.QLabel()
        self.k.setFont(font8)
        self.k.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal conductivity") + ", W/mK")
        self.k.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.k, 3, 3)
        self.first_ionization = QtWidgets.QLabel()
        self.first_ionization.setFont(font8)
        self.first_ionization.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "First ionization energy") + ", kJ/mol")
        self.first_ionization.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.first_ionization, 4, 3)

        self.oxidation = QtWidgets.QLabel()
        self.oxidation.setFont(font8)
        self.oxidation.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Oxidation states"))
        self.oxidation.setAlignment(
            QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.oxidation, 6, 1, 1, 3)

        elemento = Elemental(1)
        self.actualizar(elemento)
        logging.info(QtWidgets.QApplication.translate(
            "pychemqt", "Starting periodic table tool"))

    def populate(self):
        CATEGORIES, PROP, COLORS, PMAX = _configValues(self.Preferences)
        for i in range(1, 119):
            element = Elemental(i)
            b = boton(element, CATEGORIES, PROP, COLORS, PMAX, self)
            if element.group == 0:
                if i < 80:
                    j = i-58
                else:
                    j = i-90
                self.layout().addWidget(b, element.period+4, j+4)
            elif i == 57 or i == 89:
                self.layout().addWidget(b, element.period+4, element.group)
            else:
                self.layout().addWidget(b, element.period, element.group)

    def actualizar(self, elemento):
        """Update botton info with data for current element"""
        self.numero_atomico.setText(str(elemento.id))
        self.nombre.setText(elemento.name)
        self.simbolo.setText(elemento.symbol)
        self.textData(self.peso_atomico, elemento.atomic_mass)
        if elemento.density_Solid:
            color = "#A52A2A"
            value = elemento.density_Solid
        elif elemento.density_Liq:
            color = "#0000FF"
            value = elemento.density_Liq
        elif elemento.density_Gas:
            color = "#0A640A"
            value = elemento.density_Gas
        else:
            color = "#888888"
            value = "N/A"
        self.densidad.setText("<font color={}>{}</font>".format(color, value))

        self.textData(self.Tf, elemento.Tf, color="#0000C8")
        self.textData(self.Heat_f, elemento.Heat_f, color="#0000C8")
        self.textData(self.Tb, elemento.Tb, color="#C80000")
        self.textData(self.Heat_b, elemento.Heat_b, color="#C80000")
        self.textData(self.configuracion, elemento.electron_configuration)

        self.textData(self.atomic_volume, elemento.atomic_volume)
        self.textData(self.atomic_radius, elemento.atomic_radius)
        self.textData(self.ionic_radii, elemento.ionic_radii)
        self.textData(self.covalent_radius, elemento.covalent_radius)
        self.textData(self.vanderWaals_radius, elemento.vanderWaals_radius)

        self.textData(self.electronegativity, elemento.electronegativity)
        self.textData(self.Cp, elemento.Cp)
        self.textData(self.k, elemento.k)
        self.textData(self.first_ionization, elemento.first_ionization)
        self.textData(self.oxidation, elemento.oxidation)

    def textData(self, widget, data, color="#000000", color_dis="#888888"):
        if not data:
            widget.setText("<font color={:}>N/A</font>".format(color_dis))
        else:
            widget.setText("<font color={:}>{:}</font>".format(color, data))

    def configure(self):
        from UI.prefElemental import Dialog
        dlg = Dialog(self.Preferences)
        if dlg.exec_():
            self.Preferences = dlg.value(self.Preferences)
            self.Preferences.write(open(conf_dir+"pychemqtrc", "w"))
            self.populate()


class boton(QtWidgets.QPushButton):
    """Button widget to define a element"""
    def __init__(self, element, CATEGORIES, PROP, COLORS, PMAX, parent=None):
        """Constructor,
        element: the atomic number, used as id to tagged button
        """
        super(boton, self).__init__(parent)
        self.setFixedSize(40, 35)
        self.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Click for view properties"))
        self.parent = parent
        self.Element = element

        if COLORS and not PMAX:
            if PROP == "group_element":
                PROP = "group"
            prop = element.__getattribute__(PROP)
            color = COLORS[CATEGORIES.index(prop)]
        elif PMAX:
            prop = cleanFloat(element.__getattribute__(PROP))
            if prop:
                if prop > CATEGORIES[-1]:
                    color = COLORS[-1]
                else:
                    index = 1
                    while True:
                        if prop <= CATEGORIES[index]:
                            color = COLORS[index-1]
                            break
                        else:
                            index += 1
            else:
                color = "#DDDDDD"

        else:
            color = self.Element.color
        self.setStyleSheet("QPushButton { background-color: %s }" % color)
        self.setText(self.Element.symbol)
        self.clicked.connect(self.press)

    def enterEvent(self, event):
        """Enter event to change boton info with new element data"""
        self.parent.actualizar(self.Element)
        event.accept()

    def press(self):
        """Press button show elementDialog with other element data"""
        dialog = ElementDialog(self.Element)
        dialog.exec_()


class ElementDialog(QtWidgets.QDialog):
    """Dialog to show all element properties"""
    def __init__(self, elemento, parent=None):
        super(ElementDialog, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate(
            "pychemqt", "Properties of "+elemento.name))
        lyt = QtWidgets.QVBoxLayout(self)
        tabWidget = QtWidgets.QTabWidget()
        lyt.addWidget(tabWidget)
        btbox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        btbox.rejected.connect(self.reject)
        lyt.addWidget(btbox)

        tabGeneral = QtWidgets.QWidget()
        layoutGeneral = QtWidgets.QGridLayout(tabGeneral)
        layoutGeneral.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Name:")), 1, 1)
        layoutGeneral.addWidget(QtWidgets.QLabel(elemento.name), 1, 2)
        layoutGeneral.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Serie:")), 2, 1)
        layoutGeneral.addWidget(QtWidgets.QLabel(elemento.serie), 2, 2)
        layoutGeneral.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Group")), 3, 1)
        layoutGeneral.addWidget(QtWidgets.QLabel(str(elemento.group)), 3, 2)
        layoutGeneral.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Period")), 4, 1)
        layoutGeneral.addWidget(QtWidgets.QLabel(str(elemento.period)), 4, 2)
        layoutGeneral.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Block")), 5, 1)
        layoutGeneral.addWidget(QtWidgets.QLabel(elemento.block), 5, 2)

        layoutGeneral.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1)
        label = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "History"))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        label.setFont(font)
        layoutGeneral.addWidget(label, 7, 1, 1, 3)
        label_8 = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Discovery") + ": " +
            elemento.country + "(" + elemento.country + ")" + os.linesep +
            QtWidgets.QApplication.translate("pychemqt", "Discovered by ") +
            elemento.discover + os.linesep +
            QtWidgets.QApplication.translate("pychemqt", "Etymology") + ": " +
            elemento.etymology)
        label_8.setMargin(5)
        label_8.setWordWrap(True)
        layoutGeneral.addWidget(label_8, 8, 1, 1, 3)
        self.botoncito = QtWidgets.QLabel()
        self.botoncito.setStyleSheet(
            "background-color: %s;" % elemento.color)
        self.botoncito.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.botoncito.setFixedSize(60, 60)
        layoutGeneral.addWidget(self.botoncito, 1, 5, 3, 1)
        label = QtWidgets.QLabel()
        label.setText(str(elemento.id))
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(label, 1, 5)
        label = QtWidgets.QLabel(elemento.symbol)
        font.setPointSize(12)
        label.setFont(font)
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(label, 2, 5)
        layoutGeneral.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 9, 3)
        tabWidget.addTab(tabGeneral, QtWidgets.QApplication.translate(
            "pychemqt", "General"))

        tabFisica = QtWidgets.QWidget()
        lytphy = QtWidgets.QGridLayout(tabFisica)
        lytphy.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Fase:")), 1, 1)
        lytphy.addWidget(QtWidgets.QLabel(elemento.phase + " a 0ºC"),
                         1, 2, 1, 1)
        if elemento.density_Solid:
            lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Solid Density:")), 2, 1)
            lytphy.addWidget(self.drawData(
                unidades.Density, elemento.density_Solid,
                "gcc", txt=" @ 20ºC"), 2, 2)
        if elemento.density_Liq:
            lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Liquid Density:")), 3, 1)
            lytphy.addWidget(self.drawData(
                unidades.Density, elemento.density_Liq, "gcc", txt=" " +
                QtWidgets.QApplication.translate(
                    "pychemqt", "at melting point")), 3, 2)
        if elemento.density_Gas:
            lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Gas Density:")), 4, 1)
            lytphy.addWidget(self.drawData(
                unidades.Density, elemento.density_Gas, "gl", txt=" @ 0ºC"),
                4, 2)
        lytphy.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Appearance:")), 5, 1)
        label = QtWidgets.QLabel(elemento.appearance)
        label.setWordWrap(True)
        lytphy.addWidget(label, 5, 2)
        lytphy.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            6, 1, 1, 3)
        label = QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal properties"))
        label.setFont(font)
        lytphy.addWidget(label, 7, 1)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Melting point:")), 8, 1)
        self.punto_fusion = self.drawData(
            unidades.Temperature, elemento.Tf)
        lytphy.addWidget(self.punto_fusion, 8, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Boiling point:")), 9, 1)
        self.punto_ebullicion = self.drawData(
            unidades.Temperature, elemento.Tb)
        lytphy.addWidget(self.punto_ebullicion, 9, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of fusion:")), 10, 1)
        self.calor_fusion = self.drawData(
            unidades.MolarEnthalpy, elemento.Heat_f, "kJkmol")
        lytphy.addWidget(self.calor_fusion, 10, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of vaporization:")), 11, 1)
        self.calor_vaporizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.Heat_b, "kJkmol")
        lytphy.addWidget(self.calor_vaporizacion, 11, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Specific heat capacity:")), 12, 1)
        self.capacidad_calorifica = self.drawData(
            unidades.SpecificHeat, elemento.Cp, "JgK")
        lytphy.addWidget(self.capacidad_calorifica, 12, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal conductivity:")), 13, 1)
        self.conductividad_termica = self.drawData(
            unidades.ThermalConductivity, elemento.k,
            txt=" @ 300K")
        lytphy.addWidget(self.conductividad_termica, 13, 2)
        lytphy.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Debye Temperature:")), 14, 1)
        self.temperatura_debye = self.drawData(
            unidades.Temperature, elemento.T_debye)
        lytphy.addWidget(self.temperatura_debye, 14, 2)
        lytphy.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 15, 1, 1, 3)
        tabWidget.addTab(tabFisica, QtWidgets.QApplication.translate(
            "pychemqt", "Physical properties"))

        tabAtom = QtWidgets.QWidget()
        lyt_A = QtWidgets.QGridLayout(tabAtom)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic mass:")), 1, 1)
        if elemento.atomic_mass:
            self.masa_atomica = QtWidgets.QLabel(
                str(elemento.atomic_mass) + " g/mol")
        else:
            self.masa_atomica = QtWidgets.QLabel(elemento.atomic_mass)
        lyt_A.addWidget(self.masa_atomica, 1, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic Volume:")), 2, 1)
        self.volumen_atomico = self.drawData(unidades.MolarVolume,
                                             elemento.atomic_volume)
        lyt_A.addWidget(self.volumen_atomico, 2, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic radius:")), 3, 1)
        if elemento.atomic_radius:
            self.radio_atomico = QtWidgets.QLabel(
                str(elemento.atomic_radius) + " pm")
        else:
            self.radio_atomico = QtWidgets.QLabel(str(elemento.atomic_radius))
        lyt_A.addWidget(self.radio_atomico, 3, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Covalent radius:")), 4, 1)
        if elemento.covalent_radius:
            self.radio_covalente = QtWidgets.QLabel(
                str(elemento.covalent_radius) + " pm")
        else:
            self.radio_covalente = QtWidgets.QLabel(
                str(elemento.covalent_radius))
        lyt_A.addWidget(self.radio_covalente, 4, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Van der Waals radius:")), 5, 1)
        if elemento.vanderWaals_radius:
            self.radio_waals = QtWidgets.QLabel(
                str(elemento.vanderWaals_radius) + " pm")
        else:
            self.radio_waals = QtWidgets.QLabel(
                str(elemento.vanderWaals_radius))
        lyt_A.addWidget(self.radio_waals, 5, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Ionic radii:")), 6, 1)
        if elemento.ionic_radii:
            self.radio_ionico = QtWidgets.QLabel(
                str(elemento.ionic_radii) + " pm")
        else:
            self.radio_ionico = QtWidgets.QLabel(str(elemento.ionic_radii))
        lyt_A.addWidget(self.radio_ionico, 6, 2)
        lyt_A.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            7, 1, 1, 3)
        label = QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electronic properties"))
        label.setFont(font)
        lyt_A.addWidget(label, 8, 1)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electronic configuration:")), 9, 1)
        lyt_A.addWidget(QtWidgets.QLabel(
            elemento.electron_configuration), 9, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Oxidation states:")), 10, 1)
        lyt_A.addWidget(QtWidgets.QLabel(
            elemento.oxidation), 10, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electronegativity:")), 11, 1)
        lyt_A.addWidget(QtWidgets.QLabel(
            str(elemento.electronegativity)), 11, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electron affinity:")), 12, 1)
        self.afinidad_electronica = self.drawData(
            unidades.MolarEnthalpy, elemento.electron_affinity, "kJkmol")
        lyt_A.addWidget(self.afinidad_electronica, 12, 2)
        lyt_A.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "1st ionization energy:")), 13, 1)
        self.energia_ionizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.first_ionization, "kJkmol")
        lyt_A.addWidget(self.energia_ionizacion, 13, 2)
        lyt_A.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 14, 1, 1, 3)
        tabWidget.addTab(tabAtom, QtWidgets.QApplication.translate(
            "pychemqt", "Atomic properties"))

        tabCristal = QtWidgets.QWidget()
        lyt_C = QtWidgets.QGridLayout(tabCristal)
        lyt_C.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Lattice type:")), 1, 1)
        lyt_C.addWidget(QtWidgets.QLabel(
            elemento.lattice_type), 1, 2)
        lyt_C.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Space group:")), 2, 1)
        lyt_C.addWidget(QtWidgets.QLabel(
            elemento.space_group), 2, 2)
        lyt_C.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate(
                "pychemqt", "Lattice edge lengths:")), 3, 1)
        self.lados = QtWidgets.QLabel()
        if elemento.lattice_edges:
            self.lados.setText("%spm, %spm, %spm" % (
                elemento.lattice_edges[0], elemento.lattice_edges[1],
                elemento.lattice_edges[2]))
        else:
            self.lados.setText(elemento.lattice_edges)
        lyt_C.addWidget(self.lados, 3, 2)
        lyt_C.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Lattice angles:")), 4, 1)
        self.angulos = QtWidgets.QLabel()
        if elemento.lattice_angles:
            self.angulos.setText("%sº, %sº, %sº" % (
                elemento.lattice_angles[0], elemento.lattice_angles[1],
                elemento.lattice_angles[2]))
        else:
            self.angulos.setText(elemento.lattice_angles)
        lyt_C.addWidget(self.angulos, 4, 2)
        lyt_C.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Lattice unit volume:")), 5, 1)
        self.volumen_celda = QtWidgets.QLabel()
        if elemento.lattice_angles:
            self.volumen_celda.setText("%0.5f mm<sup>3</sup>"
                                       % elemento.lattice_volume)
        else:
            self.volumen_celda.setText(elemento.lattice_volume)
        lyt_C.addWidget(self.volumen_celda, 5, 2)
        lyt_C.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 6, 1, 1, 3)

        label = QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Isotopes"))
        label.setFont(font)
        lyt_C.addWidget(label, 8, 1)
        title = [QtWidgets.QApplication.translate("pychemqt", "Mass Number"),
                 QtWidgets.QApplication.translate("pychemqt", "Mass"),
                 QtWidgets.QApplication.translate("pychemqt", "Abundance")]
        self.isotopes = Tabla(3, horizontalHeader=title, readOnly=True,
                              stretch=True, verticalHeader=False)
        self.isotopes.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectRows)
        self.isotopes.setColumn(0, [iso[0] for iso in elemento.isotopes],
                                decimales=0)
        self.isotopes.setColumn(1, [iso[1] for iso in elemento.isotopes],
                                format=1, decimales=10)
        self.isotopes.setColumn(2, [iso[2] for iso in elemento.isotopes],
                                format=1, decimales=10)
        lyt_C.addWidget(self.isotopes, 9, 1, 1, 2)
        lyt_C.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding), 10, 1, 1, 3)

        tabWidget.addTab(
            tabCristal,
            QtWidgets.QApplication.translate("pychemqt", "Crystallographic"))

    def drawData(self, unit, data, unidad="", txt=""):
        """Return a widget with the data inprint
            unit: unidad subclass
            data: value to set
            unidad: unit of value to show
            txt: opcional txt to show for unit"""
        if data and unidad:
            value = unit(data, unidad)
            widget = Entrada_con_unidades(
                unit, readOnly=True, value=value, textounidad=txt)
        elif data:
            widget = Entrada_con_unidades(
                unit, readOnly=True, value=data, textounidad=txt)
        else:
            widget = QtWidgets.QLabel(str(data))
        return widget


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = qtelemental()
    Form.show()
    sys.exit(app.exec_())
