#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Basic periodic table with properties dialog
#   - qtelemental: Periodic table
#   - boton: Element button in periodic table
#   - ElementDialog: Dialog to show properties of components
#
# Use gelemental python binding for library data, that library now is loose so
# maybe I need add that library in the future to pychemqt
###############################################################################

import os
from configparser import ConfigParser

from PyQt5 import QtCore, QtGui, QtWidgets

from lib import unidades
from lib.config import conf_dir
from lib.elemental import Elemental, color_serie, color_block, color_phase
from UI.widgets import Entrada_con_unidades

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

Preferences = ConfigParser()
Preferences.read(conf_dir+"pychemqtrc")

class qtelemental(QtWidgets.QDialog):
    """Periodic table graph"""
    def __init__(self, parent=None):
        super(qtelemental, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/PeriodicTableIcon.png")))
        self.setWindowTitle(
            QtWidgets.QApplication.translate("pychemqt", "Periodic Table"))
        layout = QtWidgets.QGridLayout(self)
        layout.setSpacing(2)
        for i in range(1, 119):
            element = Elemental(i)
            b = boton(element, self)
            if element.group == 0:
                if i < 80:
                    j = i-58
                else:
                    j = i-90
                layout.addWidget(b, element.period+4, j+4, 1, 1)
            elif i == 57 or i == 89:
                layout.addWidget(b, element.period+4, element.group, 1, 1)
            else:
                layout.addWidget(b, element.period, element.group, 1, 1)
        layout.addItem(QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Fixed,
                                         QtWidgets.QSizePolicy.Fixed), 8, 0, 1, 20)
        layout.addItem(QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Expanding,
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
        self.densidad.setToolTip(QtWidgets.QApplication.translate("pychemqt",
            "Density:\nBrown: Solid, kg/l\nBlue: Liquid, kg/l\nGreen: Gas, g/l"))
        layoutInfo.addWidget(self.densidad, 3, 3)
        self.Tfusion = QtWidgets.QLabel()
        self.Tfusion.setFont(font8)
        self.Tfusion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Melting Point, K"))
        layoutInfo.addWidget(self.Tfusion, 4, 1)
        self.Calorfusion = QtWidgets.QLabel()
        self.Calorfusion.setFont(font8)
        self.Calorfusion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of fusion, kJmol"))
        self.Calorfusion.setAlignment(alignment)
        layoutInfo.addWidget(self.Calorfusion, 4, 3)
        self.Tebullicion = QtWidgets.QLabel()
        self.Tebullicion.setFont(font8)
        self.Tebullicion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Boiling Point, K"))
        layoutInfo.addWidget(self.Tebullicion, 5, 1)
        self.Calorebullicion = QtWidgets.QLabel()
        self.Calorebullicion.setFont(font8)
        self.Calorebullicion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of vaporization, kJmol"))
        self.Calorebullicion.setAlignment(alignment)
        layoutInfo.addWidget(self.Calorebullicion, 5, 3)

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
        self.VolumenAtomico = QtWidgets.QLabel()
        self.VolumenAtomico.setFont(font8)
        self.VolumenAtomico.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic volume")+", cm³/mol")
        layoutInfo2.addWidget(self.VolumenAtomico, 1, 1)
        self.RadioAtomico = QtWidgets.QLabel()
        self.RadioAtomico.setFont(font8)
        self.RadioAtomico.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic radius") + ", pm")
        layoutInfo2.addWidget(self.RadioAtomico, 2, 1)
        self.RadioCovalente = QtWidgets.QLabel()
        self.RadioCovalente.setFont(font8)
        self.RadioCovalente.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Covalent radius") + ", pm")
        layoutInfo2.addWidget(self.RadioCovalente, 3, 1)
        self.RadioVanderWaals = QtWidgets.QLabel()
        self.RadioVanderWaals.setFont(font8)
        self.RadioVanderWaals.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Van der Waals radius")+", pm")
        layoutInfo2.addWidget(self.RadioVanderWaals, 4, 1)
        self.RadioIonico = QtWidgets.QLabel()
        self.RadioIonico.setFont(font7)
        self.RadioIonico.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Ionic radii")+", pm")
        layoutInfo2.addWidget(self.RadioIonico, 5, 1, 1, 3)
        self.Electronegatividad = QtWidgets.QLabel()
        self.Electronegatividad.setFont(font8)
        self.Electronegatividad.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Electronegativity, Pauling scale"))
        self.Electronegatividad.setAlignment(QtCore.Qt.AlignRight |
                                             QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Electronegatividad, 1, 3)
        self.Cp = QtWidgets.QLabel()
        self.Cp.setFont(font8)
        self.Cp.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Specific heat capacitiy") + ", kJ/kgK")
        self.Cp.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Cp, 2, 3)
        self.Conductividad = QtWidgets.QLabel()
        self.Conductividad.setFont(font8)
        self.Conductividad.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal conductivity") + ", W/mK")
        self.Conductividad.setAlignment(QtCore.Qt.AlignRight |
                                        QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Conductividad, 3, 3)
        self.Ionizacion = QtWidgets.QLabel()
        self.Ionizacion.setFont(font8)
        self.Ionizacion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "First ionization energy") + ", kJ/mol")
        self.Ionizacion.setAlignment(QtCore.Qt.AlignRight |
                                     QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Ionizacion, 4, 3)

        self.EstadoOxidacion = QtWidgets.QLabel()
        self.EstadoOxidacion.setFont(font8)
        self.EstadoOxidacion.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Oxidation states"))
        self.EstadoOxidacion.setAlignment(QtCore.Qt.AlignCenter |
                                          QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.EstadoOxidacion, 6, 1, 1, 3)

        elemento = Elemental(1)
        self.actualizar(elemento)

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
        self.densidad.setText("<font color={:}>{:}</font>".format(color, value))

        self.textData(self.Tfusion, elemento.Tf, color="#0000C8")
        self.textData(self.Calorfusion, elemento.Heat_f, color="#0000C8")
        self.textData(self.Tebullicion, elemento.Tb, color="#C80000")
        self.textData(self.Calorebullicion, elemento.Heat_b, color="#C80000")
        self.textData(self.configuracion, elemento.electron_configuration)

        self.textData(self.VolumenAtomico, elemento.atomic_volume)
        self.textData(self.RadioAtomico, elemento.atomic_radius)
        self.textData(self.RadioIonico, elemento.ionic_radii)
        self.textData(self.RadioCovalente, elemento.covalent_radius)
        self.textData(self.RadioVanderWaals, elemento.vanderWaals_radius)

        self.textData(self.Electronegatividad, elemento.electronegativity)
        self.textData(self.Cp, elemento.Cp)
        self.textData(self.Conductividad, elemento.k)
        self.textData(self.Ionizacion, elemento.first_ionization)
        self.textData(self.EstadoOxidacion, elemento.oxidation)

    def textData(self, widget, data, color="#000000", color_disabled="#888888"):
        if not data:
            widget.setText("<font color={:}>N/A</font>".format(color_disabled))
        elif isinstance(data, str) or not data.code:
            widget.setText("<font color={:}>{:}</font>".format(color, data))
        else:
            widget.setText("<font color={:}>N/A</font>".format(color_disabled))


class boton(QtWidgets.QPushButton):
    """Button widget to define a element"""
    def __init__(self, element, parent=None):
        """Constructor, i parameter is the atomic number and useful ad id to
        tagged button"""
        super(boton, self).__init__(parent)
        self.setFixedSize(40, 35)
        self.setToolTip(QtWidgets.QApplication.translate(
            "pychemqt", "Click for view properties"))
        self.parent = parent
        self.Element = element
        
        if Preferences.getint("Applications", "elementalColorby") == 1:
            if self.Element.serie:
                color = color_serie[self.Element.serie]
            else:
                color = "#DDDDDD"
        elif Preferences.getint("Applications", "elementalColorby") == 2:
            if self.Element.block:
                color = color_block[self.Element.block]
            else:
                color = "#DDDDDD"
        elif Preferences.getint("Applications", "elementalColorby") == 3:
            if self.Element.phase:
                color = color_phase[self.Element.phase]
            else:
                color = "#DDDDDD"
        else:
            color = self.Element.color
        self.setPalette(QtGui.QPalette(QtGui.QColor(color)))
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
        self.setFixedSize(505, 432)
        self.verticallayout = QtWidgets.QVBoxLayout(self)
        self.tabWidget = QtWidgets.QTabWidget()
        self.verticallayout.addWidget(self.tabWidget)
        self.Aceptar = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        self.Aceptar.rejected.connect(self.reject)
        self.verticallayout.addWidget(self.Aceptar)

        self.tabGeneral = QtWidgets.QWidget()
        layoutGeneral = QtWidgets.QGridLayout(self.tabGeneral)
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
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed), 6, 1)
        label_6 = QtWidgets.QLabel(QtWidgets.QApplication.translate("pychemqt", "History"))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        label_6.setFont(font)
        layoutGeneral.addWidget(label_6, 7, 1, 1, 3)
        label_8 = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Discovery") + ": " +
            elemento.country + "(" + elemento.country + ")" +os.linesep +
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
        self.label_10 = QtWidgets.QLabel()
        self.label_10.setText(str(elemento.id))
        self.label_10.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(self.label_10, 1, 5)
        self.label_38 = QtWidgets.QLabel(elemento.symbol)
        font.setPointSize(12)
        self.label_38.setFont(font)
        self.label_38.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(self.label_38, 2, 5)
        layoutGeneral.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding), 9, 3)
        self.tabWidget.addTab(self.tabGeneral, QtWidgets.QApplication.translate(
            "pychemqt", "General"))

        self.tab_Fisica = QtWidgets.QWidget()
        layoutFisica = QtWidgets.QGridLayout(self.tab_Fisica)
        layoutFisica.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Fase:")), 1, 1)
        layoutFisica.addWidget(QtWidgets.QLabel(elemento.phase +
                                            " a 0ºC"), 1, 2, 1, 1)
        if elemento.density_Solid:
            layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Solid Density:")), 2, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_Solid, "gcc", txt=" @ 20ºC"), 2, 2)
        if elemento.density_Liq:
            layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Liquid Density:")), 3, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_Liq, "gcc", txt=" " +
                QtWidgets.QApplication.translate("pychemqt", "at melting point")), 3, 2)
        if elemento.density_Gas:
            layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
                "pychemqt", "Gas Density:")), 4, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_Gas, "gl", txt=" @ 0ºC"), 4, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Appearance:")), 5, 1)
        label = QtWidgets.QLabel(elemento.appearance)
        label.setWordWrap(True)
        layoutFisica.addWidget(label, 5, 2)
        layoutFisica.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed), 6, 1, 1, 3)
        self.label_7 = QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal properties"))
        self.label_7.setFont(font)
        layoutFisica.addWidget(self.label_7, 7, 1)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Melting point:")), 8, 1)
        self.punto_fusion = self.drawData(
            unidades.Temperature, elemento.Tf)
        layoutFisica.addWidget(self.punto_fusion, 8, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Boiling point:")), 9, 1)
        self.punto_ebullicion = self.drawData(
            unidades.Temperature, elemento.Tb)
        layoutFisica.addWidget(self.punto_ebullicion, 9, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of fusion:")), 10, 1)
        self.calor_fusion = self.drawData(
            unidades.MolarEnthalpy, elemento.Heat_f, "kJkmol")
        layoutFisica.addWidget(self.calor_fusion, 10, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Heat of vaporization:")), 11, 1)
        self.calor_vaporizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.Heat_b, "kJkmol")
        layoutFisica.addWidget(self.calor_vaporizacion, 11, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Specific heat capacity:")), 12, 1)
        self.capacidad_calorifica = self.drawData(
            unidades.SpecificHeat, elemento.Cp, "JgK")
        layoutFisica.addWidget(self.capacidad_calorifica, 12, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Thermal conductivity:")), 13, 1)
        self.conductividad_termica = self.drawData(
            unidades.ThermalConductivity, elemento.k,
            txt=" @ 300K")
        layoutFisica.addWidget(self.conductividad_termica, 13, 2)
        layoutFisica.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Debye Temperature:")), 14, 1)
        self.temperatura_debye = self.drawData(
            unidades.Temperature, elemento.T_debye)
        layoutFisica.addWidget(self.temperatura_debye, 14, 2)
        layoutFisica.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            15, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Fisica, QtWidgets.QApplication.translate(
            "pychemqt", "Physical properties"))

        self.tab_Atomo = QtWidgets.QWidget()
        layout_Atomo = QtWidgets.QGridLayout(self.tab_Atomo)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Atomic mass:")), 1, 1)
        if elemento.atomic_mass:
            self.masa_atomica = QtWidgets.QLabel(
                str(elemento.atomic_mass) +" g/mol")
        else:
            self.masa_atomica = QtWidgets.QLabel(elemento.atomic_mass)
        layout_Atomo.addWidget(self.masa_atomica, 1, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Atomic Volume:")), 2, 1)
        self.volumen_atomico = self.drawData(unidades.MolarVolume,
                                             elemento.atomic_volume)
        layout_Atomo.addWidget(self.volumen_atomico, 2, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Atomic radius:")), 3, 1)
        if elemento.atomic_radius:
            self.radio_atomico = QtWidgets.QLabel(
                str(elemento.atomic_radius) + " pm")
        else:
            self.radio_atomico = QtWidgets.QLabel(elemento.atomic_radius)
        layout_Atomo.addWidget(self.radio_atomico, 3, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Covalent radius:")), 4, 1)
        if elemento.covalent_radius:
            self.radio_covalente = QtWidgets.QLabel(
                str(elemento.covalent_radius) + " pm")
        else:
            self.radio_covalente = QtWidgets.QLabel(
                str(elemento.covalent_radius))
        layout_Atomo.addWidget(self.radio_covalente, 4, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Van der Waals radius:")), 5, 1)
        if elemento.vanderWaals_radius:
            self.radio_waals = QtWidgets.QLabel(
                str(elemento.vanderWaals_radius) + " pm")
        else:
            self.radio_waals = QtWidgets.QLabel(
                str(elemento.vanderWaals_radius))
        layout_Atomo.addWidget(self.radio_waals, 5, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Ionic radii:")), 6, 1)
        if elemento.ionic_radii:
            self.radio_ionico = QtWidgets.QLabel(
                str(elemento.ionic_radii) + " pm")
        else:
            self.radio_ionico = QtWidgets.QLabel(str(elemento.ionic_radii))
        layout_Atomo.addWidget(self.radio_ionico, 6, 2)
        layout_Atomo.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed),
            7, 1, 1, 3)
        self.label_9 = QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Electronic properties"))
        self.label_9.setFont(font)
        layout_Atomo.addWidget(self.label_9, 8, 1)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electronic configuration:")), 9, 1)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            elemento.electron_configuration), 9, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Oxidation states:")), 10, 1)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            elemento.oxidation), 10, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electronegativity:")), 11, 1)
        layout_Atomo.addWidget(QtWidgets.QLabel(
            str(elemento.electronegativity)), 11, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "Electron affinity:")), 12, 1)
        self.afinidad_electronica = self.drawData(
            unidades.MolarEnthalpy, elemento.electron_affinity, "kJkmol")
        layout_Atomo.addWidget(self.afinidad_electronica, 12, 2)
        layout_Atomo.addWidget(QtWidgets.QLabel(QtWidgets.QApplication.translate(
            "pychemqt", "1st ionization energy:")), 13, 1)
        self.energia_ionizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.first_ionization, "kJkmol")
        layout_Atomo.addWidget(self.energia_ionizacion, 13, 2)
        layout_Atomo.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            14, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Atomo, QtWidgets.QApplication.translate(
            "pychemqt", "Atomic properties"))

        self.tab_Cristalografia = QtWidgets.QWidget()
        layout_Cristalografia = QtWidgets.QGridLayout(self.tab_Cristalografia)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Lattice type:")), 1, 1)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            elemento.lattice_type), 1, 2)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Space group:")), 2, 1)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            elemento.space_group), 2, 2)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Lattice edge lengths:")), 3, 1)
        self.lados = QtWidgets.QLabel()
        if elemento.lattice_edges:
            self.lados.setText(
                "%spm, %spm, %spm" % (elemento.lattice_edges[0], 
                                       elemento.lattice_edges[1], 
                                       elemento.lattice_edges[2]))
        else:
            self.lados.setText(elemento.lattice_edges)
        layout_Cristalografia.addWidget(self.lados, 3, 2)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Lattice angles:")), 4, 1)
        self.angulos = QtWidgets.QLabel()
        if elemento.lattice_angles:
            self.angulos.setText(
                "%sº, %sº, %sº" % (elemento.lattice_angles[0], 
                                   elemento.lattice_angles[1], 
                                   elemento.lattice_angles[2]))
        else:
            self.angulos.setText(elemento.lattice_angles)
        layout_Cristalografia.addWidget(self.angulos, 4, 2)
        layout_Cristalografia.addWidget(QtWidgets.QLabel(
            QtWidgets.QApplication.translate("pychemqt", "Lattice unit volume:")), 5, 1)
        self.volumen_celda = QtWidgets.QLabel()
        if elemento.lattice_angles:
            vol = elemento.lattice_edges[0]*elemento.lattice_edges[1]*elemento.lattice_edges[2]/1e9
            self.volumen_celda.setText("%0.5f mm<sup>3</sup>" % vol)
        else:
            self.volumen_celda.setText(elemento.lattice_volume)
        layout_Cristalografia.addWidget(self.volumen_celda, 5, 2)
        layout_Cristalografia.addItem(QtWidgets.QSpacerItem(
            20, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding),
            6, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Cristalografia,
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
