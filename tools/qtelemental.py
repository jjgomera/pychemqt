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

from PyQt4 import QtCore, QtGui
import Elemental

from lib import unidades
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


class qtelemental(QtGui.QDialog):
    """Periodic table graph"""
    def __init__(self, parent=None):
        super(qtelemental, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(
            os.environ["pychemqt"]+"/images/button/PeriodicTableIcon.png")))
        self.setWindowTitle(
            QtGui.QApplication.translate("pychemqt", "Periodic Table"))
        layout = QtGui.QGridLayout(self)
        layout.setSpacing(2)
        for i in range(1, 119):
            b = boton(i, self)
            if Elemental.get_element(i).group.value == 0:
                if i < 80:
                    j = i-58
                else:
                    j = i-90
                layout.addWidget(b, Elemental.get_element(i).period.value+4,
                                 j+4, 1, 1)
            elif i == 57 or i == 89:
                layout.addWidget(b, Elemental.get_element(i).period.value+4,
                                 Elemental.get_element(i).group.value, 1, 1)
            else:
                layout.addWidget(b, Elemental.get_element(i).period.value,
                                 Elemental.get_element(i).group.value, 1, 1)
        layout.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Fixed,
                                         QtGui.QSizePolicy.Fixed), 8, 0, 1, 20)
        layout.addItem(QtGui.QSpacerItem(10, 10, QtGui.QSizePolicy.Expanding,
                                         QtGui.QSizePolicy.Expanding), 12, 0, 1, 20)
        asterisco = QtGui.QLabel("*")
        asterisco.setFont(font20)
        asterisco.setAlignment(alignment)
        layout.addWidget(asterisco, 6, 3)
        asterisco2 = QtGui.QLabel("**")
        asterisco2.setFont(font20)
        asterisco2.setAlignment(alignment)
        layout.addWidget(asterisco2, 7, 3)
        asterisco_ = QtGui.QLabel("*")
        asterisco_.setFont(font20)
        asterisco_.setAlignment(alignment)
        layout.addWidget(asterisco_, 10, 2)
        asterisco2_ = QtGui.QLabel("**")
        asterisco2_.setFont(font20)
        asterisco2_.setAlignment(alignment)
        layout.addWidget(asterisco2_, 11, 2)

        self.Info = QtGui.QFrame()
        layout.addWidget(self.Info, 0, 5, 3, 3)
        layoutInfo = QtGui.QGridLayout(self.Info)
        layoutInfo.setSpacing(1)
        layoutInfo.setContentsMargins(2, 0, 2, 0)
        self.Info.setFrameShape(QtGui.QFrame.StyledPanel)
        self.Info.setFrameShadow(QtGui.QFrame.Raised)
        self.Info.setAutoFillBackground(True)
        self.Info.setPalette(palette)
        self.numero_atomico = QtGui.QLabel()
        self.numero_atomico.setToolTip(
            QtGui.QApplication.translate("pychemqt", "Atomic number"))
        layoutInfo.addWidget(self.numero_atomico, 1, 1)
        self.simbolo = QtGui.QLabel()
        self.simbolo.setAlignment(alignment)
        self.simbolo.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Symbol"))
        self.simbolo.setFont(font11)
        layoutInfo.addWidget(self.simbolo, 1, 3)
        self.nombre = QtGui.QLabel()
        self.nombre.setAlignment(QtCore.Qt.AlignCenter)
        self.nombre.setFont(font_title)
        layoutInfo.addWidget(self.nombre, 2, 1, 1, 3)
        font8 = QtGui.QFont()
        font8.setPointSize(8)
        self.peso_atomico = QtGui.QLabel()
        self.peso_atomico.setFont(font8)
        self.peso_atomico.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Atomic mass, g/mol"))
        layoutInfo.addWidget(self.peso_atomico, 3, 1)
        self.densidad = QtGui.QLabel()
        self.densidad.setFont(font8)
        self.densidad.setAlignment(alignment)
        self.densidad.setToolTip(QtGui.QApplication.translate("pychemqt",
            "Density:\nBrown: Solid, kg/l\nBlue: Liquid, kg/l\nGreen: Gas, g/l"))
        layoutInfo.addWidget(self.densidad, 3, 3)
        self.Tfusion = QtGui.QLabel()
        self.Tfusion.setFont(font8)
        self.Tfusion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Melting Point, K"))
        layoutInfo.addWidget(self.Tfusion, 4, 1)
        self.Calorfusion = QtGui.QLabel()
        self.Calorfusion.setFont(font8)
        self.Calorfusion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Heat of fusion, kJmol"))
        self.Calorfusion.setAlignment(alignment)
        layoutInfo.addWidget(self.Calorfusion, 4, 3)
        self.Tebullicion = QtGui.QLabel()
        self.Tebullicion.setFont(font8)
        self.Tebullicion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Boiling Point, K"))
        layoutInfo.addWidget(self.Tebullicion, 5, 1)
        self.Calorebullicion = QtGui.QLabel()
        self.Calorebullicion.setFont(font8)
        self.Calorebullicion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Heat of vaporization, kJmol"))
        self.Calorebullicion.setAlignment(alignment)
        layoutInfo.addWidget(self.Calorebullicion, 5, 3)

        self.configuracion = QtGui.QLabel()
        self.configuracion.setFont(font7)
        self.configuracion.setAlignment(QtCore.Qt.AlignCenter)
        self.configuracion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Electronic configuration"))
        layoutInfo.addWidget(self.configuracion, 6, 1, 1, 3)

        self.Info2 = QtGui.QFrame()
        layout.addWidget(self.Info2, 0, 8, 3, 3)
        layoutInfo2 = QtGui.QGridLayout(self.Info2)
        layoutInfo2.setSpacing(1)
        layoutInfo2.setContentsMargins(2, 0, 2, 0)
        self.Info2.setFrameShape(QtGui.QFrame.StyledPanel)
        self.Info2.setFrameShadow(QtGui.QFrame.Raised)
        self.Info2.setAutoFillBackground(True)
        self.Info2.setPalette(palette)
        self.VolumenAtomico = QtGui.QLabel()
        self.VolumenAtomico.setFont(font8)
        self.VolumenAtomico.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Atomic volume")+u", cm³/mol")
        layoutInfo2.addWidget(self.VolumenAtomico, 1, 1)
        self.RadioAtomico = QtGui.QLabel()
        self.RadioAtomico.setFont(font8)
        self.RadioAtomico.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Atomic radius") + ", pm")
        layoutInfo2.addWidget(self.RadioAtomico, 2, 1)
        self.RadioCovalente = QtGui.QLabel()
        self.RadioCovalente.setFont(font8)
        self.RadioCovalente.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Covalent radius") + ", pm")
        layoutInfo2.addWidget(self.RadioCovalente, 3, 1)
        self.RadioVanderWaals = QtGui.QLabel()
        self.RadioVanderWaals.setFont(font8)
        self.RadioVanderWaals.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Van der Waals radius")+", pm")
        layoutInfo2.addWidget(self.RadioVanderWaals, 4, 1)
        self.RadioIonico = QtGui.QLabel()
        self.RadioIonico.setFont(font7)
        self.RadioIonico.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Ionic radii")+", pm")
        layoutInfo2.addWidget(self.RadioIonico, 5, 1, 1, 3)
        self.Electronegatividad = QtGui.QLabel()
        self.Electronegatividad.setFont(font8)
        self.Electronegatividad.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Electronegativity, Pauling scale"))
        self.Electronegatividad.setAlignment(QtCore.Qt.AlignRight |
                                             QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Electronegatividad, 1, 3)
        self.Cp = QtGui.QLabel()
        self.Cp.setFont(font8)
        self.Cp.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Specific heat capacitiy") + ", kJ/kgK")
        self.Cp.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Cp, 2, 3)
        self.Conductividad = QtGui.QLabel()
        self.Conductividad.setFont(font8)
        self.Conductividad.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Thermal conductivity") + ", W/mK")
        self.Conductividad.setAlignment(QtCore.Qt.AlignRight |
                                        QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Conductividad, 3, 3)
        self.Ionizacion = QtGui.QLabel()
        self.Ionizacion.setFont(font8)
        self.Ionizacion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "First ionization energy") + ", kJ/mol")
        self.Ionizacion.setAlignment(QtCore.Qt.AlignRight |
                                     QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.Ionizacion, 4, 3)

        self.EstadoOxidacion = QtGui.QLabel()
        self.EstadoOxidacion.setFont(font8)
        self.EstadoOxidacion.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Oxidation states"))
        self.EstadoOxidacion.setAlignment(QtCore.Qt.AlignCenter |
                                          QtCore.Qt.AlignVCenter)
        layoutInfo2.addWidget(self.EstadoOxidacion, 6, 1, 1, 3)

        self.actualizar(1)

    def actualizar(self, i):
        """Update botton info with data for current element"""
        elemento = Elemental.get_element(i)
        self.numero_atomico.setText(str(elemento.number))
        self.nombre.setText(elemento.name.get_string())
        self.simbolo.setText(elemento.symbol)
        self.textData(self.peso_atomico, elemento.atomic_mass)
        if elemento.density_solid.has_value:
            color = "#A52A2A"
            value = elemento.density_solid.get_string()
        elif elemento.density_liquid.has_value:
            color = "#0000FF"
            value = elemento.density_liquid.get_string()
        elif elemento.density_gas.has_value:
            color = "#0A640A"
            value = elemento.density_gas.get_string()
        else:
            color = "#888888"
            value = "N/A"
        self.densidad.setText("<font color={:}>{:}</font>".format(color, value))

        self.textData(self.Tfusion, elemento.melting_point, color="#0000C8")
        self.textData(self.Calorfusion, elemento.fusion_heat, color="#0000C8")
        self.textData(self.Tebullicion, elemento.boiling_point, color="#C80000")
        self.textData(self.Calorebullicion, elemento.vaporization_heat,
                      color="#C80000")
        self.textData(self.configuracion, elemento.configuration, string=True)

        self.textData(self.VolumenAtomico, elemento.atomic_volume, string=True)
        self.textData(self.RadioAtomico, elemento.atomic_radius)
        self.textData(self.RadioIonico, elemento.ionic_radii, string=True)
        self.textData(self.RadioCovalente, elemento.covalent_radius)
        self.textData(self.RadioVanderWaals, elemento.specific_heat)

        self.textData(self.Electronegatividad, elemento.van_der_waals_radius)
        self.textData(self.Cp, elemento.specific_heat)
        self.textData(self.Conductividad, elemento.thermal_conductivity)
        self.textData(self.Ionizacion, elemento.first_energy)
        self.textData(self.EstadoOxidacion, elemento.oxidation_states, string=True)

    def textData(self, widget, data, color="#000000", color_disabled="#888888",
                 string=False):
        if data.has_value:
            if string:
                widget.setText("<font color={:}>{:}</font>".format(
                    color, data.get_string()))
            else:
                widget.setText("<font color={:}>{:g}</font>".format(
                    color, data.value))
        else:
            widget.setText("<font color={:}>N/A</font>".format(color_disabled))


class boton(QtGui.QPushButton):
    """Button widget to define a element"""
    def __init__(self, i, parent=None):
        """Constructor, i parameter is the atomic number and useful ad id to
        tagged button"""
        super(boton, self).__init__(parent)
        self.setFixedSize(40, 35)
        self.setToolTip(QtGui.QApplication.translate(
            "pychemqt", "Click for view properties"))
        self.parent = parent
        self.setPalette(QtGui.QPalette(QtGui.QColor(
            Elemental.get_element(i).series.color.red*255,
            Elemental.get_element(i).series.color.green*255,
            Elemental.get_element(i).series.color.blue*255)))
        self.id = i
        self.setText(Elemental.get_element(i).symbol)
        self.clicked.connect(self.press)

    def enterEvent(self, event):
        """Enter event to change boton info with new element data"""
        self.parent.actualizar(self.id)
        event.accept()

    def press(self):
        """Press button show elementDialog with other element data"""
        dialog = ElementDialog(self.id)
        dialog.exec_()


class ElementDialog(QtGui.QDialog):
    """Dialog to show all element properties"""
    def __init__(self, i, parent=None):
        super(ElementDialog, self).__init__(parent)
        elemento = Elemental.get_element(i)
        self.setWindowTitle(QtGui.QApplication.translate(
            "pychemqt", "Properties of "+elemento.name.__str__().encode("utf-8")))
        self.setFixedSize(505, 432)
        self.verticallayout = QtGui.QVBoxLayout(self)
        self.tabWidget = QtGui.QTabWidget()
        self.verticallayout.addWidget(self.tabWidget)
        self.Aceptar = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.Aceptar.rejected.connect(self.reject)
        self.verticallayout.addWidget(self.Aceptar)

        self.tabGeneral = QtGui.QWidget()
        layoutGeneral = QtGui.QGridLayout(self.tabGeneral)
        layoutGeneral.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Name:")), 1, 1)
        layoutGeneral.addWidget(QtGui.QLabel(elemento.name.get_string()), 1, 2)
        layoutGeneral.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Serie:")), 2, 1)
        layoutGeneral.addWidget(QtGui.QLabel(elemento.series.get_string()), 2, 2)
        layoutGeneral.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Group")), 3, 1)
        layoutGeneral.addWidget(QtGui.QLabel(elemento.group.get_string()), 3, 2)
        layoutGeneral.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Period")), 4, 1)
        layoutGeneral.addWidget(QtGui.QLabel(elemento.period.get_string()), 4, 2)
        layoutGeneral.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Block")), 5, 1)
        layoutGeneral.addWidget(QtGui.QLabel(elemento.block.get_string()), 5, 2)

        layoutGeneral.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed), 6, 1)
        label_6 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "History"))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        label_6.setFont(font)
        layoutGeneral.addWidget(label_6, 7, 1, 1, 3)
        label_8 = QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Discovery") + ": " +
            elemento.discovery.get_string() + os.linesep +
            QtGui.QApplication.translate("pychemqt", "Discovered by ") +
            elemento.discovered_by.get_string() + os.linesep +
            QtGui.QApplication.translate("pychemqt", "Etymology") + ": " +
            elemento.etymology.get_string())
        label_8.setMargin(5)
        label_8.setWordWrap(True)
        layoutGeneral.addWidget(label_8, 8, 1, 1, 3)
        self.botoncito = QtGui.QLabel()
        self.botoncito.setStyleSheet(
            "background-color: %s;" % elemento.color.value.hex_spec)
        self.botoncito.setFrameShape(QtGui.QFrame.StyledPanel)
        self.botoncito.setFixedSize(60, 60)
        layoutGeneral.addWidget(self.botoncito, 1, 5, 3, 1)
        self.label_10 = QtGui.QLabel()
        self.label_10.setText(str(elemento.number))
        self.label_10.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(self.label_10, 1, 5)
        self.label_38 = QtGui.QLabel(elemento.symbol)
        font.setPointSize(12)
        self.label_38.setFont(font)
        self.label_38.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignBottom)
        layoutGeneral.addWidget(self.label_38, 2, 5)
        layoutGeneral.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding), 9, 3)
        self.tabWidget.addTab(self.tabGeneral, QtGui.QApplication.translate(
            "pychemqt", "General"))

        self.tab_Fisica = QtGui.QWidget()
        layoutFisica = QtGui.QGridLayout(self.tab_Fisica)
        layoutFisica.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Fase:")), 1, 1)
        layoutFisica.addWidget(QtGui.QLabel(elemento.phase.get_string() +
                                            u" a 0ºC"), 1, 2, 1, 1)
        if elemento.density_solid.has_value:
            layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
                "pychemqt", "Solid Density:")), 2, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_solid, "gcc", txt=u" @ 20ºC"), 2, 2)
        if elemento.density_liquid.has_value:
            layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
                "pychemqt", "Liquid Density:")), 3, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_liquid, "gcc", txt=" " +
                QtGui.QApplication.translate("pychemqt", "at melting point")), 3, 2)
        if elemento.density_gas.has_value:
            layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
                "pychemqt", "Gas Density:")), 4, 1)
            layoutFisica.addWidget(self.drawData(
                unidades.Density, elemento.density_gas, "gl", txt=u" @ 0ºC"), 4, 2)
        layoutFisica.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Appearance:")), 5, 1)
        label = QtGui.QLabel(elemento.appearance.get_string())
        label.setWordWrap(True)
        layoutFisica.addWidget(label, 5, 2)
        layoutFisica.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed), 6, 1, 1, 3)
        self.label_7 = QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Thermal properties"))
        self.label_7.setFont(font)
        layoutFisica.addWidget(self.label_7, 7, 1)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Melting point:")), 8, 1)
        self.punto_fusion = self.drawData(
            unidades.Temperature, elemento.melting_point)
        layoutFisica.addWidget(self.punto_fusion, 8, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Boiling point:")), 9, 1)
        self.punto_ebullicion = self.drawData(
            unidades.Temperature, elemento.boiling_point)
        layoutFisica.addWidget(self.punto_ebullicion, 9, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Heat of fusion:")), 10, 1)
        self.calor_fusion = self.drawData(
            unidades.MolarEnthalpy, elemento.fusion_heat, "kJkmol")
        layoutFisica.addWidget(self.calor_fusion, 10, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Heat of vaporization:")), 11, 1)
        self.calor_vaporizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.vaporization_heat, "kJkmol")
        layoutFisica.addWidget(self.calor_vaporizacion, 11, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Specific heat capacity:")), 12, 1)
        self.capacidad_calorifica = self.drawData(
            unidades.SpecificHeat, elemento.specific_heat, "JgK")
        layoutFisica.addWidget(self.capacidad_calorifica, 12, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Thermal conductivity:")), 13, 1)
        self.conductividad_termica = self.drawData(
            unidades.ThermalConductivity, elemento.thermal_conductivity,
            txt=" @ 300K")
        layoutFisica.addWidget(self.conductividad_termica, 13, 2)
        layoutFisica.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Debye Temperature:")), 14, 1)
        self.temperatura_debye = self.drawData(
            unidades.Temperature, elemento.debye_temperature)
        layoutFisica.addWidget(self.temperatura_debye, 14, 2)
        layoutFisica.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            15, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Fisica, QtGui.QApplication.translate(
            "pychemqt", "Physical properties"))

        self.tab_Atomo = QtGui.QWidget()
        layout_Atomo = QtGui.QGridLayout(self.tab_Atomo)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Atomic mass:")), 1, 1)
        if elemento.atomic_mass.has_value:
            self.masa_atomica = QtGui.QLabel(
                elemento.atomic_mass.get_string()+" g/mol")
        else:
            self.masa_atomica = QtGui.QLabel(elemento.atomic_mass.get_string())
        layout_Atomo.addWidget(self.masa_atomica, 1, 2)
        layout_Atomo.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Atomic Volume:")), 2, 1)
        self.volumen_atomico = self.drawData(unidades.MolarVolume,
                                             elemento.atomic_volume)
        layout_Atomo.addWidget(self.volumen_atomico, 2, 2)
        layout_Atomo.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Atomic radius:")), 3, 1)
        if elemento.atomic_radius.has_value:
            self.radio_atomico = QtGui.QLabel(
                elemento.atomic_radius.get_string() + " pm")
        else:
            self.radio_atomico = QtGui.QLabel(elemento.atomic_radius.get_string())
        layout_Atomo.addWidget(self.radio_atomico, 3, 2)
        layout_Atomo.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Covalent radius:")), 4, 1)
        if elemento.covalent_radius.has_value:
            self.radio_covalente = QtGui.QLabel(
                elemento.covalent_radius.get_string() + " pm")
        else:
            self.radio_covalente = QtGui.QLabel(
                elemento.covalent_radius.get_string())
        layout_Atomo.addWidget(self.radio_covalente, 4, 2)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Van der Waals radius:")), 5, 1)
        if elemento.van_der_waals_radius.has_value:
            self.radio_waals = QtGui.QLabel(
                elemento.van_der_waals_radius.get_string() + " pm")
        else:
            self.radio_waals = QtGui.QLabel(
                elemento.van_der_waals_radius.get_string())
        layout_Atomo.addWidget(self.radio_waals, 5, 2)
        layout_Atomo.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Ionic radii:")), 6, 1)
        if elemento.ionic_radii.has_value:
            self.radio_ionico = QtGui.QLabel(
                elemento.ionic_radii.get_string() + " pm")
        else:
            self.radio_ionico = QtGui.QLabel(elemento.ionic_radii.get_string())
        layout_Atomo.addWidget(self.radio_ionico, 6, 2)
        layout_Atomo.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed),
            7, 1, 1, 3)
        self.label_9 = QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Electronic properties"))
        self.label_9.setFont(font)
        layout_Atomo.addWidget(self.label_9, 8, 1)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Electronic configuration:")), 9, 1)
        layout_Atomo.addWidget(QtGui.QLabel(
            elemento.configuration.get_string()), 9, 2)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Oxidation states:")), 10, 1)
        layout_Atomo.addWidget(QtGui.QLabel(
            elemento.oxidation_states.get_string()), 10, 2)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Electronegativity:")), 11, 1)
        layout_Atomo.addWidget(QtGui.QLabel(
            elemento.electronegativity.get_string()), 11, 2)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "Electron affinity:")), 12, 1)
        self.afinidad_electronica = self.drawData(
            unidades.MolarEnthalpy, elemento.electron_affinity, "kJkmol")
        layout_Atomo.addWidget(self.afinidad_electronica, 12, 2)
        layout_Atomo.addWidget(QtGui.QLabel(QtGui.QApplication.translate(
            "pychemqt", "1st ionization energy:")), 13, 1)
        self.energia_ionizacion = self.drawData(
            unidades.MolarEnthalpy, elemento.first_energy, "kJkmol")
        layout_Atomo.addWidget(self.energia_ionizacion, 13, 2)
        layout_Atomo.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            14, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Atomo, QtGui.QApplication.translate(
            "pychemqt", "Atomic properties"))

        self.tab_Cristalografia = QtGui.QWidget()
        layout_Cristalografia = QtGui.QGridLayout(self.tab_Cristalografia)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Lattice type:")), 1, 1)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            elemento.lattice_type.get_string()), 1, 2)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Space group:")), 2, 1)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            elemento.space_group.get_string()), 2, 2)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Lattice edge lengths:")), 3, 1)
        self.lados = QtGui.QLabel()
        if elemento.lattice_edges.has_value:
            self.lados.setText(
                elemento.lattice_edges.values[0]*10 + u" Å" +
                elemento.lattice_edges.values[1]*10 + u" Å" +
                elemento.lattice_edges.values[2]*10 + u" Å")
        else:
            self.lados.setText(elemento.lattice_edges.get_string())
        layout_Cristalografia.addWidget(self.lados, 3, 2)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Lattice angles:")), 4, 1)
        self.angulos = QtGui.QLabel()
        if elemento.lattice_angles.has_value:
            self.angulos.setText(
                str(elemento.lattice_angles.values[0]) + u"º, " +
                str(elemento.lattice_angles.values[1]) + u"º, " +
                str(elemento.lattice_angles.values[2]) + u"º")
        else:
            self.angulos.setText(elemento.lattice_angles.get_string())
        layout_Cristalografia.addWidget(self.angulos, 4, 2)
        layout_Cristalografia.addWidget(QtGui.QLabel(
            QtGui.QApplication.translate("pychemqt", "Lattice unit volume:")), 5, 1)
        self.volumen_celda = QtGui.QLabel()
        if elemento.lattice_angles.has_value:
            self.volumen_celda.setText(
                elemento.lattice_volume.__str__().encode("utf-8") +
                "nm<sup>3</sup>")
        else:
            self.volumen_celda.setText(elemento.lattice_volume.get_string())
        layout_Cristalografia.addWidget(self.volumen_celda, 5, 2)
        layout_Cristalografia.addItem(QtGui.QSpacerItem(
            20, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding),
            6, 1, 1, 3)
        self.tabWidget.addTab(self.tab_Cristalografia,
            QtGui.QApplication.translate("pychemqt", "Crystallographic"))

    def drawData(self, unit, data, unidad="", txt=""):
        """Return a widget with the data inprint
            unit: unidad subclass
            data: value to set
            unidad: unit of value to show
            txt: opcional txt to show for unit"""
        if data.has_value and unidad:
            value = unit(data.value, unidad)
            widget = Entrada_con_unidades(
                unit, readOnly=True, value=value, textounidad=txt)
        elif data.has_value:
            widget = Entrada_con_unidades(
                unit, readOnly=True, value=data.value, textounidad=txt)
        else:
            widget = QtGui.QLabel(data.get_string())
        return widget


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = qtelemental()
    Form.show()
    sys.exit(app.exec_())
