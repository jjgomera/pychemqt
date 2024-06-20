#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Library to configure MEOS tool
#
#   - ColorMapCombo: Custom QComboBox to choose a matplotlib colormap
#   - Isolinea: widget to configure isolines for mEoS
#   - Widget: mEoS parameter configuration dialog
#   - Dialog: Dialog tool for standalone use
###############################################################################


import os

from matplotlib import colormaps
from numpy import arange
from tools.qt import QtCore, QtGui, QtWidgets

from lib import unidades
from lib.utilities import representacion
from UI.widgets import Entrada_con_unidades, LineConfig


class ColorMapCombo(QtWidgets.QComboBox):
    """Custom QComboBox to choose a matplotlib colormap"""

    COLORMAP = [
        # Perceptually Uniform Sequential
        'viridis', 'plasma', 'inferno', 'magma', 'cividis',
        # Sequential
        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
        # Sequential (2)
        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'hot', 'afmhot', 'gist_heat', 'copper',
        # Diverging
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
        # Miscellaneous
        'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
        'gist_rainbow', 'rainbow', 'jet', 'turbo']

    def __init__(self, *args, **kwargs):
        """Autofill of widget with accepted values"""
        super().__init__(*args, **kwargs)

        for cm in self.COLORMAP:
            maps = colormaps.get(cm)

            gradient = QtGui.QLinearGradient(0, 0, 50, 0)
            for value in arange(0, 1, 0.01):
                r, g, b, a = maps(value)
                gradient.setColorAt(
                    value, QtGui.QColor(int(r*255), int(g*255), int(b*255)))

            pix = QtGui.QPixmap(50, 50)
            pix.fill(QtGui.QColorConstants.White)

            painter = QtGui.QPainter()
            painter.begin(pix)
            brush = QtGui.QBrush(gradient)
            painter.setBrush(brush)
            painter.drawRect(0, 0, 50, 50)
            icon = QtGui.QIcon(pix)
            painter.end()
            self.addItem(icon, cm)

    def paintEvent(self, event):
        """Paint the widget with a image with the selected colormap"""
        # Paint default element
        painter = QtWidgets.QStylePainter(self)
        opt = QtWidgets.QStyleOptionComboBox()
        self.initStyleOption(opt)
        painter.drawComplexControl(
            QtWidgets.QStyle.ComplexControl.CC_ComboBox, opt)

        # Apply selected brush
        maps = colormaps.get(self.COLORMAP[self.currentIndex()])
        gradient = QtGui.QLinearGradient(0, 0, event.rect().right()-24, 0)
        for value in arange(0, 1, 0.01):
            r, g, b, a = maps(value)
            gradient.setColorAt(
                value, QtGui.QColor(int(r*255), int(g*255), int(b*255)))

        brush = QtGui.QBrush(gradient)
        painter.setBrush(brush)

        # Don't paint the arrow region in QComboBox
        rect = event.rect()
        rect.setRight(rect.right()-24)
        painter.drawRect(rect)


class Isolinea(QtWidgets.QDialog):
    """Widget for isoline configuration for mEoS plot tools"""
    def __init__(self, unit, ConfSection, config, section="MEOS", parent=None):
        """Constructor
            unit: subclass of unidad to define the isoline type
            ConfSection: title of isoline
            config: config of pychemqt project
        """
        super().__init__(parent)
        self.ConfSection = ConfSection
        self.magnitud = unit.__name__
        self.unidad = unit
        self.section = section
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel(self.tr("Start")), 1, 1)
        self.inicio = Entrada_con_unidades(unit)
        layout.addWidget(self.inicio, 1, 2, 1, 3)
        layout.addWidget(QtWidgets.QLabel(self.tr("Fin")), 2, 1)
        self.fin = Entrada_con_unidades(unit)
        layout.addWidget(self.fin, 2, 2, 1, 3)
        layout.addWidget(QtWidgets.QLabel(self.tr("Intervalo")), 3, 1)
        if unit.__name__ == "Temperature":
            self.intervalo = Entrada_con_unidades(unidades.DeltaT)
        elif unit.__name__ == "Pressure":
            self.intervalo = Entrada_con_unidades(unidades.DeltaP)
        else:
            self.intervalo = Entrada_con_unidades(unit)
        layout.addWidget(self.intervalo, 3, 2, 1, 3)
        self.Personalizar = QtWidgets.QCheckBox(self.tr("Customize"))
        layout.addWidget(self.Personalizar, 4, 1, 1, 4)
        self.Lista = QtWidgets.QLineEdit()
        layout.addWidget(self.Lista, 5, 1, 1, 4)
        self.Personalizar.toggled.connect(self.inicio.setDisabled)
        self.Personalizar.toggled.connect(self.fin.setDisabled)
        self.Personalizar.toggled.connect(self.intervalo.setDisabled)
        self.Personalizar.toggled.connect(self.Lista.setEnabled)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 6, 1, 1, 4)
        if unit.__name__ != "float" and section != "Psychr":
            self.Critica = QtWidgets.QCheckBox(
                self.tr("Include critic point line"))
            layout.addWidget(self.Critica, 7, 1, 1, 4)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 8, 1, 1, 4)

        self.lineconfig = LineConfig(ConfSection, self.tr("Line Style"))
        layout.addWidget(self.lineconfig, 9, 1, 1, 4)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 10, 1)
        self.label = QtWidgets.QCheckBox(self.tr("Label"))
        layout.addWidget(self.label, 11, 1)
        self.variable = QtWidgets.QCheckBox(self.tr("Variable in Label"))
        layout.addWidget(self.variable, 12, 1, 1, 4)
        self.unit = QtWidgets.QCheckBox(self.tr("Units in Label"))
        layout.addWidget(self.unit, 13, 1, 1, 4)
        layout.addWidget(QtWidgets.QLabel(self.tr("Position")), 14, 1)
        self.label5 = Entrada_con_unidades(
            int, value=0, width=25, frame=False, readOnly=True)
        self.label5.setFixedWidth(30)
        layout.addWidget(self.label5, 14, 2)
        self.Posicion = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.Posicion.valueChanged.connect(self.label5.setValue)
        layout.addWidget(self.Posicion, 14, 3, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 4)

        if config.has_section(section):
            self.inicio.setValue(config.getfloat(section, ConfSection+'Start'))
            self.fin.setValue(config.getfloat(section, ConfSection+'End'))
            self.intervalo.setValue(
                config.getfloat(section, ConfSection+'Step'))
            self.Personalizar.setChecked(
                config.getboolean(section, ConfSection+'Custom'))
            if config.get(section, ConfSection+'List') != "":
                T = []
                for i in config.get(section, ConfSection+'List').split(','):
                    if unit.__name__ == "float":
                        T.append(representacion(float(i)))
                    else:
                        T.append(representacion(unit(float(i)).config()))
                self.Lista.setText(",".join(T))
            if unit.__name__ != "float" and section != "Psychr":
                self.Critica.setChecked(
                    config.getboolean(section, ConfSection+'Critic'))
            self.inicio.setDisabled(self.Personalizar.isChecked())
            self.fin.setDisabled(self.Personalizar.isChecked())
            self.intervalo.setDisabled(self.Personalizar.isChecked())
            self.Lista.setEnabled(self.Personalizar.isChecked())
            self.label.setChecked(
                config.getboolean(section, ConfSection+'Label'))
            self.variable.setChecked(
                config.getboolean(section, ConfSection+'Variable'))
            self.unit.setChecked(
                config.getboolean(section, ConfSection+'Units'))
            self.Posicion.setValue(
                config.getint(section, ConfSection+'Position'))
            self.lineconfig.setConfig(config, section)

    def value(self, config):
        """Update ConfigParser instance with the config"""
        config.set(self.section, self.ConfSection+"Start",
                   str(self.inicio.value))
        config.set(self.section, self.ConfSection+"End", str(self.fin.value))
        config.set(self.section, self.ConfSection+"Step",
                   str(self.intervalo.value))
        config.set(self.section, self.ConfSection+"Custom",
                   str(self.Personalizar.isChecked()))
        T = []
        if self.Lista.text():
            T1 = self.Lista.text().split(',')
            for i in T1:
                if self.unidad.__name__ == "float":
                    T.append(str(float(i)))
                else:
                    T.append(str(self.unidad(float(i), "conf")))
        config.set(self.section, self.ConfSection+"List", ", ".join(T))
        if self.unidad.__name__ != "float" and self.section != "Psychr":
            config.set(self.section, self.ConfSection+"Critic",
                       str(self.Critica.isChecked()))
        config = self.lineconfig.value(config, self.section)

        config.set(self.section, self.ConfSection+"Label",
                   str(self.label.isChecked()))
        config.set(self.section, self.ConfSection+"Variable",
                   str(self.variable.isChecked()))
        config.set(self.section, self.ConfSection+"Units",
                   str(self.unit.isChecked()))
        config.set(self.section, self.ConfSection+"Position",
                   str(self.Posicion.value()))
        return config


class Widget(QtWidgets.QDialog):
    """Config mEoS parameter dialog"""
    lineas = [
        ("Isotherm", unidades.Temperature, QtWidgets.QApplication.translate("UI_Tables", "Isotherm")),
        ("Isobar", unidades.Pressure, QtWidgets.QApplication.translate("UI_Tables", "Isobar")),
        ("Isoenthalpic", unidades.Enthalpy, QtWidgets.QApplication.translate("UI_Tables", "Isoenthalpic")),
        ("Isoentropic", unidades.SpecificHeat, QtWidgets.QApplication.translate("UI_Tables", "Isoentropic")),
        ("Isochor", unidades.SpecificVolume, QtWidgets.QApplication.translate("UI_Tables", "Isochor")),
        ("Isoquality", float, QtWidgets.QApplication.translate("UI_Tables", "Isoquality"))]

    def __init__(self, config, parent=None):
        """constructor, config optional parameter to input project config"""
        super().__init__(parent)

        lyt = QtWidgets.QGridLayout(self)
        lyt.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setFrameStyle(QtWidgets.QFrame.Shape.NoFrame)
        lyt.addWidget(scroll)

        dlg = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout(dlg)

        self.coolProp = QtWidgets.QCheckBox(
            self.tr("Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp, 3, 1, 1, 2)
        self.refprop = QtWidgets.QCheckBox(
            self.tr("Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop, 4, 1, 1, 2)
        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 4, 1)
        self.lineconfig = LineConfig(
            "saturation", self.tr("Saturation Line Style"))
        layout.addWidget(self.lineconfig, 5, 1, 1, 2)
        group = QtWidgets.QGroupBox(self.tr("Isolines"))
        layout.addWidget(group, 6, 1, 1, 2)
        layoutgroup = QtWidgets.QGridLayout(group)
        self.comboIsolineas = QtWidgets.QComboBox()
        layoutgroup.addWidget(self.comboIsolineas, 1, 1)
        self.Isolineas = QtWidgets.QStackedWidget()
        self.comboIsolineas.currentIndexChanged.connect(
            self.Isolineas.setCurrentIndex)
        layoutgroup.addWidget(self.Isolineas, 2, 1)
        for nombre, unidad, text in self.lineas:
            self.comboIsolineas.addItem(text)
            self.Isolineas.addWidget(Isolinea(unidad, nombre, config))
        layout.addWidget(QtWidgets.QLabel(self.tr("Plot Definition")), 7, 1)
        quality = [self.tr("Very Low"),
                   self.tr("Low"),
                   self.tr("Medium"),
                   self.tr("High"),
                   self.tr("Ultra High")]
        self.definition = QtWidgets.QComboBox()
        for q in quality:
            self.definition.addItem(q)
        layout.addWidget(self.definition, 7, 2)
        self.grid = QtWidgets.QCheckBox(self.tr("Draw grid"))
        layout.addWidget(self.grid, 9, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed), 10, 2)
        group2 = QtWidgets.QGroupBox(self.tr("3D plot"))
        layout.addWidget(group2, 11, 1, 1, 2)
        layoutgroup2 = QtWidgets.QGridLayout(group2)
        self.checkMesh = QtWidgets.QCheckBox(self.tr("Draw 3D mesh"))
        layoutgroup2.addWidget(self.checkMesh, 1, 1, 1, 2)
        layoutgroup2.addWidget(QtWidgets.QLabel(self.tr("Mesh type")), 2, 1)
        self.typeMesh = QtWidgets.QComboBox()
        self.typeMesh.addItem("Surface")
        self.typeMesh.addItem("Wireframe")
        layoutgroup2.addWidget(self.typeMesh, 2, 2)

        self.widgetSurface = QtWidgets.QWidget()
        layoutgroup2.addWidget(self.widgetSurface, 3, 1, 1, 2)

        lytSurface = QtWidgets.QGridLayout(self.widgetSurface)
        lytSurface.addWidget(QtWidgets.QLabel(self.tr("Colormap")), 1, 1)
        self.colormapMesh = ColorMapCombo()
        lytSurface.addWidget(self.colormapMesh, 1, 2)
        lytSurface.addWidget(QtWidgets.QLabel(self.tr("Alpha")), 2, 1)
        self.alphaSurface = QtWidgets.QSpinBox()
        self.alphaSurface.setRange(0, 255)
        lytSurface.addWidget(self.alphaSurface, 2, 2)

        self.wireframeConfig = LineConfig(
            "3D", self.tr("Wireframe line configuration"))
        self.wireframeConfig.Marker.setEnabled(False)
        layoutgroup2.addWidget(self.wireframeConfig, 4, 1, 1, 2)

        layout.addItem(QtWidgets.QSpacerItem(
            10, 10, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding), 15, 2)

        scroll.setWidget(dlg)

        if os.environ["CoolProp"]:
            self.coolProp.setEnabled(True)
        if os.environ["refprop"]:
            self.refprop.setEnabled(True)

        if config.has_section("MEOS"):
            self.coolProp.setChecked(config.getboolean("MEOS", 'coolprop'))
            self.refprop.setChecked(config.getboolean("MEOS", 'refprop'))
            self.lineconfig.setConfig(config)
            self.definition.setCurrentIndex(
                config.getint("MEOS", 'definition'))
            self.grid.setChecked(config.getboolean("MEOS", 'grid'))

            self.checkMesh.setChecked(config.getboolean("MEOS", '3Dmesh'))
            self.typeMesh.setCurrentIndex(config.getint("MEOS", '3Dtype'))
            self.colormapMesh.setCurrentText(config.get("MEOS", '3Dcolormap'))
            self.alphaSurface.setValue(config.getint("MEOS", '3Dalphasurface'))
            self.wireframeConfig.setConfig(config)

    def value(self, config):
        """Return value for main dialog"""
        if not config.has_section("MEOS"):
            config.add_section("MEOS")

        config.set("MEOS", "coolprop", str(self.coolProp.isChecked()))
        config.set("MEOS", "refprop", str(self.refprop.isChecked()))
        config = self.lineconfig.value(config)

        for indice in range(self.Isolineas.count()):
            config = self.Isolineas.widget(indice).value(config)

        config.set("MEOS", "definition", str(self.definition.currentIndex()))
        config.set("MEOS", "grid", str(self.grid.isChecked()))

        config.set("MEOS", "3Dmesh", str(self.checkMesh.isChecked()))
        config.set("MEOS", "3Dtype", str(self.typeMesh.currentIndex()))
        config.set("MEOS", "3Dcolormap", self.colormapMesh.currentText())
        config.set("MEOS", "3Dalphasurface", str(self.alphaSurface.value()))
        config = self.wireframeConfig.value(config)
        return config


class Dialog(QtWidgets.QDialog):
    """Dialog to config thermal method calculations"""
    def __init__(self, config=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(
            self.tr("Multiparameter equation of state configuration"))
        layout = QtWidgets.QVBoxLayout(self)
        self.widget = Widget(config)
        layout.addWidget(self.widget)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
            | QtWidgets.QDialogButtonBox.StandardButton.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)

    def value(self, config):
        """Function result for wizard"""
        config = self.widget.value(config)
        return config


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    app = QtWidgets.QApplication(sys.argv)

    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    conf = ConfigParser()
    conf.read(conf_dir+"pychemqtrc")

    Dialog = Dialog(conf)
    Dialog.show()
    sys.exit(app.exec())
