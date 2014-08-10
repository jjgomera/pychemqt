#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to show/configure pychemqt general preferences
#
#   - Preferences: Preferences main dialog
###############################################################################

import os
import sys
from functools import partial
from math import pi
from ConfigParser import ConfigParser

from PyQt4 import QtCore, QtGui

from UI.widgets import Entrada_con_unidades, ColorSelector, LineConfig, PathConfig
from UI.delegate import CheckEditor, comboLine
from tools import UI_confResolution
from lib import unidades, corriente
from equipment import equipments
from lib.firstrun import calculator, editor, shell
from lib.config import representacion


def format2txt(formato):
    """Function to convert dict format config in a string equivalent"""
    if formato["signo"]:
        txt="+"
    else:
        txt=""
    if formato["format"]==0:
        txt+="{total}.{decimales} fixed".format(**formato)
    elif formato["format"]==1:
        txt+="{decimales} sign".format(**formato)
    elif formato["format"]==2:
        txt+="{decimales} exp".format(**formato)
    if formato.get("exp", False):
        txt+=" ({tol} exp)".format(**formato)
    return txt


class ConfLine(QtGui.QWidget):
    """Composite widget con las herramientas de definición del formato de líneas"""
    def __init__(self, pen=None, parent=None):
        super(ConfLine, self).__init__(parent)
        lyt = QtGui.QGridLayout(self)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Color:")),0,1)
        self.ColorButtonLine = ColorSelector()
        lyt.addWidget(self.ColorButtonLine,0,2,1,4)

        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Width:")),1,1)
        self.width=Entrada_con_unidades(float, width=50, decimales=1, spinbox=True, step=0.1, textounidad="px")
        lyt.addWidget(self.width,1,2,1,4)

        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Join:")),2,1)
        toolJoinMitter=QtGui.QToolButton()
        toolJoinMitter.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-join-miter.png")))
        toolJoinMitter.setIconSize(QtCore.QSize(24, 24))
        toolJoinMitter.setCheckable(True)
        toolJoinMitter.setToolTip(QtGui.QApplication.translate("pychemqt", "Join bevel: The triangular notch between the two lines is filled"))
        lyt.addWidget(toolJoinMitter,2,2)
        toolJoinBevel=QtGui.QToolButton()
        toolJoinBevel.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-join-bevel.png")))
        toolJoinBevel.setIconSize(QtCore.QSize(24, 24))
        toolJoinBevel.setCheckable(True)
        toolJoinBevel.setToolTip(QtGui.QApplication.translate("pychemqt", "Join bevel: The triangular notch between the two lines is filled"))
        lyt.addWidget(toolJoinBevel,2,3)
        toolJoinRound=QtGui.QToolButton()
        toolJoinRound.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-join-round.png")))
        toolJoinRound.setIconSize(QtCore.QSize(24, 24))
        toolJoinRound.setCheckable(True)
        toolJoinRound.setToolTip(QtGui.QApplication.translate("pychemqt", "Join round: A circular arc between the two lines is filled"))
        lyt.addWidget(toolJoinRound,2,4)

        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mitter Limit:")),3,1)
        self.mitterLimit=Entrada_con_unidades(float, width=50, decimales=1, spinbox=True, step=0.1)
        lyt.addWidget(self.mitterLimit,3,2,1,4)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1,1,4)
        self.groupJoint=QtGui.QButtonGroup()
        self.groupJoint.addButton(toolJoinMitter)
        self.groupJoint.addButton(toolJoinBevel)
        self.groupJoint.addButton(toolJoinRound)
        self.groupJoint.buttonClicked["int"].connect(self.mitterlimitEnabled)

        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Cap:")),5,1)
        toolCapFlat=QtGui.QToolButton()
        toolCapFlat.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-cap-butt.png")))
        toolCapFlat.setIconSize(QtCore.QSize(24, 24))
        toolCapFlat.setCheckable(True)
        toolCapFlat.setToolTip(QtGui.QApplication.translate("pychemqt", "Flat Cap: A square line end that does not cover the end point of the line"))
        lyt.addWidget(toolCapFlat,5,2)
        toolCapRound=QtGui.QToolButton()
        toolCapRound.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-cap-round.png")))
        toolCapRound.setIconSize(QtCore.QSize(24, 24))
        toolCapRound.setCheckable(True)
        toolCapRound.setToolTip(QtGui.QApplication.translate("pychemqt", "Round Cap: A rounded line end"))
        lyt.addWidget(toolCapRound,5,3)
        toolCapSquare=QtGui.QToolButton()
        toolCapSquare.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/stroke-cap-square.png")))
        toolCapSquare.setIconSize(QtCore.QSize(24, 24))
        toolCapSquare.setCheckable(True)
        toolCapSquare.setToolTip(QtGui.QApplication.translate("pychemqt", "Square Cap: A square line end that covers the end point and extends beyond it by half the line width"))
        lyt.addWidget(toolCapSquare,5,4)
        self.groupCap=QtGui.QButtonGroup()
        self.groupCap.addButton(toolCapFlat)
        self.groupCap.addButton(toolCapRound)
        self.groupCap.addButton(toolCapSquare)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,4)

        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Dashes:")),7,1)
        self.guion=comboLine()
        lyt.addWidget(self.guion,7,2,1,3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Dash offset:")),8,1)
        self.dashOffset=Entrada_con_unidades(float, width=50, decimales=1, spinbox=True, step=0.1)
        lyt.addWidget(self.dashOffset,8,2,1,4)
        lyt.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),9,5,1,1)

        if pen:
            self.ColorButtonLine.setColor(pen.color().name())
            self.groupJoint.button(([QtCore.Qt.MiterJoin, QtCore.Qt.BevelJoin, QtCore.Qt.RoundJoin].index(pen.joinStyle())+2)*-1).setChecked(True)
            self.mitterLimit.setValue(pen.miterLimit())
            self.groupCap.button(([QtCore.Qt.FlatCap, QtCore.Qt.RoundCap, QtCore.Qt.SquareCap].index(pen.capStyle())+2)*-1).setChecked(True)
            self.guion.setCurrentIndex([QtCore.Qt.SolidLine, QtCore.Qt.DashLine, QtCore.Qt.DotLine, QtCore.Qt.DashDotLine, QtCore.Qt.DashDotDotLine].index(pen.style()))
            self.dashOffset.setValue(pen.dashOffset())
            self.width.setValue(pen.widthF())


    def mitterlimitEnabled(self, id):
        self.mitterLimit.setEnabled(id==-2)

    def pen(self):
        """Return a QPen with the live configuration"""
        pen=QtGui.QPen(QtGui.QColor(self.ColorButtonLine.color.name()))
        pen.setWidthF(self.width.value)
        pen.setJoinStyle([QtCore.Qt.MiterJoin, QtCore.Qt.BevelJoin, QtCore.Qt.RoundJoin][abs(self.groupJoint.checkedId())-2])
        pen.setMiterLimit(self.mitterLimit.value)
        pen.setCapStyle([QtCore.Qt.FlatCap, QtCore.Qt.RoundCap, QtCore.Qt.SquareCap][abs(self.groupCap.checkedId())-2])
        pen.setStyle([QtCore.Qt.SolidLine, QtCore.Qt.DashLine, QtCore.Qt.DotLine, QtCore.Qt.DashDotLine, QtCore.Qt.DashDotDotLine][self.guion.currentIndex()])
        pen.setDashOffset(self.dashOffset.value)
        return pen



class ConfGeneral(QtGui.QDialog):
    """Clase con el widget de las características generales"""
    def __init__(self, config=None, parent=None):
        super(ConfGeneral, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        layout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed), 0, 1, 1, 4)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Highlight color:")), 1, 1, 1, 1)
        self.ColorButtonResaltado = ColorSelector()
        layout.addWidget(self.ColorButtonResaltado, 1, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Readonly color:")), 2, 1, 1, 1)
        self.ColorButtonReadOnly = ColorSelector()
        layout.addWidget(self.ColorButtonReadOnly, 2, 2, 1, 1)
        layout.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1)
        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Recent Files"))
        layout.addWidget(group,4,1,1,4)
        lyt=QtGui.QHBoxLayout(group)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Number of recent files:")))
        self.recentFiles=QtGui.QSpinBox()
        self.recentFiles.setRange(1, 20)
        lyt.addWidget(self.recentFiles)
        lyt.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed))
        self.loadLastProject=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Load last session project"))
        layout.addWidget(self.loadLastProject,5,1)
        self.showTrayIcon=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show tray icon"))
        layout.addWidget(self.showTrayIcon,6,1)

        layout.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("General"):
            self.ColorButtonResaltado.setColor(config.get("General", 'Color_Resaltado'))
            self.ColorButtonReadOnly.setColor(config.get("General", 'Color_ReadOnly'))
            self.recentFiles.setValue(config.getint("General", 'Recent_Files'))
            self.loadLastProject.setChecked(config.getboolean("General", 'Load_Last_Project'))
            self.showTrayIcon.setChecked(config.getboolean("General", 'Tray'))


    def value(self, config):
        if not config.has_section("General"):
            config.add_section("General")
        config.set("General", "Color_Resaltado", self.ColorButtonResaltado.color.name())
        config.set("General", "Color_ReadOnly", self.ColorButtonReadOnly.color.name())
        config.set("General", "Recent_Files", str(self.recentFiles.value()))
        config.set("General", "Load_Last_Project", str(self.loadLastProject.isChecked()))
        config.set("General", "Tray", str(self.showTrayIcon.isChecked()))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("General")
        config.set("General", "Color_Resaltado", "#ffff00")
        config.set("General", "Color_ReadOnly", "#eaeaea")
        config.set("General", "Recent_Files", "10")
        config.set("General", "Load_Last_Project", "True")
        config.set("General", "Tray", "False")
        return config


class ConfPFD(QtGui.QDialog):
    """Clase con el widget de las características del diagrama de flujo"""
    def __init__(self, config=None, parent=None):
        super(ConfPFD, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        layout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed), 0, 1, 1, 4)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Input color:")),1,1)
        self.ColorButtonEntrada = ColorSelector()
        layout.addWidget(self.ColorButtonEntrada,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output color:")),2,1)
        self.ColorButtonSalida = ColorSelector()
        layout.addWidget(self.ColorButtonSalida,2,2)

        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Line format"))
        layout.addWidget(group,3,1,1,3)
        lyt = QtGui.QHBoxLayout(group)
        self.lineFormat=ConfLine()
        lyt.addWidget(self.lineFormat)

        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "PFD resolution"))
        layout.addWidget(group,4,1,1,3)
        lyt = QtGui.QHBoxLayout(group)
        self.resolution=UI_confResolution.UI_confResolution_widget(config)
        lyt.addWidget(self.resolution)

        layout.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding), 14, 1, 1, 4)

        if config and config.has_section("PFD"):
            self.ColorButtonEntrada.setColor(config.get("PFD", 'Color_Entrada'))
            self.ColorButtonSalida.setColor(config.get("PFD", 'Color_Salida'))
            self.lineFormat.ColorButtonLine.setColor(config.get("PFD", 'Color_Stream'))
            self.lineFormat.groupJoint.button((config.getint("PFD", 'Union')+2)*-1).setChecked(True)
            self.lineFormat.mitterLimit.setValue(config.getfloat("PFD", 'Miter_limit'))
            self.lineFormat.groupCap.button((config.getint("PFD", 'Punta')+2)*-1).setChecked(True)
            self.lineFormat.guion.setCurrentIndex(config.getint("PFD", 'Guion'))
            self.lineFormat.dashOffset.setValue(config.getfloat("PFD", 'Dash_offset'))
            self.lineFormat.width.setValue(config.getfloat("PFD", 'Width'))


    def value(self, config):
        if not config.has_section("PFD"):
            config.add_section("PFD")
        config=self.resolution.value(config)
        config.set("PFD", "Color_Entrada", self.ColorButtonEntrada.color.name())
        config.set("PFD", "Color_Salida", self.ColorButtonSalida.color.name())
        config.set("PFD", "Color_Stream", self.lineFormat.ColorButtonLine.color.name())
        config.set("PFD", "Width", str(self.lineFormat.width.value))
        config.set("PFD", "Union", str(abs(self.lineFormat.groupJoint.checkedId())-2))
        config.set("PFD", "Miter_limit", str(self.lineFormat.mitterLimit.value))
        config.set("PFD", "Punta", str(abs(self.lineFormat.groupCap.checkedId())-2))
        config.set("PFD", "Guion", str(self.lineFormat.guion.currentIndex()))
        config.set("PFD", "Dash_offset", str(self.lineFormat.dashOffset.value))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("PFD")
        config.set("PFD", "x", "800")
        config.set("PFD", "y", "600")
        config.set("PFD", "Color_Entrada", "#c80000")
        config.set("PFD", "Color_Salida", "#0000c8")
        config.set("PFD", "Color_Stream", "#000000")
        config.set("PFD", "Width", "1.0")
        config.set("PFD", "Union", "0")
        config.set("PFD", "Miter_limit", "2.0")
        config.set("PFD", "Punta", "0")
        config.set("PFD", "Guion", "0")
        config.set("PFD", "Dash_offset", "0.0")
        return config


class ConfTooltipUnit(QtGui.QDialog):
    def __init__(self, config, parent=None):
        super(ConfTooltipUnit, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        self.checkShow=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show Tool Tips"))
        self.checkShow.toggled.connect(self.checkShow_Toggled)
        layout.addWidget(self.checkShow,1,1)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),2,1)

        self.groupsystems = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Systems of measurement"))
        self.groupsystems.setFixedHeight(50)
        layout.addWidget(self.groupsystems,3,1)
        lytSystems = QtGui.QHBoxLayout(self.groupsystems)
        self.SI = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "SI"))
        self.SI.toggled.connect(partial(self.systems, "si"))
        lytSystems.addWidget(self.SI)
        self.AltSI = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Alt SI"))
        self.AltSI.toggled.connect(partial(self.systems, "altsi"))
        lytSystems.addWidget(self.AltSI)
        self.English = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "English"))
        self.English.toggled.connect(partial(self.systems, "english"))
        lytSystems.addWidget(self.English)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1)

        self.eleccion=QtGui.QComboBox()
        layout.addWidget(self.eleccion,8,1)
        self.stacked = QtGui.QStackedWidget()
        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        layout.addWidget(self.stacked,9,1)

        self.tabla=[]
        for i, magnitud in enumerate(unidades._magnitudes[:-1]):
            textos=magnitud[2].__text__
            self.tabla.append(QtGui.QTableWidget())
            self.stacked.addWidget(self.tabla[i])

            self.tabla[i].setRowCount(len(textos))
            self.tabla[i].setColumnCount(1)
            self.tabla[i].setColumnWidth(0, 16)
            self.tabla[i].setItemDelegateForColumn(0, CheckEditor(self))
            self.tabla[i].horizontalHeader().setVisible(False)
            for j in range(len(textos)):
                self.tabla[i].setVerticalHeaderItem(j, QtGui.QTableWidgetItem(textos[j]))
                self.tabla[i].setRowHeight(j,24)
                self.tabla[i].setItem(j, 0, QtGui.QTableWidgetItem(""))
                self.tabla[i].item(j, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                self.tabla[i].openPersistentEditor(self.tabla[i].item(j, 0))
            self.rellenar(magnitud[0], i, config)
            self.eleccion.addItem(magnitud[1])

        if config.has_section("Tooltip"):
            self.checkShow.setChecked(config.getboolean("Tooltip", "Show"))


    def rellenar(self, magnitud, tabla, config):
        if config.has_section("Tooltip"):
            lista=eval(config.get("Tooltip",magnitud))
            for i in lista:
                self.tabla[tabla].item(i, 0).setText("true")

    def checkShow_Toggled(self, bool):
        self.eleccion.setEnabled(bool)
        self.groupsystems.setEnabled(bool)
        for tabla in self.tabla:
            tabla.setEnabled(bool)

    def systems(self, set, bool):
        if bool: txt="true"
        else: txt="false"
        for tabla, value in enumerate(unidades.units_set[set][:-1]):
            self.tabla[tabla].item(value, 0).setText(txt)


    def value(self, config):
        if not config.has_section("Tooltip"):
            config.add_section("Tooltip")
        config.set("Tooltip", "Show", self.checkShow.isChecked())
        for i, tabla in enumerate(self.tabla):
            lista=[]
            for j in range(tabla.rowCount()):
                if tabla.item(j, 0).text()=="true":
                    lista.append(j)
            config.set("Tooltip", unidades._magnitudes[i][0], str(lista))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Tooltip")
        config.set("Tooltip", "Show", "True")
        for i, magnitud in enumerate(unidades._magnitudes[:-1]):
            lista=[]
            if magnitud[2]._magnitudes:
                lista.append(magnitud[2].__units__.index(magnitud[2].__units_set__[magnitud[0]]["si"]))
                lista.append(magnitud[2].__units__.index(magnitud[2].__units_set__[magnitud[0]]["english"]))
            else:
                lista.append(magnitud[2].__units__.index(magnitud[2].__units_set__["si"]))
                lista.append(magnitud[2].__units__.index(magnitud[2].__units_set__["english"]))
            config.set("Tooltip", magnitud[0], str(lista))
        return config


class ConfTooltipEntity(QtGui.QDialog):
    def __init__(self, config, parent=None):
        super(ConfTooltipEntity, self).__init__(parent)
        layout = QtGui.QVBoxLayout(self)

        self.eleccion=QtGui.QComboBox()
        layout.addWidget(self.eleccion)
        self.stacked = QtGui.QStackedWidget()
        self.eleccion.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        layout.addWidget(self.stacked)
        layout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Preferred))


        self.tabla=[QtGui.QTableWidget()]
        self.tabla[0].setRowCount(len(corriente.Corriente.propertiesNames()))
        self.tabla[0].setColumnCount(1)
        self.tabla[0].setColumnWidth(0, 16)
        self.tabla[0].setItemDelegateForColumn(0, CheckEditor(self))
        self.tabla[0].horizontalHeader().setVisible(False)
        self.stacked.addWidget(self.tabla[0])
        self.eleccion.addItem(QtGui.QApplication.translate("pychemqt", "Stream"))

        for i, propiedad in enumerate(corriente.Corriente.propertiesNames()):
            self.tabla[0].setVerticalHeaderItem(i, QtGui.QTableWidgetItem(propiedad[0]))
            self.tabla[0].setRowHeight(i,24)
            self.tabla[0].setItem(i, 0, QtGui.QTableWidgetItem(""))
            self.tabla[0].item(i, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.tabla[0].openPersistentEditor(self.tabla[0].item(i, 0))

        if config.has_option("TooltipEntity", "Corriente"):
            lista=eval(config.get("TooltipEntity","Corriente"))
            for i in lista:
                self.tabla[0].item(i, 0).setText("true")

        for i, equipo in enumerate(equipments):
            propiedades=[propiedad[0] for propiedad in equipo.propertiesNames()]
            self.tabla.append(QtGui.QTableWidget())
            self.stacked.addWidget(self.tabla[-1])
            self.tabla[-1].setRowCount(len(propiedades))
            self.tabla[-1].setColumnCount(1)
            self.tabla[-1].setColumnWidth(0, 16)
            self.tabla[-1].setItemDelegateForColumn(0, CheckEditor(self))
            self.tabla[-1].horizontalHeader().setVisible(False)
            for j, propiedad in enumerate(propiedades):
                self.tabla[-1].setVerticalHeaderItem(j, QtGui.QTableWidgetItem(propiedad))
                self.tabla[-1].setRowHeight(j,24)
                self.tabla[-1].setItem(j, 0, QtGui.QTableWidgetItem(""))
                self.tabla[-1].item(j, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                self.tabla[-1].openPersistentEditor(self.tabla[-1].item(j, 0))
            self.rellenar(equipo.__name__, i+1, config)
            self.eleccion.addItem(equipo.title)


    def rellenar(self, equipo, tabla, config):
        if config.has_section("TooltipEntity"):
            lista=eval(config.get("TooltipEntity",equipo))
            for i in lista:
                self.tabla[tabla].item(i, 0).setText("true")

    def value(self, config):
        if not config.has_section("TooltipEntity"):
            config.add_section("TooltipEntity")
        lista=[]
        for j in range(self.tabla[0].rowCount()):
            if self.tabla[0].item(j, 0).text()=="true":
                lista.append(j)
        config.set("TooltipEntity", "Corriente", str(lista))

        for i, tabla in enumerate(self.tabla[1:]):
            lista=[]
            for j in range(tabla.rowCount()):
                if tabla.item(j, 0).text()=="true":
                    lista.append(j)
            config.set("TooltipEntity", equipments[i].__name__, str(lista))

        return config

    @classmethod
    def default(cls, config):
        config.add_section("TooltipEntity")
        config.set("TooltipEntity", "Corriente", str(range(len(corriente.Corriente.propertiesNames()))))
        for i, equipo in enumerate(equipments):
            config.set("TooltipEntity", equipments[i].__name__, str(range(len(equipo.propertiesNames()))))
        return config


class NumericFactor(QtGui.QDialog):
    """Clase que define un dialogo con las opciones de formato numérico"""
    def __init__(self, config, parent=None):
        super(NumericFactor, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Format"))
        layout=QtGui.QGridLayout(self)
        self.checkFixed=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Fixed decimal point"))
        layout.addWidget(self.checkFixed,1,1,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Total digits")),2,2)
        self.TotalDigits=Entrada_con_unidades(int, width=45, value=0, boton=False, spinbox=True, min=0, max=12, showNull=True)
        layout.addWidget(self.TotalDigits,2,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Decimal digits")),3,2)
        self.DecimalDigits=Entrada_con_unidades(int, width=45, value=4, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.DecimalDigits,3,3)
        self.checkSignificant=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Significant figures"))
        layout.addWidget(self.checkSignificant,4,1,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Figures")),5,2)
        self.FiguresSignificatives=Entrada_con_unidades(int, width=45, value=5, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.FiguresSignificatives,5,3)
        self.checkExp=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Exponential preferred"))
        layout.addWidget(self.checkExp,6,1,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Figures")),7,2)
        self.FiguresExponential=Entrada_con_unidades(int, width=45, value=5, boton=False, spinbox=True, min=1, max=12)
        layout.addWidget(self.FiguresExponential,7,3)
        layout.addItem(QtGui.QSpacerItem(30,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),8,1)
        self.checkExpVariable=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Exponential for big/small values"))
        layout.addWidget(self.checkExpVariable,9,1,1,3)
        self.labelTolerancia=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tolerance"))
        layout.addWidget(self.labelTolerancia,10,2)
        self.Tolerance=Entrada_con_unidades(int, width=45, value=4, boton=False, spinbox=True, min=0, max=12)
        layout.addWidget(self.Tolerance,10,3)
        self.checkSign=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show sign in positive values"))
        layout.addWidget(self.checkSign,11,1,1,3)
        self.checkThousand=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show thousand separator"))
        layout.addWidget(self.checkThousand,12,1,1,3)

        self.checkFixed.toggled.connect(self.TotalDigits.setNotReadOnly)
        self.checkFixed.toggled.connect(self.DecimalDigits.setNotReadOnly)
        self.checkSignificant.toggled.connect(self.FiguresSignificatives.setNotReadOnly)
        self.checkExp.toggled.connect(self.ExpToggled)
        self.checkExp.toggled.connect(self.FiguresExponential.setNotReadOnly)
        self.checkExpVariable.toggled.connect(self.Tolerance.setNotReadOnly)
        self.checkExpVariable.toggled.connect(self.labelTolerancia.setEnabled)

        layout.addItem(QtGui.QSpacerItem(20,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),13,1)
        self.muestra=QtGui.QLabel()
        layout.addWidget(self.muestra,14,1,1,3)

        buttonBox = QtGui.QDialogButtonBox()
        buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        layout.addWidget(buttonBox,15,1,1,3)

        self.checkFixed.setChecked(config["format"]==0)
        self.TotalDigits.setNotReadOnly(config["format"]==0)
        self.DecimalDigits.setNotReadOnly(config["format"]==0)
        self.checkSignificant.setChecked(config["format"]==1)
        self.FiguresSignificatives.setNotReadOnly(config["format"]==1)
        self.checkExp.setChecked(config["format"]==2)
        self.FiguresExponential.setNotReadOnly(config["format"]==2)
        if config["format"]==0:
            self.DecimalDigits.setValue(config["decimales"])
        elif config["format"]==1:
            self.FiguresSignificatives.setValue(config["decimales"])
        else:
            self.FiguresExponential.setValue(config["decimales"])
        if config.has_key("total"):
            self.TotalDigits.setValue(config["total"])
        if config.has_key("exp"):
            self.checkExpVariable.setChecked(config["exp"])
        if config.has_key("tol"):
            self.Tolerance.setValue(config["tol"])
        self.Tolerance.setNotReadOnly(config.get("exp", False))
        if config.has_key("signo"):
            self.checkSign.setChecked(config["signo"])
        if config.has_key("thousand"):
            self.checkThousand.setChecked(config["thousand"])

        self.updateMuestra()
        self.checkFixed.toggled.connect(self.updateMuestra)
        self.checkSignificant.toggled.connect(self.updateMuestra)
        self.checkExp.toggled.connect(self.updateMuestra)
        self.checkExpVariable.toggled.connect(self.updateMuestra)
        self.TotalDigits.valueChanged.connect(self.updateMuestra)
        self.DecimalDigits.valueChanged.connect(self.updateMuestra)
        self.FiguresSignificatives.valueChanged.connect(self.updateMuestra)
        self.FiguresExponential.valueChanged.connect(self.updateMuestra)
        self.Tolerance.valueChanged.connect(self.updateMuestra)
        self.checkSign.toggled.connect(self.updateMuestra)
        self.checkThousand.toggled.connect(self.updateMuestra)


    def ExpToggled(self, bool):
        self.FiguresExponential.setNotReadOnly(bool)
        self.checkExpVariable.setDisabled(bool)
        if self.checkExpVariable.isChecked():
            self.labelTolerancia.setDisabled(False)
            self.Tolerance.setReadOnly(True)

    def updateMuestra(self):
        kwargs=self.args()
        txt=QtGui.QApplication.translate("pychemqt", "Sample")+": "+representacion(pi*1e4, **kwargs)
        self.muestra.setText(txt)

    def args(self):
        kwarg={}
        if self.checkFixed.isChecked():
            kwarg["format"]=0
            kwarg["total"]=self.TotalDigits.value
            kwarg["decimales"]=self.DecimalDigits.value
        elif self.checkSignificant.isChecked():
            kwarg["format"]=1
            kwarg["decimales"]=self.FiguresSignificatives.value
        else:
            kwarg["format"]=2
            kwarg["decimales"]=self.FiguresExponential.value
        kwarg["exp"]=self.checkExpVariable.isEnabled() and self.checkExpVariable.isChecked()
        kwarg["tol"]=self.Tolerance.value
        kwarg["signo"]=self.checkSign.isChecked()
        kwarg["thousand"]=self.checkThousand.isChecked()
        return kwarg


class ConfFormat(QtGui.QTableWidget):
    """Clase con el widget de las características generales"""
    def __init__(self, config=None, parent=None):
        super(ConfFormat, self).__init__(parent)
        self.setColumnCount(2)
        self.setRowCount(len(unidades._magnitudes))
        self.setHorizontalHeaderLabels([QtGui.QApplication.translate("pychemqt", "Format"), QtGui.QApplication.translate("pychemqt", "Sample")])
        self.config=[]
        for i in range(self.rowCount()):
            self.setVerticalHeaderItem(i, QtGui.QTableWidgetItem(unidades._magnitudes[i][1]))
            self.setRowHeight(i, 22)
            self.setItem(i, 0, QtGui.QTableWidgetItem(""))
            self.item(i, 0).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.setItem(i, 1, QtGui.QTableWidgetItem(""))
            self.item(i, 1).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)

        if config.has_section("NumericFormat"):
            for i, magnitud in enumerate(unidades._magnitudes):
                formato=eval(config.get("NumericFormat", magnitud[0]))
                self.config.append(formato)
                self.item(i, 0).setText(self.txt(formato))
                self.item(i, 1).setText(representacion(pi, **formato))

        self.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.cellDoubleClicked.connect(self.showConfDialog)

    def showConfDialog(self, fila, columna):
        dialog=NumericFactor(self.config[fila], self)
        if dialog.exec_():
            config=dialog.args()
            self.config[fila]=config
            self.item(fila, 0).setText(self.txt(config))
            self.item(fila, 1).setText(representacion(pi, **config))

    def txt(self, formato):
        if formato["signo"]:
            txt="+"
        else:
            txt=""
        if formato["format"]==0:
            txt+="{total}.{decimales} fixed".format(**formato)
        elif formato["format"]==1:
            txt+="{decimales} sign".format(**formato)
        elif formato["format"]==2:
            txt+="{decimales} exp".format(**formato)
        if formato.get("exp", False):
            txt+=" ({tol} exp)".format(**formato)
        return txt


    def value(self, config):
        if not config.has_section("NumericFormat"):
            config.add_section("NumericFormat")
        for i, magnitud in enumerate(unidades._magnitudes):
            config.set("NumericFormat", magnitud[0], str(self.config[i]))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("NumericFormat")
        for i, magnitud in enumerate(unidades._magnitudes):
            kwarg={'total': 0, 'signo': False, 'decimales': 4, 'format': 0}
            config.set("NumericFormat", magnitud[0], str(kwarg))
        return config


class ConfPetro(QtGui.QDialog):
    def __init__(self, config=None, parent=None):
        super(ConfPetro, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Molecular weight:")), 1, 1, 1, 1)
        self.Peso_molecular = QtGui.QComboBox()
        self.Peso_molecular.addItem("Riazi Daubert")
        self.Peso_molecular.addItem("Riazi Daubert extended")
        self.Peso_molecular.addItem("Lee Kesler")
        self.Peso_molecular.addItem("Sim Daubert")
        self.Peso_molecular.addItem("API")
        self.Peso_molecular.addItem("ASTM")
        self.Peso_molecular.addItem("Goossens")
        self.Peso_molecular.addItem("TWu")
        layout.addWidget(self.Peso_molecular, 1, 2, 1, 1)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Critics properties:")), 2, 1, 1, 1)
        self.critical = QtGui.QComboBox()
        self.critical.addItem("Riazi Daubert")
        self.critical.addItem("Riazi Daubert extended")
        self.critical.addItem("Riazi Adwani")
        self.critical.addItem("Lee Kesler")
        self.critical.addItem("Cavett")
        self.critical.addItem("Sim Daubert")
        self.critical.addItem("Watansiri Owens Starling")
        self.critical.addItem("Edmister")
        self.critical.addItem("Magoulas")
        self.critical.addItem("Twu")
        self.critical.addItem("Tsonopoulos")
        layout.addWidget(self.critical, 2, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Critic volume:")), 3, 1, 1, 1)
        self.vc = QtGui.QComboBox()
        self.vc.addItem("Riazi Daubert")
        self.vc.addItem("Riazi Daubert extended")
        self.vc.addItem("Riazi Adwani")
        self.vc.addItem("Watansiri Owens Starling")
        self.vc.addItem("Twu")
        self.vc.addItem("Hall Yarborough")
        self.vc.addItem("API")
        layout.addWidget(self.vc, 3, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Acentric factor:")), 4, 1, 1, 1)
        self.factor_acentrico = QtGui.QComboBox()
        self.factor_acentrico.addItem("Edmister")
        self.factor_acentrico.addItem("Lee Kesler")
        self.factor_acentrico.addItem("Watansiri Owens Starling")
        self.factor_acentrico.addItem("Magoulas")
        layout.addWidget(self.factor_acentrico, 4, 2, 1, 1)
        layout.addWidget(QtGui.QLabel("Z<sub>c</sub>:"), 5, 1, 1, 1)
        self.Zc = QtGui.QComboBox()
        self.Zc.addItem("Lee Kesler")
        self.Zc.addItem("Haugen")
        self.Zc.addItem("Reid")
        self.Zc.addItem("Salerno")
        self.Zc.addItem("Nath")
        layout.addWidget(self.Zc, 5, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T boiling:")), 6, 1, 1, 1)
        self.t_ebull = QtGui.QComboBox()
        self.t_ebull.addItem("Riazi Daubert extended")
        self.t_ebull.addItem("Riazi Adwani")
        self.t_ebull.addItem("Edmister")
        self.t_ebull.addItem("Soreide")
        layout.addWidget(self.t_ebull, 6, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "PNA descomposition:")), 7, 1, 1, 1)
        self.PNA = QtGui.QComboBox()
        self.PNA.addItem("Peng Robinson")
        self.PNA.addItem("Bergman")
        self.PNA.addItem("Riazi")
        self.PNA.addItem("van Nes van Westen")
        layout.addWidget(self.PNA, 7, 2, 1, 1)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Destilate curve conversion:")), 12, 1, 1, 1)
        self.Curvas = QtGui.QComboBox()
        self.Curvas.addItem("Riazi")
        self.Curvas.addItem("Daubert")
        layout.addWidget(self.Curvas, 12, 2, 1, 1)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "% Hydrogen:")), 13, 1, 1, 1)
        self.Hidrogeno = QtGui.QComboBox()
        self.Hidrogeno.addItem("Riazi")
        self.Hidrogeno.addItem("Goossens")
        self.Hidrogeno.addItem("ASTM")
        self.Hidrogeno.addItem("Jenkins Walsh")
        layout.addWidget(self.Hidrogeno, 13, 2, 1, 1)
        layout.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding), 15, 0, 1, 3)

        if config.has_section("petro"):
            self.Peso_molecular.setCurrentIndex(config.getint("petro","molecular_weight"))
            self.critical.setCurrentIndex(config.getint("petro","critical"))
            self.vc.setCurrentIndex(config.getint("petro","vc"))
            self.factor_acentrico.setCurrentIndex(config.getint("petro","f_acent"))
            self.t_ebull.setCurrentIndex(config.getint("petro","t_ebull"))
            self.Zc.setCurrentIndex(config.getint("petro","Zc"))
            self.PNA.setCurrentIndex(config.getint("petro","PNA"))
            self.Hidrogeno.setCurrentIndex(config.getint("petro","H"))
            self.Curvas.setCurrentIndex(config.getint("petro","curva"))

    def value(self, config):
        if not config.has_section("petro"):
            config.add_section("petro")
        config.set("petro", "molecular_weight", str(self.Peso_molecular.currentIndex()))
        config.set("petro", "critical", str(self.critical.currentIndex()))
        config.set("petro", "vc", str(self.vc.currentIndex()))
        config.set("petro", "f_acent", str(self.factor_acentrico.currentIndex()))
        config.set("petro", "t_ebull", str(self.t_ebull.currentIndex()))
        config.set("petro", "Zc", str(self.Zc.currentIndex()))
        config.set("petro", "PNA", str(self.PNA.currentIndex()))
        config.set("petro", "H", str(self.Hidrogeno.currentIndex()))
        config.set("petro", "curva", str(self.Curvas.currentIndex()))
        return config

    @classmethod
    def default(cls, config):
        config.add_section("petro")
        config.set("petro", "molecular_weight", "0")
        config.set("petro", "critical", "0")
        config.set("petro", "vc", "0")
        config.set("petro", "f_acent", "0")
        config.set("petro", "t_ebull", "0")
        config.set("petro", "Zc", "0")
        config.set("petro", "PNA", "0")
        config.set("petro", "H", "0")
        config.set("petro", "curva", "0")
        return config

class ConfApplications(QtGui.QDialog):
    """Clase con el widget de las aplicaciones preferidas"""
    def __init__(self, config=None, parent=None):
        super(ConfApplications, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        label=QtGui.QApplication.translate("pychemqt", "External Calculator")+":"
        msg=QtGui.QApplication.translate("pychemqt", "Select External Calculator")
        self.calculadora=PathConfig(label, msg=msg, patron="exe")
        layout.addWidget(self.calculadora,1,1)
        label=QtGui.QApplication.translate("pychemqt", "Text viewer")+":"
        msg=QtGui.QApplication.translate("pychemqt", "Select External Text Viewer")
        self.textViewer=PathConfig(label, msg=msg, patron="exe")
        layout.addWidget(self.textViewer,2,1)

        terminal=QtGui.QGroupBox()
        layout.addWidget(terminal,3,1)
        layoutTerminal=QtGui.QGridLayout(terminal)
        msg=QtGui.QApplication.translate("pychemqt", "Select External shell")
        self.terminal=PathConfig("", msg=msg, patron="exe")
        layoutTerminal.addWidget(self.terminal,1,1,1,3)
        layoutTerminal.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Foreground color:")),2,1)
        self.ForegroundColor = ColorSelector()
        layoutTerminal.addWidget(self.ForegroundColor,2,2,1,2)
        layoutTerminal.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Background color:")),3,1)
        self.BackgroundColor = ColorSelector()
        layoutTerminal.addWidget(self.BackgroundColor,3,2,1,2)
        self.ipython=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use ipython if its available"))
        layoutTerminal.addWidget(self.ipython,4,1,1,3)
        self.maximized=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Show maximized"))
        layoutTerminal.addWidget(self.maximized,5,1,1,3)

        layout.addItem(QtGui.QSpacerItem(10,0,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1)

        terminalTitle=QtGui.QApplication.translate("pychemqt", "Shell")
        if sys.platform=="win32":
            terminal.setEnabled(False)
            terminal.setTitle(terminalTitle+" ("+QtGui.QApplication.translate("pychemqt", "Only Available on linux")+")")
        else:
            terminal.setTitle(terminalTitle)

        if config.has_section("Applications"):
            self.calculadora.setText(config.get("Applications", 'Calculator'))
            self.textViewer.setText(config.get("Applications", 'TextViewer'))
            self.terminal.setText(config.get("Applications", 'Shell'))
            self.ipython.setChecked(config.getboolean("Applications", 'ipython'))
            self.maximized.setChecked(config.getboolean("Applications", 'maximized'))
            self.ForegroundColor.setColor(config.get("Applications", 'foregroundColor'))
            self.BackgroundColor.setColor(config.get("Applications", 'backgroundColor'))

        #Habilitar cuando añada soporte para otras terminales
        self.terminal.setEnabled(False)

    def value(self, config):
        if not config.has_section("Applications"):
            config.add_section("Applications")
        config.set("Applications", "Calculator", self.calculadora.text())
        config.set("Applications", "TextViewer", self.textViewer.text())
        config.set("Applications", "Shell",  self.terminal.text())
        config.set("Applications", "ipython", self.ipython.isChecked())
        config.set("Applications", "maximized", self.maximized.isChecked())
        config.set("Applications", "foregroundColor", self.ForegroundColor.color.name())
        config.set("Applications", "backgroundColor", self.BackgroundColor.color.name())
        return config

    @classmethod
    def default(cls, config):
        config.add_section("Applications")
        config.set("Applications", "Calculator", calculator)
        config.set("Applications", "TextViewer", editor)
        config.set("Applications", "Shell", shell)
        config.set("Applications", "ipython", False)
        config.set("Applications", "maximized", False)
        config.set("Applications", "foregroundColor", "#ffffff")
        config.set("Applications", "backgroundColor", "#000000")

        return config


class Isolinea(QtGui.QDialog):
    """Clase generica de widget con la configuración de isolineas"""
    def __init__(self, unidad, ConfSection, config, parent=None):
        """unidad: la clase unidad que define la magnitud de las variables"""
        super(Isolinea, self).__init__(parent)
        self.ConfSection=ConfSection
        self.magnitud=unidad.__name__
        self.unidad=unidad
        layout = QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Start")),1,1)
        self.inicio =Entrada_con_unidades(unidad)
        layout.addWidget(self.inicio,1,2,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Fin")),2,1)
        self.fin =Entrada_con_unidades(unidad)
        layout.addWidget(self.fin,2,2,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Intervalo")),3,1)
        if unidad.__name__=="Temperature":
            self.intervalo =Entrada_con_unidades(unidades.DeltaT)
        elif unidad.__name__=="Pressure":
            self.intervalo =Entrada_con_unidades(unidades.DeltaP)
        else:
            self.intervalo =Entrada_con_unidades(unidad)
        layout.addWidget(self.intervalo,3,2,1,3)
        self.Personalizar = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Customize"))
        layout.addWidget(self.Personalizar,4,1,1,4)
        self.Lista = QtGui.QLineEdit()
        layout.addWidget(self.Lista,5,1,1,4)
        self.Personalizar.toggled.connect(self.inicio.setDisabled)
        self.Personalizar.toggled.connect(self.fin.setDisabled)
        self.Personalizar.toggled.connect(self.intervalo.setDisabled)
        self.Personalizar.toggled.connect(self.Lista.setEnabled)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,1,1,4)
        if unidad.__name__!="float":
            self.Critica = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Include critic point line"))
            layout.addWidget(self.Critica,7,1,1,4)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),8,1,1,4)

        self.lineconfig = LineConfig(ConfSection, QtGui.QApplication.translate("pychemqt", "Line Style"))
        layout.addWidget(self.lineconfig,9,1,1,4)

        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),10,1)
        self.label=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Label"))
        layout.addWidget(self.label,11,1)
        self.unit=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Units in Label"))
        layout.addWidget(self.unit,12,1,1,3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Position")),13,1)
        self.label5=Entrada_con_unidades(int, value=0, width=25, frame=False, readOnly=True)
        self.label5.setFixedWidth(30)
        layout.addWidget(self.label5,13,2)
        self.Posicion=QtGui.QSlider(QtCore.Qt.Horizontal)
        self.Posicion.valueChanged.connect(self.label5.setValue)
        layout.addWidget(self.Posicion,13,3,1,2)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),14,4)

        if config.has_section("MEOS"):
            self.inicio.setValue(config.getfloat("MEOS", ConfSection+'Start'))
            self.fin.setValue(config.getfloat("MEOS", ConfSection+'End'))
            self.intervalo.setValue(config.getfloat("MEOS", ConfSection+'Step'))
            self.Personalizar.setChecked(config.getboolean("MEOS", ConfSection+'Custom'))
            if config.get("MEOS", ConfSection+'List')<>"":
                T=[]
                for i in config.get("MEOS", ConfSection+'List').split(','):
                    if unidad.__name__=="float":
                        T.append(representacion(float(i)))
                    else:
                        T.append(representacion(unidad(float(i)).config()))
                self.Lista.setText(",".join(T))
            if self.unidad.__name__!="float":
                self.Critica.setChecked(config.getboolean("MEOS", ConfSection+'Critic'))
            self.inicio.setDisabled(self.Personalizar.isChecked())
            self.fin.setDisabled(self.Personalizar.isChecked())
            self.intervalo.setDisabled(self.Personalizar.isChecked())
            self.Lista.setEnabled(self.Personalizar.isChecked())
            self.label.setChecked(config.getboolean("MEOS", ConfSection+'Label'))
            self.unit.setChecked(config.getboolean("MEOS", ConfSection+'Units'))
            self.Posicion.setValue(config.getint("MEOS", ConfSection+'Position'))
            self.lineconfig.setConfig(config)


    def value(self, config):
        config.set("MEOS", self.ConfSection+"Start", self.inicio.value)
        config.set("MEOS", self.ConfSection+"End", self.fin.value)
        config.set("MEOS", self.ConfSection+"Step", self.intervalo.value)
        config.set("MEOS", self.ConfSection+"Custom", self.Personalizar.isChecked())
        T=[]
        if not self.Lista.text().isEmpty():
            T1=self.Lista.text().split(',')
            for i in T1:
                if self.unidad.__name__=="float":
                    T.append(str(float(i)))
                else:
                    T.append(str(self.unidad(float(i), "conf")))
        config.set("MEOS", self.ConfSection+"List", ", ".join(T))
        if self.unidad.__name__!="float":
            config.set("MEOS", self.ConfSection+"Critic", self.Critica.isChecked())
        config=self.lineconfig.value(config)

        config.set("MEOS", self.ConfSection+"Label", self.label.isChecked())
        config.set("MEOS", self.ConfSection+"Units", self.unit.isChecked())
        config.set("MEOS", self.ConfSection+"Position", self.Posicion.value())
        return config

    @classmethod
    def default(cls, config, ConfSection):
        config.set("MEOS", ConfSection+"Start", "0")
        config.set("MEOS", ConfSection+"End", "0")
        config.set("MEOS", ConfSection+"Step", "0")
        config.set("MEOS", ConfSection+"Custom", "True")
        config.set("MEOS", ConfSection+"List", "")
        if ConfSection!="Isoquality":
            config.set("MEOS", ConfSection+"Critic", "True")
        config=LineConfig.default(config, ConfSection)
        config.set("MEOS", ConfSection+"Label", "False")
        config.set("MEOS", ConfSection+"Units", "False")
        config.set("MEOS", ConfSection+"Position", "50")
        return config


class ConfmEoS(QtGui.QDialog):
    """Clase que define el widget de configuración de la herramienta de ecuaciones multiparametro"""

    lineas=[("Isotherm", QtGui.QApplication.translate("pychemqt", "Isotherm"), unidades.Temperature),
                ("Isobar", QtGui.QApplication.translate("pychemqt", "Isobar"), unidades.Pressure),
                ("Isoenthalpic", QtGui.QApplication.translate("pychemqt", "Isoenthalpic"), unidades.Enthalpy),
                ("Isoentropic", QtGui.QApplication.translate("pychemqt", "Isoentropic"), unidades.SpecificHeat),
                ("Isochor", QtGui.QApplication.translate("pychemqt", "Isochor"), unidades.SpecificVolume),
                ("Isoquality", QtGui.QApplication.translate("pychemqt", "Isoquality"), float)]

    def __init__(self, config, parent=None):
        super(ConfmEoS, self).__init__(parent)
        layout = QtGui.QGridLayout(self)

        self.iapws = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use iapws97 for water"))
        layout.addWidget(self.iapws,1,1,1,2)
        self.freesteam = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use freesteam library (faster)"))
        self.freesteam.setEnabled(False)
        layout.addWidget(self.freesteam,2,2)
        self.coolProp = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use external library coolProp (faster)"))
        self.coolProp.setEnabled(False)
        layout.addWidget(self.coolProp,3,1,1,2)
        self.refprop = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use external library refprop (fastest)"))
        self.refprop.setEnabled(False)
        layout.addWidget(self.refprop,4,1,1,2)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1)
        self.lineconfig = LineConfig("saturation", QtGui.QApplication.translate("pychemqt", "Saturation Line Style"))
        layout.addWidget(self.lineconfig,5,1,1,2)
        group=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Isolines"))
        layout.addWidget(group,7,1,1,2)
        layoutgroup=QtGui.QGridLayout(group)
        self.comboIsolineas=QtGui.QComboBox()
        layoutgroup.addWidget(self.comboIsolineas,1,1)
        self.Isolineas = QtGui.QStackedWidget()
        self.comboIsolineas.currentIndexChanged.connect(self.Isolineas.setCurrentIndex)
        layoutgroup.addWidget(self.Isolineas,2,1)
        for nombre, texto, unidad in self.lineas:
            self.comboIsolineas.addItem(texto)
            self.Isolineas.addWidget(Isolinea(unidad, nombre, config))
        self.surface = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use surface in plot (slower)"))
        layout.addWidget(self.surface,8,1,1,2)
        self.grid = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Draw grid"))
        layout.addWidget(self.grid,9,1,1,2)

        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,2)

        if os.environ["freesteam"]:
            self.iapws.toggled.connect(self.freesteam.setEnabled)
        if os.environ["CoolProp"]:
            self.coolProp.setEnabled(True)
        if os.environ["refprop"]:
            self.refprop.setEnabled(True)

        if config.has_section("MEOS"):
            self.iapws.setChecked(config.getboolean("MEOS",'iapws'))
            self.freesteam.setChecked(config.getboolean("MEOS",'freesteam'))
            self.coolProp.setChecked(config.getboolean("MEOS",'coolprop'))
            self.refprop.setChecked(config.getboolean("MEOS",'refprop'))
            self.surface.setChecked(config.getboolean("MEOS",'surface'))
            self.grid.setChecked(config.getboolean("MEOS",'grid'))
            self.lineconfig.setConfig(config)


    def value(self, config):
        if not config.has_section("MEOS"):
            config.add_section("MEOS")

        config.set("MEOS", "iapws", self.iapws.isChecked())
        config.set("MEOS", "freesteam", self.freesteam.isChecked())
        config.set("MEOS", "coolprop", self.coolProp.isChecked())
        config.set("MEOS", "refprop", self.refprop.isChecked())
        config=self.lineconfig.value(config)
        config.set("MEOS", "surface", self.surface.isChecked())
        config.set("MEOS", "grid", self.grid.isChecked())

        for indice in range(self.Isolineas.count()):
            config=self.Isolineas.widget(indice).value(config)
        return config

    @classmethod
    def default(cls, config):
        config.add_section("MEOS")
        config.set("MEOS", "iapws", "False")
        config.set("MEOS", "freesteam", "False")
        config.set("MEOS", "coolprop", "False")
        config.set("MEOS", "refprop", "False")
        config=LineConfig.default(config, "saturation")
        config.set("MEOS", "surface", "False")
        config.set("MEOS", "grid", "False")
        for nombre, texto, unidad in cls.lineas:
            config=Isolinea.default(config, nombre)
        return config


class Preferences(QtGui.QDialog):
    """Preferences main dialog"""

    classes = [
        (os.environ["pychemqt"]+"/images/pychemqt.png",
         QtGui.QApplication.translate("pychemqt", "General"), ConfGeneral),
        (os.environ["pychemqt"]+"/images/button/PFD.png",
         QtGui.QApplication.translate("pychemqt", "PFD"), ConfPFD),
        (os.environ["pychemqt"]+"/images/button/tooltip.png",
         QtGui.QApplication.translate("pychemqt", "Tooltips in units"), ConfTooltipUnit),
        (os.environ["pychemqt"]+"/images/button/format_numeric.png",
         QtGui.QApplication.translate("pychemqt", "Numeric format"), ConfFormat),
        (os.environ["pychemqt"]+"/images/button/oil.png",
         QtGui.QApplication.translate("pychemqt", "Pseudocomponents"), ConfPetro),
        (os.environ["pychemqt"]+"/images/button/applications.png",
         QtGui.QApplication.translate("pychemqt", "Applications"), ConfApplications),
        (os.environ["pychemqt"]+"/images/button/tooltip.png",
         QtGui.QApplication.translate("pychemqt", "Tooltips in PFD"), ConfTooltipEntity),
        (os.environ["pychemqt"]+"/images/button/steamTables.png",
         QtGui.QApplication.translate("pychemqt", "mEoS"), ConfmEoS)]

    def __init__(self, config=None, parent=None):
        super(Preferences, self).__init__(parent)
        if not config:
            config = self.default()
        self.config = config
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Preferences"))
        layout = QtGui.QGridLayout(self)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1, 1, 1, 1)
        self.stacked = QtGui.QStackedWidget()
        layout.addWidget(self.stacked, 1, 2, 1, 1)
        self.lista=QtGui.QListWidget()
        self.lista.setFixedWidth(175)
        self.lista.setIconSize(QtCore.QSize(30, 30))
        layout.addWidget(self.lista, 1, 0, 1, 1)
        for icon, title, dialog in self.classes:
            self.stacked.addWidget(dialog(config))
            self.lista.addItem(QtGui.QListWidgetItem(QtGui.QIcon(QtGui.QPixmap(icon)), title))

        self.lista.currentRowChanged.connect(self.stacked.setCurrentIndex)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 0, 1, 3)

    def value(self):
        config=self.config
        for indice in range(self.stacked.count()):
            config=self.stacked.widget(indice).value(config)
        return config

    @classmethod
    def default(cls):
        config=ConfigParser()
        for icon, title, dialog in cls.classes:
            config=dialog.default(config)
        return config


if __name__ == "__main__":
    import sys
    from ConfigParser import ConfigParser
    import os
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    app = QtGui.QApplication(sys.argv)

#    config={"format": 0, "decimales": 4, "signo": False}
#    dialogo=NumericFactor(config)

    config=ConfigParser()
    config.read(conf_dir+"pychemqtrc")
    dialogo = Preferences(config)

    dialogo.show()
    sys.exit(app.exec_())
