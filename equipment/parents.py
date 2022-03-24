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
# Library for equipment common functionality
#   * equipment: Base class of equipment library
#   * UI_equip: Base class of equipment UI functionality
###############################################################################


from functools import partial
import logging
import os

from PyQt5 import QtCore, QtGui, QtWidgets

from lib.config import Entity, IMAGE_PATH
from lib.thread import Evaluate
from tools.costIndex import indiceBase, indiceActual
from UI.texteditor import TextEditor
from UI.UI_corriente import Ui_corriente
from UI.widgets import Status


class equipment(Entity):
    """General structure for equipment, each child class must define the
    properties and procedures
        status: equipment estatus
            0: undefined
            1: solved
            2: error
            3: solved but with warning
        msg: status message of equipment
        title: string title of equipment
        help: help file its available
        kwargs: dict with incomming variables and a undefined (nul) value. Its
         important define fine the initial value because it used to detect
         type to save to file
            - Streams: None
            - List options and integer: 0
            - values: 0.0
        kwargsInput: Inputs of stream type
        kwargsValue: Inputs of value string
        kwargsList: Inputs of combobox (option of a list)
        kwargsCheck: Inputs of checkbox (True/False)
        calculateValue: Results values
        statusCoste: Cost section status
        indiceCostos: Index of equipment in costindex
        calculateCostos: Results values in cost section
        __doi__: Bibliography in array of dict format

        isCalculable: procedure to define solvable status of equipment, use too
         to define cost status and other definition options
        calculo: procedure to calculation
        coste: procedure to cost calculation
        cleanOldValues: procedure to update kwargs values its neccessary
        propTxt: procedure to define text format for simple text report
        propertiesEquipment: procedure to define output values in a list
         with format (Name, kwargs name, units)
    """
    status = 0
    msg = ""
    title = ""
    help = ""
    kwargs = {}
    kwargsInput = ()
    kwargsValue = ()
    kwargsList = ()
    kwargsCheck = ()
    calculateValue = ()
    statusCoste = False
    indiceCostos = None
    calculateCostos = ()
    __doi__ = []
    salida = [None]

    def __init__(self, **kwargs):
        """Class constructor, copy kwargs for child class to instance and do
        several works with that"""
        self.kwargs = self.__class__.kwargs.copy()
        self.kwargs["notas"] = ""
        self.kwargs_forbidden = self.kwargs_forbidden + list(self.kwargsInput)

        # Values defined as integer Entrada_con_unidades return as float
        # and it must be corrected
        self.kwargsInteger = []
        for key, value in list(self.kwargs.items()):
            if isinstance(value, int):
                self.kwargsInteger.append(key)

        self.cleanOldValues(**equipment.kwargs)
        if self.indiceCostos:
            self.kwargs["Base_index"] = indiceBase[self.indiceCostos]
            self.kwargs["Current_index"] = indiceActual[self.indiceCostos]

        # Calculate the maximum properties name length for best report
        # formating
        letter = 0
        for name in self.propertiesTitle():
            if len(name) > letter:
                letter = len(name)
        self.TEXT_FORMATING_LENG = (letter//4+1)*4

        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        """All equipment are callables, so we can instance or add/change
        input value with flexibility"""
        Entity.__call__(self, **kwargs)
        input = False
        for key in kwargs:
            if key in self.kwargsInput:
                input = True
                break
        if self.isCalculable and (self._oldkwargs != self.kwargs or input):
            logging.info('Calculate EQUIPMENT: %s' % self.__class__.__name__)
            kw_new = {}
            for key, value in list(kwargs.items()):
                if self.kwargs[key] != value:
                    kw_new[key] = value
            logging.debug('kwarg; %s' % kw_new)
            QtWidgets.QApplication.processEvents()
            self.calculo()
            if self.statusCoste:
                self.coste()

    @property
    def isCalculable(self):
        """Each child class must define if its calculable for input kwargs"""
        pass

    def calculo(self):
        """Procedure to calcute equipment, defined in child class"""
        pass

    def cleanOldValues(self, **kwargs):
        """Update kwargs with new input kwargs, defined in child class,
        here can be implemented kwarg incompatibiity input and more"""
        self.kwargs.update(kwargs)

    def clear(self):
        """clear equipment instance to a fresh new instance"""
        self.kwargs = self.__class__.kwargs
        self.kwargs.update(equipment.kwargs)
        self.__dict__.clear()
        self._bool = False

    def txt(self):
        """Return plain text to report with input properties of equipment"""
        txt = str(self.notasPlain)+os.linesep+os.linesep
        txt += "#---------------"
        txt += QtWidgets.QApplication.translate("pychemqt", "Input properties")
        txt += "-----------------#"+os.linesep
        mask = "%s-%is%ss" % ("%", self.TEXT_FORMATING_LENG + 1, "%")
        for key, val in list(self.kwargs.items()):
            if val and key not in ["f_install", "Base_index", "Current_index"]:
                txt += mask % (key, val) + os.linesep
        txt += os.linesep

        txt += self.propTxt()
        return txt

    def propTxt(self):
        """txt equivalent to output properties of equipment"""
        pass

    @classmethod
    def propertiesNames(cls):
        p = cls.propertiesEquipment()
        p.append((QtWidgets.QApplication.translate("pychemqt", "Notes"),
                  "notasPlain", str))
        p.append((QtWidgets.QApplication.translate("pychemqt", "Object Type"),
                  "className", str))
        return p

    @classmethod
    def propertiesEquipment(cls):
        """
        procedure to define output values in a list with format:
        (Name, kwargs name, units), if kwargs name if a combobox element the
        index isn't useful so use a tuple (Txt_Values kwargs_name)
        """
        return []


class UI_equip(QtWidgets.QDialog):
    """UI general for equipments, each child class must define specifics"""
    def __init__(self, equipment, entrada=True, salida=True, calculo=True,
                 parent=None):
        """
        equipment: name of equipment to model
        entrada: boolean to create or not the input tab
        salida: boolean to create or not the input tab
            - True para equipos con varias entradas/salidas, create de tab,
              the child must define the UI_corriente
            - False para equipos con una, create UI_corriente
            - None: Not create nothing
        calculo: boolean to create or not the calcule tab
        """
        super(UI_equip, self).__init__(parent)
        self.setWindowTitle(equipment.title)
        icono = os.path.join(IMAGE_PATH, "equipment",
                             "%s.png" % equipment.__name__.lower())
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(icono)))
        self.evaluate = Evaluate()
        self.evaluate.finished.connect(self.rellenar)

        layout = QtWidgets.QGridLayout(self)
        self.tabWidget = QtWidgets.QTabWidget()
        layout.addWidget(self.tabWidget, 0, 0, 1, 3)
        self.status = Status()
        layout.addWidget(self.status, 1, 0, 1, 1)
        self.checkIgnorar = QtWidgets.QCheckBox()
        self.checkIgnorar.setText(
            QtWidgets.QApplication.translate("pychemqt", "Ignore"))
        self.checkIgnorar.toggled.connect(self.ignorar)
        layout.addWidget(self.checkIgnorar, 1, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok |
            QtWidgets.QDialogButtonBox.Help)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.helpRequested.connect(self.ayuda)
        layout.addWidget(self.buttonBox, 1, 2, 1, 1)

        if not equipment.help:
            button = self.buttonBox.button(QtWidgets.QDialogButtonBox.Help)
            button.setVisible(False)

        # Input tab
        if entrada:
            self.Entrada = QtWidgets.QTabWidget()
            self.tabWidget.addTab(
                self.Entrada,
                QtGui.QIcon(os.path.join(IMAGE_PATH, "equipment", "in.svg")),
                QtWidgets.QApplication.translate("pychemqt", "Input"))
        elif entrada is None:
            pass
        else:
            self.Entrada = Ui_corriente()
            self.Entrada.Changed.connect(partial(self.changeParams, "entrada"))
            self.tabWidget.addTab(
                self.Entrada,
                QtGui.QIcon(os.path.join(IMAGE_PATH, "equipment", "in.svg")),
                QtWidgets.QApplication.translate("pychemqt", "Input"))

        # Calcule tab
        if calculo:
            self.tabCalculo = QtWidgets.QWidget()
            self.tabWidget.addTab(
                self.tabCalculo,
                QtGui.QIcon(os.path.join(
                    IMAGE_PATH, "button", "calculator.png")),
                QtWidgets.QApplication.translate("pychemqt", "Calculation"))

        # Cost tab
        if equipment.indiceCostos is not None:
            self.tabCostos = QtWidgets.QWidget()
            self.tabWidget.addTab(
                self.tabCostos,
                QtGui.QIcon(os.path.join(
                    IMAGE_PATH, "button", "currency.png")),
                QtWidgets.QApplication.translate("pychemqt", "Cost"))

        # Output tab
        if salida:
            self.Salida = QtWidgets.QTabWidget()
            self.tabWidget.addTab(
                self.Salida,
                QtGui.QIcon(os.path.join(IMAGE_PATH, "equipment", "out.svg")),
                QtWidgets.QApplication.translate("pychemqt", "Output"))
        elif salida is None:
            pass
        else:
            self.Salida = Ui_corriente(readOnly=True)
            self.tabWidget.addTab(
                self.Salida,
                QtGui.QIcon(os.path.join(IMAGE_PATH, "equipment", "out.svg")),
                QtWidgets.QApplication.translate("pychemqt", "Output"))

        # Notes tab
        self.tabNotas = TextEditor()
        self.tabWidget.addTab(
            self.tabNotas,
            QtGui.QIcon(os.path.join(IMAGE_PATH, "button", "editor.png")),
            QtWidgets.QApplication.translate("pychemqt", "Notes"))
        self.tabNotas.notas.textChanged.connect(self.cambiar_notas)

    def addSalida(self, title, **kw):
        widget = Ui_corriente(readOnly=True, **kw)
        self.Salida.addTab(widget, title)

    def addEntrada(self, title, key, **kw):
        widget = Ui_corriente(**kw)
        widget.Changed.connect(partial(self.changeParams, key))
        self.Entrada.addTab(widget, title)

    def ignorar(self, bool):
        """Ignore the equipment"""
        if bool:
            self.status.setState(2)
        else:
            self.status.restaurar()
        self.tabWidget.setEnabled(not bool)

    def cambiar_notas(self):
        """Change notes properties"""
        htm = self.tabNotas.notas.toHtml()
        txt = self.tabNotas.notas.toPlainText()
        self.Equipment.setNotas(htm, txt)

    def ayuda(self):
        """Show help page"""
        url = QtCore.QUrl(self.Equipment.help)
        QtGui.QDesktopServices.openUrl(url)

    def setEquipment(self, equipment):
        self.Equipment = equipment
        self.rellenar()

    def changeParams(self, key, value):
        """Change any kwargs value"""
        self.calculo(**{key: value})

    def changeParamsCoste(self, parametro, valor):
        """Change any cost kwarg value,
        separate of normal calcule to improve performance"""
        self.Equipment.cleanOldValues(**{str(parametro): valor})
        if self.Equipment.status:
            self.Equipment.coste()
            self.rellenar()

    def calculo(self, **kwargs):
        """Start equipment calcule
        use a different thread to improve UI response"""
        self.status.setState(4)
        self.evaluate.start(self.Equipment, kwargs)

    def rellenar(self):
        """Fill widget with equipment values"""
        self.rellenarInput()
        if self.Equipment.status in [1, 3]:
            self.tabNotas.setText(self.Equipment.notas)
            for variable in self.Equipment.calculateValue:
                self.__getattribute__(variable).setValue(
                    self.Equipment.__getattribute__(variable))
            if len(self.Equipment.salida) == 1:
                self.Salida.setCorriente(self.Equipment.salida[0])
            else:
                for i, salida in enumerate(self.Equipment.salida):
                    self.Salida.widget(i).setCorriente(salida)

            if self.Equipment.indiceCostos is not None and \
                    self.Equipment.statusCoste:
                for variable in self.Equipment.calculateCostos:
                    self.__getattribute__(variable).setValue(
                        self.Equipment.__getattribute__(variable))
        self.status.setState(self.Equipment.status, self.Equipment.msg)

    def rellenarInput(self):
        """Fill widget with input value of equipment"""
        self.blockSignals(True)
        if len(self.Equipment.kwargsInput) == 1:
            self.Entrada.blockSignals(True)
            entrada = self.Equipment.kwargsInput[0]
            self.Entrada.setCorriente(self.Equipment.kwargs[entrada])
            self.Entrada.blockSignals(False)
        else:
            for i, entrada in enumerate(self.Equipment.kwargsInput):
                widget = self.Entrada.widget(i)
                widget.blockSignals(True)
                widget.setCorriente(self.Equipment.kwargs[entrada])
                widget.blockSignals(False)
        for variable in self.Equipment.kwargsValue:
            self.__getattribute__(variable).setValue(
                self.Equipment.kwargs[variable])
        for combo in self.Equipment.kwargsList:
            self.__getattribute__(combo).setCurrentIndex(
                self.Equipment.kwargs[combo])
        for chck in self.Equipment.kwargsCheck:
            self.__getattribute__(chck).setChecked(self.Equipment.kwargs[chck])
        if self.Equipment.indiceCostos is not None:
            self.Costos.setFactor(self.Equipment.kwargs["f_install"])
            self.Costos.setBase(self.Equipment.kwargs["Base_index"])
            self.Costos.setActual(self.Equipment.kwargs["Current_index"])
        self.blockSignals(False)
#        self.status.setState(self.Equipment.status, self.Equipment.msg)
