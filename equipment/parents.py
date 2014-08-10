#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library for equipment common functionality
#   * equipment: libreria de funcionamiento
#   * UI_equip: Gui
###############################################################################

from functools import partial
import os

from PyQt4 import QtCore, QtGui

from lib.config import Entity
from lib.thread import Evaluate
from UI import texteditor
from UI.widgets import Status
from UI import UI_corriente
from tools.HelpView import HelpView
from tools.costIndex import indiceBase, indiceActual


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

    def __init__(self, **kwargs):
        """Class constructor, copy kwargs for child class to instance and do
        several works with that"""
        self.kwargs = self.__class__.kwargs.copy()
        self.kwargs_forbidden = self.kwargs_forbidden + list(self.kwargsInput)

        # Values defined as integer Entrada_con_unidades return as float
        # and it must be corrected
        self.kwargsInteger = []
        for key, value in self.kwargs.iteritems():
            if isinstance(value, int):
                self.kwargsInteger.append(key)

        self.cleanOldValues(**equipment.kwargs)
        if self.indiceCostos:
            self.kwargs["Base_index"] = indiceBase[self.indiceCostos]
            self.kwargs["Current_index"] = indiceActual[self.indiceCostos]

        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        """All equipment are callables, so we can instance or add/change
        input value with flexibility"""
        oldkwargs = self.kwargs.copy()
        self.cleanOldValues(**kwargs)
        self._bool = True
        txt = kwargs.get("notas", "")
        if txt:
            self.notas = txt
            self.notasPlain = txt

        for key in self.kwargsInteger:
            if key not in self.kwargs_forbidden:
                self.kwargs[key] = int(self.kwargs[key])

        if oldkwargs != self.kwargs and self.isCalculable:
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
        txt += QtGui.QApplication.translate("pychemqt", "Input properties")
        txt += "-----------------#"+os.linesep
        for key, value in self.kwargs.iteritems():
            if value and key not in ["f_install", "Base_index", "Current_index"]:
                txt += key+": "+str(value)+os.linesep
        txt += os.linesep
        txt += self.propTxt()
        return txt

    def propTxt(self):
        """txt equivalent to output properties of equipment"""
        pass

    @classmethod
    def propertiesNames(cls):
        prop = cls.propertiesEquipment()
        prop.append((QtGui.QApplication.translate("pychemqt", "Notes"),
                     "notasPlain", str))
        prop.append((QtGui.QApplication.translate("pychemqt", "Object Type"),
                     "className", str))
        return prop

    @classmethod
    def propertiesEquipment(cls):
        """
        procedure to define output values in a list with format:
        (Name, kwargs name, units), if kwargs name if a combobox element the
        index isn't useful so use a tuple (Txt_Values kwargs_name)
        """
        return []


class UI_equip(QtGui.QDialog):
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
        icono = os.environ["pychemqt"] + \
            "/images/equipment/%s.png" % equipment.__name__.lower()
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(icono)))
        self.evaluate = Evaluate()
        self.evaluate.finished.connect(self.rellenar)

        layout = QtGui.QGridLayout(self)
        self.tabWidget = QtGui.QTabWidget()
        layout.addWidget(self.tabWidget, 0, 0, 1, 3)
        self.status = Status()
        layout.addWidget(self.status, 1, 0, 1, 1)
        self.checkIgnorar = QtGui.QCheckBox()
        self.checkIgnorar.setText(QtGui.QApplication.translate("pychemqt",
                                                               "Ignore"))
        self.checkIgnorar.toggled.connect(self.ignorar)
        layout.addWidget(self.checkIgnorar, 1, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel |
                                                QtGui.QDialogButtonBox.Ok |
                                                QtGui.QDialogButtonBox.Help)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.helpRequested.connect(self.ayuda)
        layout.addWidget(self.buttonBox, 1, 2, 1, 1)

        if not equipment.help:
            self.buttonBox.button(QtGui.QDialogButtonBox.Help).setVisible(False)

        # Input tab
        if entrada:
            self.entrada = QtGui.QTabWidget()
            self.tabWidget.addTab(
                self.entrada, QtGui.QApplication.translate("pychemqt", "Input"))
        elif entrada is None:
            pass
        else:
            self.entrada = UI_corriente.Ui_corriente()
            self.entrada.Changed.connect(partial(self.changeParams, "entrada"))
            self.tabWidget.addTab(
                self.entrada, QtGui.QApplication.translate("pychemqt", "Input"))

        # Calcule tab
        if calculo:
            self.tabCalculo = QtGui.QWidget()
            self.tabWidget.addTab(self.tabCalculo,
                QtGui.QApplication.translate("pychemqt", "Calculation"))

        # Cost tab
        if equipment.indiceCostos is not None:
            self.tabCostos = QtGui.QWidget()
            self.tabWidget.addTab(
                self.tabCostos, QtGui.QApplication.translate("pychemqt", "Cost"))

        # Output tab
        if salida:
            self.Salida = QtGui.QTabWidget()
            self.tabWidget.addTab(
                self.Salida, QtGui.QApplication.translate("pychemqt", "Output"))
        elif salida is None:
            pass
        else:
            self.Salida = UI_corriente.Ui_corriente(readOnly=True)
            self.tabWidget.addTab(
                self.Salida, QtGui.QApplication.translate("pychemqt", "Output"))

        # Notes tab
        self.tabNotas = texteditor.TextEditor()
        self.tabWidget.addTab(
            self.tabNotas, QtGui.QApplication.translate("pychemqt", "Notes"))
        self.tabNotas.notas.textChanged.connect(self.cambiar_notas)

    def addSalida(self, title):
        widget = UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(widget, title)

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
        Dialog = HelpView(self.windowTitle(), QtCore.QUrl(self.Equipment.help))
        Dialog.exec_()

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
        for entrada in self.Equipment.kwargsInput:
            self.__getattribute__(entrada).blockSignals(True)
            self.__getattribute__(entrada).setCorriente(
                self.Equipment.kwargs[entrada])
            self.__getattribute__(entrada).blockSignals(False)
        for variable in self.Equipment.kwargsValue:
            self.__getattribute__(variable).setValue(
                self.Equipment.kwargs[variable])
        for combo in self.Equipment.kwargsList:
            self.__getattribute__(combo).setCurrentIndex(
                self.Equipment.kwargs[combo])
        for check in self.Equipment.kwargsCheck:
            self.__getattribute__(check).setChecked(self.Equipment.kwargs[check])
        if self.Equipment.indiceCostos is not None:
            self.Costos.setFactor(self.Equipment.kwargs["f_install"])
            self.Costos.setBase(self.Equipment.kwargs["Base_index"])
            self.Costos.setActual(self.Equipment.kwargs["Current_index"])
        self.blockSignals(False)
        self.status.setState(self.Equipment.status, self.Equipment.msg)
