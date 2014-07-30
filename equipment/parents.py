#!/usr/bin/python
# -*- coding: utf-8 -*-


from functools import partial
import os

from PyQt4 import QtCore, QtGui

from lib.config import Entity
from lib.thread import Evaluate
from lib.unidades import unidad, Fouling, Length, ThermalConductivity
from lib.table import finnedTube_database, thermalConductivity

from UI import texteditor
from UI.widgets import Status, Entrada_con_unidades
from UI import UI_corriente
from tools.HelpView import HelpView
from tools.costIndex import indiceBase, indiceActual



class equipment(Entity):
    """Esquema general de un equipo, define las propiedades comunes a todos ellos"""
    status=0
    msg=""
    statusCoste=False
    kwargs={}
    kwargsInput=()
    kwargsValue=()
    kwargsList=()
    kwargsCheck=()
    calculateValue=()
    calculateCostos=()

    indiceCostos=None

    def __init__(self, **kwargs):
        """Lee la documentación de cada equpo para conocer los kwargs aceptados"""
        self.kwargs=self.__class__.kwargs.copy()
        self.kwargs_forbidden=self.kwargs_forbidden+list(self.kwargsInput)
        
        #Los valores que estan definidos como integer en la Entrada_con_unidades los devuelve como float y hay que rectificarlos
        self.kwargsInteger=[]
        for key, value in self.kwargs.iteritems():
            if isinstance(value, int):
                self.kwargsInteger.append(key)
        
        self.cleanOldValues(**equipment.kwargs)
        if self.indiceCostos:
            self.kwargs["Base_index"]=indiceBase[self.indiceCostos]
            self.kwargs["Current_index"]=indiceActual[self.indiceCostos]
        
        if kwargs:
            self.__call__(**kwargs)
            
    def __call__(self, **kwargs):
        """Todos los equipos son calables, al invocarlos sin argumentos solo se instancia, si se invocan con argumentos se invoca el método de cálculo."""
        oldkwargs=self.kwargs.copy()
        self.cleanOldValues(**kwargs)
        self._bool=True
        txt=kwargs.get("notas", "")
        if txt:
            self.notas=txt
            self.notasPlain=txt
        
        for key in self.kwargsInteger:
            if key not in self.kwargs_forbidden:
                self.kwargs[key]=int(self.kwargs[key])
            
        if oldkwargs!=self.kwargs and self.isCalculable:
            self.calculo()
            if self.statusCoste:
                self.coste()

        
    @property
    def isCalculable(self):
        """Método que estima si el equipo es calculable en función de los datos disponibles, definido por cada equipo"""
        pass
    def calculo(self):
        """Procedimiento de cálculo propiamente dicho, definido por cada equipo"""
        pass
    def cleanOldValues(self, **kwargs):
        """Actualización de los kwargs con los nuevos introducidos si es necesario para cada equipo"""
        self.kwargs.update(kwargs)
        
    def clear(self):
        self.kwargs=self.__class__.kwargs
        self.kwargs.update(equipment.kwargs)
        self.__dict__.clear()
        self._bool=False       
        
    def txt(self):
        """Devuelve el texto para report en formato texto"""
        txt=str(self.notasPlain)+os.linesep+os.linesep
        txt+="#---------------"+QtGui.QApplication.translate("pychemqt", "Input properties")+"-----------------#"+os.linesep
        for key, value in self.kwargs.iteritems():
            if value and key not in ["f_install", "Base_index", "Current_index"]:
                txt+=key+": "+str(value)+os.linesep
        txt+=os.linesep
        txt+=self.propTxt()
        return txt
    
    def propTxt(self):
        """texto especifico de cada equipo"""
        pass

    @classmethod
    def propertiesNames(cls):
        propiedades=cls.propertiesEquipment()
        propiedades.append((QtGui.QApplication.translate("pychemqt", "Notes"), "notasPlain", str))
        propiedades.append((QtGui.QApplication.translate("pychemqt", "Object Type"), "className", str)) 
        return propiedades

    @classmethod
    def propertiesEquipment(cls):
        """Implementar en cada equipo"""
        return []


class UI_equip(QtGui.QDialog):
    """Diálogo general de definición de equipos con la estructura general común a todos ellos"""
    def __init__(self, tipo, entrada=True, salida=True, costos=True, calculo=True, parent=None):
        """
        tipo: nombre de la clase que modela el equipo
        costos: boolean que indica si se crea la pestaña de costos
        entrada: boolean que indica si se crea la pestaña de entrada
        salida: boolean que indica si se crea la pestaña de salida
            (generalmente True para equipos con una entrada/salida, False para equipos con varias)
        """
        super(UI_equip, self).__init__(parent)
        self.setWindowTitle(tipo.title)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equipment/%s.png" % tipo.__name__.lower())))
        self.evaluate=Evaluate()
        self.evaluate.finished.connect(self.rellenar)

        layout = QtGui.QGridLayout(self)
        self.tabWidget = QtGui.QTabWidget()
        layout.addWidget(self.tabWidget, 0, 0, 1, 3)
        self.status=Status()
        layout.addWidget(self.status, 1, 0, 1, 1)
        self.checkIgnorar=QtGui.QCheckBox()
        self.checkIgnorar.setText(QtGui.QApplication.translate("pychemqt", "Ignore"))
        self.checkIgnorar.toggled.connect(self.ignorar)
        layout.addWidget(self.checkIgnorar, 1, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Help)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.helpRequested.connect(self.ayuda)
        layout.addWidget(self.buttonBox, 1, 2, 1, 1)

        if not tipo.help:
            self.buttonBox.button(QtGui.QDialogButtonBox.Help).setVisible(False)
        
        #Pestaña entrada
        if entrada:
            self.entrada = QtGui.QTabWidget()
            self.tabWidget.addTab(self.entrada,QtGui.QApplication.translate("pychemqt", "Input"))
        elif entrada==None:
            pass
        else:
            self.entrada= UI_corriente.Ui_corriente()
            self.entrada.Changed.connect(partial(self.changeParams, "entrada"))
            self.tabWidget.addTab(self.entrada, QtGui.QApplication.translate("pychemqt", "Input"))

        
        #Pestaña cálculo
        if calculo:
            self.tabCalculo = QtGui.QWidget()
            self.tabWidget.addTab(self.tabCalculo,QtGui.QApplication.translate("pychemqt", "Calculation"))

        #Pestaña costos
        if costos:
            self.tabCostos = QtGui.QWidget()
            self.tabWidget.addTab(self.tabCostos,QtGui.QApplication.translate("pychemqt", "Cost"))
        
        #Pestaña salida
        if salida:
            self.Salida = QtGui.QTabWidget()
            self.tabWidget.addTab(self.Salida,QtGui.QApplication.translate("pychemqt", "Output"))
        elif salida==None:
            pass
        else:
            self.Salida= UI_corriente.Ui_corriente(readOnly=True)
            self.tabWidget.addTab(self.Salida,QtGui.QApplication.translate("pychemqt", "Output"))


        #Notas
        self.tabNotas = texteditor.TextEditor()
        self.tabWidget.addTab(self.tabNotas,QtGui.QApplication.translate("pychemqt", "Notes"))
        self.tabNotas.notas.textChanged.connect(self.cambiar_notas)


    def addSalida(self, title):
        widget=UI_corriente.Ui_corriente(readOnly=True)
        self.Salida.addTab(widget, title)

    def ignorar(self, bool):
        if bool:
            self.status.setState(2)
        else:
            self.status.restaurar()
        self.tabWidget.setEnabled(not bool)
        
    def cambiar_notas(self):
        self.Equipment.setNotas(self.tabNotas.notas.toHtml(), self.tabNotas.notas.toPlainText())

    def ayuda(self):
        Dialog = HelpView(self.windowTitle(), QtCore.QUrl(self.Equipment.help))
        Dialog.exec_()     

    def setEquipment(self, equipment):
        self.Equipment=equipment
        self.rellenar()
        
    def changeParams(self, parametro, valor):
        self.calculo(**{parametro: valor})
        
    def changeParamsCoste(self, parametro, valor):
        self.Equipment.cleanOldValues(**{str(parametro): valor})
        if self.Equipment.status:
            self.Equipment.coste()
            self.rellenar()

    def calculo(self, **kwargs):
        self.status.setState(4)
        self.evaluate.start(self.Equipment, kwargs)

    def rellenar(self):
        """Método que rellena los widgets con los valores del equipo"""
        self.rellenarInput()
        if self.Equipment.status in [1, 3]:
            self.tabNotas.setText(self.Equipment.notas)
            for variable in self.Equipment.calculateValue:
                self.__getattribute__(variable).setValue(self.Equipment.__getattribute__(variable))
            if len(self.Equipment.salida)==1:
                self.Salida.setCorriente(self.Equipment.salida[0])
            else:
                for i, salida in enumerate(self.Equipment.salida):
                    self.Salida.widget(i).setCorriente(salida)
                
            if self.Equipment.indiceCostos!=None and self.Equipment.statusCoste:
                for variable in self.Equipment.calculateCostos:
                    self.__getattribute__(variable).setValue(self.Equipment.__getattribute__(variable))
        self.status.setState(self.Equipment.status, self.Equipment.msg)
        
    def rellenarInput(self):
        """Método que rellena los widgets de entrada de datos del equipo"""
        self.blockSignals(True)
        for entrada in self.Equipment.kwargsInput:
            self.__getattribute__(entrada).blockSignals(True)
            self.__getattribute__(entrada).setCorriente(self.Equipment.kwargs[entrada])
            self.__getattribute__(entrada).blockSignals(False)
        for variable in self.Equipment.kwargsValue:
            self.__getattribute__(variable).setValue(self.Equipment.kwargs[variable])
        for combo in self.Equipment.kwargsList:
            self.__getattribute__(combo).setCurrentIndex(self.Equipment.kwargs[combo])
        for check in self.Equipment.kwargsCheck:
            self.__getattribute__(check).setChecked(self.Equipment.kwargs[check])
        if self.Equipment.indiceCostos!=None:
            self.Costos.setFactor(self.Equipment.kwargs["f_install"])
            self.Costos.setBase(self.Equipment.kwargs["Base_index"])
            self.Costos.setActual(self.Equipment.kwargs["Current_index"])
        self.blockSignals(False)
        self.status.setState(self.Equipment.status, self.Equipment.msg)

class UI_equipment(UI_equip):
    pass

class FoulingWidget(QtGui.QWidget):
    """Widget con los parametros de fouling de tuberias"""
    valueChanged = QtCore.pyqtSignal(float)
    Fouling_Factor={
        "Industrial": {
            "Fuel oil no.2": 0.000352, 
            "Fuel oil no.6": 0.000881, 
            "Transformer oil": 0.000173, 
            "Engine Lube oil": 0.000173, 
            "Quench oil": 0.000705, 
            "Manufactured gas": 0.001761, 
            "Engine exhaust gas": 0.001761, 
            "Steam (nonoil bearing)": 0.000088, 
            "Exhaust steam (oil bearing)": 0.0003, 
            "Refrigerant vapors (Oil bearing)": 0.000352, 
            "Compressed air": 0.000176, 
            "Ammonia vapor": 0.000176, 
            "CO2 vapor": 0.000176, 
            "Chlorine vapor": 0.000352, 
            "Coal flue gas": 0.001761, 
            "Natural gas flue gas": 0.000881, 
            "Molten heat transfer salts": 0.000088, 
            "Refrigerant liquids": 0.000176, 
            "Hydraulic fluid": 0.000176, 
            "Industrial organic heat transfer media": 0.000352, 
            "Ammonia liquid": 0.000176, 
            "Ammonia liquid (oil bearing)": 0.000528, 
            "Calcium chloride solutions": 0.000528, 
            "Sodium chloride solutions": 0.000528, 
            "CO2 liquid": 0.000176, 
            "Chlorine liquid": 0.000352, 
            "Methanol solutions": 0.000352, 
            "Ethanol solutions": 0.000352, 
            "Ethilene glycol solutions": 0.000352}, 
        "Chemical": {
            "Acid gases": 0.00044, 
            "Solvent vapors": 0.000176, 
            "Stable overhead products": 0.000176, 
            "MEA and DEA solutions": 0.000352, 
            "DEG and TEG solutions": 0.000352, 
            "Stable side draw and bottom product": 0.00026, 
            "Caustic solutions": 0.000352, 
            "Vegetable oils": 0.000528}, 
        "Natural Gas-Gasoline": {
            "Natural gas": 0.00026, 
            "Overhead products": 0.00026, 
            "Lean oil": 0.000352, 
            "Rich oil": 0.00026, 
            "Natural gasoline": 0.00026, 
            "Liquified petroleum gases": 0.00026}, 
        u"Water (T<50C, v<0.9)": {
            "Seawater": 0.000088, 
            "Brackish water": 0.000352, 
            "Cooling tower (treated)": 0.000176, 
            "Cooling tower (untreaterd)": 0.000528, 
            "City or well water": 0.000176, 
            "River water minimum": 0.000352, 
            "River water average": 0.000528, 
            "Muddy or silty": 0.000528, 
            "Hard (>15 grains/gal": 0.000528, 
            "Engine jacket": 0.000176, 
            "Distilled, condensate": 0.000088, 
            "Distilled, boiler blowdown": 0.000352, 
            "Distilled, treated boiler feedwater": 0.000176}, 
        u"Water (T<50C, v>0.9)": {
            "Seawater": 0.000088, 
            "Brackish water": 0.000176, 
            "Cooling tower (treated)": 0.000176, 
            "Cooling tower (untreaterd)": 0.000528, 
            "City or well water": 0.000176, 
            "River water minimum": 0.000176, 
            "River water average": 0.000352, 
            "Muddy or silty": 0.000352, 
            "Hard (>15 grains/gal": 0.000528, 
            "Engine jacket": 0.000176, 
            "Distilled, condensate": 0.000088, 
            "Distilled, boiler blowdown": 0.000352, 
            "Distilled, treated boiler feedwater": 0.000088}, 
        u"Water (T>50C, v<0.9)": {
            "Seawater": 0.000176, 
            "Brackish water": 0.000528, 
            "Cooling tower (treated)": 0.000352, 
            "Cooling tower (untreaterd)": 0.000881, 
            "City or well water": 0.000352, 
            "River water minimum": 0.000528, 
            "River water average": 0.000705, 
            "Muddy or silty": 0.000705, 
            "Hard (>15 grains/gal": 0.000881, 
            "Engine jacket": 0.000176, 
            "Distilled, condensate": 0.000088, 
            "Distilled, boiler blowdown": 0.000352, 
            "Distilled, treated boiler feedwater": 0.000176}, 
        u"Water (T>50C, v>0.9)": {
            "Seawater": 0.000176, 
            "Brackish water": 0.000352, 
            "Cooling tower (treated)": 0.000352, 
            "Cooling tower (untreaterd)": 0.000705, 
            "City or well water": 0.000352, 
            "River water minimum": 0.000352, 
            "River water average": 0.000528, 
            "Muddy or silty": 0.000528, 
            "Hard (>15 grains/gal": 0.000881, 
            "Engine jacket": 0.000176, 
            "Distilled, condensate": 0.000088, 
            "Distilled, boiler blowdown": 0.000352, 
            "Distilled, treated boiler feedwater": 0.000176}, 
        u"Refinery vapors": {
            "Atmospheric tower overhead vapors": 0.000176, 
            "Light naphthas": 0.000176, 
            "Vacuum overhead vapors": 0.000352}, 
        u"Refinery liq.": {
            "Crude oil dry T<120C, v<0.6": 0.000528, 
            "Crude oil salt T<120C, v<0.6": 0.000528, 
            "Crude oil dry T<120C, 0.6<v<1.2": 0.000352, 
            "Crude oil salt T<120C, 0.6<v<1.2": 0.000352, 
            "Crude oil dry T<120C, 1.2<v": 0.000352, 
            "Crude oil salt T<120C, 1.2<v": 0.000352, 
            "Crude oil dry 120C<T<175C, v<0.6": 0.000528, 
            "Crude oil salt 120C<T<175C, v<0.6": 0.000881, 
            "Crude oil dry 120C<T<175C, 0.6<v<1.2": 0.000352, 
            "Crude oil salt 120C<T<175C, 0.6<v<1.2": 0.000705, 
            "Crude oil dry 120C<T<175C, 1.2<v": 0.000352, 
            "Crude oil salt 120C<T<175C, 1.2<v": 0.000705, 
            "Crude oil dry 175C<T<230C, v<0.6": 0.000705, 
            "Crude oil salt 175C<T<230C, v<0.6": 0.001057, 
            "Crude oil dry 175C<T<230C, 0.6<v<1.2": 0.000528, 
            "Crude oil salt 175C<T<230C, 0.6<v<1.2": 0.000881, 
            "Crude oil dry 175C<T<230C, 1.2<v": 0.000528, 
            "Crude oil salt 175C<T<230C, 1.2<v": 0.000881, 
            "Crude oil dry T>230C, v<0.6": 0.000881, 
            "Crude oil salt T>230C, v<0.6": 0.001233, 
            "Crude oil dry T>230C, 0.6<v<1.2": 0.000705, 
            "Crude oil salt T>230C, 0.6<v<1.2": 0.001057, 
            "Crude oil dry T>230C, 1.2<v": 0.000705, 
            "Crude oil salt T>230C, 1.2<v": 0.001057, 
            "Gasoline": 0.000352, 
            "Naphtha and light distillates": 0.00044, 
            "Kerosene": 0.00044, 
            "Light gas oil": 0.00044, 
            "Heavy gas oil": 0.00067, 
            "Heavy fuel oils": 0.00105}, 
        u"Refinery Asphalt": {
            "Vacuum tower bottoms": 0.001761, 
            "Atmosphere tower bottoms": 0.001233}, 
        u"Refinery Cracking and caking": {
            "Overhead vapors": 0.000352, 
            "Light cycle oil": 0.00044, 
            "Heavy cycle oil": 0.00061, 
            "Light coker gas oil": 0.00061, 
            "Heavy coker gas oil": 0.00079, 
            "Bottoms slurry oil": 0.000528, 
            "Light liquid products": 0.000176}, 
        u"Refinery Reforming": {
            "Reformer charge": 0.000264, 
            "Reformer effluent": 0.000264, 
            "Hydrocracker charge and effluent": 0.000352, 
            "Recycle gas": 0.000176, 
            "Overhead vapors": 0.000176, 
            "Liquid product >50 API": 0.000176, 
            "Liquid product 30-50 API": 0.000352}, 
        u"Refinery Light Ends": {
            "Overhead vapors and gases": 0.000176, 
            "Liquid products": 0.000176, 
            "Absorption oils": 0.00044, 
            "Alkylation trace acid streams": 0.000352, 
            "Reboiler streams": 0.00044}, 
        u"Refinery Lube oil": {
            "Feed stock": 0.000352, 
            "Solvent feed mix": 0.000352, 
            "Solvent": 0.000176, 
            "Extract": 0.000528, 
            "Rafftnate": 0.000176, 
            "Asphalt": 0.000881, 
            "Wax slurries": 0.000528, 
            "Refined lube oil": 0.000176}, 
        u"Refinery Visbreaker": {
            "Overhead vapor": 0.000528, 
            "Visbreaker bottoms": 0.001761}, 
        u"Refinery Naphtha Hydrotreater": {
            "Feed": 0.000528, 
            "Effluent": 0.000352, 
            "Naphfthas": 0.000352, 
            "Overhead vapors": 0.000264}, 
        u"Refinery Catalytic": {
            "Charge": 0.00079, 
            "Effluent": 0.000352, 
            "H.T. separator": 0.000352, 
            "Stripper charge": 0.000528, 
            "Liquid products": 0.000352}, 
        u"Refinery HF Alky": {
            "Alkylate, deprop. bottons, main fract": 0.000528, 
            "Other": 0.000352}}


    def __init__(self, parent=None):
        super(FoulingWidget, self).__init__(parent)
        layout = QtGui.QHBoxLayout(self)
        margins=layout.contentsMargins()
        margins.setTop(0)
        margins.setBottom(0)
        layout.setContentsMargins(margins)
        self.list=QtGui.QComboBox()
        self.list.addItem("")
        layout.addWidget(self.list)     
        self.value=Entrada_con_unidades(Fouling, decimales=6)
        self.value.valueChanged.connect(self.valueChanged.emit)
        layout.addWidget(self.value)     

        for tipo in sorted(self.Fouling_Factor):
            self.list.insertSeparator(self.list.count()+1)
            for componente in sorted(self.Fouling_Factor[tipo]):
                self.list.addItem(" - ".join([tipo, componente]))
        self.list.currentIndexChanged["QString"].connect(self.rellenar)

    def setValue(self, value):
        self.value.setValue(value)
        
    def rellenar(self, txt):
        if txt:
            tipo, componente=txt.split(" - ")
            value=self.Fouling_Factor[str(tipo)][str(componente)]
            self.value.setReadOnly(True)
            self.value.setValue(value)
            self.valueChanged.emit(value)
        else:
            self.value.setReadOnly(False)


class Dialog_Finned(QtGui.QDialog):
    def __init__(self, kwarg=None, parent=None):
        super(Dialog_Finned, self).__init__(parent=parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Specify tube finned characteristics"))
        layout = QtGui.QGridLayout(self)
        self.listTube=QtGui.QComboBox()
        self.listTube.addItem("")
        layout.addWidget(self.listTube,0,1,1,2)

        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,1,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material")),2,1)
        self.listMaterial=QtGui.QComboBox()
        self.listMaterial.addItem("")
        self.listMaterial.addItem(QtGui.QApplication.translate("pychemqt", "Carbon Steel"))
        layout.addWidget(self.listMaterial,2,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermal Conductivity")),3,1)
        self.kFin=Entrada_con_unidades(ThermalConductivity)
        layout.addWidget(self.kFin,3,2)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1,1,2)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Root diameter")),5,1)
        self.RootD=Entrada_con_unidades(Length, "PipeDiameter")
        layout.addWidget(self.RootD,5,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Fin Height")),6,1)
        self.hFin=Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.hFin,6,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Base Fin Thickness")),7,1)
        self.BaseThickness=Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.BaseThickness,7,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Top Fin Thickness")),8,1)
        self.TopThickness=Entrada_con_unidades(Length, "Thickness")
        layout.addWidget(self.TopThickness,8,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Number of fins")),9,1)
        self.Nfin=Entrada_con_unidades(float, textounidad="fins/m")
        layout.addWidget(self.Nfin,9,2)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,1,1,2)
        
        for tuberia in finnedTube_database:
            self.listTube.addItem("%s %s" % (tuberia[0], tuberia[1]))
        self.listTube.currentIndexChanged.connect(self.rellenarData)
        self.listTube.currentIndexChanged.connect(self.setDisabled)
        
        if kwarg:
            self.hFin.setValue(kwarg["hFin"])
            self.BaseThickness.setValue(kwarg["thicknessBaseFin"])
            self.TopThickness.setValue(kwarg["thicknessTopFin"])
            self.kFin.setValue(kwarg["kFin"])
            self.Nfin.setValue(kwarg["nFin"])
            self.RootD.setValue(kwarg["rootDoFin"])
            

    def rellenarData(self, ind):
        tuberia=finnedTube_database[ind-1]
        if tuberia[0]=="HPT":
            self.Nfin.setValue(int(tuberia[1][:2]))
            self.BaseThickness.setValue(tuberia[12]/1000.)
            self.TopThickness.setValue(tuberia[12]/1000.)
            self.RootD.setValue(tuberia[6]/1000.)
            self.hFin.setValue(tuberia[13]/1000.)

    def setDisabled(self, bool):
        self.RootD.setReadOnly(bool)
        self.BaseThickness.setReadOnly(bool)
        self.TopThickness.setReadOnly(bool)
        self.Nfin.setReadOnly(bool)
        self.hFin.setReadOnly(bool)
        
    def kwarg(self):
        kwarg={"hFin": self.hFin.value, 
                    "thicknessBaseFin": self.BaseThickness.value, 
                    "thicknessTopFin": self.TopThickness.value, 
                    "kFin": self.kFin.value, 
                    "nFin": self.Nfin.value, 
                    "rootDoFin": self.RootD.value}
        return kwarg


if __name__ == "__main__":
    import sys        
    from lib.corriente import Mezcla, Solid, Corriente
    from equipment.gas_solid import GravityChamber
    app = QtGui.QApplication(sys.argv)
#    distribucion=[[96.5, 0.02],
#                        [105, 0.05], 
#                        [110,  0.1],  
#                        [118, 0.15],
#                        [125, 0.25], 
#                        [130, 0.2], 
#                        [140, 0.15], 
#                        [150, 0.05], 
#                        [170, 0.03]]
#
#    solido=Solid([638], [5000], distribucion)
#    aire=Corriente(300, 1, 12100,  Mezcla([475], [1.]), solido)
#    dialogo = UI_equip(GravityChamber, aire)
    dialogo=Dialog_Finned()
    dialogo.show()
    sys.exit(app.exec_())
    
