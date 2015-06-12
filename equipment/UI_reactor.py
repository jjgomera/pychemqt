#!/usr/bin/python
# -*- coding: utf-8 -*-

###################################
###   Diálogo de definición de reactores UI_reactor  ###
###################################

import cPickle
from functools import partial
import os

from PyQt4 import QtCore, QtGui
from scipy.optimize import leastsq
from scipy import array, exp, log

from lib import unidades, reaction, plot
from lib.thread import Evaluate
from lib.config import getComponents
from UI.widgets import Status
from equipment.parents import UI_equip
from equipment.reactor import Reactor
from UI import UI_corriente, entrada_datos
from UI.widgets import Entrada_con_unidades, Tabla


class widgetReacciones(QtGui.QWidget):
    """Widget con la tabla de reacciones y los botones para modificar la lista de reacciones"""
    changed = QtCore.pyqtSignal()
    reacciones=[]
    reaccion=None
    activo=None
    ajuste=None
    
    def __init__(self, parent=None):
        super(widgetReacciones, self).__init__(parent)
        self.indices, self.nombres, M=getComponents()
        gridLayout = QtGui.QGridLayout(self)        

        self.TablaReacciones=Tabla(5, horizontalHeader=[QtGui.QApplication.translate("pychemqt", "Reaction"), u"ΔHr, %s" %unidades.MolarEnthalpy(None).text(), QtGui.QApplication.translate("pychemqt", "Type"), QtGui.QApplication.translate("pychemqt", "Phase"), QtGui.QApplication.translate("pychemqt", "Description")], dinamica=False, verticalHeader=True, orientacion=QtCore.Qt.AlignLeft)
        self.TablaReacciones.setMinimumWidth(500)
        self.TablaReacciones.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.TablaReacciones.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.TablaReacciones.horizontalHeader().setStretchLastSection(True)
        self.TablaReacciones.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.TablaReacciones.itemSelectionChanged.connect(self.actualizarBotones)
        gridLayout.addWidget(self.TablaReacciones,1,1,6,4)
        
        self.botonAbrir=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileOpen.png")), QtGui.QApplication.translate("pychemqt", "Open"))
        self.botonAbrir.clicked.connect(self.botonAbrirClicked)
        gridLayout.addWidget(self.botonAbrir,1,5)
        self.botonGuardar=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileSave.png")), QtGui.QApplication.translate("pychemqt", "Save"))
        self.botonGuardar.clicked.connect(self.botonGuardarClicked)
        self.botonGuardar.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
        self.botonGuardar.setEnabled(False)
        gridLayout.addWidget(self.botonGuardar,2,5)

        self.botonNew=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileNew.png")), QtGui.QApplication.translate("pychemqt", "New"))
        self.botonNew.clicked.connect(self.botonNewClicked)
        gridLayout.addWidget(self.botonNew,3,5)
        self.botonEdit=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editor.png")), QtGui.QApplication.translate("pychemqt", "Edit"))
        self.botonEdit.setEnabled(False)
        self.botonEdit.setCheckable(True)
        self.botonEdit.clicked.connect(self.botonEditClicked)
        gridLayout.addWidget(self.botonEdit,4,5)
        self.botonDelete=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editDelete.png")), QtGui.QApplication.translate("pychemqt", "Delete"))
        self.botonDelete.setEnabled(False)
        self.botonDelete.clicked.connect(self.botonDeleteClicked)
        gridLayout.addWidget(self.botonDelete,5,5)
        self.botonClear=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtGui.QApplication.translate("pychemqt", "Clear"))
        self.botonClear.clicked.connect(self.botonClearClicked)
        gridLayout.addWidget(self.botonClear,6,5)
        gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,1)


    def actualizarBotones(self, bool=True):
        self.botonEdit.setEnabled(bool)
        self.botonDelete.setEnabled(bool)

    def botonAbrirClicked(self):
        fname = unicode(QtGui.QFileDialog.getOpenFileName(self, QtGui.QApplication.translate("pychemqt", "Open reaction file"), "./", QtGui.QApplication.translate("pychemqt", "reaction file")+" (*.rec);;"+QtGui.QApplication.translate("pychemqt", "All files")+" (*.*)"))
        if fname:
            with open(fname, "r") as archivo:
                reacciones=cPickle.load(archivo)
            print reacciones
            self.reacciones=reacciones
            self.botonGuardar.setEnabled(True)
            for fila, reaccion in enumerate(reacciones):
                self.TablaReacciones.addRow()
                self.TablaReacciones.setValue(fila, 0, reaccion.text)
                self.TablaReacciones.setValue(fila, 1, "%0.4e" %reaccion.Hr.config(), QtCore.Qt.AlignRight)
                self.TablaReacciones.setValue(fila, 2, str(reaccion.tipo+1)+" - "+reaction.Reaction.TEXT_TYPE[reaccion.tipo])
                self.TablaReacciones.setValue(fila, 3, reaction.Reaction.TEXT_PHASE[reaccion.fase])
                self.TablaReacciones.item(fila, 4).setFlags(QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            for i in range(4):
                self.TablaReacciones.resizeColumnToContents(i)
        self.changed.emit()
    
    def botonGuardarClicked(self):
        fname = unicode(QtGui.QFileDialog.getSaveFileName(self, QtGui.QApplication.translate("pychemqt", "Save reaction to file"), "./", QtGui.QApplication.translate("pychemqt", "reaction file")+" (*.rec)"))
        if fname:
            if fname.split(".")[-1]!="rec":
                fname+=".rec"
            cPickle.dump(self.reacciones, open(fname, "w"))

    def botonNewClicked(self):
        dialog=UI_reacciones(parent=self)
        if dialog.exec_():
            pass

        
    def botonEditClicked(self, bool):
        if bool:
            indice=self.TablaReacciones.currentRow()
            reaccion=self.reacciones[indice]
            dialogo=UI_reacciones(reaccion, self)
            dialogo.exec_()
#            self.rellenar(self.reaccion)
#            self.activo=indice
        else:
            self.botonAddClicked(self.activo, False)
            self.reacciones[self.activo]=self.reaccion
            self.TablaReacciones.setCurrentCell(self.activo, 0)
            self.activo=-1
            self.changed.emit()

        self.botonNew.setEnabled(not bool)
        self.botonDelete.setEnabled(not bool)
        self.botonClear.setEnabled(not bool)
        self.botonAdd.setEnabled(not bool)
        self.botonAbrir.setEnabled(not bool)
        self.botonGuardar.setEnabled(not bool) 
            
        
    def botonDeleteClicked(self):
        indice=self.TablaReacciones.currentRow()
        self.TablaReacciones.removeRow(indice)
        del self.reacciones[indice]
        self.TablaReacciones.clearSelection()
        self.actualizarBotones(False)
        self.changed.emit()

    def botonClearClicked(self):
        if self.reacciones:
            self.reacciones=[]
            self.TablaReacciones.setRowCount(0)
            self.botonGuardar.setEnabled(False)

    def botonAddClicked(self, fila, add=True):
        if add:
            fila=self.TablaReacciones.rowCount()
            self.TablaReacciones.addRow()
        self.TablaReacciones.setValue(fila, 0, self.Formula.text())
        self.TablaReacciones.setValue(fila, 1, "%0.4e" %self.Hr.value.config(), QtCore.Qt.AlignRight)
        self.TablaReacciones.setValue(fila, 2, str(self.tipo.currentIndex()+1)+" - "+self.tipo.currentText())
        self.TablaReacciones.setValue(fila, 3, self.Fase.currentText())
        self.TablaReacciones.item(fila, 4).setFlags(QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
        for i in range(4):
            self.TablaReacciones.resizeColumnToContents(i)
        self.reacciones.insert(fila, self.reaccion)
        self.botonGuardar.setEnabled(True)
        self.changed.emit()
    
class UI_reacciones(QtGui.QDialog):
    reaction=reaction.Reaction()
    def __init__(self, reaccion=None, parent=None):
        super(UI_reacciones, self).__init__(parent)
        self.evaluate=Evaluate()
        self.evaluate.finished.connect(self.rellenar)
        self.indices, self.nombres, M=getComponents()
        gridLayout = QtGui.QGridLayout(self)        

        lyt=QtGui.QHBoxLayout()
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Key component")))        
        self.key=QtGui.QComboBox()
        for i, nombre in enumerate(self.nombres):
            self.key.addItem("%i - %s" %(i+1, nombre))
        self.key.currentIndexChanged.connect(partial(self.changeParams, "key"))
        lyt.addWidget(self.key)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed))
        gridLayout.addLayout(lyt,1,1,1,5)

        lyt=QtGui.QHBoxLayout()
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Phase")))        
        self.fase=QtGui.QComboBox()
        for txt in reaction.Reaction.TEXT_PHASE:
            self.fase.addItem(txt)
        self.fase.currentIndexChanged.connect(partial(self.changeParams, "fase"))
        lyt.addWidget(self.fase)
        self.Formula=QtGui.QLabel()
        self.Formula.setAlignment(QtCore.Qt.AlignCenter)
        self.Formula.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed)
        lyt.addWidget(self.Formula)
        gridLayout.addLayout(lyt,2,1,1,5)
        
        lyt=QtGui.QVBoxLayout()
        title=self.nombres[:]
        title.append("")
        self.Estequiometria=Tabla(1, verticalHeaderLabels=title, horizontalHeader=[QtGui.QApplication.translate("pychemqt", "Coefficients")], filas=len(self.indices))
        self.Estequiometria.setFixedHeight(22*len(self.indices)+22+4+22)
        lyt.addWidget(self.Estequiometria)
        self.Estequiometria.addRow()
        brush=QtGui.QBrush(QtGui.QColor("#eaeaea"))
        self.Estequiometria.item(len(self.indices), 0).setBackground(brush)
        self.Estequiometria.item(len(self.indices), 0).setFlags(QtCore.Qt.NoItemFlags)
        self.Estequiometria.cellChanged.connect(self.reaccionCambiada)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Expanding))
        gridLayout.addLayout(lyt,3,1,1,2)
        
        lyt=QtGui.QGridLayout()
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,1)
        self.formula=QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Use name in formula"))
        self.formula.toggled.connect(partial(self.changeParams, "formula"))
        lyt.addWidget(self.formula,1,2,1,2)
        self.customHr=QtGui.QCheckBox(u"ΔHr "+QtGui.QApplication.translate("pychemqt", "user specified"))
        self.customHr.toggled.connect(self.changeHr)
        lyt.addWidget(self.customHr,2,2,1,2)
        lyt.addWidget(QtGui.QLabel(u"ΔHr<sup>o</sup>"),3,2)
        self.Hr=Entrada_con_unidades(unidades.MolarEnthalpy, readOnly=True)
        self.Hr.valueChanged.connect(partial(self.changeParams, "Hr"))
        lyt.addWidget(self.Hr,3,3)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Expanding))
        gridLayout.addLayout(lyt,3,3,1,2)
        
        gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,2)

        lyt=QtGui.QHBoxLayout()
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Type")))
        self.tipo=QtGui.QComboBox()
        for txt in reaction.Reaction.TEXT_TYPE:
            self.tipo.addItem(txt)
        self.tipo.currentIndexChanged.connect(partial(self.changeParams, "tipo"))
        lyt.addWidget(self.tipo)
        lyt.addItem(QtGui.QSpacerItem(20,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed))
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Concentration")))
        self.base=QtGui.QComboBox()
        for txt in reaction.Reaction.TEXT_BASE:
            self.base.addItem(txt)
        self.base.currentIndexChanged.connect(partial(self.changeParams, "base"))
        lyt.addWidget(self.base)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed))
        gridLayout.addLayout(lyt,5,1,1,5)

        self.stacked = QtGui.QStackedWidget()
        self.tipo.currentIndexChanged.connect(self.stacked.setCurrentIndex)
        gridLayout.addWidget(self.stacked,6,1,1,5)
        gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),7,1,1,5)

        widget=QtGui.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtGui.QGridLayout(widget)
        lyt.addWidget(QtGui.QLabel("<h3>"+QtGui.QApplication.translate("pychemqt", "Estequiometric reaction")+"</h3>"),1,1,1,4)
        self.Conversion=Tabla(1, verticalHeaderModel="C", filas=3)
        self.Conversion.setConnected()
        self.Conversion.setFixedWidth(100)
        lyt.addWidget(self.Conversion,2,1,3,1)
        label=QtGui.QLabel()
        label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/reaction_conversion.png"))
        lyt.addWidget(label,2,2,1,3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Temperature unit")),3,2)
        self.unidadesTemperatura=QtGui.QComboBox()
        for i in unidades.Temperature.__text__:
            self.unidadesTemperatura.addItem(i)
        lyt.addWidget(self.unidadesTemperatura,3,3)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),4,4)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),5,1,1,5)


        widget=QtGui.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtGui.QGridLayout(widget)
        self.check_KFijo=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Fixed"))
        self.check_KFijo.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KFijo,1,1,1,2)
        lyt.addWidget(QtGui.QLabel("K<sub>eq</sub>"),1,3)
        self.Keq=Entrada_con_unidades(float)
        lyt.addWidget(self.Keq,1,4)
        label=QtGui.QLabel()
        label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/reaction_equilibrium.png"))
        lyt.addWidget(label,1,5,1,4)

        self.check_KEq=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Equation"))
        self.check_KEq.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KEq,2,1,1,2)
        self.check_KTabla=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Table"))
        self.check_KTabla.toggled.connect(self.KeqChanged)
        lyt.addWidget(self.check_KTabla,2,5,1,2)
        self.KEq_Dat=Tabla(1, verticalHeaderLabels=["A", "B", "C", "D", "E", "F", "G", "H"], filas=8)
        self.KEq_Dat.setFixedHeight(22*8+4)
        self.KEq_Dat.setFixedWidth(120)
        lyt.addWidget(self.KEq_Dat,3,3,1,2)
        self.KEq_Tab=Tabla(4, horizontalHeader=["T, K", "Keq", "Kcalc", "%Error"], verticalHeader=False, columnReadOnly=[False, False, True, True])
        self.KEq_Tab.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        self.KEq_Tab.setFixedWidth(400)
        self.KEq_Tab.setConnected()
        self.KEq_Tab.rowFinished.connect(self.Regresion)
        self.KEq_Tab.setAlternatingRowColors(False)
        lyt.addWidget(self.KEq_Tab,3,5,1,4)
        lyt.addWidget(QtGui.QLabel(u"r²"),4,5)
        self.r2=Entrada_con_unidades(float, readOnly=True)
        lyt.addWidget(self.r2,4,6)
        self.botonTablaPlot=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/plot.png")), QtGui.QApplication.translate("pychemqt", "Plot"))
        self.botonTablaPlot.clicked.connect(self.Plot)
        lyt.addWidget(self.botonTablaPlot,4,7)
        self.botonTablaClear=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtGui.QApplication.translate("pychemqt", "Clear"))
        self.botonTablaClear.clicked.connect(self.KEq_Tab.clear)
        lyt.addWidget(self.botonTablaClear,4,8)
        label=QtGui.QLabel()
        label.setPixmap(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equation/reaction_equilibrium2.png"))
        label.setAlignment(QtCore.Qt.AlignCenter)
        lyt.addWidget(label,5,1,1,8)

        
        self.checkGibbs=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "From Gibbs free energy minimization"))
        lyt.addWidget(self.checkGibbs,6,1,1,4)
        
        self.check_KFijo.setChecked(True)
        

        widget=QtGui.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtGui.QGridLayout(widget)

        widget=QtGui.QWidget()
        self.stacked.addWidget(widget)
        lyt=QtGui.QGridLayout(widget)

        self.status=Status()
        gridLayout.addWidget(self.status, 10,1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox,10,2,1,4)
        
        if reaccion:
            self.setReaction(reaccion)
            
    def changeParams(self, parametro, valor):
        self.calculo(**{parametro: valor})

    def calculo(self, **kwargs):
        self.status.setState(4)
        self.evaluate.start(self.reaction, kwargs)

    def changeHr(self, bool):
        self.Hr.setReadOnly(not bool)
        self.changeParams("customHr", bool)

    def reaccionCambiada(self):
        kwargs={"componentes": self.indices, 
                    "coeficientes": self.Estequiometria.getColumn(0)[:-1]}
        self.calculo(**kwargs)
        
    def setReaction(self, reaction):
        self.reaction=reaction
        self.rellenar()

#        if self.Estequiometria.getValue(0, self.Base.currentIndex()):
#            reaccion=reaction.Reaction(self.indices, self.Estequiometria.getColumn(0), base=self.Base.currentIndex(), estequiometria=[0, 0, 0.5], formulas=self.checkFormula.isChecked(), calor=self.checkCalorEspecificado.isChecked(), Hr=self.Hr.value, tipo=self.tipo.currentIndex(), conversion=self.Conversion.getColumn(0)[-1::-1])
#            self.Balance.setValue(reaccion.error)
#            if reaccion.state:
#                self.Formula.setText(reaccion._txt(self.checkFormula.isChecked()))
#                self.Hr.setValue(reaccion.Hr)
#            else:
#                self.Formula.clear()
#                self.Hr.clear()
#            self.botonAdd.setEnabled(reaccion.state and not self.botonEdit.isChecked())
#            self.reaccion=reaccion

    def rellenar(self):
        self.blockSignals(True)
        for variable in self.reaction.kwargsValue:
            self.__getattribute__(variable).setValue(self.reaction.kwargs[variable])
        for combo in self.reaction.kwargsList:
            self.__getattribute__(combo).setCurrentIndex(self.reaction.kwargs[combo])
        for check in self.reaction.kwargsCheck:
            self.__getattribute__(check).setChecked(self.reaction.kwargs[check])
            
        self.Estequiometria.setColumn(0, self.reaction.kwargs["coeficientes"])
#        self.Conversion.setColumn(0, self.reaction.estequiometria[-1::-1])            
        self.blockSignals(False)

        self.status.setState(self.reaction.status, self.reaction.msg)
        self.Estequiometria.item(len(self.indices), 0).setText(str(self.reaction.error))
        if self.reaction.status:
            self.Formula.setText(self.reaction._txt())
            self.Hr.setValue(self.reaction.Hr)

    def KeqChanged(self):
        self.Keq.setReadOnly(not self.check_KFijo.isChecked())
        self.KEq_Dat.setEnabled(self.check_KEq.isChecked())
        self.KEq_Tab.setEnabled(self.check_KTabla.isChecked())
        self.botonTablaClear.setEnabled(self.check_KTabla.isChecked())
        self.botonTablaPlot.setEnabled(self.check_KTabla.isChecked())

    def Regresion(self):
        t=array(self.KEq_Tab.getColumn(0)[:-1])
        k=array(self.KEq_Tab.getColumn(1)[:-1])
        if len(t)>=4:
            if 4<=len(t)<8:
                inicio=r_[0, 0, 0, 0]
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
                resto=lambda par, T, k: k-f(par, T)
            else:
                inicio=r_[0, 0, 0, 0, 0, 0, 0, 0]
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)
                resto=lambda par, T, k: k-f(par, T)
        
            ajuste=leastsq(resto,inicio,args=(t, k))
            kcalc=f(ajuste[0], t)
            error=(k-kcalc)/k*100
            self.KEq_Dat.setColumn(0, ajuste[0])
            self.KEq_Tab.setColumn(2, kcalc)
            self.KEq_Tab.setColumn(3, error)
            
            if ajuste[1] in [1, 2, 3, 4]:
                self.ajuste=ajuste[0]
            
    def Plot(self):
        if self.ajuste!=None:
            t=array(self.KEq_Tab.getColumn(0)[:-1])
            k=array(self.KEq_Tab.getColumn(1)[:-1])
            if 4<=len(t)<8:
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T)
            else:
                f=lambda par, T: exp(par[0]+par[1]/T+par[2]*log(T)+par[3]*T+par[4]*T**2+par[5]*T**3+par[6]*T**4+par[7]*T**5)
                
            grafico=plot.Plot()
            grafico.data(t, k, 'ro', t, f(self.ajuste, t))
            grafico.exec_()



class UI_equipment(UI_equip):
    """Diálogo de definición de reactores"""
    profile_T=None
    Equipment=Reactor()
    def __init__(self, equipment=None, parent=None):
        """entrada: Parametro opcional de clase corriente que indica la corriente de entrada en el reactor"""
        super(UI_equipment, self).__init__(Reactor, entrada=False, salida=False, parent=parent)
        
        #Pestaña reacciones
        self.Reacciones= widgetReacciones()
        self.Reacciones.changed.connect(self.calculo)
        self.tabWidget.insertTab(1, self.Reacciones, QtGui.QApplication.translate("pychemqt", "Reactions"))
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure")),2,0,1,1)
        self.P=Entrada_con_unidades(unidades.Pressure)
        self.P.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.P,2,1,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure drop")),3,0,1,1)
        self.DeltaP=Entrada_con_unidades(unidades.Pressure)
        self.DeltaP.valueChanged.connect(self.calculo)
        gridLayout_Calculo.addWidget(self.DeltaP,3,1,1,1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,0,1,5)
        lyt=QtGui.QHBoxLayout()
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Type")))
        self.tipo=QtGui.QComboBox()
        self.tipo.addItem(QtGui.QApplication.translate("pychemqt", "CSTR, continuous stirred-tank"))
        self.tipo.addItem(QtGui.QApplication.translate("pychemqt", "PFR, plug flow"))
        self.tipo.currentIndexChanged.connect(self.tipoCambiado)
        lyt.addWidget(self.tipo)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed))
        gridLayout_Calculo.addLayout(lyt,5,0,1,5)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),6,0,1,5)

        groupbox=QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Thermal mode"))
        layout=QtGui.QGridLayout(groupbox)
        self.checkAdiabatico=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Adiabatic"))
        self.checkAdiabatico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkAdiabatico, 1, 1, 1, 1)
        self.checkIsotermico=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Isothermal"))
        self.checkIsotermico.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIsotermico, 2, 1, 1, 1)
        self.checkFlux=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Heat duty"))
        self.checkFlux.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkFlux, 3, 1, 1, 1)
        self.checkIntercambio=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Heat transfer"))
        self.checkIntercambio.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkIntercambio, 4, 1, 1, 1)
        self.checkPerfil=QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "PFR temperature profile"))
        self.checkPerfil.toggled.connect(self.heatChanged)
        layout.addWidget(self.checkPerfil, 5, 1, 1, 1)
        self.T=Entrada_con_unidades(unidades.Temperature)
        self.T.valueChanged.connect(self.calculo)
        layout.addWidget(self.T, 2, 2, 1, 2)
        self.Q=Entrada_con_unidades(unidades.Power)
        self.Q.valueChanged.connect(self.calculo)
        layout.addWidget(self.Q, 3, 2, 1, 2)
        self.T_ext=Entrada_con_unidades(unidades.Temperature)
        self.T_ext.valueChanged.connect(self.calculo)
        layout.addWidget(self.T_ext,4,2,1,2)
        self.U=Entrada_con_unidades(unidades.HeatTransfCoef)
        self.U.valueChanged.connect(self.calculo)
        layout.addWidget(self.U,4,4)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow")),4,5)
        self.direccion=QtGui.QComboBox()
        self.direccion.addItem(QtGui.QApplication.translate("pychemqt", "Countercurrent"))
        self.direccion.addItem(QtGui.QApplication.translate("pychemqt", "Cocurrent"))
        layout.addWidget(self.direccion,4,6)
        self.botonPerfil=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/table.png")), QtGui.QApplication.translate("pychemqt", "Add Profile"))
        self.botonPerfil.clicked.connect(self.editorPerfil)
        layout.addWidget(self.botonPerfil,5,2,1,1)
        gridLayout_Calculo.addWidget(groupbox, 7, 0, 1, 5)
        
        self.groupBox_Diseno= QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Design"))
        gridLayout_Calculo.addWidget(self.groupBox_Diseno,8,0,1,5)
        lyt = QtGui.QGridLayout(self.groupBox_Diseno)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Mode")),1,1)
        self.modo=QtGui.QComboBox()
        self.modo.addItem(QtGui.QApplication.translate("pychemqt", "Rating: calculate conversión"))
        self.modo.addItem(QtGui.QApplication.translate("pychemqt", "Design, calculate volumen"))
        self.modo.currentIndexChanged.connect(self.calculo)
        lyt.addWidget(self.modo,1,2,1,3)
        lyt.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,5)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reactor Volume")),2,1)
        self.V=Entrada_con_unidades(unidades.Volume, "VolLiq")
        lyt.addWidget(self.V,2,2)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Key Component")),3,1)
        self.key=QtGui.QComboBox()
#        for i, nombre in enumerate(self.nombres):
#            self.key.addItem("%i - %s" %(i+1, nombre))
        lyt.addWidget(self.key,3,2)
        lyt.addItem(QtGui.QSpacerItem(20,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,3)
        lyt.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Conversion")),3,4)
        self.conversion=Entrada_con_unidades(float, max=1)
        lyt.addWidget(self.conversion,3,5)
        
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),10,0,1,5)
        
        groupBox_Calculo = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(groupBox_Calculo,11,0,1,5)
        gridLayout_1 = QtGui.QGridLayout(groupBox_Calculo)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T output")),0,1,1,1)
        self.TCalc=Entrada_con_unidades(unidades.Temperature, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.TCalc,0,2,1,1)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heat")),1,1,1,1)
        self.HeatCalc=Entrada_con_unidades(unidades.Power, retornar=False, readOnly=True)
        gridLayout_1.addWidget(self.HeatCalc,1,2,1,1)

        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Conversion")),0,4)
        self.conversionCalc=Entrada_con_unidades(float, readOnly=True)
        gridLayout_1.addWidget(self.conversionCalc,0,5)
        gridLayout_1.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reactor Volume")),1,4)
        self.VCalc=Entrada_con_unidades(unidades.Volume, "VolLiq", readOnly=True)
        gridLayout_1.addWidget(self.VCalc,1,5)

        self.checkAdiabatico.setChecked(True)
        self.tipoCambiado(0)
        

    def heatChanged(self):
        self.T.setReadOnly(not self.checkIsotermico.isChecked())
        self.Q.setReadOnly(not self.checkFlux.isChecked())
        self.T_ext.setReadOnly(not self.checkIntercambio.isChecked())
        self.U.setReadOnly(not self.checkIntercambio.isChecked())
        self.direccion.setEnabled(self.checkIntercambio.isChecked() and self.tipo.currentIndex())
        self.botonPerfil.setEnabled(self.checkPerfil.isChecked())
    
    def tipoCambiado(self, ind):
        self.checkPerfil.setEnabled(ind)
        self.direccion.setEnabled(self.checkIntercambio.isChecked() and ind)
        self.calculo()
            
        
    def editorPerfil(self):
        dialog=entrada_datos.Entrada_Datos(data=self.profile_T, title=QtGui.QApplication.translate("pychemqt", "Temperature profile"), horizontalHeader=["x", "T, "+unidades.Temperature(None).text()])
        if dialog.exec_():
            self.profile_T=dialog.data
        
#    def cambiar_entrada(self, corriente):
#        self.entrada=corriente
#        self.calculo()
#        
#    def todos_datos(self):
#        if self.checkAdiabatico.isChecked():
#            heat=True
#        elif self.checkIsotermico.isChecked():
#            heat=self.T.value
#        elif self.checkFlux.isChecked():
#            heat=self.Q.value
#        else: heat=self.T_ext.value and self.U.value
#        return  self.Reacciones.reacciones and heat

#    def calculo(self):
#        if self.todos_datos():
#            self.status.setState(4)
#            if self.checkAdiabatico.isChecked():
#                thermal=0
#            elif self.checkIsotermico.isChecked():
#                thermal=1
#            elif self.checkFlux.isChecked():
#                thermal=2
#            else: thermal=3
#        
#            if self.P.value:
#                P=self.P.value.atm
#            else:
#                P=None
#
#            self.Equipment(entrada=self.entrada, thermal=thermal, reaccion=self.Reacciones.reacciones, P=P, T=self.T.value, Q=self.Q.value, Text=self.T_ext.value, U=self.U.value)
#
#            self.rellenoSalida()
#            self.status.setState(1)
#
#
#    def rellenoSalida(self):
#        self.TCalc.setValue(self.Equipment.Salida.T)
#        self.HeatCalc.setValue(self.Equipment.Heat)
#        self.Salida.rellenar(self.Equipment.Salida)
        

if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente, Mezcla, Solid
    from numpy import r_
    app = QtGui.QApplication(sys.argv)
#    dialogo = UI_reacciones()
    dialogo = UI_equipment()
    dialogo.show()
    sys.exit(app.exec_())
