#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################
###    Diálogo de definición de tuberías, UI_pipe    ###
##################################

from functools import partial
import os

from PyQt4 import QtCore, QtGui
from scipy import pi

from lib.config import Preferences
from lib.utilities import representacion
from lib.unidades import Length, Temperature, HeatTransfCoef, Pressure, Power, Speed, Currency
from lib.table import Material_Tuberia, Accesorios_Tuberia, Fitting_icon
from lib.physics import K_contraction, K_enlargement, K_flush, K_MitreBend, Ft, K_longBend
from tools.costIndex import CostData
from equipment.parents import UI_equip
from equipment.pipe import Pipe
from UI import UI_corriente
from UI.delegate import SpinEditor, CellEditor
from UI.widgets import Entrada_con_unidades


class Dialog(QtGui.QDialog):
    """Diálogo genérico de definición de elementos especiales"""
    def __init__(self, tipo=0, titulo="", icon=None, textos=["", ""], parent=None):
        """tipo: define el tipo de dialogo:
            0 - Contracción, expansión
            1 - Codo segmentado de angulo arbitrario
            2 - Codo largo 
            3 - Entrada redondeada"""
        super(Dialog, self).__init__(parent)
        self.tipo=tipo
        self.setWindowTitle(titulo)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(icon)))
        gridLayout = QtGui.QGridLayout(self)
        self.label = QtGui.QLabel(textos[0])
        self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.D1=Entrada_con_unidades(Length)
        self.D1.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.D1, 0, 1, 1, 2)
        self.label_2 = QtGui.QLabel(textos[1])
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.D2=Entrada_con_unidades(Length)
        self.D2.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.D2, 1, 1, 1, 2)
        self.checkBox = QtGui.QCheckBox()
        self.checkBox.setText(QtGui.QApplication.translate("pychemqt", "Gradual"))
        self.checkBox.toggled.connect(self.check_toggled)
        self.checkBox.toggled.connect(self.calculaK)
        gridLayout.addWidget(self.checkBox, 2, 0, 1, 1)
        self.label_3 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Angle"))
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        gridLayout.addWidget(self.label_3, 2, 1, 1, 1)
        self.angulo = QtGui.QSpinBox()
        self.angulo.setFixedWidth(50)
        self.angulo.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.angulo.setRange(1, 90)
        gridLayout.addWidget(self.angulo, 2, 2, 1, 1)
        self.regla = QtGui.QSlider()
        self.regla.setMaximum(90)
        self.regla.setOrientation(QtCore.Qt.Horizontal)
        self.regla.setTickPosition(QtGui.QSlider.TicksAbove)
        self.regla.setTickInterval(5)
        self.regla.valueChanged.connect(self.angulo.setValue)
        self.regla.valueChanged.connect(self.calculaK)
        self.angulo.valueChanged.connect(self.regla.setValue)
        self.angulo.valueChanged.connect(self.calculaK)
        gridLayout.addWidget(self.regla, 3, 0, 1, 3)
        self.K_text=QtGui.QLabel("K=")
        self.K_text.setAlignment(QtCore.Qt.AlignCenter|QtCore.Qt.AlignVCenter)
        gridLayout.addWidget(self.K_text, 4, 0, 1, 3)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.setCenterButtons(True)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        gridLayout.addWidget(self.buttonBox, 5, 0, 1, 3)
        
        self.check_toggled(False)
        
        if tipo==1:
            self.label_2.setVisible(False)
            self.D2.setVisible(False)
            self.checkBox.setVisible(False)
            self.check_toggled(True)
        
        if tipo>=2:
            self.checkBox.setVisible(False)
            self.label_3.setVisible(False)
            self.angulo.setVisible(False)
            self.regla.setVisible(False)

    def check_toggled(self, bool):
        self.angulo.setEnabled(bool)
        self.regla.setEnabled(bool)
    
    def todos_datos(self):
        if self.tipo==1:
            return True
        else:
            return self.D1.value and self.D2.value
            
    def calculaK(self):
        if self.tipo==1:
            tita=self.angulo.value()
            f=Ft(self.D1.value.mm)
            K=K_MitreBend(tita)*f
            self.K_text.setText("K = %0.3f" %K)
            self.K=K
        else:
            if self.D1.value and self.D2.value:
                beta=self.D2.value/self.D1.value
                if self.checkBox.isChecked():
                    tita=self.angulo.value()
                else: tita=180
                if self.tipo==0:
                    if beta>1:
                        K=K_enlargement(tita, beta)
                    else:
                        K=K_contraction(tita, beta)
                elif self.tipo==2:
                    f=Ft(self.D2.value.mm)
                    rD=self.D1.value/self.D2.value
                    K=K_longBend(rD)*f
                else:
                    K=K_flush(1./beta)
                self.K_text.setText("K = %0.3f" %K)
                self.K=K
        
class Catalogo_Materiales(QtGui.QWidget):
    """Diálogo del catalogo de materiales y tamaños de tuberías"""
    valueChanged = QtCore.pyqtSignal(list)
    def __init__(self, parent=None):
        super(Catalogo_Materiales, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        
        self.TablaMaterial=QtGui.QTableWidget()
        self.TablaMaterial.verticalHeader().hide()
        self.TablaMaterial.setRowCount(1)
        self.TablaMaterial.setColumnCount(3)
        self.TablaMaterial.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Material")))
        self.TablaMaterial.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Class")))
        self.TablaMaterial.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Roughness") +", mm"))
        self.TablaMaterial.setSelectionMode(QtGui.QAbstractItemView.NoSelection)
        self.TablaMaterial.setColumnWidth(0, 200)
        self.TablaMaterial.setColumnWidth(1, 80)
        self.TablaMaterial.setColumnWidth(2, 100)
        self.TablaMaterial.setRowHeight(0, 20)
        self.TablaMaterial.setFixedHeight(47)
        self.TablaMaterial.setFixedWidth(384)
        self.TablaMaterial.setItemDelegateForColumn(2, CellEditor(self))
        self.TablaMaterial.itemChanged.connect(self.emitirSignal)
        gridLayout.addWidget(self.TablaMaterial, 1, 0, 1, 2)

        self.TablaDiametro=QtGui.QTableWidget()
        self.TablaDiametro.verticalHeader().hide()
        self.TablaDiametro.setRowCount(1)
        self.TablaDiametro.setColumnCount(7)
        self.TablaDiametro.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Nº", None, QtGui.QApplication.UnicodeUTF8)))
        self.TablaDiametro.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "D int")+", mm" ))
        self.TablaDiametro.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Thickness") +", mm"))
        self.TablaDiametro.setHorizontalHeaderItem(3, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "D ext")+", mm" ))
        self.TablaDiametro.setHorizontalHeaderItem(4, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Weight")+", kg/m" ))
        self.TablaDiametro.setHorizontalHeaderItem(5, QtGui.QTableWidgetItem("V"+u", m³/100m" ))
        self.TablaDiametro.setHorizontalHeaderItem(6, QtGui.QTableWidgetItem("S"+u", m²/100m" ))
        self.TablaDiametro.setSelectionMode(QtGui.QAbstractItemView.NoSelection)
        self.TablaDiametro.horizontalHeader().setStretchLastSection(True)
        self.TablaDiametro.setRowHeight(0, 20)
        self.TablaDiametro.setFixedHeight(47)
        self.TablaDiametro.setFixedWidth(704)
        self.TablaDiametro.setItemDelegateForColumn(1, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(2, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(3, CellEditor(self))
        self.TablaDiametro.setItemDelegateForColumn(4, CellEditor(self))
        self.TablaDiametro.itemChanged.connect(self.emitirSignal)
        gridLayout.addWidget(self.TablaDiametro, 2, 0, 1, 2)
                
        self.botonClear=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtGui.QApplication.translate("pychemqt", "Reset"))
        self.botonClear.clicked.connect(self.clear)
        gridLayout.addWidget(self.botonClear, 1, 2, 1, 1)
        self.botonEdit=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/filenew.png")), QtGui.QApplication.translate("pychemqt", "Edit"))
        self.botonEdit.clicked.connect(self.edit)
        self.botonEdit.setCheckable(True)
        gridLayout.addWidget(self.botonEdit, 2, 2, 1, 1)
        
        self.line = QtGui.QFrame()
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        gridLayout.addWidget(self.line, 3, 0, 1, 3)
        
        label=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Material:"))
        label.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        gridLayout.addWidget(label, 4, 0, 1, 1)
        self.ComboMaterial=QtGui.QComboBox()
        gridLayout.addWidget(self.ComboMaterial, 4, 1, 1, 1)
        
        self.TablaDiametro2=QtGui.QTableWidget()
        self.TablaDiametro2.verticalHeader().hide()
        self.TablaDiametro2.setColumnCount(7)
        self.TablaDiametro2.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Nº", None, QtGui.QApplication.UnicodeUTF8)))
        self.TablaDiametro2.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "D int")+", mm" ))
        self.TablaDiametro2.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Thickness") +", mm"))
        self.TablaDiametro2.setHorizontalHeaderItem(3, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "D ext")+", mm" ))
        self.TablaDiametro2.setHorizontalHeaderItem(4, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Weight")+", kg/m" ))
        self.TablaDiametro2.setHorizontalHeaderItem(5, QtGui.QTableWidgetItem("V"+u", m³/100m" ))
        self.TablaDiametro2.setHorizontalHeaderItem(6, QtGui.QTableWidgetItem("S"+u", m²/100m" ))
        self.TablaDiametro2.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.TablaDiametro2.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.TablaDiametro2.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.TablaDiametro2.setFixedWidth(720)
        gridLayout.addWidget(self.TablaDiametro2, 5, 0, 1, 2)
        
        self.botonSelect=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/arrow-up-double.png")), QtGui.QApplication.translate("pychemqt", "Select"))
        self.botonSelect.setDisabled(True)
        self.TablaDiametro2.currentCellChanged.connect(self.selectedChanged)
        self.TablaDiametro2.cellDoubleClicked.connect(self.add)
        self.botonSelect.clicked.connect(self.add)
        gridLayout.addWidget(self.botonSelect, 4, 2, 1, 1)
        gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Preferred),6,0,1,3)

        materiales=[]
        diametros=[]
        for material in Material_Tuberia:
            txt=material[0:2]
            if txt not in materiales:
                materiales.append(material[0:2])
                diametros.append([])
            indice=materiales.index(txt)
            diametros[indice].append(material[2:])
        for material in materiales:
            self.ComboMaterial.addItem(material[0]+" "+material[1])
        self.ComboMaterial.currentIndexChanged.connect(self.rellenarDiametros)
        self.diametros=diametros
        self.materiales=materiales
        self.rellenarDiametros(0)
        
        self.material=None
        self.diametro=None
        self.rellenarActivo()
        self.edit(False)
        
    def rellenarActivo(self):
        self.TablaDiametro.setItem(0, 0, QtGui.QTableWidgetItem(""))
        for i in range(1, 7):
            self.TablaDiametro.setItem(0, i, QtGui.QTableWidgetItem(""))
            self.TablaDiametro.item(0, i).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
        for i in range(3):
            self.TablaMaterial.setItem(0, i, QtGui.QTableWidgetItem(""))
        self.TablaMaterial.item(0, 2).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)

    def selectedChanged(self):
        self.botonSelect.setEnabled(True)

    def clear(self):
        self.material=None
        self.diametro=None
        for i in range(7):
            self.TablaDiametro.item(0, i).setText("")
        for i in range(3):
            self.TablaMaterial.item(0, i).setText("")

    def add(self, i=None, j=None, load=True):
        if load:
            i=self.ComboMaterial.currentIndex()
            j=self.TablaDiametro2.currentRow()
        self.blockSignals(True)
        self.TablaMaterial.item(0, 0).setText(self.materiales[i][0])
        self.TablaMaterial.item(0, 1).setText(self.materiales[i][1])
        self.TablaMaterial.item(0, 2).setText(str(self.diametros[i][j][0]))
        for columna in range(7):
            self.TablaDiametro.item(0, columna).setText(self.TablaDiametro2.item(j, columna).text())
        self.material=self.materiales[i]
        self.diametro=self.diametros[j]
        self.blockSignals(False)
        self.emitirSignal()

    def addCustom(self, lista):
        self.blockSignals(True)
        self.TablaMaterial.item(0, 0).setText(lista[0])
        self.TablaMaterial.item(0, 1).setText(lista[1])
        self.TablaMaterial.item(0, 2).setText(str(lista[2]))
        for columna in range(7):
            self.TablaDiametro.item(0, columna).setText(str(lista[columna+3]))
        self.blockSignals(False)

    def edit(self, bool):
        self.blockSignals(True)
        inactivo=QtGui.QColor(Preferences.get("General", 'Color_ReadOnly'))
        activo=QtGui.QColor(Preferences.get("General", 'Color_Resaltado'))
        if bool:
            flags=QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable
            color=activo            
        else:
            flags=QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable
            color=inactivo
            
        for i in [0, 1, 2]:
            self.TablaMaterial.item(0, i).setFlags(flags)
            self.TablaMaterial.item(0, i).setBackground(color)
        for i in range(0, 5):
            self.TablaDiametro.item(0, i).setFlags(flags)
            self.TablaDiametro.item(0, i).setBackground(color)
            
        self.TablaDiametro.item(0, 5).setBackground(inactivo)
        self.TablaDiametro.item(0, 6).setBackground(inactivo)
        self.TablaDiametro.item(0, 5).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
        self.TablaDiametro.item(0, 6).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
        self.blockSignals(False)
        
    def rellenarDiametros(self, indice):
        filas=len(self.diametros[indice])
        self.TablaDiametro2.setRowCount(filas)
        for i in range(filas):
            di=self.diametros[indice][i][4]-self.diametros[indice][i][3]*2
            V=pi/4*di**2*100/1e6
            S=pi*self.diametros[indice][i][4]*100/1e3
            self.TablaDiametro2.setRowHeight(i, 20)
            self.TablaDiametro2.setItem(i, 0, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 0).setText(self.diametros[indice][i][2])
            self.TablaDiametro2.setItem(i, 1, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 1).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 1).setText(str(round(di, 2)))
            self.TablaDiametro2.setItem(i, 2, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 2).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 2).setText(str(self.diametros[indice][i][3]))
            self.TablaDiametro2.setItem(i, 3, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 3).setText(str(self.diametros[indice][i][4]))
            self.TablaDiametro2.setItem(i, 4, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 4).setText(str(round(self.diametros[indice][i][5], 3)))
            self.TablaDiametro2.setItem(i, 5, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 5).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 5).setText(str(round(V, 3)))
            self.TablaDiametro2.setItem(i, 6, QtGui.QTableWidgetItem(""))
            self.TablaDiametro2.item(i, 6).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.TablaDiametro2.item(i, 6).setText(str(round(S, 3)))
                
    def getMaterial(self):
        try:
            lista=[]
            lista.append(self.TablaMaterial.item(0, 0).text())
            lista.append(self.TablaMaterial.item(0, 1).text())
            lista.append(float(self.TablaMaterial.item(0, 2).text())/1000.)
            lista.append(self.TablaDiametro.item(0, 0).text())
            for i in range(1, 7):
                lista.append(float(self.TablaDiametro.item(0, i).text())/1000.)
            if self.botonEdit.isChecked():
                lista.append(self.ComboMaterial.currentIndex())
                lista.append(self.TablaDiametro2.currentRow())
            else:
                lista.append(-1)
                lista.append(-1)
                
        except: 
            return False
        else: 
            return lista
        
    def emitirSignal(self):
        material=self.getMaterial()
        if material:
            self.valueChanged.emit(material)

class Catalogo_Materiales_Dialog(QtGui.QDialog, Catalogo_Materiales):
    valueChanged = QtCore.pyqtSignal(list)
    def __init__(self, parent=None):
        super(Catalogo_Materiales_Dialog, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Select Pipe from Database"))
        buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.layout().addWidget(buttonBox,10,0,1,3)


class Catalogo_Accesorios(QtGui.QWidget):
    """Diálogo del catalogo de accesorios de tuberías"""
    valueChanged = QtCore.pyqtSignal(list)
    def __init__(self, parent=None):
        super(Catalogo_Accesorios, self).__init__(parent)
        self.semaforo=QtCore.QSemaphore(1)
        self.accesorios=[]
        self.K=0
        gridLayout = QtGui.QGridLayout(self)
        self.Accesorios = QtGui.QTableWidget()
        self.Accesorios.verticalHeader().hide()
        self.Accesorios.setRowCount(0)
        self.Accesorios.setColumnCount(6)
        self.Accesorios.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Type")))
        self.Accesorios.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem("D,mm"))
        self.Accesorios.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem("D,in"))
        self.Accesorios.setHorizontalHeaderItem(3, QtGui.QTableWidgetItem("K"))
        self.Accesorios.setHorizontalHeaderItem(4, QtGui.QTableWidgetItem(u"Nº"))
        self.Accesorios.setItemDelegateForColumn(4, SpinEditor(self))
        self.Accesorios.setHorizontalHeaderItem(5, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Description")))
        self.Accesorios.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.Accesorios.setSortingEnabled(True)
        self.Accesorios.horizontalHeader().setStretchLastSection(True)
        self.Accesorios.setColumnWidth(0, 80)
        self.Accesorios.setColumnWidth(1, 60)
        self.Accesorios.setColumnWidth(2, 60)
        self.Accesorios.setColumnWidth(3, 50)
        self.Accesorios.setColumnWidth(4, 50)
        self.Accesorios.cellChanged.connect(self.CalcularK)
        self.Accesorios.currentCellChanged.connect(self.selectedChangedAccesorios)
        gridLayout.addWidget(self.Accesorios, 0, 0, 3, 1)
        
        self.k_total=QtGui.QLabel()
        self.k_total.setText("K = 0")
        self.k_total.setAlignment(QtCore.Qt.AlignCenter|QtCore.Qt.AlignTop)
        self.k_total.setFixedSize(100, 25)
        gridLayout.addWidget(self.k_total, 0, 1, 1, 1)
        self.botonBorrar=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editDelete.png")), QtGui.QApplication.translate("pychemqt", "Delete"))
        self.botonBorrar.clicked.connect(self.borrar)
        gridLayout.addWidget(self.botonBorrar, 1, 1, 1, 1)
        self.botonClear=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/clear.png")), QtGui.QApplication.translate("pychemqt", "Clear"))
        self.botonClear.clicked.connect(self.clear)
        gridLayout.addWidget(self.botonClear, 2, 1, 1, 1)
        
        self.line = QtGui.QFrame()
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        gridLayout.addWidget(self.line, 3, 0, 1, 2)
        
        self.TablaAccesorios=QtGui.QTableWidget()
        self.TablaAccesorios.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.TablaAccesorios.horizontalHeader().setStretchLastSection(True)
        self.TablaAccesorios.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.TablaAccesorios.currentCellChanged.connect(self.selectedChanged)
        self.TablaAccesorios.cellDoubleClicked.connect(self.add)
        gridLayout.addWidget(self.TablaAccesorios, 4, 0, 7, 1)
        
        self.botonAdd=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/add.png")), QtGui.QApplication.translate("pychemqt", "Add"))
        self.botonAdd.setDisabled(True)
        self.botonAdd.clicked.connect(self.add)
        gridLayout.addWidget(self.botonAdd, 4, 1, 1, 1)
        
        self.botonEntrada=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/ER.png")), QtGui.QApplication.translate("pychemqt", "Entrance"))
        self.botonEntrada.clicked.connect(self.EntradaRedondeada)
        gridLayout.addWidget(self.botonEntrada, 6, 1, 1, 1)
        self.botonExpansion=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/GE.png")), QtGui.QApplication.translate("pychemqt", "Enlargement"))
        self.botonExpansion.clicked.connect(self.expansion)
        gridLayout.addWidget(self.botonExpansion, 7, 1, 1, 1)
        self.botonContraccion=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/GC.png")), QtGui.QApplication.translate("pychemqt", "Contraction"))
        self.botonContraccion.clicked.connect(self.contraccion)
        gridLayout.addWidget(self.botonContraccion, 8, 1, 1, 1)
        self.botonCodoLargo=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/LB.png")), QtGui.QApplication.translate("pychemqt", "Long bend"))
        self.botonCodoLargo.clicked.connect(self.CodoLargo)
        gridLayout.addWidget(self.botonCodoLargo, 9, 1, 1, 1)
        self.botonCodoSegmentado=QtGui.QPushButton(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/MB45.png")), QtGui.QApplication.translate("pychemqt", "Mitre bend"))
        self.botonCodoSegmentado.clicked.connect(self.CodoSegmentado)
        gridLayout.addWidget(self.botonCodoSegmentado, 10, 1, 1, 1)
        
        self.rellenarTabla()
        
    def rellenarTabla(self):
        self.TablaAccesorios.verticalHeader().hide()
        self.TablaAccesorios.setRowCount(0)
        self.TablaAccesorios.setColumnCount(5)
        self.TablaAccesorios.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Type")))
        self.TablaAccesorios.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem("D,mm"))
        self.TablaAccesorios.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem("D,in"))
        self.TablaAccesorios.setHorizontalHeaderItem(3, QtGui.QTableWidgetItem("K"))
        self.TablaAccesorios.setHorizontalHeaderItem(4, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Description")))
        self.ListaAccesorios=[]
        self.diametros=[]
        self.pulgadas=[]
        self.k=[]
        self.comboxPulgadas=[]
        self.comboxmm=[]
        for i, fitting in enumerate(Accesorios_Tuberia):
            if fitting[2] in self.ListaAccesorios:
                indice=self.ListaAccesorios.index(fitting[2])
                self.diametros[indice].append(fitting[0])
                self.pulgadas[indice].append(fitting[3])
                self.k[indice].append(fitting[5])
                self.comboxPulgadas[indice].addItem(fitting[3])
                self.comboxmm[indice].addItem(str(int(fitting[0])))
            else:
                indice=self.TablaAccesorios.rowCount()
                self.ListaAccesorios.append(fitting[2])
                self.diametros.append([fitting[0]])
                self.pulgadas.append([fitting[3]])
                self.k.append([fitting[5]])
                self.TablaAccesorios.setRowCount(indice+1)
                self.TablaAccesorios.setItem(indice, 0, QtGui.QTableWidgetItem(QtGui.QIcon(QtGui.QPixmap(Fitting_icon[int(fitting[1])])), fitting[2]))
                self.comboxmm.append(QtGui.QComboBox())
                self.TablaAccesorios.setCellWidget(indice, 1, self.comboxmm[-1])
                self.comboxmm[-1].addItem(str(int(fitting[0])))
                self.comboxmm[-1].currentIndexChanged.connect(self.combommChanged)
                self.comboxPulgadas.append(QtGui.QComboBox())
                self.TablaAccesorios.setCellWidget(indice, 2, self.comboxPulgadas[-1])
                self.comboxPulgadas[-1].addItem(fitting[3])
                self.comboxPulgadas[-1].currentIndexChanged.connect(self.comboinchChanged)
                self.TablaAccesorios.setItem(indice, 3, QtGui.QTableWidgetItem("%0.3f" %fitting[5]))
                self.TablaAccesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                self.TablaAccesorios.setItem(indice, 4, QtGui.QTableWidgetItem(fitting[4]))
                self.TablaAccesorios.setRowHeight(self.TablaAccesorios.rowCount()-1, 20)
                
        self.TablaAccesorios.resizeColumnToContents(0)
        self.TablaAccesorios.setColumnWidth(1, 60)
        self.TablaAccesorios.setColumnWidth(2, 60)
        self.TablaAccesorios.setColumnWidth(3, 50)

        
    def expansion(self):
        dialog=Dialog(0, QtGui.QApplication.translate("pychemqt", "Expansion"), os.environ["pychemqt"]+"/images/equip/GE.png", [QtGui.QApplication.translate("pychemqt", "Input diameter"), QtGui.QApplication.translate("pychemqt", "Output diameter")])
        if dialog.exec_():
            indice=self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            if dialog.checkBox.isChecked():
                icon=QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/GE.png"))
                type="GE"
                txt=QtGui.QApplication.translate("pychemqt", "Gradual enlargement")
                dia=u"%0.3f->%0.3f θ=%iº" %(dialog.D1.value, dialog.D2.value, dialog.angulo.value())
            else:
                icon=QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/SE.png"))
                type="SE"
                txt=QtGui.QApplication.translate("pychemqt", "Sudden enlargement")
                dia= "%0.3f->%0.3f" %(dialog.D1.value, dialog.D2.value)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(icon, type))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem(dia))
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(txt))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()
    
    def contraccion(self):
        dialog=Dialog(0, QtGui.QApplication.translate("pychemqt", "Contraction"), os.environ["pychemqt"]+"/images/equip/GC.png", [QtGui.QApplication.translate("pychemqt", "Input diameter"), QtGui.QApplication.translate("pychemqt", "Output diameter")])
        if dialog.exec_():
            indice=self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            if dialog.checkBox.isChecked():
                icon=QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/GC.png"))
                type="GC"
                txt=QtGui.QApplication.translate("pychemqt", "Gradual contraction")
                dia=u"%0.3f->%0.3f θ=%iº" %(dialog.D1.value, dialog.D2.value, dialog.angulo.value())
            else:
                icon=QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/SC.png"))
                type="SC"
                txt=QtGui.QApplication.translate("pychemqt", "Sudden contraction")
                dia= "%0.3f->%0.3f" %(dialog.D1.value, dialog.D2.value)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(icon, type))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem(dia))
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(txt))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()
            
    def EntradaRedondeada(self):
        dialog=Dialog(3, QtGui.QApplication.translate("pychemqt", "Rounded entrance"), os.environ["pychemqt"]+"/images/equip/ER.png", [QtGui.QApplication.translate("pychemqt", "Exit radio"), QtGui.QApplication.translate("pychemqt", "Pipe Diameter")])
        if dialog.exec_():
            indice=self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/ER.png")), "ER"))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem("r=%0.3f, D=%0.3f" %(dialog.D1.value.mm, dialog.D2.value.mm)))
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Pipe entrance rounded intake")))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()
            
    def CodoLargo(self):
        dialog=Dialog(2, QtGui.QApplication.translate("pychemqt", "Long Pipe Bend"), os.environ["pychemqt"]+"/images/equip/LB.png", [QtGui.QApplication.translate("pychemqt", "Bend radio"), QtGui.QApplication.translate("pychemqt", "Pipe diameter")])
        if dialog.exec_():
            indice=self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/LB.png")), "LB"))
            self.Accesorios.setSpan(indice, 1, 1, 2)
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem("r=%0.3f, D=%0.3f" %(dialog.D1.value.mm, dialog.D2.value.mm)))
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Long Pipe Bend")))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()
            
    def CodoSegmentado(self):
        dialog=Dialog(1, QtGui.QApplication.translate("pychemqt", "Mitre bend with custom angle"), os.environ["pychemqt"]+"/images/equip/MB45.png", [QtGui.QApplication.translate("pychemqt", "Pipe diameter"), ""])
        if dialog.exec_():
            indice=self.Accesorios.rowCount()
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/MB45.png")), "MBx"))
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem("%0.3f" %dialog.D1.value.mm))
            self.Accesorios.setItem(indice, 2, QtGui.QTableWidgetItem(u"θ=%iº" %dialog.angulo.value()))
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(representacion(dialog.K, 3)))
            self.Accesorios.item(indice, 3).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem(str(1)))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Mitre Bend")))
            self.Accesorios.setRowHeight(indice, 20)
            self.CalcularK()

    def combommChanged(self, i):
        if self.semaforo.available()>0:
            self.semaforo.acquire(1)
            indice=self.comboxmm.index(self.focusWidget())
            self.comboxPulgadas[indice].setCurrentIndex(i)
            self.TablaAccesorios.item(indice, 3).setText("%0.3f" %self.k[indice][i])
            self.TablaAccesorios.selectRow(indice)
            self.semaforo.release(1)

    def comboinchChanged(self, i):
        if self.semaforo.available()>0:
            self.semaforo.acquire(1)
            indice=self.comboxPulgadas.index(self.focusWidget())
            self.comboxmm[indice].setCurrentIndex(i)
            self.TablaAccesorios.item(indice, 3).setText("%0.3f" %self.k[indice][i])
            self.TablaAccesorios.selectRow(indice)
            self.semaforo.release(1)

    def borrar(self):
        del self.accesorios[self.Accesorios.currentRow()]
        self.Accesorios.removeRow(self.Accesorios.currentRow())
        if self.accesorios:
            self.CalcularK()
        else:
            self.clear()
    
    def clear(self):
        self.Accesorios.clearContents()
        self.Accesorios.setRowCount(0)
        self.k_total.setText("K = 0")
        self.K=0
        self.accesorios=[]
    
    def CalcularK(self):
        k=0
        for i in range(self.Accesorios.rowCount()):
            k+=float(self.Accesorios.item(i, 3).text())*float(self.Accesorios.item(i, 4).text())
        self.k_total.setText("K = %0.3f" %k)
        self.K=k
    
    def selectedChanged(self, i):
        self.botonAdd.setEnabled(i>-1)
        
    def selectedChangedAccesorios(self, newi, newj, oldi, oldj):
        self.Accesorios.openPersistentEditor(self.Accesorios.item(newi, 4))
        self.Accesorios.closePersistentEditor(self.Accesorios.item(oldi, 4))

        
    def add(self, indicetabla=None, indicecombo=None, emit=True):
#        self.Accesorios.blockSignals(True)
        if indicetabla==None:
            indicetabla=self.TablaAccesorios.currentRow()
            indicecombo=self.comboxmm[indicetabla].currentIndex()
        print indicetabla, indicecombo, self.accesorios
        if [indicetabla, indicecombo] not in self.accesorios:
            indice=self.Accesorios.rowCount()
            self.accesorios.append([indicetabla, indicecombo])
            self.Accesorios.setRowCount(indice+1)
            self.Accesorios.setItem(indice, 0, QtGui.QTableWidgetItem(self.TablaAccesorios.item(indicetabla, 0)))
            self.Accesorios.item(indice, 0).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            self.Accesorios.setItem(indice, 1, QtGui.QTableWidgetItem(self.comboxmm[indicetabla].currentText()))
            self.Accesorios.item(indice, 1).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.item(indice, 1).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            self.Accesorios.setItem(indice, 2, QtGui.QTableWidgetItem(self.comboxPulgadas[indicetabla].currentText()))
            self.Accesorios.item(indice, 2).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.item(indice, 2).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            self.Accesorios.setItem(indice, 3, QtGui.QTableWidgetItem(self.TablaAccesorios.item(indicetabla, 3)))
            self.Accesorios.item(indice, 3).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            self.Accesorios.setItem(indice, 4, QtGui.QTableWidgetItem("1"))
            self.Accesorios.item(indice, 4).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            self.Accesorios.setItem(indice, 5, QtGui.QTableWidgetItem(self.TablaAccesorios.item(indicetabla, 4)))
            self.Accesorios.item(indice, 5).setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable)
            self.Accesorios.setRowHeight(indice, 20)
        else:
            indice=self.accesorios.index([indicetabla, indicecombo])
            self.Accesorios.item(indice, 4).setText(str(int(self.Accesorios.item(indice, 4).text())+1))
        self.CalcularK()
#        self.Accesorios.blockSignals(False)
            
        if emit:
            array=[]
            for i, [indicetabla, indicecombo] in enumerate(self.accesorios):
                lista=[indicetabla, indicecombo]
                lista.append(float(self.Accesorios.item(i, 3).text()))
                lista.append(int(self.Accesorios.item(i, 4).text()))
                lista.append(str(self.Accesorios.item(i, 0).text().toUtf8()))
                lista.append(str(self.Accesorios.item(i, 1).text().toUtf8()))
                lista.append(str(self.Accesorios.item(i, 2).text().toUtf8()))
                lista.append(str(self.Accesorios.item(i, 5).text().toUtf8()))
                array.append(lista)
            self.valueChanged.emit(array)

        
class UI_equipment(UI_equip):
    """Diálogo de definición de tuberías"""
    Equipment=Pipe()
    def __init__(self, equipment=None, parent=None):
        """
        equipment: instancia de equipo inicial
        """
        super(UI_equipment, self).__init__(Pipe, entrada=False, salida=False, parent=parent)
        
        #Pestaña Catálogo
        self.tabCatalogo = Catalogo_Materiales()
        self.tabCatalogo.valueChanged.connect(partial(self.changeParams, "material"))
        self.tabWidget.insertTab(1, self.tabCatalogo,QtGui.QApplication.translate("pychemqt", "Catalog"))
        
        #Pestaña Accesorios
        self.tabAccesorios = Catalogo_Accesorios()
        self.tabAccesorios.valueChanged.connect(partial(self.changeParams, "accesorios"))
        self.tabWidget.insertTab(2, self.tabAccesorios,QtGui.QApplication.translate("pychemqt", "Fittings"))
        
        #Pestaña calculo
        gridLayout_Calculo = QtGui.QGridLayout(self.tabCalculo)
#        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Modo:")),0,0)
#        self.Modo=QtGui.QComboBox()
#        gridLayout_Calculo.addWidget(self.Modo,0,1,1,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method")),1,0)
        self.metodo=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_METODO:
            self.metodo.addItem(txt)
        self.metodo.currentIndexChanged.connect(self.metodoCambiado)
        gridLayout_Calculo.addWidget(self.metodo,1,1,1,5)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,0,1,6)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Length")),2,0)
        self.l=Entrada_con_unidades(Length, resaltado=True)
        self.l.valueChanged.connect(partial(self.changeParams, "l"))
        gridLayout_Calculo.addWidget(self.l,2,1)
        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Elevation")),3,0)
        self.h=Entrada_con_unidades(Length, min=float("-inf"))
        self.h.valueChanged.connect(partial(self.changeParams, "h"))
        gridLayout_Calculo.addWidget(self.h,3,1)
        self.labelC=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "C Factor"))
        gridLayout_Calculo.addWidget(self.labelC,2,3)
        self.C=Entrada_con_unidades(float)
        self.C.valueChanged.connect(partial(self.changeParams, "C"))
        gridLayout_Calculo.addWidget(self.C,2,4)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,0,1,6)

        gridLayout_Calculo.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Thermal Mode")),11,0)
        self.thermal=QtGui.QComboBox()
        for txt in self.Equipment.TEXT_THERMAL:
            self.thermal.addItem(txt)
        self.thermal.currentIndexChanged.connect(self.Heat_visible)
        gridLayout_Calculo.addWidget(self.thermal,11,1,1,5)
        
        self.groupBox_Heat = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Heat transfer to surroundings"))
        gridLayout_Calculo.addWidget(self.groupBox_Heat,12,0,3,2)        
        layout = QtGui.QGridLayout(self.groupBox_Heat)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "External temperature")),1,1)
        self.T_ext=Entrada_con_unidades(Temperature)
        self.T_ext.valueChanged.connect(partial(self.changeParams, "T_ext"))
        layout.addWidget(self.T_ext,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heat transfer coefficient")),2,1)
        self.U=Entrada_con_unidades(HeatTransfCoef)
        self.U.valueChanged.connect(partial(self.changeParams, "U"))
        layout.addWidget(self.U,2,2)
 
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),12,0,1,6) 
        self.labelQ=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Heat Flux"))
        gridLayout_Calculo.addWidget(self.labelQ,13,0)
        self.Calor=Entrada_con_unidades(Power)
        self.Calor.valueChanged.connect(partial(self.changeParams, "Q"))
        gridLayout_Calculo.addWidget(self.Calor,13,1)
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),14,0,1,6) 
        
        gridLayout_Calculo.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),15,0,1,6)

        self.groupBox_Resultados = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Results"))
        gridLayout_Calculo.addWidget(self.groupBox_Resultados,16,0,1,6)
        gridLayout_3 = QtGui.QGridLayout(self.groupBox_Resultados)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP total", None, QtGui.QApplication.UnicodeUTF8)),0,0)
        self.DeltaP=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP,0,1)        
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP friction", None, QtGui.QApplication.UnicodeUTF8)),1,0)
        self.DeltaP_f=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP_f,1,1)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP fittings", None, QtGui.QApplication.UnicodeUTF8)),2,0)
        self.DeltaP_ac=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP_ac,2,1)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP elevation", None, QtGui.QApplication.UnicodeUTF8)),3,0)
        self.DeltaP_h=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP_h,3,1)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "ΔP acceleration", None, QtGui.QApplication.UnicodeUTF8)),4,0)
        self.DeltaP_v=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP_v,4,1)

        gridLayout_3.addWidget(QtGui.QLabel(u"ΔP/100ft"),0,4)
        self.DeltaP_100ft=Entrada_con_unidades(Pressure, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.DeltaP_100ft,0,5)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Speed")),1,4)
        self.V=Entrada_con_unidades(Speed, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.V,1,5)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Friction factor")),2,4)
        self.f=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.f,2,5)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reynolds")),3,4)
        self.Re=Entrada_con_unidades(float, tolerancia=5, decimales=1, readOnly=True)
        gridLayout_3.addWidget(self.Re,3,5)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Output Temperature")),4,4)
        self.Tout=Entrada_con_unidades(Temperature, decimales=2, retornar=False, readOnly=True)
        gridLayout_3.addWidget(self.Tout,4,5)

        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Flow region")),0,7)
        self.regimen=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.regimen,0,8)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reynolds liquid")),1,7)
        self.Re_liq=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.Re_liq,1,8)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Reynolds gas")),2,7)
        self.Re_gas=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.Re_gas,2,8)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Friction liquid")),3,7)
        self.f_liq=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.f_liq,3,8)
        gridLayout_3.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Friction gas")),4,7)
        self.f_gas=Entrada_con_unidades(float, readOnly=True)
        gridLayout_3.addWidget(self.f_gas,4,8)

        gridLayout_3.addItem(QtGui.QSpacerItem(40,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Minimum),0,3)
        gridLayout_3.addItem(QtGui.QSpacerItem(40,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Minimum),0,6)
#        gridLayout_Calculo.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),17,1,1,6)

        
        #Pestaña costos
        layout_costos = QtGui.QGridLayout(self.tabCostos)
        self.labelAvisoCostos=QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Costs only available for steel pipes"))
        self.labelAvisoCostos.setVisible(False)
        self.labelAvisoCostos.setStyleSheet("color: #FF0000;")
        layout_costos.addWidget(self.labelAvisoCostos,0,0,1,2)

        self.Costos=CostData(self.Equipment)
        self.Costos.valueChanged.connect(self.changeParamsCoste)
        layout_costos.addWidget(self.Costos,1,0,1,2)
        
        layout_costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),2,0,1,2)
        groupBox_Costos = QtGui.QGroupBox(QtGui.QApplication.translate("pychemqt", "Stimated Cost"))
        layout_costos.addWidget(groupBox_Costos,3,0,1,2)
        gridLayout_5 = QtGui.QGridLayout(groupBox_Costos)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Purchase Costs")),0,1)
        self.C_adq=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_adq,0,2)
        gridLayout_5.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Installed Costs")),1,1)
        self.C_inst=Entrada_con_unidades(Currency, retornar=False, readOnly=True)
        gridLayout_5.addWidget(self.C_inst,1,2)
        layout_costos.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),4,0,1,2)

        self.Heat_visible(0)
        self.metodoCambiado(0)
        if equipment:
            self.setEquipment(equipment)


    def Heat_visible(self, indice):
        if indice==1:
            self.groupBox_Heat.setVisible(False)
            self.labelQ.setVisible(True)
            self.Calor.setVisible(True)
        elif indice==2:
            self.groupBox_Heat.setVisible(True)
            self.labelQ.setVisible(False)
            self.Calor.setVisible(False)
        else:
            self.groupBox_Heat.setVisible(False)
            self.labelQ.setVisible(False)
            self.Calor.setVisible(False)
        self.changeParams("thermal", indice)
        
    def metodoCambiado(self, indice):
        self.C.setEnabled(indice==1)
        self.labelC.setEnabled(indice==1)
        self.changeParams("metodo", indice)
        
        
    def rellenar(self):
        UI_equip.rellenar(self)
        if self.Equipment.status==1:
            if self.Equipment.statusCoste:
                self.labelAvisoCostos.setVisible(False)
            else:
                self.labelAvisoCostos.setVisible(True)
        
    def rellenarInput(self):
        self.blockSignals(True)
        UI_equip.rellenarInput(self)
        if self.Equipment.kwargs["material"]:
            if self.Equipment.kwargs["material"][10]>=0:
                self.tabCatalogo.add(self.Equipment.kwargs["material"][10], self.Equipment.kwargs["material"][11], False)
            else:
                self.tabCatalogo.addCustom(self.Equipment.kwargs["material"]) 
                
        if self.Equipment.kwargs["accesorios"]:
            self.tabAccesorios.clear()
            for accesorio in self.Equipment.kwargs["accesorios"]:
                for numero in range(accesorio[3]):
                    self.tabAccesorios.add(accesorio[0], accesorio[1], False)
        self.blockSignals(False)


if __name__ == "__main__":
    import sys        
    from lib.corriente import Corriente
    app = QtGui.QApplication(sys.argv)
    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Steel (ANSI)", "Sch. 40", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, -1, 2], 
                 notas="Tuberia")
    dialogo = UI_equipment(tuberia)
    dialogo.show()
    sys.exit(app.exec_())
