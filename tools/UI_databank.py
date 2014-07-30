#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3

from PyQt4 import QtCore, QtGui, QtSql

from lib import sql
from UI import viewComponents

class UI_databank_widget(QtGui.QWidget):
    def __init__(self, parent=None):
        super(UI_databank_widget, self).__init__(parent)
        gridLayout = QtGui.QGridLayout(self)
        
        layoutTitle=QtGui.QHBoxLayout()
        layoutTitle.setSpacing(5)
        self.buttonNew=QtGui.QToolButton(self)
        self.buttonNew.setToolTip(QtGui.QApplication.translate("pychemqt", "Create new element"))
        self.buttonNew.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/fileNew.png")))
        self.buttonNew.clicked.connect(self.newComponent)
        layoutTitle.addWidget(self.buttonNew)
        self.buttonCopy=QtGui.QToolButton(self)
        self.buttonCopy.setToolTip(QtGui.QApplication.translate("pychemqt", "Clone this element"))
        self.buttonCopy.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editCopy.png")))
        self.buttonCopy.clicked.connect(self.copyComponent)
        layoutTitle.addWidget(self.buttonCopy)
        self.buttonDelete=QtGui.QToolButton(self)
        self.buttonDelete.setToolTip(QtGui.QApplication.translate("pychemqt", "Delete element"))
        self.buttonDelete.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/editDelete.png")))
        self.buttonDelete.clicked.connect(self.deleteComponent)
        self.buttonDelete.setEnabled(False)
        layoutTitle.addWidget(self.buttonDelete)
        gridLayout.addItem(layoutTitle,1,1)

        gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Find")+": "),1,2)
        self.Busqueda = QtGui.QLineEdit()
        self.Busqueda.textChanged.connect(self.buscar)
        gridLayout.addWidget(self.Busqueda,1,3)
        self.siguiente = QtGui.QPushButton(QtGui.QApplication.translate("pychemqt", "Next"))
        self.siguiente.clicked.connect(self.Next)
        gridLayout.addWidget(self.siguiente,1,4)

        self.BaseDatos = QtGui.QTableWidget()
        self.BaseDatos.setMinimumWidth(375)
        self.BaseDatos.verticalHeader().hide()
        self.BaseDatos.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.BaseDatos.setShowGrid(False)
        self.BaseDatos.setRowCount(0)
        self.BaseDatos.setColumnCount(3)
        self.BaseDatos.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem("Id"))
        self.BaseDatos.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Name")))
        self.BaseDatos.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem(QtGui.QApplication.translate("pychemqt", "Formula")))
        self.BaseDatos.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.BaseDatos.horizontalHeader().setStretchLastSection(True)
        self.BaseDatos.currentCellChanged.connect(self.checkButton)
        self.BaseDatos.doubleClicked.connect(self.mostrarPropiedades)
        gridLayout.addWidget(self.BaseDatos,2,1,1,4)
        self.correctos=[]
        self.indice=0
        self.rellenar()
        
    def rellenar(self):
        self.BaseDatos.setRowCount(0)
        sql.databank.execute("select * from compuestos")
        for i in sql.databank:
            self.BaseDatos.setRowCount(self.BaseDatos.rowCount()+1)
            self.BaseDatos.setItem(i[0]-1, 0, QtGui.QTableWidgetItem(str(i[0])))
            self.BaseDatos.setItem(i[0]-1, 1, QtGui.QTableWidgetItem(i[2]))
            self.BaseDatos.setItem(i[0]-1, 2, QtGui.QTableWidgetItem(i[1]))
            self.BaseDatos.setRowHeight(self.BaseDatos.rowCount()-1, 20)
            
        sql.databank_Custom.execute("select * from compuestos")
        for i in sql.databank_Custom:
            filas=self.BaseDatos.rowCount()
            self.BaseDatos.setRowCount(filas+1)
            self.BaseDatos.setItem(filas, 0, QtGui.QTableWidgetItem(str(i[0])))
            self.BaseDatos.setItem(filas, 1, QtGui.QTableWidgetItem(i[2]))
            self.BaseDatos.setItem(filas, 2, QtGui.QTableWidgetItem(i[1]))
            self.BaseDatos.setRowHeight(self.BaseDatos.rowCount()-1, 20)
        
        self.BaseDatos.resizeColumnsToContents()
        
    def buscar(self):
        self.indice=0
        texto="%"+self.Busqueda.text()+"%"
        sql.databank.execute("select * from compuestos where nombre LIKE '%s' or formula LIKE '%s'" % (texto, texto))
        self.correctos=[]
        for i in sql.databank:
            self.correctos.append(i[0])
        self.BaseDatos.setCurrentCell(self.correctos[self.indice]-1, 0)
    
    def Next(self):
        if self.indice < len(self.correctos)-1:
            self.indice+=1
        else:
            self.indice=0
        self.BaseDatos.setCurrentCell(self.correctos[self.indice]-1, 0)
            
    def checkButton(self, indice):
        if indice>=sql.N_comp:
            self.buttonDelete.setEnabled(True)
        else:
            self.buttonDelete.setEnabled(False)
        
    def mostrarPropiedades(self):
        indice=self.currentIndice
        if indice >0:
            Dialog = viewComponents.View_Component(indice)
            Dialog.exec_()
            
    @property
    def currentIndice(self):
        return float(self.BaseDatos.item(self.BaseDatos.currentRow(), 0).text())
        
    def currentRow(self):
        return self.BaseDatos.currentRow()
        
    def newComponent(self):
        Dialog = viewComponents.View_Component(0)
        if Dialog.exec_():
            self.rellenar()
        
    def copyComponent(self):
        sql.copyElement(self.currentIndice)
        self.rellenar()
        
    def deleteComponent(self):
        sql.deleteElement(self.currentIndice)
        self.rellenar()


class UI_databank(QtGui.QDialog):
    def __init__(self, parent=None):
        super(UI_databank, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Components database"))
        layout = QtGui.QVBoxLayout(self)
        self.databank=UI_databank_widget()
        layout.addWidget(self.databank)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = UI_databank()
    Dialog.show()
    sys.exit(app.exec_())
