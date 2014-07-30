#!/usr/bin/python
# -*- coding: utf-8 -*-

import cPickle
from functools import partial

from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
#import matplotlib.gridspec as gridspec #necesita matplotlib >=1.0
from pylab import Figure
from numpy import transpose

from lib.unidades import Length, VolFlow, Power
from lib import config
from widgets import Entrada_con_unidades, Tabla


class Plot(FigureCanvasQTAgg):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

        
class Ui_bombaCurva(QtGui.QDialog):
    def __init__(self, curva=[], parent=None):
        """curva: Parametro opcional indicando la curva de la bomba"""
        super(Ui_bombaCurva, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Pump curves dialog"))
        self.showMaximized()

        self.gridLayout = QtGui.QGridLayout(self)
        self.botones = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Apply|QtGui.QDialogButtonBox.Open|QtGui.QDialogButtonBox.Save|QtGui.QDialogButtonBox.Reset, QtCore.Qt.Vertical)
        self.botones.clicked.connect(self.botones_clicked)
        self.gridLayout.addWidget(self.botones,1,1,3,1)
        self.gridLayout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,2,3,1)
        self.gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Curves")),1,3)
        self.lista=QtGui.QComboBox()
        self.lista.currentIndexChanged.connect(self.cambiarCurvaVista)
        self.gridLayout.addWidget(self.lista,1,4)
        self.gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Diameter")),2,3)
        self.diametro=Entrada_con_unidades(int, width=60, textounidad='"')
        self.gridLayout.addWidget(self.diametro,2,4)
        self.gridLayout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "RPM")),3,3)
        self.rpm=Entrada_con_unidades(int, width=60, textounidad="rpm")
        self.gridLayout.addWidget(self.rpm,3,4)
        self.gridLayout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),4,1,1,4)

        self.Tabla=Tabla(4, horizontalHeader=[QtGui.QApplication.translate("pychemqt", "Flowrate"), QtGui.QApplication.translate("pychemqt", "Head"), QtGui.QApplication.translate("pychemqt", "Power"), QtGui.QApplication.translate("pychemqt", "NPSH")], verticalOffset=1, filas=1, stretch=False)
        self.Tabla.setColumnWidth(0, 100)
        self.unidadesCaudal = QtGui.QComboBox()
        self.Tabla.setCellWidget(0, 0, self.unidadesCaudal)
        self.Tabla.setColumnWidth(1, 80)
        self.unidadesCarga = QtGui.QComboBox()
        self.Tabla.setCellWidget(0, 1, self.unidadesCarga)
        self.Tabla.setColumnWidth(2, 80)
        self.unidadesPotencia = QtGui.QComboBox()
        self.Tabla.setCellWidget(0, 2, self.unidadesPotencia)
        self.Tabla.setColumnWidth(3, 80)
        self.unidadesNPSH = QtGui.QComboBox()
        self.Tabla.setCellWidget(0, 3, self.unidadesNPSH)
        self.Tabla.setRowHeight(0, 24)
        self.Tabla.setFixedWidth(360)
        self.Tabla.setConnected()
        self.Tabla.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.gridLayout.addWidget(self.Tabla,5,1,1,4)
        self.gridLayout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),1,5,5,1)

        self.checkCarga = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Heat"))
        self.gridLayout.addWidget(self.checkCarga,1,6)
        self.checkPotencia = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Power"))
        self.gridLayout.addWidget(self.checkPotencia,1,7)
        self.checkNPSH = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "NPSH"))
        self.gridLayout.addWidget(self.checkNPSH,1,8)
        self.rejilla = QtGui.QCheckBox(QtGui.QApplication.translate("pychemqt", "Grid"))
        self.rejilla.toggled.connect(self.rejilla_toggled)
        self.gridLayout.addWidget(self.rejilla,1,9)
        self.gridLayout.addItem(QtGui.QSpacerItem(1000,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1,10)
        self.Plot = Plot(self, width=5, height=1, dpi=100)
        self.gridLayout.addWidget(self.Plot,2,6,4,5)
        

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.gridLayout.addWidget(self.buttonBox,6,9,1,2)
        
        if curva:
            self.curvas=curva
            for i in curva:
                self.lista.addItem(str(i[0])+'", '+str(i[1])+" rpm")
            self.lista.setCurrentIndex(self.lista.count()-1)
            self.cambiarCurvaVista(self.lista.count()-1)

#            self.curva=curva[-1]
#            self.diametro.setValue(self.curva[0])
#            self.rpm.setValue(self.curva[1])
#            self.Tabla.setMatrix(self.curva[2:])
#            self.Tabla.addRow()
            self.actualizarPlot()
        else:
            self.curva=[]
            
        for i in Length.__text__:
            self.unidadesCarga.addItem(i)
            self.unidadesNPSH.addItem(i)
        for i in VolFlow.__text__:
            self.unidadesCaudal.addItem(i)
        for i in Power.__text__:
            self.unidadesPotencia.addItem(i)
            
        self.oldIndices=[0, 0, 0, 0]
        self.unidadesCaudal.currentIndexChanged.connect(partial(self.cambiar_unidades, 0, VolFlow, "VolFlow"))
        self.unidadesCarga.currentIndexChanged.connect(partial(self.cambiar_unidades, 1, Length,"Head", ))
        self.unidadesPotencia.currentIndexChanged.connect(partial(self.cambiar_unidades, 2, Power, "Power"))
        self.unidadesNPSH.currentIndexChanged.connect(partial(self.cambiar_unidades, 3, Length,"Head", ))
        
        Config=config.getMainWindowConfig()
        self.unidadesCaudal.setCurrentIndex(Config.getint("Units","QLiq"))
        self.unidadesCarga.setCurrentIndex(Config.getint("Units","Head"))
        self.unidadesPotencia.setCurrentIndex(Config.getint("Units","Power"))
        self.unidadesNPSH.setCurrentIndex(Config.getint("Units","Head"))
        
        self.checkCarga.toggled.connect(self.actualizarPlot)
        self.checkPotencia.toggled.connect(self.actualizarPlot)
        self.checkNPSH.toggled.connect(self.actualizarPlot)


    def botones_clicked(self, boton):
        if boton == self.botones.button(QtGui.QDialogButtonBox.Reset):
            self.diametro.clear()
            self.rpm.clear()
            self.Tabla.clear()
            self.checkCarga.setChecked(False)
            self.checkPotencia.setChecked(False)
            self.checkNPSH.setChecked(False)
            self.rejilla.setChecked(False)
            
        elif boton == self.botones.button(QtGui.QDialogButtonBox.Open):
            fname = unicode(QtGui.QFileDialog.getOpenFileName(self, QtGui.QApplication.translate("pychemqt", "Open curve file"), "./", "cpickle file (*.pkl);;All files (*.*)"))
            if fname:
                with open(fname, "r") as archivo:
                    curvas=cPickle.load(archivo)
                self.curvas=curvas
                self.curva=curvas[-1]
                self.lista.clear()
                self.lista.blockSignals(True)
                for i in curvas:
                    self.lista.addItem(str(i[0])+'", '+str(i[1])+" rpm")
                self.lista.blockSignals(False)
                self.lista.setCurrentIndex(self.lista.count()-1)
                self.cambiarCurvaVista(self.lista.count()-1)
                self.actualizarPlot()
                
        elif boton == self.botones.button(QtGui.QDialogButtonBox.Apply):
            txt=str(self.diametro.value)+'", '+str(self.rpm.value)+" rpm"
            indice=self.lista.findText(txt)
            if indice ==-1:
                self.lista.addItem(txt)
                self.curva=self.CalcularCurva()
                self.curvas.append(self.curva)
                self.lista.setCurrentIndex(self.lista.count()-1)
                self.actualizarPlot()
            else:
                self.curva=self.CalcularCurva()
                self.curvas[indice]=self.curva
                self.lista.setCurrentIndex(indice)
                self.actualizarPlot()
                
        else:
            fname = unicode(QtGui.QFileDialog.getSaveFileName(self, QtGui.QApplication.translate("pychemqt", "Save curve to file"), "./", "cpickle file (*.pkl)"))
            if fname:
                if fname.split(".")[-1]!="pkl":
                    fname+=".pkl"
                cPickle.dump(self.curvas, open(fname, "w"))

    def CalcularCurva(self):
        caudal=[]
        carga=[]
        potencia=[]
        NPSH=[]
        for i in range(1, self.Tabla.rowCount()-1):
            q=VolFlow(float(self.Tabla.item(i, 0).text()), VolFlow.__units__[self.unidadesCaudal.currentIndex()])
            caudal.append(q)
            h=Length(float(self.Tabla.item(i, 1).text()), Length.__units__[self.unidadesCarga.currentIndex()])
            carga.append(h)
            P=Power(float(self.Tabla.item(i, 2).text()),Power.__units__[self.unidadesPotencia.currentIndex()])
            potencia.append(P)
            n=Length(float(self.Tabla.item(i, 3).text()), Length.__units__[self.unidadesNPSH.currentIndex()])
            NPSH.append(n)
        return [self.diametro.value, self.rpm.value, caudal, carga, potencia, NPSH]

    def cambiar_unidades(self, col, unit, magn, i):
        text=self.Tabla.getColumn(col)
        old_i=self.oldIndices[col]
        self.Tabla.blockSignals(True)
        for fila in range(0, len(text)):
            a=unit(text[fila], units[magn][old_i]).__getattribute__(units[magn][i])
            self.Tabla.item(fila+1, col).setText(config.representacion(a))
        self.Tabla.blockSignals(False)
        self.oldIndices[col]=i
        
    def cambiarCurvaVista(self, ind):
        self.curva=self.curvas[ind]
        self.diametro.setValue(self.curva[0])
        self.rpm.setValue(self.curva[1])
        q=[VolFlow(i).__getattribute__(VolFlow.__units__[self.unidadesCaudal.currentIndex()]) for i in self.curva[2]]
        h=[Length(i).__getattribute__(Length.__units__[self.unidadesCarga.currentIndex()]) for i in self.curva[3]]
        P=[Power(i).__getattribute__(Power.__units__[self.unidadesPotencia.currentIndex()]) for i in self.curva[4]]
        n=[Length(i).__getattribute__(Length.__units__[self.unidadesNPSH.currentIndex()]) for i in self.curva[5]]
        self.Tabla.setMatrix(transpose([q, h, P, n]))
        self.Tabla.addRow()

    def actualizarPlot(self):
        #FIXME: mejorar el manejo de gráficas, básicamente empollarse matplotlib, opcional para profesionalizar el gráfico
        Carga=self.checkCarga.isChecked()
        Potencia=self.checkPotencia.isChecked()
        NPSH=self.checkNPSH.isChecked()
        self.Plot.fig.clear()
        share=None
        self.Plot.carga=[]
        self.Plot.potencia=[]
        self.Plot.npsh=[]
        total=Carga+Potencia+NPSH
        actual=1
        if Carga:
            self.Plot.ax1 = self.Plot.fig.add_subplot(total, 1, actual)
            self.Plot.ax1.grid(self.rejilla.isChecked())
            self.Plot.ax1.set_title("Curva H/Q", size='14')
            self.Plot.ax1.set_ylabel("H, m", size='12')
            share=self.Plot.ax1
            actual+=1
            for curva in self.curvas:
                Q=curva[2]
                h=curva[3]
                self.Plot.carga+=self.Plot.ax1.plot(Q, h, label=str(curva[0])+'", '+str(curva[1])+" rpm")
            self.Plot.ax1.legend()

        if Potencia:
            self.Plot.ax2 = self.Plot.fig.add_subplot(total, 1, actual, sharex=share)
            self.Plot.ax2.grid(self.rejilla.isChecked())
            self.Plot.ax2.set_ylabel("P, kW", size='12')
            share=self.Plot.ax2
            actual+=1
            for curva in self.curvas:
                Q=curva[2]
                Pot=curva[4]
                self.Plot.potencia+=self.Plot.ax2.plot(Q, Pot, label=str(curva[0])+'", '+str(curva[1])+" rpm")
            self.Plot.ax2.legend(loc='upper left')
                
        if NPSH:
            self.Plot.ax3 = self.Plot.fig.add_subplot(total, 1, actual, sharex=share)
            self.Plot.ax3.grid(self.rejilla.isChecked())
            self.Plot.ax3.set_ylabel("NPSH, m", size='12')
            for curva in self.curvas:
                Q=curva[2]
                npsh=curva[5]
                self.Plot.npsh+=self.Plot.ax3.plot(Q, npsh, label=str(curva[0])+'", '+str(curva[1])+" rpm")
            self.Plot.ax3.legend(loc='upper left')

        if share:
            share.set_xlabel("Q, $m^3/s$", horizontalalignment='left', size='12')


        
#        self.Plot.ax2.set_title("Curva P/Q", size='14')
#        self.Plot.ax3.set_title("Curva NPSH/Q", size='14')
#        self.Plot.ax1.set_xlabel("Q, $m^3/s$", horizontalalignment='left', size='12')
#        self.Plot.ax2.set_xlabel("Q, $m^3/s$", horizontalalignment='left', size='12')
        #otras rpm
#        for rpm in [1200., 1500., 1800., 2100.]:
#            Q2=[]
#            h2=[]
#            for i in range(len(Q)):
#                Q2.append(Q[i]*rpm/1400)
#                h2.append(h[i]*rpm**2/1400**2)
#            self.Plot.ax1.plot(Q2, h2)


#        self.Plot.tight_layout()       #Necesita matploblib >=1.1
        self.Plot.draw()

    def rejilla_toggled(self, bool):
        self.Plot.ax1.grid(bool)
        self.Plot.ax2.grid(bool)
        self.Plot.ax3.grid(bool)
        self.Plot.draw()


if __name__ == "__main__":
    import sys
    from numpy import r_, zeros, arange
    app = QtGui.QApplication(sys.argv)
    Q=arange(0, 21, 1.)
    h=[15.5, 15.468115, 15.40372, 15.319705, 15.23896, 15.154375, 15.06884, 14.945245, 14.84648, 14.725435, 14.605, 14.418065, 14.18752, 13.906255, 13.56716, 13.163125, 12.68704, 12.131795, 11.49028, 10.755385, 9.92]
#    h=r_[15.5, 15.49, 15.45, 15.4, 15.3, 15.2, 15.1, 15., 14.8, 14.65, 14.505, 14.35, 14.15, 13.9, 13.55, 13.2, 12.75, 12.3, 11.6, 10.85, 9.7]
    Pot=r_[0.5, 0.51, 0.53, 0.55, 0.575, 0.595, 0.605, 0.625, 0.66, 0.685, 0.705, 0.73, 0.765, 0.795, 0.82, 0.855, 0.895, 0.92, 0.975, 1.005, 1.06]*1000.
    NHPS=zeros(21)
    curva=[[212, 1400, Q/3600., h, Pot, NHPS], ]


    bombaCurva = Ui_bombaCurva(curva)
    bombaCurva.show()
    sys.exit(app.exec_())
