#!/usr/bin/python
# -*- coding: utf-8 -*-

###Módulo que define los dialogos de definición de gráficos

from PyQt4 import QtCore, QtGui
from scipy import arange
from numpy import transpose

from UI.widgets import Tabla
from lib.plot import mpl
from lib import config, unidades
from lib.corriente import Corriente, Mezcla

class Binary_distillation(QtGui.QDialog):
    title=QtGui.QApplication.translate("pychemqt", "x-y Distillation")
    def __init__(self, indices=None, nombres=None, x=None, y=None, parent=None):
        super(Binary_distillation, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Component 1:")),1,1)
        self.Comp1=QtGui.QComboBox()
        layout.addWidget(self.Comp1,1,2)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("equipment", "Component 2:")),1,4)
        self.Comp2=QtGui.QComboBox()
        layout.addWidget(self.Comp2,1,5)
        
        self.indices=indices
        self.nombres=nombres
        for i, nombre in enumerate(nombres):
            self.Comp1.addItem("%i - %s" %(i+1, nombre))
            self.Comp2.addItem("%i - %s" %(i+1, nombre))
        self.Comp2.setCurrentIndex(1)
        tab=QtGui.QTabWidget()
        layout.addWidget(tab,2,1,1,5)
        self.plot=mpl()
        tab.addTab(self.plot, QtGui.QApplication.translate("equipment", "Plot"))
        self.tabla=Tabla(2, horizontalHeader=["x", "y"], stretch=False, readOnly=True)
        tab.addTab(self.tabla, QtGui.QApplication.translate("equipment", "Table"))
        
        self.Comp1.currentIndexChanged.connect(self.calculo)
        self.Comp2.currentIndexChanged.connect(self.calculo)
        
        if x and y:
            self.rellenar(x, y)
        else:
            self.calculo()
        
    def rellenar(self, x, y):
        self.x=x
        self.y=y
        self.plot.axes2D.clear()
        self.plot.data([0, 1], [0, 1], x, y, 'ro')

        self.tabla.setMatrix(transpose([x, y]))

    def calculo(self):
        ind1=self.Comp1.currentIndex()
        ind2=self.Comp2.currentIndex()
        if ind1!=ind2:
            zi=arange(0.025, 1., 0.025)
            id1=self.indices[ind1]
            id2=self.indices[ind2]
            x=[0]
            y=[0]
            for z in zi:
                try:
                    fraccion=[0.]*len(self.indices)
                    fraccion[ind1]=z
                    fraccion[ind2]=1-z
                    mez=Mezcla(tipo=3, fraccionMolar=fraccion, caudalMasico=1.)
                    tb=mez.componente[0].Tb
                    corr=Corriente(T=tb, P=101325., mezcla=mez)
                    T=corr.eos._Dew_T()
                    corr=Corriente(T=T, P=101325., mezcla=mez)
                    while corr.Liquido.fraccion[0]==corr.Gas.fraccion[0] and corr.T<corr.mezcla.componente[1].Tb:
                        corr=Corriente(T=corr.T-0.1, P=101325., mezcla=mez)
                    x.append(corr.Liquido.fraccion[0])
                    y.append(corr.Gas.fraccion[0])
                except:
                    pass
            x.append(1)
            y.append(1)
            self.rellenar(x, y)

    def writeToStream(self, stream):
        stream.writeInt32(self.widget().Comp1.currentIndex())
        stream.writeInt32(self.widget().Comp2.currentIndex())
        stream.writeInt32(len(self.widget().x))
        for i in self.widget().x:
            stream.writeFloat(i)
        for i in self.widget().y:
            stream.writeFloat(i)
    
    @classmethod
    def readToStream(cls, stream):
        id1=stream.readInt32()
        id2=stream.readInt32()
        len=stream.readInt32()
        x=[]
        for i in range(len):
            x.append(stream.readFloat())
        y=[]
        for i in range(len):
            y.append(stream.readFloat())
        self.plot(0, x, y)


class Plot_Distribucion(mpl):
    title=QtGui.QApplication.translate("pychemqt", "Solid Distribution")
    def __init__(self, id, solido, parent=None):
        super(Plot_Distribucion, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Stream")+" "+str(id)+" - "+self.title)
        self.fill(solido)
        self.axes2D.legend(loc=4)
        self.axes2D.set_ylim(0, 1)
        self.axes2D.set_xlabel("Dp, "+unidades.Length(None).text("ParticleDiameter"), horizontalalignment='right', size='12')
        self.axes2D.set_ylabel(QtGui.QApplication.translate("pychemqt", "Accumulated fraction"), horizontalalignment='right', size='12')


    def fill(self, solido):
        self.data(solido.diametros, solido.fracciones_acumuladas, 'ro')
        for i in range(6):
            x, y, leyenda=solido.ajustar_distribucion(i)
            self.data(x, y, label=leyenda)
        
        
        
__all__=Binary_distillation, 

if __name__ == "__main__":
    import sys
    from ConfigParser import ConfigParser
    configuracion=ConfigParser()
    configuracion.read(config.conf_dir+"pychemqtrc")

    app = QtGui.QApplication(sys.argv)
    Dialog = Binary_distillation(configuracion)
    Dialog.show()
    sys.exit(app.exec_())

