#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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



###Módulo que define los dialogos de definición de gráficos


from numpy import arange, transpose

from tools.qt import QtWidgets, translate

from lib import config, unidades
from lib.corriente import Corriente
from lib.mezcla import Mezcla
from lib.plot import PlotWidget
from UI.widgets import Tabla


class Binary_distillation(QtWidgets.QDialog):
    title = translate("Binary_distillation", "x-y Distillation")
    def __init__(self, indices=None, nombres=None, x=None, y=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.title)
        lyt = QtWidgets.QGridLayout(self)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Component 1:")),1,1)
        self.Comp1 = QtWidgets.QComboBox()
        lyt.addWidget(self.Comp1,1,2)
        lyt.addWidget(QtWidgets.QLabel(self.tr("Component 2:")),1,4)
        self.Comp2 = QtWidgets.QComboBox()
        lyt.addWidget(self.Comp2,1,5)
        lyt.addItem(QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Fixed), 1, 6)

        self.indices = indices
        self.nombres = nombres
        for i, nombre in enumerate(nombres):
            self.Comp1.addItem("%i - %s" %(i+1, nombre))
            self.Comp2.addItem("%i - %s" %(i+1, nombre))
        self.Comp2.setCurrentIndex(1)
        tab = QtWidgets.QTabWidget()
        lyt.addWidget(tab,2,1,1,6)
        self.plot = PlotWidget(parent=self)
        tab.addTab(self.plot, self.tr("Plot"))
        self.tabla = Tabla(2, horizontalHeader=["x", "y"], stretch=False, readOnly=True)
        tab.addTab(self.tabla, self.tr("Table"))

        self.Comp1.currentIndexChanged.connect(self.calculo)
        self.Comp2.currentIndexChanged.connect(self.calculo)

        if x and y:
            self.rellenar(x, y)
        # else:
            # self.calculo()

    def rellenar(self, x, y):
        self.x = x
        self.y = y
        # self.plot.axes2D.clear()
        # self.plot.data([0, 1], [0, 1], x, y, 'ro')

        # self.tabla.setData(transpose([x, y]))

    def calculo(self):
        ind1 = self.Comp1.currentIndex()
        ind2 = self.Comp2.currentIndex()
        print(ind1, ind2)
        if ind1 != ind2:
            zi = arange(0.025, 1., 0.025)
            id1 = self.indices[ind1]
            id2 = self.indices[ind2]
            x = [0]
            y = [0]
            for z in zi:
                # try:
                fraccion = [0.]*len(self.indices)
                fraccion[ind1] = z
                fraccion[ind2] = 1-z
                mez = Mezcla(tipo=3, fraccionMolar=fraccion, caudalMasico=1.)
                tb = mez.componente[0].Tb
                # corr = Corriente(T=tb, P=101325., mezcla=mez)
                corr = Corriente(T=tb, P=101325., ids=self.indices,
                                 fraccionMolar=fraccion, caudalMasico=1)
                print(corr.status, corr.msg, corr.eos)
                T = corr.eos._Dew_T(101325)
                corr = Corriente(T=T, P=101325., mezcla=mez)
                while corr.Liquido.fraccion[0] == corr.Gas.fraccion[0] and \
                        corr.T < corr.mezcla.componente[1].Tb:
                    corr = Corriente(T=corr.T-0.1, P=101325., mezcla=mez)
                x.append(corr.Liquido.fraccion[0])
                y.append(corr.Gas.fraccion[0])
                # except:
                    # pass
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

    def readToStream(self, stream):
        id1 = stream.readInt32()
        id2 = stream.readInt32()
        len = stream.readInt32()
        x = []
        for i in range(len):
            x.append(stream.readFloat())
        y = []
        for i in range(len):
            y.append(stream.readFloat())
        self.plot(0, x, y)


class Plot_Distribucion(PlotWidget):
    title = translate("Plot_Distribucion", "Solid Distribution")
    def __init__(self, id, solido, parent=None):
        super().__init__(parent=parent)
        self.setWindowTitle(self.tr("Stream")+" "+str(id)+" - "+self.title)
        self.fill(solido)
        self.ax.legend(loc=4)
        self.ax.set_ylim(0, 1)
        self.ax.set_xlabel("Dp, "+unidades.Length(None).text("ParticleDiameter"),
                               horizontalalignment='right', size='12')
        self.ax.set_ylabel(self.tr("Accumulated fraction"),
                               horizontalalignment='right', size='12')


    def fill(self, solido):
        self.plot(solido.diametros, solido.fracciones_acumuladas, 'ro')
        for i in range(6):
            x, y, leyenda = solido.ajustar_distribucion(i)
            self.plot(x, y, label=leyenda)


_all = Binary_distillation, Plot_Distribucion


if __name__ == "__main__":
    import sys
    from lib.solids import Solid
    # from configparser import ConfigParser
    # configuracion=ConfigParser()
    # configuracion.read(config.conf_dir+"pychemqtrc_temporal")
    # print(dict(configuracion))

    app = QtWidgets.QApplication(sys.argv)
    # kw = {"indices": [5, 6, 7, 8, 10],
          # "nombres": ["i-butane", "n-butane", "i-pentane", "n-pentane", "n-hexane"]}
    # Dialog = Binary_distillation(**kw)

    distribucion = [[17.5, 0.02],
                    [22.4, 0.03],
                    [26.2, 0.05],
                    [31.8, 0.1],
                    [37, 0.1],
                    [42.4, 0.1],
                    [48, 0.1],
                    [54, 0.1],
                    [60, 0.1],
                    [69, 0.1],
                    [81.3, 0.1],
                    [96.5, 0.05],
                    [109, 0.03],
                    [127, 0.02]]
    diametros = []
    fracciones = []
    for diametro, fraccion in distribucion:
        diametros.append(diametro)
        fracciones.append(fraccion)

    solido = Solid(caudalSolido=[5], distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones)
    Dialog = Plot_Distribucion(638, solido)
    Dialog.show()
    sys.exit(app.exec())
