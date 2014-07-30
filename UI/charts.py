#!/usr/bin/python
# -*- coding: utf-8 -*-

###Módulo que define los dialogos de definición de gráficos

from PyQt4 import QtCore, QtGui
from scipy import arange, logspace, log10, pi, arctan, concatenate, r_
from scipy.optimize import fsolve
from matplotlib.patches import ConnectionPatch

from UI.widgets import Entrada_con_unidades
from lib.petro import Z_list
from lib.physics import f_list
from lib.plot import mpl
from lib import config

from equipment.UI_heatExchanger import chart as chartHE


class Standing_Katz(QtGui.QDialog):
    title=QtGui.QApplication.translate("pychemqt", "Standing Katz chart")
    def __init__(self, parent=None):
        super(Standing_Katz, self).__init__(parent)
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method:")),1,1)
        self.metodos=QtGui.QComboBox()
        self.metodos.addItem("Hall Yarborough")
        self.metodos.addItem("Dranchuk Abu-Kassem")
        self.metodos.addItem("Dranchuk Purvis Robinson")
        self.metodos.addItem("Beggs Brill")
        self.metodos.addItem("Sarem")
        self.metodos.addItem("Gopal")
        self.metodos.addItem("Papay")
        self.metodos.currentIndexChanged.connect(self.plot_Z)
        layout.addWidget(self.metodos,1,2)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Fixed),1, 3)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pr<sub>min</sub>")),1,4)
        self.Prmin=Entrada_con_unidades(float, spinbox=True, value=0.0, width=60, decimales=1, step=0.1)
        layout.addWidget(self.Prmin,1,5)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pr<sub>max</sub>")),1,6)
        self.Prmax=Entrada_con_unidades(float, spinbox=True, value=8.0, width=60, decimales=1, step=0.1)
        layout.addWidget(self.Prmax,1,7)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Tr")),1, 8)
        self.Tr=QtGui.QLineEdit(", ".join(str(i) for i in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.]))
        self.Tr.setMinimumWidth(400)
        layout.addWidget(self.Tr,1, 9)
        self.diagrama = mpl(self, dpi=90)
        layout.addWidget(self.diagrama,2,1,1,9)
        
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 5, 1, 1, 6)

        self.plot_Z(0)
        
    def plot_Z(self, indice):
        Prmin=self.Prmin.value
        Prmax=self.Prmax.value
        self.diagrama.config(self.Prmin.value, self.Prmax.value)
        
        try:
            Tr=self.Tr.text().split()
        except:
            Tr=[1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.]
        Z=Z_list[indice]
        
        P=arange(Prmin, Prmax, 0.1)
        for Tr in Tr:
            self.diagrama.plot(P, [Z(Tr, Pr) for Pr in P], "k")
        self.diagrama.axes2D.set_title(QtGui.QApplication.translate("pychemqt", "Standing and Katz compressivitity factors chart for natural gas"), size='12')
        self.diagrama.draw()
   

class Moody(QtGui.QDialog):
    title=QtGui.QApplication.translate("pychemqt", "Moody Diagram")
    def __init__(self, parent=None):
        super(Moody, self).__init__(parent)
        self.showMaximized()
        self.setWindowTitle(self.title)
        layout=QtGui.QGridLayout(self)
        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Method:")),1,1)
        self.metodos=QtGui.QComboBox()
        self.metodos.addItem("Colebrook")
        self.metodos.addItem("Chen (1979")
        self.metodos.addItem("Romeo (2002)")
        self.metodos.addItem("Goudar-Sonnad")
        self.metodos.addItem("Manadilli (1997)")
        self.metodos.addItem("Serghides")
        self.metodos.addItem("Churchill (1977)")
        self.metodos.addItem("Zigrang-Sylvester (1982)")
        self.metodos.addItem("Swamee-Jain (1976)")        
        
        self.metodos.currentIndexChanged.connect(self.cambiar)
        layout.addWidget(self.metodos,1,2)
        layout.setColumnStretch(3, 1)
        self.diagrama = mpl(self, dpi=90)
        layout.addWidget(self.diagrama,2,1,1,4)
        self.diagrama.fig.text(0.95, 0.4, QtGui.QApplication.translate("pychemqt", "Relative roughness")+", "+r"$r=\frac{\epsilon}{D}$", rotation=90, size='14')
        
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 3, 1, 1, 4)
        
        self.cambiar(0)

        
    def cambiar(self, int):
        self.diagrama.axes2D.clear()
        self.diagrama.axes2D.set_autoscale_on(False)
        self.diagrama.axes2D.set_xlabel(QtGui.QApplication.translate("pychemqt", "Reynolds number")+ ",  " +r"$Re=\frac{V\rho D}{\mu}$" , horizontalalignment='center', size='16')
        self.diagrama.axes2D.set_ylabel(QtGui.QApplication.translate("pychemqt", "Friction factor")+",  "+r"$f=\frac{2hDg}{LV^2}$", size='14')
        self.diagrama.axes2D.grid(True)
        self.diagrama.axes2D.set_xlim(600, 1e8)
        self.diagrama.axes2D.set_ylim(0.008, 0.1)
        self.diagrama.axes2D.set_xscale("log")
        self.diagrama.axes2D.set_yscale("log")
        self.diagrama.axes2D.set_title(QtGui.QApplication.translate("pychemqt", "Moody Diagram"), size='12')
        self.diagrama.axes2D.set_xticks([7e2, 8e2, 9e2, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7], minor=False)
        self.diagrama.axes2D.set_yticks([9e-3, 1e-2, 1.1e-2, 1.2e-2, 1.3e-2, 1.4e-2, 1.5e-2, 1.6e-2, 1.7e-2, 1.8e-2, 1.9e-2, 2e-2, 2.1e-2, 2.2e-2, 2.3e-2, 2.4e-2, 2.5e-2, 2.6e-2, 2.7e-2, 2.8e-2, 2.9e-2, 3e-2, 3.2e-2, 3.4e-2, 3.6e-2, 3.8e-2, 4e-2, 4.2e-2, 4.4e-2, 4.6e-2, 4.8e-2, 5e-2, 5.2e-2, 5.4e-2, 5.6e-2, 5.8e-2, 6e-2, 6.2e-2, 6.4e-2, 6.6e-2, 6.8e-2, 7e-2, 7.5e-2, 8e-2, 8.5e-2, 9e-2, 9.5e-2, 1e-1])
        self.diagrama.axes2D.set_yticks([], minor=True)
        self.diagrama.axes2D.set_yticklabels([9, r"$10^{-2}$", "", 1.2, "", 1.4, "", 1.6, "", 1.8, "", 2, "", "", "", "", 2.5, "", "", "", "", 3, "", "", "", "", 4, "", "", "", "", 5, "", "", "", "", 6, "", "", "", "", 7, "", 8, "", 9, "", r"$10^{-1}$"])
        self.__plot(int)
        self.diagrama.draw()


    def __plot(self, metodo=0, eD=[]):
        """Plot the Moody chart using the indicate method
        método de cálculo:
            0   -   Colebrook
            1   -   Chen (1979)
            2   -   Romeo (2002)
            3   -   Goudar-Sonnad
            4   -   Manadilli (1997)
            5   -   Serghides
            6   -   Churchill (1977)
            7   -   Zigrang-Sylvester (1982)
            8   -   Swamee-Jain (1976)")      
            
        eD: lista con las líneas de rugosidades relativas a dibujar
        Prmin: escala del eje x, minimo valor de Pr a representar
        Prmax: escala del eje y, maximo valor de Pr a representar
        """
        if not eD:
            eD=[0, 1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07]
        F=f_list[metodo]
        
        #laminar
        Re=[600, 2400]
        f=[64./R for R in Re]
        self.diagrama.axes2D.plot(Re, f, "k")
        #turbulento
        Re=logspace(log10(2400), 8, 50)
        for e in eD:
            self.diagrama.axes2D.plot(Re, [F(Rei, e) for Rei in Re], "k")
            self.diagrama.axes2D.annotate(config.representacion(e, tol=4.5), (Re[45], F(Re[45], e)), size="small", horizontalalignment="center", verticalalignment="bottom", rotation=arctan((log10(F(Re[47], e))-log10(F(Re[35], e)))/(log10(Re[47])-log10(Re[35])))*360/2/pi)

        #Transición
        f=[(1/(1.14-2*log10(3500/R)))**2 for R in Re]
        self.diagrama.axes2D.plot(Re, f, "k", lw=0.5, linestyle=":")
        
        self.diagrama.axes2D.add_artist(ConnectionPatch((600, 0.009), (2400, 0.009), "data", "data", arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.axes2D.add_artist(ConnectionPatch((2400, 0.009), (6000, 0.009), "data", "data", arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.axes2D.add_artist(ConnectionPatch((6000, 0.095), (40000, 0.095), "data", "data", arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.axes2D.add_artist(ConnectionPatch((40000, 0.095), (9.9e7, 0.095), "data", "data", arrowstyle="<|-|>", mutation_scale=20, fc="w"))
        self.diagrama.axes2D.text(15000, 0.094, QtGui.QApplication.translate("pychemqt", "Transition Zone"), size="small", verticalalignment="top", horizontalalignment="center")
        self.diagrama.axes2D.text(2e6, 0.094, QtGui.QApplication.translate("pychemqt", "Turbulent flux fully desarrolled"), size="small", verticalalignment="top", horizontalalignment="center")
        self.diagrama.axes2D.text(4000, 0.0091, QtGui.QApplication.translate("pychemqt", "Critic\nzone"), size="small", verticalalignment="bottom", horizontalalignment="center")
        self.diagrama.axes2D.text(1200, 0.0091, QtGui.QApplication.translate("pychemqt", "Laminar flux"), size="small", verticalalignment="bottom", horizontalalignment="center")






__all__={QtGui.QApplication.translate("pychemqt", "Petro"): (Standing_Katz, ), 
             QtGui.QApplication.translate("pychemqt", "Fluid Flow"): (Moody, ), 
             QtGui.QApplication.translate("pychemqt", "Heat Exchanger"): chartHE}


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = chartHE[2]()
    Dialog.show()
    sys.exit(app.exec_())

