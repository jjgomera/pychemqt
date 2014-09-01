#!/usr/bin/python
# -*- coding: utf-8 -*-

###Módulo en el que encuentran metodos relacionados con la psicrometria

import os

from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib import pyplot
from scipy import arange, concatenate, pi, arctan, log, exp, linspace
from scipy.optimize import fsolve
from scipy.constants import lb

from lib.compuestos import Componente
from lib.psycrometry import Psychrometry
from lib import unidades, config
from lib.utilities import representacion
from lib.physics import R_atml
from UI.widgets import Entrada_con_unidades


class Plot(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = pyplot.Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        self.axes2D = self.fig.add_subplot(111)
        self.axes2D.figure.subplots_adjust(left=0.01, right=0.9, bottom=0.08, top=0.98)
        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    
    def config(self):
        self.axes2D.set_autoscale_on(False)
        self.axes2D.set_xlabel("Tdb"+unidades.Temperature(None).text())
        self.axes2D.set_ylabel(QtGui.QApplication.translate("pychemqt", "Absolute humidity")+", "+unidades.Mass(None).text()+"/"+unidades.Mass(None).text())
        self.axes2D.yaxis.set_ticks_position("right")
        self.axes2D.yaxis.set_label_position("right")
        
        self.lx = self.axes2D.axhline(color='b')  # the horiz line
        self.ly = self.axes2D.axvline(color='b')  # the vert line
        
        tmin=unidades.Temperature(274).config()
        tmax=unidades.Temperature(329).config()

        self.axes2D.set_xlim(tmin, tmax)
        self.axes2D.set_ylim(0, 0.04)

    def plot(self, x, y, color="#000000", grosor=1, linestyle="-"):
        self.axes2D.plot(x, y, color=color, lw=grosor, ls=linestyle)


        

class UI_Psychrometry(QtGui.QMainWindow):
    """Aplicación grafica de psychrometria"""
    def __init__(self, parent=None):
        super(UI_Psychrometry, self).__init__(parent)
        self.resize(800, 600)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/psychrometric.png")))
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Psychrometric chart"))
        self.AireHumedo=Psychrometry(1)
        
        #menu
        self.menubar = QtGui.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0,0,700,30))
        self.menuArchivo = QtGui.QMenu(self.menubar)
        self.menuArchivo.setTitle(QtGui.QApplication.translate("pychemqt", "File"))
        self.setMenuBar(self.menubar)
        self.actionSalir = QtGui.QAction(self)
        self.actionSalir.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/exit.png")))
        self.actionSalir.setShortcut("Alt+F4")
        self.actionSalir.setText(QtGui.QApplication.translate("pychemqt", "Close"))
        self.actionSalir.triggered.connect(self.close)
        self.menuArchivo.addAction(self.actionSalir)
        self.menubar.addAction(self.menuArchivo.menuAction())
        self.statusbar = QtGui.QStatusBar()
        self.statusbar.setSizeGripEnabled(False)
        self.setStatusBar(self.statusbar)
        self.labelPresion=QtGui.QLabel()
        self.labelPresion.setFrameShadow(QtGui.QFrame.Sunken)
        self.labelPresion.setFrameShape(QtGui.QFrame.WinPanel)
        self.labelPresion.setText(representacion(self.AireHumedo.P.config())+" "+unidades.Pressure(None).text())
        self.statusbar.addPermanentWidget(self.labelPresion)

        self.centralWidget=QtGui.QWidget()
        self.setCentralWidget(self.centralWidget)
        self.layout = QtGui.QGridLayout(self.centralWidget)
        self.diagrama2D = Plot(self, dpi=90)
        self.diagrama2D.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        self.diagrama2D.fig.canvas.mpl_connect('button_press_event', self.click)
        self.layout.addWidget(self.diagrama2D,0,1,1,13)
        self.checkMouse=QtGui.QRadioButton()
        self.checkMouse.setChecked(True)
        self.checkMouse.setText(QtGui.QApplication.translate("pychemqt", "Float"))
        self.layout.addWidget(self.checkMouse,1,1,1,1)
        self.checkPoint=QtGui.QRadioButton()
        self.checkPoint.setText(QtGui.QApplication.translate("pychemqt", "Click"))
        self.layout.addWidget(self.checkPoint,2,1,1,1)
        self.layout.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),1,2,2,1)
        
        self.layout.addWidget(QtGui.QLabel("T<sub>db</sub>"),1,3,1,1)
        self.tdb=Entrada_con_unidades(unidades.Temperature, readOnly=True, boton=False, width=50, decimales=2, frame=False)
        self.layout.addWidget(self.tdb,1,4,1,1)
        self.layout.addWidget(QtGui.QLabel("T<sub>wb</sub>"),2,3,1,1)
        self.twb=Entrada_con_unidades(unidades.Temperature, readOnly=True, boton=False, width=50, decimales=2, frame=False)
        self.layout.addWidget(self.twb,2,4,1,1)
        
        self.layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "T<sub>dew</sub>")),1,6,1,1)
        self.trocio=Entrada_con_unidades(unidades.Temperature, readOnly=True, boton=False, width=50, decimales=2, frame=False)
        self.layout.addWidget(self.trocio,1,7,1,1)
        self.layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Pressure")),2,6,1,1)
        self.presion=Entrada_con_unidades(unidades.Pressure, readOnly=True, boton=False, width=50, frame=False)
        self.layout.addWidget(self.presion,2,7,1,1)
        
        self.layout.addWidget(QtGui.QLabel("H<sub>abs</sub>:"),1,9,1,1)
        self.H=Entrada_con_unidades(float, readOnly=True, width=50, frame=False)
        self.layout.addWidget(self.H,1,10,1,1)
        self.layout.addWidget(QtGui.QLabel("H<sub>rel</sub>"),2,9,1,1)
        self.HR=Entrada_con_unidades(float, readOnly=True, width=50, frame=False)
        self.layout.addWidget(self.HR,2,10,1,1)
        
        self.layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Volume")),1,12,1,1)
        self.volumen=Entrada_con_unidades(unidades.SpecificVolume, boton=False, readOnly=True, width=50, frame=False)
        self.layout.addWidget(self.volumen,1,13,1,1)
        self.layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Enthalpy")),2,12,1,1)
        self.entalpia=Entrada_con_unidades(unidades.Enthalpy, boton=False, readOnly=True, width=50, frame=False)
        self.layout.addWidget(self.entalpia,2,13,1,1)
        
       #Toolbox
        self.dockWidget_Mouse = QtGui.QDockWidget()
        self.controles = QtGui.QWidget()
        self.gridLayout_1 = QtGui.QGridLayout(self.controles)
        self.checkPresion = QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Pressure"))
        self.gridLayout_1.addWidget(self.checkPresion,1,1,1,1)
        self.PresionDiagrama=Entrada_con_unidades(unidades.Pressure, value=101325, width=60)
        self.PresionDiagrama.valueChanged.connect(self.cambiarPresion)
        self.gridLayout_1.addWidget(self.PresionDiagrama,1,2,1,1)
        self.checkAltitud = QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Altitude"))
        self.gridLayout_1.addWidget(self.checkAltitud,2,1,1,1)
        self.Altitud=Entrada_con_unidades(unidades.Length, value=0.0, min=-1e6, width=60, decimales=2)
        self.checkPresion.toggled.connect(self.PresionDiagrama.setEnabled)
        self.checkAltitud.toggled.connect(self.Altitud.setEnabled)
        self.Altitud.valueChanged.connect(self.cambiarAltitud)
        self.checkPresion.setChecked(True)
        self.Altitud.setEnabled(False)
        self.Altitud.setValue(0)
        self.gridLayout_1.addWidget(self.Altitud,2,2,1,1)
        self.gridLayout_1.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)
        self.label_11 = QtGui.QLabel()
        self.label_11.setText(QtGui.QApplication.translate("pychemqt", "Select point"))
        self.gridLayout_1.addWidget(self.label_11,4,1,1,2)
        self.variables=QtGui.QComboBox()
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, H absolute"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, H relative"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, T dew point"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, Tª wet bulb"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dew point, H relative"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, Enthalpy"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo seco, Densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, H absoluta"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, H relativa"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Entalpia"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Tª rocio"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, densidad"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, entalpía"))
#        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, densidad"))
        self.gridLayout_1.addWidget(self.variables,5,1,1,2)
        
        self.label_12 = QtGui.QLabel("T db:")
        self.gridLayout_1.addWidget(self.label_12,6,1,1,1)
        self.EntradaTdb=Entrada_con_unidades(unidades.Temperature, width=60)
        self.EntradaTdb.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaTdb,6,2,1,1)
        self.label_13 = QtGui.QLabel("T wb:")
        self.gridLayout_1.addWidget(self.label_13,7,1,1,1)
        self.EntradaTwb=Entrada_con_unidades(unidades.Temperature, width=60)
        self.EntradaTwb.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaTwb,7,2,1,1)
        self.label_14 = QtGui.QLabel("T dp:")
        self.gridLayout_1.addWidget(self.label_14,8,1,1,1)
        self.EntradaTdp=Entrada_con_unidades(unidades.Temperature, width=60)
        self.EntradaTdp.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaTdp,8,2,1,1)
        self.label_15 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "H absolute:"))
        self.gridLayout_1.addWidget(self.label_15,9,1,1,1)
        self.EntradaH=Entrada_con_unidades(float)
        self.EntradaH.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaH,9,2,1,1)
        self.label_16 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "H relative:"))
        self.gridLayout_1.addWidget(self.label_16,10,1,1,1)
        self.EntradaHR=Entrada_con_unidades(float)
        self.EntradaHR.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaHR,10,2,1,1)
        self.label_17 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Volume"))
        self.gridLayout_1.addWidget(self.label_17,11,1,1,1)
        self.EntradaVolumen=Entrada_con_unidades(unidades.SpecificVolume, width=60)
        self.EntradaVolumen.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaVolumen,11,2,1,1)
        self.label_18 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Enthalpy"))
        self.gridLayout_1.addWidget(self.label_18,12,1,1,1)
        self.EntradaEntalpia=Entrada_con_unidades(unidades.Enthalpy, width=60)
        self.EntradaEntalpia.valueChanged.connect(self.NuevoPunto)
        self.gridLayout_1.addWidget(self.EntradaEntalpia,12,2,1,1)
        
        self.gridLayout_1.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),13,1,1,2)

        self.dockWidget_Mouse.setWidget(self.controles)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.dockWidget_Mouse)
        self.dockWidget_Mouse.setFeatures(QtGui.QDockWidget.NoDockWidgetFeatures)
        self.dockWidget_Mouse.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)

        self.plot()
        
    def plot(self):
        self.diagrama2D.axes2D.clear()
        self.diagrama2D.config()
        tmin=round(unidades.Temperature(274).config())
        tmax=round(unidades.Temperature(330).config())
        t=arange(tmin, tmax, 1.)
        T=[unidades.Temperature(i, unidades.Temperature.func()) for i in t]
        Hs=[self.AireHumedo.Humedad_Absoluta(i) for i in T]
        self.diagrama2D.plot(t, Hs, color="#000000")
        for i in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            H0=[H*i/100 for H in Hs]
            self.diagrama2D.plot(t, H0, grosor=0.5, linestyle="--", color="#000000")
            self.diagrama2D.axes2D.annotate(str(i)+"%", (t[35], H0[35]), rotation=arctan((H0[35]-H0[34])/0.04/(t[35]-t[34])*56)*360/2/pi, size="small", horizontalalignment="center", verticalalignment="center")
        tredondeada=arange
        
        for i, T in enumerate(t):
            self.diagrama2D.plot([T, T], [0, Hs[i]], grosor=0.5, linestyle=":", color="#000000")
            
        H=arange(0, 0.04, 0.001)
        t=[self.AireHumedo.Tdp(h).config() for h in H]
        for i, H in enumerate(H):
            self.diagrama2D.plot([t[i], tmax], [H, H], grosor=0.5, linestyle=":", color="#000000")
            
        for T in arange(250, 320, 1.):
            H=concatenate((arange(self.AireHumedo.Humedad_Absoluta(T), 0, -0.001), [0.]))
            Tw=[self.AireHumedo.Isoentalpica(T, h).config() for h in H]
            self.diagrama2D.plot(Tw, H, linestyle=":", color="#652A00")
  

        for v in arange(0.8, 1, 0.01):
            def f(Ts):
                return v-R_atml*Ts/self.AireHumedo.P.atm/self.AireHumedo.aire.M*(1+self.AireHumedo.Humedad_Absoluta(Ts)*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            ts=fsolve(f, 300)
            T=linspace(ts, self.AireHumedo.P.atm*v*self.AireHumedo.aire.M/R_atml, 50)
            Td=[unidades.Temperature(i).config() for i in T]
            H=[self.AireHumedo.Isocora_H(i, v) for i in T]
            self.diagrama2D.plot(Td, H, grosor=0.8, linestyle=":", color="green")
#            self.diagrama2D.axes2D.annotate(representacion(unidades.SpecificVolume(v).config())+config.Configuracion("SpecificVolume").text(), (Td[5], H[5]), rotation=arctan((H[5]-H[4])/0.04/(Td[5]-Td[4])*56)*360/2/pi, size="small", horizontalalignment="center", verticalalignment="center")

    def mouse_move(self, event):
        if event.xdata and event.ydata and self.checkMouse.isChecked():
            self.ActualizarDatos(event.xdata, event.ydata)
            
    def click(self, event):
        if event.xdata and event.ydata:
            self.ActualizarDatos(event.xdata, event.ydata)
            self.diagrama2D.lx.set_ydata(event.ydata)
            self.diagrama2D.ly.set_xdata(event.xdata)
            self.diagrama2D.draw()

    def ActualizarDatos(self, x, y):
            tdb=unidades.Temperature(x, unidades.Temperature(None).func())
            H=y
            if x and y:
                punto=self.AireHumedo.definirPunto(0, tdb=tdb, H=y)
                if H<punto.Hs:
                    self.tdb.setValue(punto.Tdb)
                    self.twb.setValue(punto.Twb)
                    self.H.setValue(punto.H)
                    self.HR.setValue(punto.HR)
                    self.trocio.setValue(punto.Tdp)
                    self.presion.setValue(punto.P)
                    self.entalpia.setValue(punto.entalpia)
                    self.volumen.setValue(punto.volumen)
                else:
                    self.tdb.clear()
                    self.twb.clear()
                    self.H.clear()
                    self.HR.clear()
                    self.trocio.clear()
                    self.presion.clear()
                    self.entalpia.clear()
                    self.volumen.clear()
                    
    def ActualizarEntradas(self, punto):
        if punto.H<punto.Hs:
            self.EntradaTdb.setValue(punto.Tdb)
            self.EntradaTwb.setValue(punto.Twb)
            self.EntradaH.setValue(punto.H)
            self.EntradaHR.setValue(punto.HR)
            self.EntradaTdp.setValue(punto.Tdp)
            self.EntradaEntalpia.setValue(punto.entalpia)
            self.EntradaVolumen.setValue(punto.volumen)
                
    def cambiarPresion(self, value):
        self.AireHumedo=Psychrometry(value.atm)
        self.Altitud.setValue(self.AireHumedo.altura())
        self.plot()
        self.labelPresion.setText(representacion(self.AireHumedo.P.config())+" "+unidades.Pressure(None).text())
        self.diagrama2D.draw()
        
    def cambiarAltitud(self, value):
        presion=self.AireHumedo.Presion(value)
        self.AireHumedo=Psychrometry(presion.atm)
        self.PresionDiagrama.setValue(presion)
        self.plot()
        self.labelPresion.setText(representacion(self.AireHumedo.P.config())+" "+unidades.Pressure(None).text())
        self.diagrama2D.draw()
        
    def NuevoPunto(self):
        modo=self.variables.currentIndex()
        if modo==0:
            calcular=self.EntradaTdb.value and self.EntradaH.value
        elif modo==1:
            calcular=self.EntradaTdb.value and self.EntradaHR.value
        elif modo==2:
            calcular=self.EntradaTdb.value and self.EntradaTdp.value
        elif modo==3:
            calcular=self.EntradaTdb.value and self.EntradaTwb.value
        elif modo==4:
            calcular=self.EntradaTdp.value and self.EntradaHR.value
        elif modo==5:
            calcular=self.EntradaTdb.value and self.EntradaEntalpia.value
            
        if calcular:
            punto=self.AireHumedo.definirPunto(self.variables.currentIndex(), tdb=self.EntradaTdb.value, twb=self.EntradaTwb.value, tdp=self.EntradaTdp.value, H=self.EntradaH.value, HR=self.EntradaHR.value, volumen=self.EntradaVolumen.value, entalpia=self.EntradaEntalpia.value)
            self.ActualizarEntradas(punto)
            self.ActualizarDatos(punto.Tdb, punto.H)
            self.diagrama2D.lx.set_ydata(punto.H)
            self.diagrama2D.ly.set_xdata(punto.Tdb)
            self.diagrama2D.draw()


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    aireHumedo=UI_Psychrometry()
    aireHumedo.show()
    sys.exit(app.exec_())
