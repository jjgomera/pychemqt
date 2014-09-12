#!/usr/bin/python
# -*- coding: utf-8 -*-

###Módulo en el que encuentran metodos relacionados con la psicrometria

import os
import cPickle
from ConfigParser import ConfigParser

from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib import pyplot
from scipy import arange, concatenate, pi, arctan, log, exp, linspace
from scipy.optimize import fsolve
from scipy.constants import lb

from lib.compuestos import Componente
from lib.psycrometry import Psychrometry, _Pbar, _height
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
        self.axes2D.figure.subplots_adjust(left=0.01, right=0.92, bottom=0.05, top=0.98)
        FigureCanvasQTAgg.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
        self.notes = []

    
    def config(self):
        self.axes2D.set_autoscale_on(False)
        self.axes2D.set_xlabel("Tdb, "+unidades.Temperature(None).text())
        self.axes2D.set_ylabel(QtGui.QApplication.translate("pychemqt", "Absolute humidity")+", "+unidades.Mass(None).text()+"/"+unidades.Mass(None).text())
        self.axes2D.yaxis.set_ticks_position("right")
        self.axes2D.yaxis.set_label_position("right")
        
        self.lx = self.axes2D.axhline(color='b')  # the horiz line
        self.ly = self.axes2D.axvline(color='b')  # the vert line
        
        tmin=unidades.Temperature(274).config()
        tmax=unidades.Temperature(329).config()

        self.axes2D.set_xlim(tmin, tmax)
        self.axes2D.set_ylim(0, 0.04)

    def plot(self, x, y, **kwargs):
        self.axes2D.plot(x, y, **kwargs)
    
    def showPointData(self, psyState):
        self.clearPointData()

        yi = 0.99
        for key in ("Tdb", "Tdp", "Twb", "HR", "w", "h", "V", "rho"):
            self.notes.append(self.axes2D.annotate(
                "%s: %s" %(key, psyState.__getattribute__(key).str), (0.01, yi),
                xycoords='axes fraction', size="small", va="center"))
            yi-=0.025
        self.draw()

    def clearPointData(self):
        while self.notes:
            anotation = self.notes.pop()
            anotation.remove()
        self.draw()


class UI_Psychrometry(QtGui.QDialog):
    """Aplicación grafica de psychrometria"""
    def __init__(self, parent=None):
        super(UI_Psychrometry, self).__init__(parent)
        self.showMaximized()
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/psychrometric.png")))
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Psychrometric chart"))
        self.AireHumedo=Psychrometry(1)
        
        layout = QtGui.QGridLayout(self)
        self.diagrama2D = Plot(self, dpi=90)
        self.diagrama2D.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        self.diagrama2D.fig.canvas.mpl_connect('button_press_event', self.click)
        layout.addWidget(self.diagrama2D,1,3,1,2)
        
       #Toolbox
        self.controles = QtGui.QWidget()
        self.controles.setSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        layoutToolbox = QtGui.QGridLayout(self.controles)
        self.checkPresion = QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Pressure"))
        layoutToolbox.addWidget(self.checkPresion,1,1,1,1)
        self.PresionDiagrama=Entrada_con_unidades(unidades.Pressure, value=101325)
        self.PresionDiagrama.valueChanged.connect(self.cambiarPresion)
        layoutToolbox.addWidget(self.PresionDiagrama,1,2,1,1)
        self.checkAltitud = QtGui.QRadioButton(QtGui.QApplication.translate("pychemqt", "Altitude"))
        layoutToolbox.addWidget(self.checkAltitud,2,1,1,1)
        self.Altitud=Entrada_con_unidades(unidades.Length, value=0)
        self.checkPresion.toggled.connect(self.PresionDiagrama.setEnabled)
        self.checkAltitud.toggled.connect(self.Altitud.setEnabled)
        self.Altitud.valueChanged.connect(self.cambiarAltitud)
        self.checkPresion.setChecked(True)
        self.Altitud.setEnabled(False)
        layoutToolbox.addWidget(self.Altitud,2,2,1,1)
        layoutToolbox.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed),3,1,1,2)
        self.label_11 = QtGui.QLabel()
        self.label_11.setText(QtGui.QApplication.translate("pychemqt", "Select point"))
        layoutToolbox.addWidget(self.label_11,4,1,1,2)
        self.variables=QtGui.QComboBox()
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, H absolute"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, H relative"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, T dew point"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, Tª wet bulb"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dew point, H relative"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "T dry bulb, Enthalpy"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo seco, Densidad"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, H absoluta"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, H relativa"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Entalpia"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Densidad"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª bulbo húmedo, Tª rocio"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, entalpía"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, entalpía"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H absoluta, densidad"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "H relativa, densidad"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, entalpía"))
        self.variables.addItem(QtGui.QApplication.translate("pychemqt", "Tª rocio, densidad"))
        layoutToolbox.addWidget(self.variables,5,1,1,2)
        
        self.label_12 = QtGui.QLabel("T db:")
        layoutToolbox.addWidget(self.label_12,6,1,1,1)
        self.EntradaTdb=Entrada_con_unidades(unidades.Temperature)
        self.EntradaTdb.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaTdb,6,2,1,1)
        self.label_13 = QtGui.QLabel("T wb:")
        layoutToolbox.addWidget(self.label_13,7,1,1,1)
        self.EntradaTwb=Entrada_con_unidades(unidades.Temperature)
        self.EntradaTwb.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaTwb,7,2,1,1)
        self.label_14 = QtGui.QLabel("T dp:")
        layoutToolbox.addWidget(self.label_14,8,1,1,1)
        self.EntradaTdp=Entrada_con_unidades(unidades.Temperature)
        self.EntradaTdp.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaTdp,8,2,1,1)
        self.label_15 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "H absolute:"))
        layoutToolbox.addWidget(self.label_15,9,1,1,1)
        self.EntradaH=Entrada_con_unidades(float)
        self.EntradaH.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaH,9,2,1,1)
        self.label_16 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "H relative:"))
        layoutToolbox.addWidget(self.label_16,10,1,1,1)
        self.EntradaHR=Entrada_con_unidades(float)
        self.EntradaHR.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaHR,10,2,1,1)
        self.label_17 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Volume"))
        layoutToolbox.addWidget(self.label_17,11,1,1,1)
        self.EntradaVolumen=Entrada_con_unidades(unidades.SpecificVolume)
        self.EntradaVolumen.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaVolumen,11,2,1,1)
        self.label_18 = QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Enthalpy"))
        layoutToolbox.addWidget(self.label_18,12,1,1,1)
        self.EntradaEntalpia=Entrada_con_unidades(unidades.Enthalpy)
        self.EntradaEntalpia.valueChanged.connect(self.NuevoPunto)
        layoutToolbox.addWidget(self.EntradaEntalpia,12,2,1,1)
        
        layoutToolbox.addItem(QtGui.QSpacerItem(20,20,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),13,1,1,2)
        layout.addWidget(self.controles,1,1,2,1)
        
        self.buttonShowToolbox = QtGui.QToolButton()
        self.buttonShowToolbox.setCheckable(True)
        self.buttonShowToolbox.toggled.connect(self.showToolBar)
        layout.addWidget(self.buttonShowToolbox,1,2,2,1)
        self.line = QtGui.QFrame()
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        layout.addWidget(self.line,1,3,2,1)
        
        self.progressBar=QtGui.QProgressBar()
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar,2,3)
        self.status = QtGui.QLabel()
        layout.addWidget(self.status,2,3)
        
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
        butonPNG = QtGui.QPushButton(QtGui.QIcon(os.environ["pychemqt"]+
            "images"+os.sep+"button"+os.sep+"image.png"), 
            QtGui.QApplication.translate("pychemqt", "Save as PNG"))
        self.buttonBox.addButton(butonPNG, QtGui.QDialogButtonBox.AcceptRole)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.accepted.connect(self.savePNG)
        layout.addWidget(self.buttonBox,2,4)
        
        self.showToolBar(False)
        self.plot()
        
    def savePNG(self):
        fname=unicode(QtGui.QFileDialog.getSaveFileName(
            self, QtGui.QApplication.translate("pychemqt", "Save chart to file"),
            "./", "Portable Network Graphics (*.png)"))
        self.diagrama2D.fig.savefig(fname, facecolor='#eeeeee')

    def showToolBar(self, checked):
        self.controles.setVisible(not checked)
        if checked:
            image = "arrow-right-double.png"
        else:
            image = "arrow-left-double.png"
        self.buttonShowToolbox.setIcon(QtGui.QIcon(os.environ["pychemqt"]+
            "images"+os.sep+"button"+os.sep+image))     
    
    def LineList(self, name, Preferences):
        if Preferences.getboolean("Psychr", name+"Custom"):
            t=[]
            for i in Preferences.get("Psychr", name+'List').split(','):
                if i:
                    t.append(float(i))
        else:
            start=Preferences.getfloat("Psychr", name+"Start")
            end=Preferences.getfloat("Psychr", name+"End")
            step=Preferences.getfloat("Psychr", name+"Step")
            t=arange(start, end, step)
        return t
        
    def calculate(self):
        Preferences=ConfigParser()
        Preferences.read(config.conf_dir+"pychemqtrc")
        
        data = {}
        
        t=self.LineList("isotdb", Preferences)
        
        # Saturation line
        Hs = []
        for ti in t:
            Hs.append(self.AireHumedo.Humedad_Absoluta(ti))
            self.progressBar.setValue(5*len(Hs)/len(t))
        data["t"] = t
        data["Hs"] = Hs
        
        # Humidity ratio lines
        hr = self.LineList("isohr", Preferences)
        Hr = {}
        cont = 0
        for i in hr:
            Hr[i] = []
            for H in Hs:
                Hr[i].append(H*i/100 for H in Hs)
                cont += 1
                self.progressBar.setValue(5+10*cont/len(hr)/len(Hs))

            Hr[i] = [H*i/100 for H in Hs]
        data["Hr"] = Hr
        
        # Twb
        lines = self.LineList("isotwb", Preferences)
        Twb = {}
        cont = 0
        for T in lines:
            H=concatenate((arange(self.AireHumedo.Humedad_Absoluta(T), 0, -0.001), [0.]))
            Tw = []
            for h in H:
                Tw.append(self.AireHumedo.Isoentalpica(T, h).config())
            cont += 1
            self.progressBar.setValue(15+75*cont/len(lines))
            Twb[T] = (H, Tw)
        data["Twb"] = Twb
        
        # v
        lines = self.LineList("isochor", Preferences)
        V = {}
        for cont, v in enumerate(lines):
            def f(Ts):
                return v-R_atml*Ts/self.AireHumedo.P.atm/self.AireHumedo.aire.M*(1+self.AireHumedo.Humedad_Absoluta(Ts)*self.AireHumedo.aire.M/self.AireHumedo.agua.M)
            ts=fsolve(f, 300)
            T=linspace(ts, self.AireHumedo.P.atm*v*self.AireHumedo.aire.M/R_atml, 50)
            Td=[unidades.Temperature(i).config() for i in T]
            H=[self.AireHumedo.Isocora_H(i, v) for i in T]
            self.progressBar.setValue(90+10*cont/len(lines))
            V[v] = (Td, H)
        data["v"] = V
        
        return data
    
    def drawlabel(self, name, Preferences, t, W, label, unit):
        if Preferences.getboolean("Psychr", name+"label"):
            tmin=unidades.Temperature(Preferences.getfloat("Psychr", "isotdbStart")).config()
            tmax=unidades.Temperature(Preferences.getfloat("Psychr", "isotdbEnd")).config()
            x = tmax-tmin
            wmin=Preferences.getfloat("Psychr", "isowStart")
            wmax=Preferences.getfloat("Psychr", "isowEnd")
            y = wmax-wmin

            i = 0
            for ti, wi in zip(t, W):
                if tmin <= ti <= tmax and wmin <= wi <= wmax:
                    i += 1
            label = str(label)
            if Preferences.getboolean("Psychr", name+"units"):
                label += unit
            pos = Preferences.getfloat("Psychr", name+"position")
            p = int(i*pos/100-1)
            rot=arctan((W[p]-W[p-1])/y/(t[p]-t[p-1])*x)*360/2/pi
            self.diagrama2D.axes2D.annotate(label, (t[p], W[p]),
                rotation=rot, size="small", ha="center", va="center")

    def plot(self):
        Preferences=ConfigParser()
        Preferences.read(config.conf_dir+"pychemqtrc")

        self.diagrama2D.axes2D.clear()
        self.diagrama2D.config()
        
        filename = config.conf_dir+"psy_%i.pkl" % self.AireHumedo.P
        if os.path.isfile(filename):
            with open(filename, "r") as archivo:
                data=cPickle.load(archivo)
                self.status.setText(QtGui.QApplication.translate("pychemqt", "Loading cached data..."))
                QtGui.QApplication.processEvents()
        else:
            self.progressBar.setVisible(True)
            self.status.setText(QtGui.QApplication.translate("pychemqt", "Calculating data, be patient..."))
            QtGui.QApplication.processEvents()
            data = self.calculate()
            cPickle.dump(data, open(filename, "w"))
            self.progressBar.setVisible(False)
        self.status.setText(QtGui.QApplication.translate("pychemqt", "Plotting..."))
        QtGui.QApplication.processEvents()

        tmin=unidades.Temperature(Preferences.getfloat("Psychr", "isotdbStart")).config()
        tmax=unidades.Temperature(Preferences.getfloat("Psychr", "isotdbEnd")).config()
        x = tmax-tmin
        
        wmin=Preferences.getfloat("Psychr", "isowStart")
        wmax=Preferences.getfloat("Psychr", "isowEnd")
        y = wmax-wmin
        
        t = [unidades.Temperature(ti).config() for ti in data["t"]]
        Hs = data["Hs"]
        format={}
        format["ls"]=Preferences.get("Psychr", "saturationlineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "saturationlineWidth")
        format["color"]=Preferences.get("Psychr", "saturationColor")
        format["marker"]=Preferences.get("Psychr", "saturationmarker")
        format["markersize"]=3
        self.diagrama2D.plot(t, Hs, **format)

        format={}
        format["ls"]=Preferences.get("Psychr", "isotdblineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "isotdblineWidth")
        format["color"]=Preferences.get("Psychr", "isotdbColor")
        format["marker"]=Preferences.get("Psychr", "isotdbmarker")
        format["markersize"]=3
        for i, T in enumerate(t):
            self.diagrama2D.plot([T, T], [0, Hs[i]], **format)

        H = self.LineList("isow", Preferences)
        th=[self.AireHumedo.Tdp(h).config() for h in H]
        format={}
        format["ls"]=Preferences.get("Psychr", "isowlineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "isowlineWidth")
        format["color"]=Preferences.get("Psychr", "isowColor")
        format["marker"]=Preferences.get("Psychr", "isowmarker")
        format["markersize"]=3
        for i, H in enumerate(H):
            self.diagrama2D.plot([th[i], tmax], [H, H], **format)

        format={}
        format["ls"]=Preferences.get("Psychr", "isohrlineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "isohrlineWidth")
        format["color"]=Preferences.get("Psychr", "isohrColor")
        format["marker"]=Preferences.get("Psychr", "isohrmarker")
        format["markersize"]=3
        for Hr, H0 in data["Hr"].iteritems():
            self.diagrama2D.plot(t, H0, **format)
            self.drawlabel("isohr", Preferences, t, H0, Hr, "%")
                    
        format={}
        format["ls"]=Preferences.get("Psychr", "isotwblineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "isotwblineWidth")
        format["color"]=Preferences.get("Psychr", "isotwbColor")
        format["marker"]=Preferences.get("Psychr", "isotwbmarker")
        format["markersize"]=3
        for T, (H, Tw) in data["Twb"].iteritems():
            self.diagrama2D.plot(Tw, H, **format)
            value = unidades.Temperature(T).config()
            txt = unidades.Temperature.text()
            self.drawlabel("isotwb", Preferences, Tw, H, value, txt)
  
        format={}
        format["ls"]=Preferences.get("Psychr", "isochorlineStyle")
        format["lw"]=Preferences.getfloat("Psychr", "isochorlineWidth")
        format["color"]=Preferences.get("Psychr", "isochorColor")
        format["marker"]=Preferences.get("Psychr", "isochormarker")
        format["markersize"]=3
        for v, (Td, H) in data["v"].iteritems():
            self.diagrama2D.plot(Td, H, **format)
            value = unidades.SpecificVolume(v).config()
            txt = unidades.SpecificVolume.text()
            self.drawlabel("isochor", Preferences, Td, H, value, txt)

        self.diagrama2D.draw()
        self.status.clear()

    def mouse_move(self, event):
        if event.xdata and event.ydata:
            self.ActualizarDatos(event.xdata, event.ydata)
            
    def click(self, event):
        if event.xdata and event.ydata:
            punto = self.createPoint(event.xdata, event.ydata)
            self.ActualizarEntradas(punto)
            self.diagrama2D.lx.set_ydata(event.ydata)
            self.diagrama2D.ly.set_xdata(event.xdata)
            self.diagrama2D.draw()

    def createPoint(self, x, y):
        tdb=unidades.Temperature(x, "conf")
        punto=self.AireHumedo.definirPunto(0, tdb=tdb, w=y)
        return punto

    def ActualizarDatos(self, x, y):
            punto = self.createPoint(x, y)
            if y<punto.Hs:
                self.diagrama2D.showPointData(punto)
            else:
                self.diagrama2D.clearPointData()
                    
    def ActualizarEntradas(self, punto):
        if punto.w<punto.Hs:
            self.EntradaTdb.setValue(punto.Tdb)
            self.EntradaTwb.setValue(punto.Twb)
            self.EntradaH.setValue(punto.w)
            self.EntradaHR.setValue(punto.HR)
            self.EntradaTdp.setValue(punto.Tdp)
            self.EntradaEntalpia.setValue(punto.h)
            self.EntradaVolumen.setValue(punto.V)
                
    def cambiarPresion(self, value):
        self.AireHumedo=Psychrometry(value/101325.)
        self.Altitud.setValue(_height(value))
        self.plot()
        self.diagrama2D.draw()
        
    def cambiarAltitud(self, value):
        presion=_Pbar(value)
        self.PresionDiagrama.setValue(presion)
        self.AireHumedo=Psychrometry(presion/101325.)
        self.plot()
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
            punto=self.AireHumedo.definirPunto(modo, tdb=self.EntradaTdb.value, 
                twb=self.EntradaTwb.value, tdp=self.EntradaTdp.value, 
                w=self.EntradaH.value, HR=self.EntradaHR.value,
                V=self.EntradaVolumen.value, h=self.EntradaEntalpia.value)
            self.ActualizarEntradas(punto)
            self.diagrama2D.lx.set_ydata(punto.w)
            self.diagrama2D.ly.set_xdata(punto.Tdb.config())
            self.diagrama2D.draw()


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    aireHumedo=UI_Psychrometry()
    aireHumedo.show()
    sys.exit(app.exec_())
