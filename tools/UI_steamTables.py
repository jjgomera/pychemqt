#!/usr/bin/python
# -*- coding: utf-8 -*-

from csv import writer
from configparser import ConfigParser
import os

from PyQt5 import QtCore, QtGui, QtWidgets

#from freesteam import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from numpy import meshgrid, zeros, arange, linspace, concatenate, max, min, transpose, logspace, log, arctan, pi

from lib import unidades, config
from lib.utilities import colors, representacion


class Plot(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100, dim=3):
        self.fig = pyplot.Figure(figsize=(width, height), dpi=dpi)
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)
        self.dim=dim
        
        if dim==2:
            self.axes2D = self.fig.add_subplot(111)
            self.axes2D.figure.subplots_adjust(left=0.09, right=0.98, bottom=0.08, top=0.98)
        else:
            self.axes3D = Axes3D(self.fig)
        
        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)


    def plot_sat(self, xsat, ysat, zsat=0):
        """Método que dibuja la línea de saturación"""
        if self.dim==3:
            self.satliq=self.axes3D.plot3D(xsat[0], ysat[0], zsat[0],'k-')
            self.satgas=self.axes3D.plot3D(xsat[1], ysat[1], zsat[1],'k-')
        else:
            self.satliq=self.axes2D.plot(xsat[0], ysat[0],'k-')
            self.satgas=self.axes2D.plot(xsat[1], ysat[1],'k-')

    def plot_3D(self, labels, xdata, ydata, zdata):
        """Método que dibuja la matriz de datos"""
        self.axes3D.clear()
        self.axes3D.plot_wireframe(xdata, ydata, zdata, rstride=1, cstride=1)
        self.axes3D.set_xlabel(labels[0])
        self.axes3D.set_ylabel(labels[1])
        self.axes3D.set_zlabel(labels[2])
        self.axes3D.mouse_init()
        
    def plot_2D(self, labels, bool):
        self.axes2D.clear()
        self.axes2D.grid(bool)        
        self.axes2D.axes.set_xlabel(labels[0], size='12')
        self.axes2D.axes.set_ylabel(labels[1], size='12')
        
    def plot_labels(self, tipo, x, y, label, angle=0):
        linea=[]
        for i in range(len(label)):
            linea.append(self.axes2D.axes.annotate(label[i], (x[i], y[i]), rotation=angle[i], size="xx-small", horizontalalignment="center", verticalalignment="center"))
        if tipo =="T":
            self.Isoterma_label=linea
        elif tipo =="P":
            self.Isobara_label=linea
        elif tipo =="V":
            self.Isocora_label=linea
        elif tipo =="S":
            self.Isoentropica_label=linea
        elif tipo =="H":
            self.Isoentalpica_label=linea
        elif tipo =="X":
            self.IsoX_label=linea

    def plot_puntos(self, x, y, z=0):
        """Método que dibuja puntos individuales"""
        colores=colors(len(x))
        self.puntos=[]
        if self.dim==3:
            for i in range(len(x)):
                self.puntos.append(self.axes3D.plot3D([x[i]], [y[i]], [z[i]], color=colores[i], marker="o"))
        else:
            for i in range(len(x)):
                self.puntos.append(self.axes2D.plot([x[i]], [y[i]], color=colores[i], marker="o"))
            
    def plot_Isolinea(self, tipo, x, y, z=0, color="#000000", grosor=1, estilo=0):
        """Método que dibuja las isolineas"""
        if estilo==0:
            linestyle='-'
        elif estilo==1:
            linestyle='--'
        elif estilo==2:
            linestyle='-.'
        elif estilo==3:
            linestyle=':'
        
        linea=[]
        if self.dim==3:
            for i in range(len(x)):
                linea.append(self.axes3D.plot3D(x[i], y[i], z[i], color, lw=grosor, ls=linestyle))
        else:
            for i in range(len(x)):
                linea.append(self.axes2D.plot(x[i], y[i], color, lw=grosor, ls=linestyle))

        if tipo =="T":
            self.Isoterma=linea
        elif tipo =="P":
            self.Isobara=linea
        elif tipo =="V":
            self.Isocora=linea
        elif tipo =="S":
            self.Isoentropica=linea
        elif tipo =="H":
            self.Isoentalpica=linea
        elif tipo =="X":
            self.IsoX=linea

class Ventana_Lista_Puntos(QtWidgets.QDialog):
    """Dialogo que muestra la lísta de puntos especificados por el usuario así como sus propiedades"""
    def __init__(self, puntos, parent=None):
        super(Ventana_Lista_Puntos, self).__init__(parent)
        self.setWindowTitle(QtWidgets.QApplication.translate("SteamTables", "Puntos individuales", None))
        self.puntos=puntos
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.listWidget = QtWidgets.QListWidget()
        self.gridLayout.addWidget(self.listWidget, 0, 0, 3, 1)
        self.tablaPropiedades=QtWidgets.QTableWidget()
        self.tablaPropiedades.setVisible(False)
        self.tablaPropiedades.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tablaPropiedades.setRowCount(16)
        self.tablaPropiedades.setVerticalHeaderItem(0, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Presión", None)+", %s" %config.Configuracion("Pressure").text()))
        self.tablaPropiedades.setVerticalHeaderItem(1, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None)+", %s" %config.Configuracion("Temperature").text()))
        self.tablaPropiedades.setVerticalHeaderItem(2, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Volumen", None)+", %s" %config.Configuracion("SpecificVolume").text()))
        self.tablaPropiedades.setVerticalHeaderItem(3, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Entalpía", None)+", %s" %config.Configuracion("Enthalpy").text()))
        self.tablaPropiedades.setVerticalHeaderItem(4, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Entropía", None)+", %s" %config.Configuracion("SpecificHeat", "Entropy").text()))
        self.tablaPropiedades.setVerticalHeaderItem(5, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None)))
        self.tablaPropiedades.setVerticalHeaderItem(6, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Energía interna", None)+", %s" %config.Configuracion("Enthalpy").text()))
        self.tablaPropiedades.setVerticalHeaderItem(7, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None)+", %s" %config.Configuracion("Enthalpy").text()))
        self.tablaPropiedades.setVerticalHeaderItem(8, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None)+", %s" %config.Configuracion("Enthalpy").text()))
        self.tablaPropiedades.setVerticalHeaderItem(9, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Densidad", None)+", %s" %config.Configuracion("Density").text()))
        self.tablaPropiedades.setVerticalHeaderItem(10, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Cp", None)+", %s" %config.Configuracion("SpecificHeat").text()))
        self.tablaPropiedades.setVerticalHeaderItem(11, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Cv", None)+", %s" %config.Configuracion("SpecificHeat").text()))
        self.tablaPropiedades.setVerticalHeaderItem(12, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Conductividad", None)+", %s" %config.Configuracion("ThermalConductivity").text()))
        self.tablaPropiedades.setVerticalHeaderItem(13, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None)+", %s" %config.Configuracion("Viscosity").text()))
        self.tablaPropiedades.setVerticalHeaderItem(14, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None)+", %s" %config.Configuracion("Speed").text()))
        self.tablaPropiedades.setVerticalHeaderItem(15, QtWidgets.QTableWidgetItem(QtWidgets.QApplication.translate("SteamTables", "Región", None)))
        self.gridLayout.addWidget(self.tablaPropiedades, 3, 0, 1, 2)
        self.botonBorrar = QtWidgets.QPushButton()
        self.botonBorrar.setText(QtWidgets.QApplication.translate("SteamTables", "Borrar", None))
        self.botonBorrar.clicked.connect(self.on_botonBorrar_clicked)
        self.gridLayout.addWidget(self.botonBorrar, 0, 1, 1, 1)
        self.botonPropiedades = QtWidgets.QPushButton()
        self.botonPropiedades.setCheckable(True)
        self.botonPropiedades.setText(QtWidgets.QApplication.translate("SteamTables", "Propiedades", None))
        self.botonPropiedades.toggled.connect(self.tablaPropiedades.setVisible)
        self.gridLayout.addWidget(self.botonPropiedades, 1, 1, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok, QtCore.Qt.Vertical)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.gridLayout.addWidget(self.buttonBox, 2, 1, 1, 1)
        self.setFixedSize(500, self.minimumHeight())
        
        self.colores=colors(len(self.puntos))
        self.rellenarLista()
        self.rellenarTabla()

    def rellenarLista(self):
        self.listWidget.clear()
        for i, punto in enumerate(self.puntos):
            self.listWidget.addItem(str(i+1)+" - "+representacion(unidades.Temperature(punto.T).config())+" "+config.Configuracion("Temperature").text()+", "+representacion(unidades.Pressure(punto.p).config())+" "+config.Configuracion("Pressure").text()+", x="+representacion(punto.x))
            pixmap=QtGui.QPixmap(10, 10)
            pixmap.fill(QtGui.QColor(self.colores[i]))
            self.listWidget.item(i).setIcon(QtGui.QIcon(pixmap))
            
    def rellenarTabla(self):
        if len(self.puntos)==0:
            self.tablaPropiedades.setFixedHeight(404)
        elif len(self.puntos)<5:
            self.tablaPropiedades.setFixedHeight(428)
        else:
            self.tablaPropiedades.setFixedHeight(444)
        self.tablaPropiedades.setColumnCount(len(self.puntos))
        for i in range(16):
            self.tablaPropiedades.setRowHeight(i,25)
        for i, punto in enumerate(self.puntos):
            pixmap=QtGui.QPixmap(10, 10)
            pixmap.fill(QtGui.QColor(self.colores[i]))
            self.tablaPropiedades.setHorizontalHeaderItem(i, QtWidgets.QTableWidgetItem(QtGui.QIcon(pixmap), str(i+1)))
            self.tablaPropiedades.setItem(0, i, QtWidgets.QTableWidgetItem(representacion(unidades.Pressure(punto.p).config())))
            self.tablaPropiedades.setItem(1, i, QtWidgets.QTableWidgetItem(representacion(unidades.Temperature(punto.T).config())))
            self.tablaPropiedades.setItem(2, i, QtWidgets.QTableWidgetItem(representacion(unidades.SpecificVolume(punto.v).config())))
            self.tablaPropiedades.setItem(3, i, QtWidgets.QTableWidgetItem(representacion(unidades.Enthalpy(punto.h).config())))
            self.tablaPropiedades.setItem(4, i, QtWidgets.QTableWidgetItem(representacion(unidades.SpecificHeat(punto.s).config("Entropy"))))
            self.tablaPropiedades.setItem(5, i, QtWidgets.QTableWidgetItem(representacion(punto.x)))
            self.tablaPropiedades.setItem(6, i, QtWidgets.QTableWidgetItem(representacion(unidades.Enthalpy(punto.u).config())))
            self.tablaPropiedades.setItem(7, i, QtWidgets.QTableWidgetItem(representacion(unidades.Enthalpy(punto.s).config())))
            self.tablaPropiedades.setItem(8, i, QtWidgets.QTableWidgetItem(representacion(unidades.Enthalpy(punto.s).config())))
            self.tablaPropiedades.setItem(9, i, QtWidgets.QTableWidgetItem(representacion(unidades.Density(punto.rho).config())))
            self.tablaPropiedades.setItem(10, i, QtWidgets.QTableWidgetItem(representacion(unidades.SpecificHeat(punto.cp).config())))
            self.tablaPropiedades.setItem(11, i, QtWidgets.QTableWidgetItem(representacion(unidades.SpecificHeat(punto.cv).config())))
            self.tablaPropiedades.setItem(12, i, QtWidgets.QTableWidgetItem(representacion(unidades.ThermalConductivity(punto.k).config())))
            self.tablaPropiedades.setItem(13, i, QtWidgets.QTableWidgetItem(representacion(unidades.Viscosity(punto.mu).config())))
            if punto.region !='\x04' and  punto.region !='\x03':
                self.tablaPropiedades.setItem(14, i, QtWidgets.QTableWidgetItem(representacion(unidades.Speed(punto.w).config())))
            else:
                self.tablaPropiedades.setItem(14, i, QtWidgets.QTableWidgetItem('nan'))
            if punto.region =='\x01':
                self.tablaPropiedades.setItem(15, i, QtWidgets.QTableWidgetItem('1'))
            elif punto.region =='\x02':
                self.tablaPropiedades.setItem(15, i, QtWidgets.QTableWidgetItem('2'))
            elif punto.region =='\x03':
                self.tablaPropiedades.setItem(15, i, QtWidgets.QTableWidgetItem('3'))
            elif punto.region =='\x04':
                self.tablaPropiedades.setItem(15, i, QtWidgets.QTableWidgetItem('4'))
            elif punto.region =='\x05':
                self.tablaPropiedades.setItem(15, i, QtWidgets.QTableWidgetItem('5'))
            for j in range(16):
                self.tablaPropiedades.item(j, i).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                
        self.tablaPropiedades.resizeColumnsToContents()

    def on_botonBorrar_clicked(self):
        """Borra el punto seleccionado de la lista"""
        if self.listWidget.currentRow()>=0:
            self.puntos.pop(self.listWidget.currentRow())
        self.rellenarLista()
        self.rellenarTabla()



class Ui_SteamTables(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(Ui_SteamTables, self).__init__(parent)
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/steamTables.png")))
        self.setWindowTitle(QtWidgets.QApplication.translate("SteamTables", "Tablas de Vapor", None))
        self.Config=ConfigParser()
        self.Config.read("UI_steamTablesrc")

        self.showMaximized()
        self.centralwidget = QtWidgets.QWidget()
        self.setCentralWidget(self.centralwidget)
        
        #menus
        self.menubar = QtWidgets.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0,0,700,30))
        self.menuArchivo = QtWidgets.QMenu(self.menubar)
        self.menuArchivo.setTitle(QtWidgets.QApplication.translate("SteamTables", "Archivo", None))
        self.menuGrafico = QtWidgets.QMenu(self.menubar)
        self.menuGrafico.setTitle(QtWidgets.QApplication.translate("SteamTables", "Gráfico", None))
        self.menuAyuda = QtWidgets.QMenu(self.menubar)
        self.menuAyuda.setTitle(QtWidgets.QApplication.translate("SteamTables", "Ayuda", None))
        self.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusbar)
        self.progresbar=QtWidgets.QProgressBar()
        self.progresbar.setVisible(False)
        self.progresbar.setFixedSize(80, 15)
        self.statusbar.addPermanentWidget(self.progresbar)
        self.Preferencias = QtWidgets.QAction(self)
        self.Preferencias.setText(QtWidgets.QApplication.translate("SteamTables", "Preferencias", None))
        self.Preferencias.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/configure.png")))
        self.Preferencias.triggered.connect(self.preferencias)
        self.actionCSV = QtWidgets.QAction(self)
        self.actionCSV.setText(QtWidgets.QApplication.translate("SteamTables", "Guardar", None))
        self.actionCSV.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/filesave.png")))
        self.actionCSV.setShortcut("Ctrl+E")
        self.actionCSV.setEnabled(False)
        self.actionCSV.triggered.connect(self.exporttoCSV)
        self.actionSalir = QtWidgets.QAction(self)
        self.actionSalir.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/button/exit.png")))
        self.actionSalir.setShortcut("Alt+F4")
        self.actionSalir.setText(QtWidgets.QApplication.translate("SteamTables", "Salir", None))
        self.actionSalir.triggered.connect(self.close)
        self.actionTipoGrafico=QtWidgets.QActionGroup(self)
        self.action2D = QtWidgets.QAction(self.actionTipoGrafico)
        self.action2D.setCheckable(True)
        self.action2D.setShortcut("Ctrl+2")
        self.action2D.setText(QtWidgets.QApplication.translate("SteamTables", "Gráfico 2D", None))
        self.action2D.toggled.connect(self.d2)
        self.action3D = QtWidgets.QAction(self.actionTipoGrafico)
        self.action3D.setCheckable(True)
        self.action3D.setChecked(True)
        self.action3D.setShortcut("Ctrl+3")
        self.action3D.setText(QtWidgets.QApplication.translate("SteamTables", "Gráfico 3D", None))
        self.actionMostrarBarra = QtWidgets.QAction(self)
        self.actionMostrarBarra.setCheckable(True)
        self.actionMostrarBarra.setChecked(False)
        self.actionMostrarBarra.setShortcut("Ctrl+T")
        self.actionMostrarBarra.triggered.connect(self.mostrarBarra)
        self.actionMostrarBarra.setText(QtWidgets.QApplication.translate("SteamTables", "Mostrar barra de herramientas", None))
        self.actionDibujarSaturacion = QtWidgets.QAction(self)
        self.actionDibujarSaturacion.setCheckable(True)
        self.actionDibujarSaturacion.setChecked(True)
        self.actionDibujarSaturacion.setShortcut("Ctrl+S")
        self.actionDibujarSaturacion.triggered.connect(self.mostrarSaturacion)
        self.actionDibujarSaturacion.setText(QtWidgets.QApplication.translate("SteamTables", "Dibujar línea de saturación", None))
        self.actionMostrarPuntos = QtWidgets.QAction(self)
        self.actionMostrarPuntos.setCheckable(True)
        self.actionMostrarPuntos.setChecked(True)
        self.actionMostrarPuntos.setShortcut("Ctrl+P")
        self.actionMostrarPuntos.triggered.connect(self.mostrarPuntos)
        self.actionMostrarPuntos.setText(QtWidgets.QApplication.translate("SteamTables", "Mostrar puntos definidos por el usuario", None))
        self.actionAcerca_de = QtWidgets.QAction(self)
        self.actionAcerca_de.setText(QtWidgets.QApplication.translate("SteamTables", "Acerca de freesteam", None))
        self.actionAcerca_de.triggered.connect(self.acerca)
        self.actionAcerca_deQt= QtWidgets.QAction(self)
        self.actionAcerca_deQt.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/AboutQt.png")))
        self.actionAcerca_deQt.setText(QtWidgets.QApplication.translate("SteamTables", "Acerca de Qt", None))
        self.actionAcerca_deQt.triggered.connect(self.acercaQt)
        self.menuArchivo.addAction(self.Preferencias)
        self.menuArchivo.addAction(self.actionCSV)
        self.menuArchivo.addSeparator()
        self.menuArchivo.addAction(self.actionSalir)
        self.menuGrafico.addAction(self.action2D)
        self.menuGrafico.addAction(self.action3D)
        self.menuGrafico.addSeparator()
        self.menuGrafico.addAction(self.actionMostrarBarra)
        self.menuGrafico.addAction(self.actionDibujarSaturacion)
        self.menuGrafico.addAction(self.actionMostrarPuntos)
        self.menuAyuda.addAction(self.actionAcerca_de)
        self.menuAyuda.addAction(self.actionAcerca_deQt)
        self.menubar.addAction(self.menuArchivo.menuAction())
        self.menubar.addAction(self.menuGrafico.menuAction())
        self.menubar.addAction(self.menuAyuda.menuAction())
                
        #Ventana principal
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.toolBox = QtWidgets.QToolBox()
        self.page_1 = QtWidgets.QWidget()
        self.page_1.setGeometry(QtCore.QRect(0,0,302,390))
        self.gridLayout_2 = QtWidgets.QGridLayout(self.page_1)
        self.tabla = QtWidgets.QTableWidget(self.page_1)
        self.gridLayout_2.addWidget(self.tabla,0, 0, 1, 1)
        self.toolBox.addItem(self.page_1,"")
        self.page_Plot = QtWidgets.QWidget()
        self.page_Plot.setGeometry(QtCore.QRect(0,0,274,406))
        self.gridLayout_3 = QtWidgets.QGridLayout(self.page_Plot)
        self.checkIsoTherm=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoTherm.setText(QtWidgets.QApplication.translate("SteamTables", "Isotérmicas", None))
        self.checkIsoTherm.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas isotérmicas", None))
        self.gridLayout_3.addWidget(self.checkIsoTherm,0,0,1,1)
        self.checkIsoBar=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoBar.setText(QtWidgets.QApplication.translate("SteamTables", "Isobáricas", None))
        self.checkIsoBar.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas isobáricas", None))
        self.gridLayout_3.addWidget(self.checkIsoBar,0,1,1,1)
        self.checkIsoEnth=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoEnth.setText(QtWidgets.QApplication.translate("SteamTables", "Isoentálpicas", None))
        self.checkIsoEnth.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas isoentálpicas", None))
        self.gridLayout_3.addWidget(self.checkIsoEnth,0,2,1,1)
        self.checkIsoEntr=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoEntr.setText(QtWidgets.QApplication.translate("SteamTables", "Isoentrópicas", None))
        self.checkIsoEntr.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas isoentrópicas", None))
        self.gridLayout_3.addWidget(self.checkIsoEntr,0,3,1,1)
        self.checkIsoVol=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoVol.setText(QtWidgets.QApplication.translate("SteamTables", "Isocóricas", None))
        self.checkIsoVol.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas isocóricas", None))
        self.gridLayout_3.addWidget(self.checkIsoVol,0,4,1,1)
        self.checkIsoX=QtWidgets.QCheckBox(self.page_Plot)
        self.checkIsoX.setText(QtWidgets.QApplication.translate("SteamTables", "Isocalidad", None))
        self.checkIsoX.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Dibujar curvas con igual fracción de vapor", None))
        self.gridLayout_3.addWidget(self.checkIsoX,0,5,1,1)
        self.diagrama2D = Plot(self.page_Plot, dpi=self.Config.getfloat("General","Dpi"), dim=2)
        self.gridLayout_3.addWidget(self.diagrama2D,1,0,1,6)
        self.diagrama3D = Plot(self.page_Plot, dpi=self.Config.getfloat("General","Dpi"), dim=3)
        self.diagrama3D.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Pinchar y arrastrar para mover la orientación del gráfico", None))
        self.gridLayout_3.addWidget(self.diagrama3D,1,0,1,6)
        self.toolbar2D=NavigationToolbar2QT(self.diagrama2D, self.diagrama2D)
        self.gridLayout_3.addWidget(self.toolbar2D,3,0,1,6)
        self.toolbar3D=NavigationToolbar2QT(self.diagrama3D, self.diagrama3D)
        self.gridLayout_3.addWidget(self.toolbar3D,3,0,1,6)
        self.toolBox.addItem(self.page_Plot,"")
        self.gridLayout.addWidget(self.toolBox,0,0,1,1)
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_1), QtWidgets.QApplication.translate("SteamTables", "Tablas", None))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_Plot), QtWidgets.QApplication.translate("SteamTables", "Gráfico", None))

        #Toolbox parámetros Tabla
        self.setTabPosition(QtCore.Qt.DockWidgetArea(1), 0)
        self.dockWidget_Tabla = QtWidgets.QDockWidget()
        self.dockWidget_Tabla.setWindowTitle(QtWidgets.QApplication.translate("SteamTables", "Tabla", None))
        self.dockWidgetContents = QtWidgets.QWidget()
        self.gridLayout_4 = QtWidgets.QGridLayout(self.dockWidgetContents)
        self.label = QtWidgets.QLabel()
        self.label.setText(QtWidgets.QApplication.translate("SteamTables", "Ejes", None))
        self.gridLayout_4.addWidget(self.label,0,0,1,1)
        self.ejesTabla = QtWidgets.QComboBox()
        self.ejesTabla.setFixedSize(QtCore.QSize(100,20))
        self.ejesTabla.setToolTip(QtWidgets.QApplication.translate("SteamTables", "p\tPresión\nT\tTemperatura\nh\tEntalpía\ns\tEntropía\nv\tVolumen específico\nx\tCalidad (cuando es vapor saturado)", None))
        self.ejesTabla.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Definir variables impuestas", None))
        self.gridLayout_4.addWidget(self.ejesTabla,0,1,1,2)
        self.label_14 = QtWidgets.QLabel()
        self.label_14.setText(QtWidgets.QApplication.translate("SteamTables", "Calcular", None))
        self.gridLayout_4.addWidget(self.label_14,1,0,1,1)
        self.variableTabla = QtWidgets.QComboBox()
        self.variableTabla.setFixedWidth(150)
        self.variableTabla.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Definir variables a calcular", None))
        self.gridLayout_4.addWidget(self.variableTabla,1,1,1,2)
        self.gridLayout_4.addItem(QtWidgets.QSpacerItem(20,20,QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum),2,0,4,1)
        self.label_26 = QtWidgets.QLabel()
        self.gridLayout_4.addWidget(self.label_26,3,1,1,1)
        self.label_26.setAlignment(QtCore.Qt.AlignCenter)
        self.label_27 = QtWidgets.QLabel()
        self.gridLayout_4.addWidget(self.label_27,3,2,1,1)
        self.label_27.setAlignment(QtCore.Qt.AlignCenter)
        self.label_17 = QtWidgets.QLabel()
        self.label_17.setText(QtWidgets.QApplication.translate("SteamTables", "Intervalo", None))
        self.gridLayout_4.addWidget(self.label_17,6,0,1,1)
        self.abscisaInicio = QtWidgets.QLineEdit()
        self.abscisaInicio.setMaximumSize(QtCore.QSize(80,16777215))
        self.abscisaInicio.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.abscisaInicio,4,1,1,1)
        self.ordenadaInicio = QtWidgets.QLineEdit()
        self.ordenadaInicio.setMaximumSize(QtCore.QSize(80,16777215))
        self.ordenadaInicio.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.ordenadaInicio,4,2,1,1)
        self.label_15 = QtWidgets.QLabel()
        self.label_15.setText(QtWidgets.QApplication.translate("SteamTables", "Inicio", None))
        self.gridLayout_4.addWidget(self.label_15,4,0,1,1)
        self.label_16 = QtWidgets.QLabel()
        self.label_16.setText(QtWidgets.QApplication.translate("SteamTables", "Fin", None))
        self.gridLayout_4.addWidget(self.label_16,5,0,1,1)
        self.abscisaFin = QtWidgets.QLineEdit()
        self.abscisaFin.setMaximumSize(QtCore.QSize(80,16777215))
        self.abscisaFin.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.abscisaFin,5,1,1,1)
        self.ordenadaFin = QtWidgets.QLineEdit()
        self.ordenadaFin.setMaximumSize(QtCore.QSize(80,16777215))
        self.ordenadaFin.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.ordenadaFin,5,2,1,1)
        self.abscisaIntervalo = QtWidgets.QLineEdit()
        self.abscisaIntervalo.setMaximumSize(QtCore.QSize(80,16777215))
        self.abscisaIntervalo.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.abscisaIntervalo,6,1,1,1)
        self.ordenadaIntervalo = QtWidgets.QLineEdit()
        self.ordenadaIntervalo.setMaximumSize(QtCore.QSize(80,16777215))
        self.ordenadaIntervalo.setAlignment(QtCore.Qt.AlignRight)        
        self.gridLayout_4.addWidget(self.ordenadaIntervalo,6,2,1,1)
        self.botonCalcular = QtWidgets.QPushButton()
        self.botonCalcular.setMaximumWidth(80)
        self.botonCalcular.setText(QtWidgets.QApplication.translate("SteamTables", "Calcular", None))
        self.botonCalcular.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/calculate.png")))
        self.botonCalcular.clicked.connect(self.botonCalcular_clicked)
        self.gridLayout_4.addWidget(self.botonCalcular,7,2,1,1)
        self.dockWidget_Tabla.setWidget(self.dockWidgetContents)
        self.dockWidget_Tabla.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dockWidget_Tabla)
        self.dockWidget_Tabla.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)


        #Toolbox graficos 2D
        self.dockWidget_2D = QtWidgets.QDockWidget()
        self.dockWidget_2D.setWindowTitle(QtWidgets.QApplication.translate("SteamTables", "Gráfico 2D", None))
        self.Dock_2D= QtWidgets.QWidget()
        self.gridLayout_13 = QtWidgets.QGridLayout(self.Dock_2D)
        self.label_28=QtWidgets.QLabel()
        self.label_28.setText(QtWidgets.QApplication.translate("SteamTables", "Variable", None))
        self.gridLayout_13.addWidget(self.label_28,0,0,1,1)
        self.ejeX = QtWidgets.QComboBox()
        self.ejeX.setFixedWidth(50)
        self.gridLayout_13.addWidget(self.ejeX,0,1,1,1)
        self.label_36=QtWidgets.QLabel()
        self.label_36.setText(QtWidgets.QApplication.translate("SteamTables", "Mínimo", None))
        self.gridLayout_13.addWidget(self.label_36,1,0,1,1)
        self.ejeX_min = Entrada_con_unidades(float)
        self.ejeX_min.entrada.editingFinished.connect(self.diagrama2D_ejeX)
        self.gridLayout_13.addWidget(self.ejeX_min,1,1,1,1)
        self.label_37=QtWidgets.QLabel()
        self.label_37.setText(QtWidgets.QApplication.translate("SteamTables", "Máximo", None))
        self.gridLayout_13.addWidget(self.label_37,1,2,1,1)
        self.ejeX_max = Entrada_con_unidades(float)
        self.ejeX_max.entrada.editingFinished.connect(self.diagrama2D_ejeX)
        self.gridLayout_13.addWidget(self.ejeX_max,1,3,1,1)
        self.ejeX_escala=QtWidgets.QCheckBox()
        self.ejeX_escala.setText(QtWidgets.QApplication.translate("SteamTables", "Escala logarítmica", None))
        self.ejeX_escala.toggled.connect(self.ejeX_log)
        self.gridLayout_13.addWidget(self.ejeX_escala,2,0,1,4)
        self.label_29=QtWidgets.QLabel()
        self.label_29.setText(QtWidgets.QApplication.translate("SteamTables", "Variable", None))
        self.gridLayout_13.addWidget(self.label_29,3,0,1,1)
        self.ejeY = QtWidgets.QComboBox()
        self.ejeY.setFixedWidth(50)
        self.gridLayout_13.addWidget(self.ejeY,3,1,1,1)
        self.label_38=QtWidgets.QLabel()
        self.label_38.setText(QtWidgets.QApplication.translate("SteamTables", "Mínimo", None))
        self.gridLayout_13.addWidget(self.label_38,4,0,1,1)
        self.ejeY_min = Entrada_con_unidades(float)
        self.ejeY_min.entrada.editingFinished.connect(self.diagrama2D_ejeY)
        self.gridLayout_13.addWidget(self.ejeY_min,4,1,1,1)
        self.label_39=QtWidgets.QLabel()
        self.label_39.setText(QtWidgets.QApplication.translate("SteamTables", "Máximo", None))
        self.gridLayout_13.addWidget(self.label_39,4,2,1,1)
        self.ejeY_max = Entrada_con_unidades(float)
        self.ejeY_max.entrada.editingFinished.connect(self.diagrama2D_ejeY)
        self.gridLayout_13.addWidget(self.ejeY_max,4,3,1,1)
        self.ejeY_escala=QtWidgets.QCheckBox()
        self.ejeY_escala.setText(QtWidgets.QApplication.translate("SteamTables", "Escala logarítmica", None))
        self.ejeY_escala.toggled.connect(self.ejeY_log)
        self.gridLayout_13.addWidget(self.ejeY_escala,5,0,1,4)
        self.rejilla=QtWidgets.QCheckBox()
        self.rejilla.setText(QtWidgets.QApplication.translate("SteamTables", "Rejilla", None))
        self.rejilla.toggled.connect(self.rejilla_toggled)
        self.gridLayout_13.addWidget(self.rejilla,6,0,1,2)
        self.botonCalcular2 = QtWidgets.QPushButton()
        self.botonCalcular2.setMaximumWidth(80)
        self.botonCalcular2.setText(QtWidgets.QApplication.translate("SteamTables", "Calcular", None))
        self.botonCalcular2.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/button/calculate.png")))
        self.botonCalcular2.clicked.connect(self.botonCalcular_clicked)
        self.gridLayout_13.addWidget(self.botonCalcular2,6,3,1,1)
        self.dockWidget_2D.setWidget(self.Dock_2D)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dockWidget_2D)
        self.dockWidget_2D.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
        self.dockWidget_2D.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)

        
        #Toolbox parámetros Puntos individuales
        self.dockWidget_Puntos = QtWidgets.QDockWidget()
        self.dockWidget_Puntos.setWindowTitle(QtWidgets.QApplication.translate("SteamTables", "Puntos específicos", None))
        self.dockWidgetContents_2 = QtWidgets.QWidget()
        self.gridLayout_1 = QtWidgets.QGridLayout(self.dockWidgetContents_2)
        self.label_1 = QtWidgets.QLabel()
        self.label_1.setText(QtWidgets.QApplication.translate("SteamTables", "Datos conocidos", None))
        self.gridLayout_1.addWidget(self.label_1,0,0,1,1)
        self.variablesCalculo = QtWidgets.QComboBox()
        self.variablesCalculo.setFixedWidth(100)
        self.variablesCalculo.setToolTip(QtWidgets.QApplication.translate("SteamTables", "p\tPresión\nT\tTemperatura\nh\tEntalpía\ns\tEntropía\nv\tVolumen específico\nx\tCalidad (cuando es vapor saturado)", None))
        self.variablesCalculo.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Definir variables impuestas", None))
        self.gridLayout_1.addWidget(self.variablesCalculo,0,1,1,2)
        self.gridLayout_1.addItem(QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed), 1, 0, 1, 3)
        self.label_2 = QtWidgets.QLabel()
        self.label_2.setText(QtWidgets.QApplication.translate("SteamTables", "Presión", None))
        self.gridLayout_1.addWidget(self.label_2,6,0,1,1)
        self.presion=Entrada_con_unidades(unidades.Pressure, color=self.Config.get("General","Color"))
        self.presion.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Presión", None))
        self.presion.entrada.editingFinished.connect(self.presion_editingFinished)
        self.gridLayout_1.addWidget(self.presion,6,1,1,2)
        self.label_3 = QtWidgets.QLabel()
        self.label_3.setText(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None))
        self.gridLayout_1.addWidget(self.label_3,7,0,1,1)
        self.temperatura=Entrada_con_unidades(unidades.Temperature, color=self.Config.get("General","Color"))
        self.temperatura.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None))
        self.temperatura.entrada.editingFinished.connect(self.temperatura_editingFinished)
        self.gridLayout_1.addWidget(self.temperatura,7,1,1,2)
        self.label_8 = QtWidgets.QLabel()
        self.label_8.setText(QtWidgets.QApplication.translate("SteamTables", "Volumen", None))
        self.gridLayout_1.addWidget(self.label_8,8,0,1,1)
        self.volumen=Entrada_con_unidades(unidades.SpecificVolume, color=self.Config.get("General","Color"))
        self.volumen.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None))
        self.volumen.entrada.editingFinished.connect(self.volumen_editingFinished)
        self.gridLayout_1.addWidget(self.volumen,8,1,1,2)
        self.label_4 = QtWidgets.QLabel()
        self.label_4.setText(QtWidgets.QApplication.translate("SteamTables", "Entalpia", None))
        self.gridLayout_1.addWidget(self.label_4,9,0,1,1)
        self.entalpia=Entrada_con_unidades(unidades.Enthalpy, color=self.Config.get("General","Color"))
        self.entalpia.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Entalpía", None))
        self.entalpia.entrada.editingFinished.connect(self.entalpia_editingFinished)
        self.gridLayout_1.addWidget(self.entalpia,9,1,1,2)
        self.label_5 = QtWidgets.QLabel()
        self.label_5.setText(QtWidgets.QApplication.translate("SteamTables", "Entropía", None))
        self.gridLayout_1.addWidget(self.label_5,10,0,1,1)
        self.entropia=Entrada_con_unidades(unidades.SpecificHeat, "Entropy", color=self.Config.get("General","Color"))
        self.entropia.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Entropía", None))
        self.entropia.entrada.editingFinished.connect(self.entropia_editingFinished)
        self.gridLayout_1.addWidget(self.entropia,10,1,1,2)
        self.label_12 = QtWidgets.QLabel()
        self.label_12.setText(QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None))
        self.gridLayout_1.addWidget(self.label_12,11,0,1,1)
        self.fraccionVapor = Entrada_con_unidades(float, color=self.Config.get("General","Color"))
        self.fraccionVapor.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None))
        self.fraccionVapor.entrada.editingFinished.connect(self.fraccionVapor_editingFinished)
        self.gridLayout_1.addWidget(self.fraccionVapor,11,1,1,1)
        self.label_13 = QtWidgets.QLabel()
        self.label_13.setText(QtWidgets.QApplication.translate("SteamTables", "Energía interna", None))
        self.gridLayout_1.addWidget(self.label_13,12,0,1,1)
        self.energiaInterna=Entrada_con_unidades(unidades.Enthalpy, readOnly=True, retornar=False)
        self.energiaInterna.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Energía Interna", None))
        self.gridLayout_1.addWidget(self.energiaInterna,12,1,1,1)
        self.label_18 = QtWidgets.QLabel()
        self.label_18.setText(QtWidgets.QApplication.translate("SteamTables", "E. Gibbs", None))
        self.gridLayout_1.addWidget(self.label_18,13,0,1,1)
        self.energiaGibbs=Entrada_con_unidades(unidades.Enthalpy, readOnly=True, retornar=False)
        self.energiaGibbs.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None))
        self.gridLayout_1.addWidget(self.energiaGibbs,13,1,1,1)
        self.label_19 = QtWidgets.QLabel()
        self.label_19.setText(QtWidgets.QApplication.translate("SteamTables", "E. Helmholtz", None))
        self.gridLayout_1.addWidget(self.label_19,14,0,1,1)
        self.energiaHelmholtz=Entrada_con_unidades(unidades.Enthalpy, readOnly=True, retornar=False)
        self.energiaHelmholtz.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None))
        self.gridLayout_1.addWidget(self.energiaHelmholtz,14,1,1,1)
        self.label_20 = QtWidgets.QLabel()
        self.label_20.setText(QtWidgets.QApplication.translate("SteamTables", "Densidad", None))
        self.gridLayout_1.addWidget(self.label_20,15,0,1,1)
        self.densidad=Entrada_con_unidades(unidades.Density, readOnly=True, retornar=False)
        self.densidad.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Densidad", None))
        self.gridLayout_1.addWidget(self.densidad,15,1,1,1)
        self.label_21 = QtWidgets.QLabel()
        self.label_21.setText(QtWidgets.QApplication.translate("SteamTables", "Cp", None))
        self.gridLayout_1.addWidget(self.label_21,16,0,1,1)
        self.cp=Entrada_con_unidades(unidades.SpecificHeat, readOnly=True, retornar=False)
        self.cp.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Capacidad calorífica a presión constante", None))
        self.gridLayout_1.addWidget(self.cp,16,1,1,1)
        self.label_22 = QtWidgets.QLabel()
        self.label_22.setText(QtWidgets.QApplication.translate("SteamTables", "Cv", None))
        self.gridLayout_1.addWidget(self.label_22,17,0,1,1)
        self.cv=Entrada_con_unidades(unidades.SpecificHeat, readOnly=True, retornar=False)
        self.cv.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Capacidad calorífica a volumen constante", None))
        self.gridLayout_1.addWidget(self.cv,17,1,1,1)
        self.label_23 = QtWidgets.QLabel()
        self.label_23.setText(QtWidgets.QApplication.translate("SteamTables", "Conductividad", None))
        self.gridLayout_1.addWidget(self.label_23,18,0,1,1)
        self.conductividad=Entrada_con_unidades(unidades.ThermalConductivity, readOnly=True, retornar=False)
        self.conductividad.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None))
        self.gridLayout_1.addWidget(self.conductividad,18,1,1,1)
        self.label_24 = QtWidgets.QLabel()
        self.label_24.setText(QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None))
        self.gridLayout_1.addWidget(self.label_24,19,0,1,1)
        self.viscosidad=Entrada_con_unidades(unidades.Viscosity, readOnly=True, retornar=False)
        self.viscosidad.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Viscosidad dinámica", None))
        self.gridLayout_1.addWidget(self.viscosidad, 19,1,1,1)
        self.label_25 = QtWidgets.QLabel()
        self.label_25.setText(QtWidgets.QApplication.translate("SteamTables", "V sonido", None))
        self.gridLayout_1.addWidget(self.label_25,20,0,1,1)
        self.velocidad=Entrada_con_unidades(unidades.Speed, readOnly=True, retornar=False)
        self.velocidad.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None))
        self.gridLayout_1.addWidget(self.velocidad,20,1,1,1)
        self.label_30 = QtWidgets.QLabel()
        self.label_30.setText(QtWidgets.QApplication.translate("SteamTables", "Región", None))
        self.gridLayout_1.addWidget(self.label_30,21,0,1,1)
        self.region = QtWidgets.QLineEdit()
        self.region.setFixedSize(QtCore.QSize(85,24))
        self.region.setReadOnly(True)
        self.region.setAlignment(QtCore.Qt.AlignRight)
        self.gridLayout_1.addWidget(self.region,21,1,1,1)
        self.botonAdd=QtWidgets.QPushButton()
        self.botonAdd.setEnabled(False)
        self.botonAdd.setText(QtWidgets.QApplication.translate("SteamTables", "Añadir", None))
        self.botonAdd.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Añadir a la lista de puntos representados en la gráfica", None))
        self.botonAdd.clicked.connect(self.botonAdd_clicked)
        self.gridLayout_1.addWidget(self.botonAdd,22,0,1,1)
        self.botonLista=QtWidgets.QPushButton()
        self.botonLista.setMaximumWidth(80)
        self.botonLista.setText(QtWidgets.QApplication.translate("SteamTables", "Lista", None))
        self.botonLista.setStatusTip(QtWidgets.QApplication.translate("SteamTables", "Mostrar lista de puntos específicos representados en la gráfica", None))
        self.botonLista.clicked.connect(self.botonLista_clicked)
        self.gridLayout_1.addWidget(self.botonLista,22,1,1,2)
        spacerItem = QtWidgets.QSpacerItem(72, 34, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_1.addItem(spacerItem, 23, 1, 1, 3)
        self.mostrarUnidades()
        self.dockWidget_Puntos.setWidget(self.dockWidgetContents_2)
        self.addDockWidget(QtCore.Qt.DockWidgetArea(1), self.dockWidget_Puntos)
        self.dockWidget_Puntos.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)
        self.dockWidget_Puntos.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
        
        #Iniciar valores de widgets
        Ejes=["p,T", "p,h", "p,s", "p,v", "T,s", "T,x"]
        for i in Ejes:
            self.ejesTabla.addItem(i)
            self.variablesCalculo.addItem(i)
        Ejes2D=["p", "T", "s", "h", "u", "v"]
        for i in Ejes2D:
            self.ejeX.addItem(i)
            self.ejeY.addItem(i)
            
        self.matriz=[]
        self.saturacion=[]
        self.isoterma=[]
        self.isobara=[]
        self.isoentropica=[]
        self.isoentalpica=[]
        self.isocora=[]
        self.isoX=[]
        self.factorx, self.factory, self.factorz, self.factorx2, self.factory2=(0, 0, 0, 0, 0)
        
        #Cargar configuración
        if self.Config.has_section("Table"):
            self.ejesTabla.setCurrentIndex(self.Config.getint("Table","Axis"))
            self.ejesTabla_currentIndexChanged(self.Config.getint("Table","Axis"))
            self.variableTabla.setCurrentIndex(self.Config.getint("Table","Calculate"))
            self.abscisaInicio.setText(self.Config.get("Table","x_start"))
            self.abscisaFin.setText(self.Config.get("Table","x_end"))
            self.abscisaIntervalo.setText(self.Config.get("Table","x_step"))
            self.ordenadaInicio.setText(self.Config.get("Table","y_start"))
            self.ordenadaFin.setText(self.Config.get("Table","y_end"))
            self.ordenadaIntervalo.setText(self.Config.get("Table","y_step"))
        if self.Config.has_section("General"):
            self.actionDibujarSaturacion.setChecked(self.Config.getboolean("General","Sat"))
            self.actionMostrarPuntos.setChecked(self.Config.getboolean("General","Points"))
            self.actionMostrarBarra.setChecked(self.Config.getboolean("General","Toolbar"))
            self.toolbar2D.setVisible(self.actionMostrarBarra.isChecked())
            self.toolbar3D.setVisible(self.actionMostrarBarra.isChecked())
            self.action2D.setChecked(self.Config.getboolean("General","2D"))
            self.d2(self.action2D.isChecked())
            if self.Config.getboolean("General","Plot"):
                self.toolBox.setCurrentIndex(1)
            self.checkIsoTherm.setChecked(self.Config.getboolean("General","Isotherm"))
            self.checkIsoBar.setChecked(self.Config.getboolean("General","Isobar"))
            self.checkIsoEnth.setChecked(self.Config.getboolean("General","Isoenthalpic"))
            self.checkIsoEntr.setChecked(self.Config.getboolean("General","Isoentropic"))
            self.checkIsoVol.setChecked(self.Config.getboolean("General","Isochor"))
            self.checkIsoX.setChecked(self.Config.getboolean("General","Isoquality"))
        if self.Config.has_section("2D"):
            self.ejeX.setCurrentIndex(self.ejeX.findText(self.Config.get("2D", "Xvariable")))
            self.rellenar_ejeY(self.ejeX.currentIndex())
            self.ejeX_escala.setChecked(self.Config.getboolean("2D", "XScale"))
            self.ejeX_min.setValue(self.Config.getfloat("2D", "XMin"))
            self.ejeX_max.setValue(self.Config.getfloat("2D", "XMax"))
            self.ejeY.setCurrentIndex(self.ejeY.findText(self.Config.get("2D", "Yvariable")))
            self.ejeY_escala.setChecked(self.Config.getboolean("2D", "YScale"))
            self.ejeY_min.setValue(self.Config.getfloat("2D", "YMin"))
            self.ejeY_max.setValue(self.Config.getfloat("2D", "YMax"))
            self.rejilla.setChecked(self.Config.getboolean("2D", "Grid"))
        if self.Config.has_section("Points"):
            self.variablesCalculo.setCurrentIndex(self.Config.getint("Points", 'Variable'))
            self.variablesCalculo_currentIndexChanged(self.Config.getint("Points", 'Variable'))
            i=0
            self.puntos=[]
            while True:
                try:
                    punto=list(map(float, self.Config.get("Points", str(i)).split(",")))
                except: break
                if punto[0]==0.0:
                    self.puntos.append(steam_Tx(punto[1], punto[2]))
                else:
                    self.puntos.append(steam_pT(punto[2], punto[1]))
                i+=1
            self.punto=self.puntos[-1]
            self.botonAdd.setEnabled(True)
            self.mostrarPropiedades()
        
        self.Isoterma=Iso()
        self.Isobara=Iso()
        self.Isoentropica=Iso()
        self.Isoentalpica=Iso()
        self.Isocora=Iso()
        self.IsoX=Iso()
        self.definirIsolineas()
    
        self.ejesTabla.currentIndexChanged.connect(self.ejesTabla_currentIndexChanged)
        self.variablesCalculo.currentIndexChanged.connect(self.variablesCalculo_currentIndexChanged)
        self.variableTabla.currentIndexChanged.connect(self.Calcular_Propiedades)
        self.ejeX.currentIndexChanged.connect(self.rellenar_ejeY)
        self.checkIsoTherm.toggled.connect(self.mostrarIsoterma)
        self.checkIsoBar.toggled.connect(self.mostrarIsobara)
        self.checkIsoEnth.toggled.connect(self.mostrarIsoentalpica)
        self.checkIsoEntr.toggled.connect(self.mostrarIsoentropica)
        self.checkIsoVol.toggled.connect(self.mostrarIsocora)
        self.checkIsoX.toggled.connect(self.mostrarIsoX)


    def acerca(self):
        QtWidgets.QMessageBox.information(self,"Acerca de" ,"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd">\n<html><head><meta name="qrichtext" content="1" /><style type="text/css">\np, li { white-space: pre-wrap; }\n</style></head><body style=" font-family:'Nimbus Sans L'; font-size:9pt; font-weight:400; font-style:normal;">\n<table border="0" style=" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px;" cellspacing="2" cellpadding="0">\n<tr>\n<td style=" vertical-align:top;">\n<p style=" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;"><span style=" font-weight:600;">freesteam</span> is an open source implementation of international-standard IAPWS-IF97 steam tables from the <a href="http://www.iapws.org"><span style=" text-decoration: underline; color:#0000ff;">International Association for the Properties of Water and Steam</span></a> (IAPWS). <span style=" font-weight:600;">freesteam</span> lets you compute water and steam properties for a wide range of pressures and temperatures: you can specify the state of the steam in terms of a variety of combinations of 'known' properties, then freesteam will solve and allow you to query to find the values of the 'unknown' properties.</p>\n<p style=" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;">Website: <a href="http://freesteam.sourceforge.net/"><span style=" text-decoration: underline; color:#0000ff;">http://freesteam.sourceforge.net/</span></a></p></td></tr></table></body></html>""")
        
    def acercaQt(self):
        QtWidgets.QMessageBox.aboutQt(self,"Acerca de Qt" )

    def exporttoCSV(self):
        """Guarda los datos de la tabla en un archivo csv"""
        fname = QtWidgets.QFileDialog.getSaveFileName(self, QtWidgets.QApplication.translate("SteamTables", "Exportar datos", None), "./", "CSV (*.csv);;All archives (*.*)")[0]
        if fname:
            texto = writer(open(fname, 'wb'), delimiter='\t')
            texto.writerow([""]+[str(i) for i in self.xdata[0]])
            for i, fila in enumerate(self.zdata):
                texto.writerow([str(self.ydata[i][0])]+[str(i) for i in fila])

    def preferencias(self):
        """Muestra el diálogo de configuración de preferencias"""
        dialog=Preferences()
        if dialog.exec_():
            self.Config.read("UI_steamTablesrc")
            self.actualizarConfiguracion()
            self.definirIsolineas()
        
    def closeEvent(self, event):
        """Guarda la configuración antes de salir"""
        Config=ConfigParser()
        Config.read("UI_steamTablesrc")
        if not Config.has_section("Table"):
            Config.add_section("Table")
        Config.set("Table", "Axis", self.ejesTabla.currentIndex())
        Config.set("Table", "Calculate", self.variableTabla.currentIndex())
        Config.set("Table", "x_start", self.abscisaInicio.text())
        Config.set("Table", "x_end", self.abscisaFin.text())
        Config.set("Table", "x_step", self.abscisaIntervalo.text())
        Config.set("Table", "y_start", self.ordenadaInicio.text())
        Config.set("Table", "y_end", self.ordenadaFin.text())
        Config.set("Table", "y_step", self.ordenadaIntervalo.text())
        if not Config.has_section("General"):
            Config.add_section("General")
        Config.set("General", "Sat", self.actionDibujarSaturacion.isChecked())
        Config.set("General", "Points", self.actionMostrarPuntos.isChecked())
        Config.set("General", "Toolbar", self.actionMostrarBarra.isChecked())
        Config.set("General", "2D", self.action2D.isChecked())
        Config.set("General", "Plot", self.page_Plot.isVisible())
        Config.set("General", "Isotherm", self.checkIsoTherm.isChecked())
        Config.set("General", "Isobar", self.checkIsoBar.isChecked())
        Config.set("General", "Isoenthalpic", self.checkIsoEnth.isChecked())
        Config.set("General", "Isoentropic", self.checkIsoEntr.isChecked())
        Config.set("General", "Isochor", self.checkIsoVol.isChecked())
        Config.set("General", "Isoquality", self.checkIsoX.isChecked())
        if not Config.has_section("2D"):
            Config.add_section("2D")
        Config.set("2D", "Xvariable", self.ejeX.currentText())
        Config.set("2D", "XScale", self.ejeX_escala.isChecked())
        Config.set("2D", "XMin", self.ejeX_min.value)
        Config.set("2D", "XMax", self.ejeX_max.value)
        Config.set("2D", "Yvariable", self.ejeY.currentText())
        Config.set("2D", "YScale", self.ejeY_escala.isChecked())
        Config.set("2D", "YMin", self.ejeY_min.value)
        Config.set("2D", "YMax", self.ejeY_max.value)
        Config.set("2D", "Grid", self.rejilla.isChecked())
        if len(self.puntos)>0:
            Config.remove_section("Points")
            Config.add_section("Points")
            for i, punto in enumerate(self.puntos):
                if punto.region=="\x04":
                    Config.set("Points", str(i), "0"+","+str(punto.T)+","+str(punto.x))
                else:
                    Config.set("Points", str(i), "1"+","+str(punto.T)+","+str(punto.p))
            Config.set("Points", "Variable", self.variablesCalculo.currentIndex())
        Config.write(open("UI_steamTablesrc", "w"))
        event.accept()
                
                
    #Controles
    def ejesTabla_currentIndexChanged(self, indice):
        """Hace los cambios pertinentes en la gui cuando se cambian los ejes de la tabla 3D:
            Actualiza unidades mostradas en la entrada de datos de la tabla
            Actualiza los checkbox de isolineas habilitados (todos menos los de los ejes x e y)"""
        if indice==0:
            self.rellenar_variableTabla(0, 1)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Presión", None)+", %s" %config.Configuracion("Pressure").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None)+", %s" %config.Configuracion("Temperature").text())
            self.checkIsoBar.setEnabled(False)
            self.checkIsoTherm.setEnabled(False)
            self.checkIsoBar.setChecked(False)
            self.checkIsoTherm.setChecked(False)
            self.checkIsoEntr.setEnabled(True)
            self.checkIsoEnth.setEnabled(True)
            self.checkIsoVol.setEnabled(True)
        elif indice==1:
            self.rellenar_variableTabla(0, 3)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Presión", None)+", %s" %config.Configuracion("Pressure").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Entalpía", None)+", %s" %config.Configuracion("Enthalpy").text())
            self.checkIsoBar.setEnabled(False)
            self.checkIsoBar.setChecked(False)
            self.checkIsoTherm.setEnabled(True)
            self.checkIsoEntr.setEnabled(True)
            self.checkIsoEnth.setEnabled(False)
            self.checkIsoEnth.setChecked(False)
            self.checkIsoVol.setEnabled(True)
        elif indice==2:
            self.rellenar_variableTabla(0, 4)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Presión", None)+", %s" %config.Configuracion("Pressure").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Entropía", None)+", %s" %config.Configuracion("SpecificHeat", "Entropy").text())
            self.checkIsoBar.setEnabled(False)
            self.checkIsoBar.setChecked(False)
            self.checkIsoTherm.setEnabled(True)
            self.checkIsoEntr.setEnabled(False)
            self.checkIsoEntr.setChecked(False)
            self.checkIsoEnth.setEnabled(True)
            self.checkIsoVol.setEnabled(True)
        elif indice==3:
            self.rellenar_variableTabla(0, 2)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Presión", None)+", %s" %config.Configuracion("Pressure").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Volumen", None)+", %s" %config.Configuracion("SpecificVolume").text())
            self.checkIsoBar.setEnabled(False)
            self.checkIsoBar.setChecked(False)
            self.checkIsoTherm.setEnabled(True)
            self.checkIsoEntr.setEnabled(True)
            self.checkIsoEnth.setEnabled(True)
            self.checkIsoVol.setEnabled(False)
            self.checkIsoVol.setChecked(False)
        elif indice==4:
            self.rellenar_variableTabla(1, 4)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None)+", %s" %config.Configuracion("Temperature").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Entropía", None)+", %s" %config.Configuracion("SpecificHeat", "Entropy").text())
            self.checkIsoBar.setEnabled(True)
            self.checkIsoTherm.setEnabled(False)
            self.checkIsoTherm.setChecked(False)
            self.checkIsoEntr.setEnabled(False)
            self.checkIsoEntr.setChecked(False)
            self.checkIsoEnth.setEnabled(True)
            self.checkIsoVol.setEnabled(True)
        elif indice==5:
            self.rellenar_variableTabla(1, 13)
            self.label_26.setText(QtWidgets.QApplication.translate("SteamTables", "Temperatura", None)+", %s" %config.Configuracion("Temperature").text())
            self.label_27.setText(QtWidgets.QApplication.translate("SteamTables", "Calidad", None))
            self.checkIsoBar.setEnabled(True)
            self.checkIsoTherm.setEnabled(False)
            self.checkIsoTherm.setChecked(False)
            self.checkIsoEntr.setEnabled(True)
            self.checkIsoEnth.setEnabled(True)
            self.checkIsoVol.setEnabled(True)

    def variablesCalculo_currentIndexChanged(self, indice):
        """Hace los cambios pertinentes en la gui cuando se cambian las variables impuestas de los puntos especificados:
            Resalta las variables a introducir
            Fija las el resto de variables para evitar su edición por el usuario"""
        if indice==0:
            self.presion.setReadOnly(False)
            self.presion.setResaltado(True)
            self.temperatura.setReadOnly(False)
            self.temperatura.setResaltado(True)
            self.entalpia.setReadOnly(True)
            self.entalpia.setResaltado(False)
            self.entropia.setReadOnly(True)
            self.entropia.setResaltado(False)
            self.volumen.setReadOnly(True)
            self.volumen.setResaltado(False)
            self.fraccionVapor.setReadOnly(True)
            self.fraccionVapor.setResaltado(False)
            self.presion.setRetornar(True)
            self.temperatura.setRetornar(True)
            self.entalpia.setRetornar(False)
            self.entropia.setRetornar(False)
            self.volumen.setRetornar(False)
        elif indice==1:
            self.entalpia.setReadOnly(False)
            self.entalpia.setResaltado(True)
            self.presion.setReadOnly(False)
            self.presion.setResaltado(True)
            self.temperatura.setReadOnly(True)
            self.temperatura.setResaltado(False)
            self.entropia.setReadOnly(True)
            self.entropia.setResaltado(False)
            self.volumen.setReadOnly(True)
            self.volumen.setResaltado(False)
            self.fraccionVapor.setReadOnly(True)
            self.fraccionVapor.setResaltado(False)
            self.presion.setRetornar(True)
            self.temperatura.setRetornar(False)
            self.entalpia.setRetornar(True)
            self.entropia.setRetornar(False)
            self.volumen.setRetornar(False)
        elif indice==2:
            self.entropia.setReadOnly(False)
            self.entropia.setResaltado(True)
            self.presion.setReadOnly(False)
            self.presion.setResaltado(True)
            self.temperatura.setReadOnly(True)
            self.temperatura.setResaltado(False)
            self.entalpia.setReadOnly(True)
            self.entalpia.setResaltado(False)
            self.volumen.setReadOnly(True)
            self.volumen.setResaltado(False)
            self.fraccionVapor.setReadOnly(True)
            self.fraccionVapor.setResaltado(False)
            self.presion.setRetornar(True)
            self.temperatura.setRetornar(False)
            self.entalpia.setRetornar(False)
            self.entropia.setRetornar(True)
            self.volumen.setRetornar(False)
        elif indice==3:
            self.volumen.setReadOnly(False)
            self.volumen.setResaltado(True)
            self.presion.setReadOnly(False)
            self.presion.setResaltado(True)
            self.temperatura.setReadOnly(True)
            self.temperatura.setResaltado(False)
            self.entropia.setReadOnly(True)
            self.entropia.setResaltado(False)
            self.entalpia.setReadOnly(True)
            self.entalpia.setResaltado(False)
            self.fraccionVapor.setReadOnly(True)
            self.fraccionVapor.setResaltado(False)
            self.presion.setRetornar(True)
            self.temperatura.setRetornar(False)
            self.entalpia.setRetornar(False)
            self.entropia.setRetornar(False)
            self.volumen.setRetornar(True)
        elif indice==4:
            self.entropia.setReadOnly(False)
            self.entropia.setResaltado(True)
            self.temperatura.setReadOnly(False)
            self.temperatura.setResaltado(True)
            self.entalpia.setReadOnly(True)
            self.entalpia.setResaltado(False)
            self.presion.setReadOnly(True)
            self.presion.setResaltado(False)
            self.volumen.setReadOnly(True)
            self.volumen.setResaltado(False)
            self.fraccionVapor.setReadOnly(True)
            self.fraccionVapor.setResaltado(False)
            self.presion.setRetornar(False)
            self.temperatura.setRetornar(True)
            self.entalpia.setRetornar(False)
            self.entropia.setRetornar(True)
            self.volumen.setRetornar(False)
        elif indice==5:
            self.fraccionVapor.setReadOnly(False)
            self.fraccionVapor.setResaltado(True)
            self.temperatura.setReadOnly(False)
            self.temperatura.setResaltado(True)
            self.entalpia.setReadOnly(True)
            self.entalpia.setResaltado(False)
            self.presion.setReadOnly(True)
            self.presion.setResaltado(False)
            self.volumen.setReadOnly(True)
            self.volumen.setResaltado(False)
            self.entropia.setReadOnly(True)
            self.entropia.setResaltado(False)
            self.presion.setRetornar(False)
            self.temperatura.setRetornar(True)
            self.entalpia.setRetornar(False)
            self.entropia.setRetornar(False)
            self.volumen.setRetornar(False)

    def presion_editingFinished(self):
        if self.variablesCalculo.currentIndex()==0:
            otro=self.temperatura
        elif self.variablesCalculo.currentIndex()==1:
            otro=self.entalpia
        elif self.variablesCalculo.currentIndex()==2:
            otro=self.entropia
        else:
            otro=self.volumen
        if otro.value and self.presion.value:
            self.calcularPropiedades()
        
    def temperatura_editingFinished(self):
        if self.variablesCalculo.currentIndex()==0:
            otro=self.presion
        elif self.variablesCalculo.currentIndex()==4:
            otro=self.entropia
        else:
            otro=self.fraccionVapor
        if otro.value and self.temperatura.value:
            self.calcularPropiedades()
        
    def entalpia_editingFinished(self):
        if self.presion.value and self.entalpia.value:
            self.calcularPropiedades()

    def entropia_editingFinished(self):
        if self.variablesCalculo.currentIndex()==2:
            otro=self.presion
        else:
            otro=self.temperatura
        if otro.value and self.entropia.value:
            self.calcularPropiedades()

    def volumen_editingFinished(self):
        if self.presion.value and self.volumen.value:
            self.calcularPropiedades()

    def fraccionVapor_editingFinished(self):
        if self.fraccionVapor.value>1.0 or self.fraccionVapor.value<0.0:
            self.fraccionVapor.clear()
            self.fraccionVapor.setFocus()
            valido=False
        else:
            valido=True
        if self.temperatura.value and valido:
            self.calcularPropiedades()

    def botonAdd_clicked(self):
        """Añade el punto especificado actual a la lista"""
        self.puntos.append(self.punto)
        self.calcularPuntos()
        self.dibujar()
    
    def botonLista_clicked(self):
        """Muestra la lista de puntos especificados por el usuario"""
        dialog=Ventana_Lista_Puntos(self.puntos)
        if dialog.exec_():
            self.puntos=dialog.puntos
            self.calcularPuntos()
            self.dibujar()

    def botonCalcular_clicked(self):
        """Método que calcula todos los datos de las tablas y gráficos"""
        self.progresbar.setVisible(True)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando matriz...", None))
        self.progresbar.setValue(0)
        QtWidgets.QApplication.processEvents()
        self.calcularMatriz(0, 30)
        self.Calcular_Propiedades(35, 5)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando línea de saturación...", None))
        QtWidgets.QApplication.processEvents()
        self.calcularSaturacion()
        if len(self.puntos)>0:
            self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando puntos personalizados...", None))
            QtWidgets.QApplication.processEvents()
            self.calcularPuntos()
        self.calcularIsoterma(40, 10)
        self.calcularIsobara(50, 10)
        self.calcularIsoentropica(60, 10)
        self.calcularIsoentalpica(70, 10)
        self.calcularIsocora(80, 10)
        self.calcularIsoX(90, 10)
        self.progresbar.setValue(100)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Dibujando...", None))
        QtWidgets.QApplication.processEvents()
        self.dibujar()
        self.progresbar.setVisible(False)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Listo", None))

    #Opciones gráficos
    def d2(self, activado):
        """Se ejecuta si se cambia el modo de gráfico a 2D"""
        self.dockWidget_Tabla.setVisible(not activado)
        self.dockWidget_2D.setVisible(activado)
        self.tabla.setEnabled(not activado)
        self.diagrama2D.setVisible(activado)
        self.toolbar2D.setVisible(self.actionMostrarBarra.isChecked() and activado)
        self.diagrama3D.setVisible(not activado)
        self.toolbar3D.setVisible(self.actionMostrarBarra.isChecked() and not activado)
        self.checkIsoBar.setEnabled(True)
        self.checkIsoTherm.setEnabled(True)
        self.checkIsoEntr.setEnabled(True)
        self.checkIsoEnth.setEnabled(True)
        self.checkIsoVol.setEnabled(True)

    def mostrarBarra(self, bool):
        """Muestra la toolbar de matplotlib"""
        self.toolbar2D.setVisible(bool and self.action2D.isChecked())
        self.toolbar3D.setVisible(bool and self.action3D.isChecked())

    def mostrarPuntos(self, bool):
        """Muestra los puntos específicos"""
        for i in self.diagrama3D.puntos:
            i[0].set_visible(bool)
        for i in self.diagrama2D.puntos:
            i[0].set_visible(bool)
        self.diagrama3D.draw()
        self.diagrama2D.draw()
        
    def mostrarSaturacion(self, bool):
        """Muestra la línea de saturación"""
        self.diagrama3D.satliq[0].set_visible(bool)
        self.diagrama3D.satgas[0].set_visible(bool)
        self.diagrama2D.satliq[0].set_visible(bool)
        self.diagrama2D.satgas[0].set_visible(bool)
        self.diagrama3D.draw()
        self.diagrama2D.draw()

    def mostrarIsoentropica(self):
        """Muestra las líneas isoentrópicas"""
        for i in self.diagrama3D.Isoentropica:
            i[0].set_visible(self.checkIsoEntr.isChecked() and self.checkIsoEntr.isEnabled())
        for i in self.diagrama2D.Isoentropica:
            i[0].set_visible(self.checkIsoEntr.isChecked())
        for i in self.diagrama2D.Isoentropica_label:
            i.set_visible(self.checkIsoEntr.isChecked() and self.Isoentropica.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()
        
    def mostrarIsoentalpica(self, bool):
        """Muestra las líneas isoentálpicas"""
        for i in self.diagrama3D.Isoentalpica:
            i[0].set_visible(bool and  self.checkIsoEnth.isEnabled())
        for i in self.diagrama2D.Isoentalpica:
            i[0].set_visible(bool)
        for i in self.diagrama2D.Isoentalpica_label:
            i.set_visible(self.checkIsoEnth.isChecked() and self.Isoentalpica.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()

    def mostrarIsobara(self):
        """Muestra las líneas isobáras"""
        for i in self.diagrama3D.Isobara:
            i[0].set_visible(self.checkIsoBar.isChecked() and self.checkIsoBar.isEnabled())
        for i in self.diagrama2D.Isobara:
            i[0].set_visible(self.checkIsoBar.isChecked())
        for i in self.diagrama2D.Isobara_label:
            i.set_visible(self.checkIsoBar.isChecked() and self.Isobara.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()
        
    def mostrarIsoterma(self):
        """Muestra las líneas isotermas"""
        for i in self.diagrama3D.Isoterma:
            i[0].set_visible(self.checkIsoTherm.isChecked() and self.checkIsoTherm.isEnabled())
        for i in self.diagrama2D.Isoterma:
            i[0].set_visible(self.checkIsoTherm.isChecked())
        for i in self.diagrama2D.Isoterma_label:
            i.set_visible(self.checkIsoTherm.isChecked() and self.Isoterma.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()
        
    def mostrarIsocora(self):
        """Muestra las líneas isocoras"""
        for i in self.diagrama3D.Isocora:
            i[0].set_visible(self.checkIsoVol.isChecked() and self.checkIsoVol.isEnabled())
        for i in self.diagrama2D.Isocora:
            i[0].set_visible(self.checkIsoVol.isChecked())
        for i in self.diagrama2D.Isocora_label:
            i.set_visible(self.checkIsoVol.isChecked() and self.Isocora.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()
        
    def mostrarIsoX(self):
        """Muestra las líneas con igual fración de vapor"""
        for i in self.diagrama3D.IsoX:
            i[0].set_visible(self.checkIsoX.isChecked())
        for i in self.diagrama2D.IsoX:
            i[0].set_visible(self.checkIsoX.isChecked())
        for i in self.diagrama2D.IsoX_label:
            i.set_visible(self.checkIsoX.isChecked() and self.IsoX.label)
        self.diagrama3D.draw()
        self.diagrama2D.draw()


    def rellenar_variableTabla(self, i, j):
        """Actualiza los elementos disponibles para el tercer eje en función de los elejidos como ejes principales i, j"""
        self.variableTabla.clear()
        variables=[QtWidgets.QApplication.translate("SteamTables", "Presión", None), QtWidgets.QApplication.translate("SteamTables", "Temperatura", None),  QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None), QtWidgets.QApplication.translate("SteamTables", "Entalpía", None), QtWidgets.QApplication.translate("SteamTables", "Entropía", None), QtWidgets.QApplication.translate("SteamTables", "Energía interna", None), QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None), QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None), QtWidgets.QApplication.translate("SteamTables", "Cp", None), QtWidgets.QApplication.translate("SteamTables", "Cv", None), QtWidgets.QApplication.translate("SteamTables", "Densidad", None), QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None), QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None),  QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None), QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None)]
        del variables[j]
        del variables[i]
        for nombre in variables:
            self.variableTabla.addItem(nombre)

    def rellenar_ejeY(self, int):
        """Rellena las variables disponibles para el ejeY en el gráfico 2D, todos menos el que este activo en el ejeX"""
        self.ejeY.clear()
        Ejes2D=["p", "T", "s", "h", "u", "v"]
        del Ejes2D[int]
        for i in Ejes2D:
            self.ejeY.addItem(i)
            
    def dibujar(self):
        """Método que dibuja todos los datos pedidos"""
        self.progresbar.setValue(100)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Dibujando...", None))
        QtWidgets.QApplication.processEvents()
        self.diagrama3D.plot_3D(self.etiquetas, self.xdata, self.ydata, self.zdata) 
        if len(self.puntos)>0:
            self.diagrama3D.plot_puntos(self.xi, self.yi, self.zi)
        self.diagrama3D.plot_sat(self.xsat, self.ysat, self.zsat)
        self.diagrama3D.plot_Isolinea("T", self.isoterma[0], self.isoterma[1], self.isoterma[2], self.Isoterma.Color, self.Isoterma.Grosor, self.Isoterma.Linea)
        self.diagrama3D.plot_Isolinea("P", self.isobara[0], self.isobara[1], self.isobara[2], self.Isobara.Color, self.Isobara.Grosor, self.Isobara.Linea)
        self.diagrama3D.plot_Isolinea("V", self.isocora[0], self.isocora[1], self.isocora[2], self.Isocora.Color, self.Isocora.Grosor, self.Isocora.Linea)
        self.diagrama3D.plot_Isolinea("H", self.isoentalpica[0], self.isoentalpica[1], self.isoentalpica[2], self.Isoentalpica.Color, self.Isoentalpica.Grosor, self.Isoentalpica.Linea)
        self.diagrama3D.plot_Isolinea("S", self.isoentropica[0], self.isoentropica[1], self.isoentropica[2], self.Isoentropica.Color, self.Isoentropica.Grosor, self.Isoentropica.Linea)
        self.diagrama3D.plot_Isolinea("X", self.isoX[0], self.isoX[1], self.isoX[2], self.IsoX.Color, self.IsoX.Grosor, self.IsoX.Linea)
        self.diagrama3D.axes3D.set_xlim3d(self.xdata[0][0], self.xdata[-1][-1])
        self.diagrama3D.axes3D.set_ylim3d(self.ydata[0][0], self.ydata[-1][-1])
        self.diagrama3D.axes3D.set_zlim3d(min(self.zdata), max(self.zdata))
        self.diagrama3D.draw()
        
        self.diagrama2D.plot_2D(self.etiquetas2, self.rejilla.isChecked())
        self.diagrama2D.plot_sat(self.xsat2, self.ysat2)
        self.diagrama2D.plot_Isolinea("T", self.isoterma2[0], self.isoterma2[1], self.isoterma[2], self.Isoterma.Color, self.Isoterma.Grosor, self.Isoterma.Linea)
        self.diagrama2D.plot_Isolinea("P", self.isobara2[0], self.isobara2[1], self.isobara[2], self.Isobara.Color, self.Isobara.Grosor, self.Isobara.Linea)
        self.diagrama2D.plot_Isolinea("V", self.isocora2[0], self.isocora2[1], self.isocora[2], self.Isocora.Color, self.Isocora.Grosor, self.Isocora.Linea)
        self.diagrama2D.plot_Isolinea("H", self.isoentalpica2[0], self.isoentalpica2[1], self.isoentalpica[2], self.Isoentalpica.Color, self.Isoentalpica.Grosor, self.Isoentalpica.Linea)
        self.diagrama2D.plot_Isolinea("S", self.isoentropica2[0], self.isoentropica2[1], self.isoentropica[2], self.Isoentropica.Color, self.Isoentropica.Grosor, self.Isoentropica.Linea)
        self.diagrama2D.plot_Isolinea("X", self.isoX2[0], self.isoX2[1], self.isoX[2], self.IsoX.Color, self.IsoX.Grosor, self.IsoX.Linea)
        if len(self.puntos)>0:
            self.diagrama2D.plot_puntos(self.xi2, self.yi2)

        self.diagrama2D.plot_labels("X", self.isoX2[2], self.isoX2[3], self.isoX2[4], self.isoX2[5])
        self.diagrama2D.plot_labels("S", self.isoentropica2[2], self.isoentropica2[3], self.isoentropica2[4], self.isoentropica2[5])
        self.diagrama2D.plot_labels("H", self.isoentalpica2[2], self.isoentalpica2[3], self.isoentalpica2[4], self.isoentalpica2[5])
        self.diagrama2D.plot_labels("V", self.isocora2[2], self.isocora2[3], self.isocora2[4], self.isocora2[5])
        self.diagrama2D.plot_labels("P", self.isobara2[2], self.isobara2[3], self.isobara2[4], self.isobara2[5])
        self.diagrama2D.plot_labels("T", self.isoterma2[2], self.isoterma2[3], self.isoterma2[4], self.isoterma2[5])

        self.mostrarSaturacion(self.actionDibujarSaturacion.isChecked())
        self.mostrarPuntos(self.actionMostrarPuntos.isChecked())
        self.mostrarIsoentropica()
        self.mostrarIsocora()
        self.mostrarIsoentalpica(self.checkIsoEnth.isChecked())
        self.mostrarIsobara()
        self.mostrarIsoterma()
        self.mostrarIsoX()
        self.mostrarPuntos(self.actionMostrarPuntos.isChecked())
        self.diagrama2D_ejeX()
        self.diagrama2D_ejeY()
        self.ejeX_log(self.ejeX_escala.isChecked())
        self.ejeY_log(self.ejeY_escala.isChecked())
        self.progresbar.setVisible(False)
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Listo...", None))

    def rejilla_toggled(self, bool):
        """Muestra o esconde la rejilla del gráfico 2D"""
        self.diagrama2D.axes2D.grid(bool)
        self.diagrama2D.draw()
        
    def diagrama2D_ejeX(self):
        """Define la orientación del eje x, creciente o decreciente"""
        if self.ejeX_min.value and self.ejeX_max.value:
            xmin=self.ejeX_min.value
            xmax=self.ejeX_max.value
            self.diagrama2D.axes2D.set_xlim(xmin, xmax)
            self.diagrama2D.axes2D.set_autoscalex_on(False)
        else:
            self.diagrama2D.axes2D.set_autoscalex_on(True)
        self.diagrama2D.draw()

    def diagrama2D_ejeY(self):
        """Define la orientación del eje y, creciente o decreciente"""
        if self.ejeY_min.value and self.ejeY_max.value:
            ymin=self.ejeY_min.value
            ymax=self.ejeY_max.value
            self.diagrama2D.axes2D.set_ylim(ymin, ymax)
            self.diagrama2D.axes2D.set_autoscaley_on(False)
        else:
            self.diagrama2D.axes2D.set_autoscaley_on(True)
        self.diagrama2D.draw()
        
    def ejeX_log(self, bool):
        """Define la escala del eje x, normal o logarítmica"""
        if bool:
            self.diagrama2D.axes2D.set_xscale("log")
        else:
            self.diagrama2D.axes2D.set_xscale("linear")
        self.diagrama2D.draw()

    def ejeY_log(self, bool):
        """Define la escala del eje y, normal o logarítmica"""
        if bool:
            self.diagrama2D.axes2D.set_yscale("log")
        else:
            self.diagrama2D.axes2D.set_yscale("linear")
        self.diagrama2D.draw()
    
    #Métodos de cálculo    
    def definirIsolineas(self):
        Isolineas=[self.Isoterma, self.Isobara, self.Isoentropica, self.Isoentalpica, self.Isocora, self.IsoX]
        iso=["Isotherm", "Isobar", "Isoentropic", "Isoenthalpic", "Isochor", "Isoquality"]
        for i, propiedad in enumerate(iso):
            if self.Config.has_section(propiedad):
                Isolineas[i].inicio=self.Config.getfloat(propiedad, 'Start')
                Isolineas[i].fin=self.Config.getfloat(propiedad, 'End')
                Isolineas[i].intervalo=self.Config.getfloat(propiedad, 'Step')
                Isolineas[i].Personalizar=self.Config.getboolean(propiedad, 'Custom')
                lista=[]
                if self.Config.get(propiedad, 'List')!="":
                    for j in self.Config.get(propiedad, 'List').split(','):
                        lista.append(float(j))
                Isolineas[i].Lista=lista
                Isolineas[i].Critica=self.Config.getboolean(propiedad, 'Critic')
                Isolineas[i].Color=self.Config.get(propiedad, 'Color')
                Isolineas[i].Grosor=self.Config.getfloat(propiedad, 'lineWidth')
                Isolineas[i].Linea=self.Config.getint(propiedad, 'lineStyle')
                Isolineas[i].label=self.Config.getboolean(propiedad, 'Label')
                Isolineas[i].units=self.Config.getboolean(propiedad, 'Units')
                Isolineas[i].posicion=self.Config.getfloat(propiedad, 'Position')
    
    def factores_conversion(self):
        """Método que calcula los factores de conversión de unidades necesarios, tambien los textos"""
        indiceT=self.Config.getint("Units", "Temperature")
        if indiceT==0:
            self.conv_T=float
            self.conv_T_inv=float
        elif indiceT==1:
            self.conv_T=config.C2K
            self.conv_T_inv=config.K2C
        elif indiceT==2:
            self.conv_T=config.F2K
            self.conv_T_inv=config.K2F
        elif indiceT==3:
            self.conv_T=config.R2K
            self.conv_T_inv=config.K2R
        elif indiceT==4:
            self.conv_T=config.Re2K
            self.conv_T_inv=config.K2Re
            
        if self.ejesTabla.currentIndex()==0:
            abcisa="P, %s" % config.Configuracion("Pressure").text()
            ordenada="T, %s" % config.Configuracion("Temperature").text()
            self.factorx=unidades.Pressure(1, config.Configuracion("Pressure").func())
            self.factory=0
        elif self.ejesTabla.currentIndex()==1:
            abcisa="P, %s" % config.Configuracion("Pressure").text()
            ordenada="H, %s" % config.Configuracion("Enthalpy").text()
            self.factorx=unidades.Pressure(1, config.Configuracion("Pressure").func())
            self.factory=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.ejesTabla.currentIndex()==2:
            abcisa="P, %s" % config.Configuracion("Pressure").text()
            ordenada="S, %s" % config.Configuracion("SpecificHeat","Entropy").text()
            self.factorx=unidades.Pressure(1, config.Configuracion("Pressure").func())
            self.factory=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat","Entropy").func())
        elif self.ejesTabla.currentIndex()==3:
            abcisa="P, %s" % config.Configuracion("Pressure").text()
            ordenada="v, %s" % config.Configuracion("SpecificVolume").text()
            self.factorx=unidades.Pressure(1, config.Configuracion("Pressure").func())
            self.factory=unidades.SpecificVolume(1, config.Configuracion("SpecificVolume").func())
        elif self.ejesTabla.currentIndex()==4:
            abcisa="T, %s" % config.Configuracion("Temperature").text()
            ordenada="S, %s" % config.Configuracion("SpecificHeat","Entropy").text()
            self.factorx=0
            self.factory=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat","Entropy").func())
        elif self.ejesTabla.currentIndex()==5:
            abcisa="T, %s" % config.Configuracion("Temperature").text()
            ordenada="x"
            self.factorx=0
            self.factory=1
            
        if self.ejeX.currentText()=="p":
            abcisa2="P, %s" % config.Configuracion("Pressure").text()
            self.factorx2=unidades.Pressure(1, config.Configuracion("Pressure").func())
        elif self.ejeX.currentText()=="T":
            abcisa2="T, %s" % config.Configuracion("Temperature").text()
            self.factorx2=0
        elif self.ejeX.currentText()=="h":
            abcisa2="H, %s" % config.Configuracion("Enthalpy").text()
            self.factorx2=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.ejeX.currentText()=="v":
            abcisa2="v, %s" % config.Configuracion("SpecificVolume").text()
            self.factorx2=unidades.SpecificVolume(1, config.Configuracion("SpecificVolume").func())
        elif self.ejeX.currentText()=="s":
            abcisa2="S, %s" % config.Configuracion("SpecificHeat","Entropy").text()
            self.factorx2=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat","Entropy").func())
        elif self.ejeX.currentText()=="u":
            abcisa2="U, %s" %config.Configuracion("Enthalpy").text()
            self.factorx2=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        if self.ejeY.currentText()=="p":
            ordenada2="P, %s" % config.Configuracion("Pressure").text()
            self.factory2=unidades.Pressure(1, config.Configuracion("Pressure").func())
        elif self.ejeY.currentText()=="T":
            ordenada2="T, %s" % config.Configuracion("Temperature").text()
            self.factory2=0
        elif self.ejeY.currentText()=="h":
            ordenada2="H, %s" % config.Configuracion("Enthalpy").text()
            self.factory2=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.ejeY.currentText()=="v":
            ordenada2="v, %s" % config.Configuracion("SpecificVolume").text()
            self.factory2=unidades.SpecificVolume(1, config.Configuracion("SpecificVolume").func())
        elif self.ejeY.currentText()=="s":
            ordenada2="S, %s" % config.Configuracion("SpecificHeat","Entropy").text()
            self.factory2=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat","Entropy").func())
        elif self.ejeY.currentText()=="u":
            ordenada2="U, %s" %config.Configuracion("Enthalpy").text()
            self.factory2=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())

        if self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Presión", None):
            texto="P, %s" %config.Configuracion("Pressure").text()
            self.factorz=unidades.Pressure(1, config.Configuracion("Pressure").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Temperatura", None):
            texto="T, %s" %config.Configuracion("Temperature").text()
            self.factorz=0
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None):
            texto="v, %s" %config.Configuracion("SpecificVolume").text()
            self.factorz=unidades.SpecificVolume(1, config.Configuracion("SpecificVolume").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entalpía", None):
            texto="H, %s" %config.Configuracion("Enthalpy").text()
            self.factorz=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entropía", None):
            texto="S, %s" %config.Configuracion("SpecificHeat", "Entropy").text()
            self.factorz=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat", "Entropy").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía interna", None):
            texto="U, %s" %config.Configuracion("Enthalpy").text()
            self.factorz=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cp", None):
            texto="Cp, %s" %config.Configuracion("SpecificHeat").text()
            self.factorz=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cv", None):
            texto="Cv, %s" %config.Configuracion("SpecificHeat").text()
            self.factorz=unidades.SpecificHeat(1, config.Configuracion("SpecificHeat").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Densidad", None):
            texto="ρ, %s" %config.Configuracion("Density").text()
            self.factorz=unidades.Density(1, config.Configuracion("Density").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None):
            texto="k, %s" %config.Configuracion("ThermalConductivity").text()
            self.factorz=unidades.ThermalConductivity(1, config.Configuracion("ThermalConductivity").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None):
            texto="μ, %s" %config.Configuracion("Viscosity").text()
            self.factorz=unidades.Viscosity(1, config.Configuracion("Viscosity").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None):
            texto="w, %s" %config.Configuracion("Speed").text()
            self.factorz=unidades.Speed(1, config.Configuracion("Speed").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None):
            texto="x"
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None):
            texto="G, %s" %config.Configuracion("Enthalpy").text()
            self.factorz=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None):
            texto="A, %s" %config.Configuracion("Enthalpy").text()
            self.factorz=unidades.Enthalpy(1, config.Configuracion("Enthalpy").func())
            
        self.etiquetas=[abcisa, ordenada, texto]
        self.etiquetas2=[abcisa2, ordenada2]


    def calcularMatriz(self, start=0, rango=40):
        """Método que actualiza los datos de matriz"""        
        xini=float(self.abscisaInicio.text())
        xfin=float(self.abscisaFin.text())
        xsalto=float(self.abscisaIntervalo.text())                
        xn=int((xfin-xini)/xsalto+1)
        yini=float(self.ordenadaInicio.text())
        yfin=float(self.ordenadaFin.text())
        ysalto=float(self.ordenadaIntervalo.text())
        yn=int((yfin-yini)/ysalto+1)
        self.tabla.setRowCount(yn)
        self.tabla.setColumnCount(xn)
        
        self.factores_conversion()
                    
        xi=arange(xini, xfin, xsalto)
        if (xfin-xini)/xsalto==float(int((xfin-xini)/xsalto)):
            xi=concatenate((xi, [xfin]))
        yi=arange(yini, yfin, ysalto)
        if (yfin-yini)/ysalto==float(int((yfin-yini)/ysalto)):
            yi=concatenate((yi, [yfin]))

        for i in range(len(xi)):
            headerItem = QtWidgets.QTableWidgetItem()
            headerItem.setText(str(xi[i]))
            self.tabla.setHorizontalHeaderItem(i,headerItem)
        for i in range(len(yi)):
            headerItem = QtWidgets.QTableWidgetItem()
            headerItem.setText(str(yi[i]))
            self.tabla.setVerticalHeaderItem(i,headerItem)
            
        xdata,ydata = meshgrid(xi, yi)
        self.matriz=[]

        for i in range(len(xi)):
            self.progresbar.setValue(start+rango*(i+1.)/len(xi))
            QtWidgets.QApplication.processEvents()
            fila=[]
            for j in range(len(yi)):
                if self.ejesTabla.currentIndex()==0:
                    vapor=steam_pT(xi[i]*self.factorx, self.conv_T(yi[j]))
                elif self.ejesTabla.currentIndex()==1:
                    vapor=steam_ph(xi[i]*self.factorx, yi[j]*self.factory)
                elif self.ejesTabla.currentIndex()==2:
                    vapor=steam_ps(xi[i]*self.factorx, yi[j]*self.factory)       
                elif self.ejesTabla.currentIndex()==3:
                    vapor=steam_pv(xi[i]*self.factorx, yi[j]*self.factory)     
                elif self.ejesTabla.currentIndex()==4:
                    if bounds_Ts(xi[i], yi[j], 0)==0:
                        vapor=steam_Ts(self.conv_T(xi[i]), yi[j]*self.factory)
                    else:
                        vapor=steam_Ts(TCRIT, steam_pT(PCRIT, TCRIT).s)
                elif self.ejesTabla.currentIndex()==5:
                    vapor=steam_Tx(self.conv_T(xi[i]), yi[j]) 
                fila.append(vapor)
                    
            self.matriz.append(fila)

        self.xdata=xdata
        self.ydata=ydata
        self.actionCSV.setEnabled(True)


    def Calcular_Propiedades(self, start, rango=5):
        """Método que actualiza los datos al cambiar la propiedad a mostrar"""
        if len(self.matriz)!=0:
            zdata = zeros(self.xdata.shape)
            
            for i, fila in enumerate(self.matriz):
                self.progresbar.setValue(start+rango*(i+1.)/len(self.matriz))
                QtWidgets.QApplication.processEvents()
                for j, vapor in enumerate(fila):
                    if self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Presión", None):
                        dato=vapor.p/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Temperatura", None):
                        dato=self.conv_T_inv(vapor.T)
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None):
                        dato=vapor.v/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entalpía", None):
                        dato=vapor.h/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entropía", None):
                        dato=vapor.s/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía interna", None):
                        dato=vapor.u/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cp", None):
                        dato=vapor.cp/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cv", None):
                        dato=vapor.cv/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Densidad", None):
                        dato=vapor.rho/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None):
                        dato=vapor.k/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None):
                        dato=vapor.mu/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None):
                        if vapor.region !='\x04' and vapor.region !='\x03':
                            dato=vapor.w/self.factorz
                        else:
                            dato=0.0
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None):
                        dato=vapor.x
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None):
                        dato=(vapor.h-vapor.T*vapor.s)/self.factorz
                    elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None):
                        dato=(vapor.u-vapor.T*vapor.s)/self.factorz
                    zdata[j, i]=dato

            for i, fila in enumerate(zdata):
                for j, dato in enumerate(fila):
                    self.tabla.setItem(i, j,QtWidgets.QTableWidgetItem(representacion(dato)))
                    self.tabla.item(i, j).setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)

            self.tabla.resizeColumnsToContents()
            self.toolBox.setItemText(self.toolBox.indexOf(self.page_1), QtWidgets.QApplication.translate("SteamTables", "Tabla", None)+"   %s - %s - %s" %(self.etiquetas[0], self.etiquetas[1], self.etiquetas[2]))
            self.zdata=zdata


    def calcularSaturacion(self):
        """Método que calcula datos de la línea de saturación"""
        TT0 = linspace(273.15, TCRIT, 200)
        psat = [psat_T(T) for T in TT0]
        
        if self.ejesTabla.currentIndex()==0:
            self.xsat=[[P/self.factorx for P in psat], [P/self.factorx for P in psat]]
            self.ysat=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
        elif self.ejesTabla.currentIndex()==1:
            self.xsat=[[P/self.factorx for P in psat], [P/self.factorx for P in psat]]
            self.ysat=[[region4_Tx(T,0).h/self.factory for T in TT0], [region4_Tx(T,1).h/self.factory for T in TT0]]
        elif self.ejesTabla.currentIndex()==2:
            self.xsat=[[P/self.factorx for P in psat], [P/self.factorx for P in psat]]
            self.ysat=[[region4_Tx(T,0).s/self.factory for T in TT0], [region4_Tx(T,1).s/self.factory for T in TT0]]
        elif self.ejesTabla.currentIndex()==3:
            self.xsat=[[P/self.factorx for P in psat], [P/self.factorx for P in psat]]
            self.ysat=[[region4_Tx(T,0).v/self.factory for T in TT0], [region4_Tx(T,1).v/self.factory for T in TT0]]
        elif self.ejesTabla.currentIndex()==4:
            self.xsat=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
            self.ysat=[[region4_Tx(T,0).s/self.factory for T in TT0], [region4_Tx(T,1).s/self.factory for T in TT0]]
        elif self.ejesTabla.currentIndex()==5:
            self.xsat=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
            self.ysat=[[0]*100, [1]*100]
        
        if self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Presión", None):
            self.zsat=[[P/self.factorz for P in psat], [P/self.factorz for P in psat]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Temperatura", None):
            self.zsat=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None):
            self.zsat=[[region4_Tx(T,0).v/self.factorz for T in TT0], [region4_Tx(T,1).v/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entalpía", None):
            self.zsat=[[region4_Tx(T,0).h/self.factorz for T in TT0], [region4_Tx(T,1).h/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entropía", None):
            self.zsat=[[region4_Tx(T,0).s/self.factorz for T in TT0], [region4_Tx(T,1).s/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía interna", None):
            self.zsat=[[region4_Tx(T,0).u/self.factorz for T in TT0], [region4_Tx(T,1).u/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cp", None):
            self.zsat=[[region4_Tx(T,0).cp/self.factorz for T in TT0], [region4_Tx(T,1).cp/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cv", None):
            self.zsat=[[region4_Tx(T,0).cv/self.factorz for T in TT0], [region4_Tx(T,1).cv/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Densidad", None):
            self.zsat=[[region4_Tx(T,0).rho/self.factorz for T in TT0], [region4_Tx(T,1).rho/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None):
            self.zsat=[[region4_Tx(T,0).k/self.factorz for T in TT0], [region4_Tx(T,1).k/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None):
            self.zsat=[[region4_Tx(T,0).mu/self.factorz for T in TT0], [region4_Tx(T,1).mu/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None):
            self.zsat=[[0 for T in TT0], [0 for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None):
            self.zsat=[[0]*len(TT0), [1]*len(TT0)]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None):
            self.zsat=[[(region4_Tx(T,0).h-T*region4_Tx(T,0).s)/self.factorz for T in TT0], [(region4_Tx(T,1).h-T*region4_Tx(T,1).s)/self.factorz for T in TT0]]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None):
            self.zsat=[[(region4_Tx(T,0).u-T*region4_Tx(T,0).s)/self.factorz for T in TT0], [(region4_Tx(T,1).u-T*region4_Tx(T,1).s)/self.factorz for T in TT0]]
        
        for i in range(len(self.xsat)):
            for j in range(len(self.xsat[i])-1, -1, -1):
                if self.xsat[i][j]<self.xdata[0][0] or self.xsat[i][j]>self.xdata[-1][-1] or self.ysat[i][j]<self.ydata[0][0] or self.ysat[i][j]>self.ydata[-1][-1]:
                    del self.xsat[i][j]
                    del self.ysat[i][j]
                    del self.zsat[i][j]
    
        if self.ejeX.currentText()=="p":
            self.xsat2=[[P/self.factorx2 for P in psat], [P/self.factorx2 for P in psat]]
        elif self.ejeX.currentText()=="T":
            self.xsat2=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
        elif self.ejeX.currentText()=="h":
            self.xsat2=[[region4_Tx(T,0).h/self.factorx2 for T in TT0], [region4_Tx(T,1).h/self.factorx2 for T in TT0]]
        elif self.ejeX.currentText()=="v":
            self.xsat2=[[region4_Tx(T,0).v/self.factorx2 for T in TT0], [region4_Tx(T,1).v/self.factorx2 for T in TT0]]
        elif self.ejeX.currentText()=="s":
            self.xsat2=[[region4_Tx(T,0).s/self.factorx2 for T in TT0], [region4_Tx(T,1).s/self.factorx2 for T in TT0]]
        elif self.ejeX.currentText()=="u":
            self.xsat2=[[region4_Tx(T,0).u/self.factorx2 for T in TT0], [region4_Tx(T,1).u/self.factorx2 for T in TT0]]
        if self.ejeY.currentText()=="p":
            self.ysat2=[[P/self.factory2 for P in psat], [P/self.factory2 for P in psat]]
        elif self.ejeY.currentText()=="T":
            self.ysat2=[[self.conv_T_inv(T) for T in TT0], [self.conv_T_inv(T) for T in TT0]]
        elif self.ejeY.currentText()=="h":
            self.ysat2=[[region4_Tx(T,0).h/self.factory2 for T in TT0], [region4_Tx(T,1).h/self.factory2 for T in TT0]]
        elif self.ejeY.currentText()=="v":
            self.ysat2=[[region4_Tx(T,0).v/self.factory2 for T in TT0], [region4_Tx(T,1).v/self.factory2 for T in TT0]]
        elif self.ejeY.currentText()=="s":
            self.ysat2=[[region4_Tx(T,0).s/self.factory2 for T in TT0], [region4_Tx(T,1).s/self.factory2 for T in TT0]]
        elif self.ejeY.currentText()=="u":
            self.ysat2=[[region4_Tx(T,0).u/self.factory2 for T in TT0], [region4_Tx(T,1).u/self.factory2 for T in TT0]]


    def calcularPuntos(self):
        """Método que actualiza los datos de puntos definidos por el usuario"""
        if self.ejesTabla.currentIndex()==0:
            self.xi=[punto.p/self.factorx for punto in self.puntos]
            self.yi=[self.conv_T_inv(punto.T) for punto in self.puntos]
        elif self.ejesTabla.currentIndex()==1:
            self.xi=[punto.p/self.factorx for punto in self.puntos]
            self.yi=[punto.h/self.factory for punto in self.puntos]
        elif self.ejesTabla.currentIndex()==2:
            self.xi=[punto.p/self.factorx for punto in self.puntos]
            self.yi=[punto.s/self.factory for punto in self.puntos]
        elif self.ejesTabla.currentIndex()==3:
            self.xi=[punto.p/self.factorx for punto in self.puntos]
            self.yi=[punto.v/self.factory for punto in self.puntos]
        elif self.ejesTabla.currentIndex()==4:
            self.xi=[self.conv_T_inv(punto.T) for punto in self.puntos]
            self.yi=[punto.s/self.factory for punto in self.puntos]
        elif self.ejesTabla.currentIndex()==5:
            self.xi=[self.conv_T_inv(punto.T) for punto in self.puntos]
            self.yi=[punto.x for punto in self.puntos]
            
        if self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Presión", None):
            self.zi=[punto.p/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Temperatura", None):
            self.zi=[self.conv_T_inv(punto.T) for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None):
            self.zi=[punto.v/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entalpía", None):
            self.zi=[punto.h/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entropía", None):
            self.zi=[punto.s/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía interna", None):
            self.zi=[punto.u/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cp", None):
            self.zi=[punto.cp/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cv", None):
            self.zi=[punto.cv/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Densidad", None):
            self.zi=[punto.rho/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None):
            self.zi=[punto.k/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None):
            self.zi=[punto.mu/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None):
            self.zi=[0 for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None):
            self.zi=[punto.x for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None):
            self.zi=[(punto.h-punto.T*punto.s)/self.factorz for punto in self.puntos]
        elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None):
            self.zi=[(punto.u-punto.T*punto.s)/self.factorz for punto in self.puntos]
                
        if self.ejeX.currentText()=="p":
            self.xi2=[punto.p/self.factorx2 for punto in self.puntos]
        elif self.ejeX.currentText()=="T":
            self.xi2=[self.conv_T_inv(punto.T) for punto in self.puntos]
        elif self.ejeX.currentText()=="h":
            self.xi2=[punto.h/self.factorx2 for punto in self.puntos]
        elif self.ejeX.currentText()=="v":
            self.xi2=[punto.v/self.factorx2 for punto in self.puntos]
        elif self.ejeX.currentText()=="s":
            self.xi2=[punto.s/self.factorx2 for punto in self.puntos]
        elif self.ejeX.currentText()=="u":
            self.xi2=[punto.u/self.factorx2 for punto in self.puntos]
        if self.ejeY.currentText()=="p":
            self.yi2=[punto.p/self.factory2 for punto in self.puntos]
        elif self.ejeY.currentText()=="T":
            self.yi2=[self.conv_T_inv(punto.T) for punto in self.puntos]
        elif self.ejeY.currentText()=="h":
            self.yi2=[punto.h/self.factory2 for punto in self.puntos]
        elif self.ejeY.currentText()=="v":
            self.yi2=[punto.v/self.factory2 for punto in self.puntos]
        elif self.ejeY.currentText()=="s":
            self.yi2=[punto.s/self.factory2 for punto in self.puntos]
        elif self.ejeY.currentText()=="u":
            self.yi2=[punto.u/self.factory2 for punto in self.puntos]

        
    def isolineas(self, S, X, funcion, start, rango):
        """Librería de cálculo de los parámetros de las isolineas"""
        x=[]
        y=[]
        z=[]
        x2=[]
        y2=[]
        for i, propiedad in enumerate(S):
            self.progresbar.setValue(start+rango*(i+1.)/len(S))
            QtWidgets.QApplication.processEvents()
            if self.ejesTabla.currentIndex()==0:
                xi=[funcion(i, propiedad).p/self.factorx for i in X]
                yi=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
            elif self.ejesTabla.currentIndex()==1:
                xi=[funcion(i, propiedad).p/self.factorx for i in X]
                yi=[funcion(i, propiedad).h/self.factory for i in X]
            elif self.ejesTabla.currentIndex()==2:
                xi=[funcion(i, propiedad).p/self.factorx for i in X]
                yi=[funcion(i, propiedad).s/self.factory for i in X]
            elif self.ejesTabla.currentIndex()==3:
                xi=[funcion(i, propiedad).p/self.factorx for i in X]
                yi=[funcion(i, propiedad).v/self.factory for i in X]
            elif self.ejesTabla.currentIndex()==4:
                xi=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
                yi=[funcion(i, propiedad).s/self.factory for i in X]
            elif self.ejesTabla.currentIndex()==5:
                xi=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
                yi=[funcion(i, propiedad).x for i in X]
                
            if self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Presión", None):
                zi=[funcion(i, propiedad).p/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Temperatura", None):
                zi=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Volumen específico", None):
                zi=[funcion(i, propiedad).v/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entalpía", None):
                zi=[funcion(i, propiedad).h/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Entropía", None):
                zi=[funcion(i, propiedad).s/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía interna", None):
                zi=[funcion(i, propiedad).u/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cp", None):
                zi=[funcion(i, propiedad).cp/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Cv", None):
                zi=[funcion(i, propiedad).cv/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Densidad", None):
                zi=[funcion(i, propiedad).rho/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Conductividad térmica", None):
                zi=[funcion(i, propiedad).k/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Viscosidad", None):
                zi=[funcion(i, propiedad).mu/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Velocidad del sonido", None):
                zi=[0 for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Fracción de vapor", None):
                zi=[funcion(i, propiedad).x for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Gibbs", None):
                zi=[(funcion(i, propiedad).h-funcion(i, propiedad).T*funcion(i, propiedad).s)/self.factorz for i in X]
            elif self.variableTabla.currentText()==QtWidgets.QApplication.translate("SteamTables", "Energía de Helmholtz", None):
                zi=[(funcion(i, propiedad).u-funcion(i, propiedad).T*funcion(i, propiedad).s)/self.factorz for i in X]
                
            if self.ejeX.currentText()=="p":
                xi2=[funcion(i, propiedad).p/self.factorx2 for i in X]
            elif self.ejeX.currentText()=="T":
                xi2=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
            elif self.ejeX.currentText()=="h":
                xi2=[funcion(i, propiedad).h/self.factorx2 for i in X]
            elif self.ejeX.currentText()=="v":
                xi2=[funcion(i, propiedad).v/self.factorx2 for i in X]
            elif self.ejeX.currentText()=="s":
                xi2=[funcion(i, propiedad).s/self.factorx2 for i in X]
            elif self.ejeX.currentText()=="u":
                xi2=[funcion(i, propiedad).u/self.factorx2 for i in X]
            if self.ejeY.currentText()=="p":
                yi2=[funcion(i, propiedad).p/self.factory2 for i in X]
            elif self.ejeY.currentText()=="T":
                yi2=[self.conv_T_inv(funcion(i, propiedad).T) for i in X]
            elif self.ejeY.currentText()=="h":
                yi2=[funcion(i, propiedad).h/self.factory2 for i in X]
            elif self.ejeY.currentText()=="v":
                yi2=[funcion(i, propiedad).v/self.factory2 for i in X]
            elif self.ejeY.currentText()=="s":
                yi2=[funcion(i, propiedad).s/self.factory2 for i in X]
            elif self.ejeY.currentText()=="u":
                yi2=[funcion(i, propiedad).u/self.factory2 for i in X]
            x.append(xi)
            y.append(yi)
            z.append(zi)
            x2.append(xi2)
            y2.append(yi2)
        return x, y, z, x2, y2
 
    def labels(self, isolineas, x, y, pos):
        x_label=[]
        y_label=[]
        label=[]
        angle=[]
        for i in range(len(isolineas)):
            j=int(pos/100.*len(x[i]))
            x_label.append(x[i][j])
            y_label.append(y[i][j])
            label.append(representacion(isolineas[i]))
            if self.ejeX_escala.isChecked():
                fraccionx=(log(x[i][j+1])-log(x[i][j]))/(log(self.ejeX_max.value)-log(self.ejeX_min.value))
            else:
                fraccionx=(x[i][j+1]-x[i][j])/(self.ejeX_max.value-self.ejeX_min.value)
            if self.ejeY_escala.isChecked():
                fracciony=(log(y[i][j+1])-log(y[i][j]))/(log(self.ejeY_max.value)-log(self.ejeY_min.value))
            else:
                fracciony=(y[i][j+1]-y[i][j])/(self.ejeY_max.value-self.ejeY_min.value)
            try:
                angle.append(arctan(fracciony/fraccionx)*360/2/pi)
            except ZeroDivisionError:
                angle.append(90)
        return x_label, y_label, label, angle
        
        
    def calcularIsoentropica(self, start=85, rango=5):
        """Método que actualiza los datos de isoentrópicas"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isoentrópicas...", None))
        if self.Isoentropica.Personalizar:
            S=self.Isoentropica.Lista
        else:
            S=arange(self.Isoentropica.inicio, self.Isoentropica.fin, self.Isoentropica.intervalo).tolist()
        if self.Isoentropica.Critica:
            SCRIT=steam_pT(PCRIT, TCRIT).s
            S.append(SCRIT)
        S2=S[:]
        X= logspace(-3, 3, 100)*1e5
        x, y, z, x2, y2=self.isolineas(S, X, steam_ps, start, rango)
        for i in range(len(S)-1, -1, -1):
            for j in range(len(X)-1, -1, -1):
                if bounds_ps(X[j], S[i], 0)!=0:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
                    del x2[i][j]
                    del y2[i][j]
            if len(x2[i])==0:
                del x2[i]
                del y2[i]
                S2.remove(S2[i])
        for i in range(len(S)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
            if len(x[i])==0:
                del x[i]
                del y[i]
                del z[i]        
                S.remove(S[i])
        x_label, y_label, label, angle=self.labels([unidades.SpecificHeat(i).config("Entropy") for i in S2], x2, y2, self.Isoentropica.posicion)
        self.isoentropica=[x, y, z]
        if self.Isoentropica.units:
            self.isoentropica2=[x2, y2, x_label, y_label, ["S="+i+config.Configuracion("SpecificHeat", "Entropy").text() for i in label], angle]
        else:
            self.isoentropica2=[x2, y2, x_label, y_label, label, angle]

    def calcularIsoentalpica(self, start=90, rango=5):
        """Método que actualiza los datos de isoentálpicas"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isoentálpicas...", None))
        if self.Isoentalpica.Personalizar:
            H=self.Isoentalpica.Lista
        else:
            H=arange(self.Isoentalpica.inicio, self.Isoentalpica.fin, self.Isoentalpica.intervalo).tolist()
        if self.Isoentalpica.Critica:
            HCRIT=steam_pT(PCRIT, TCRIT).h
            H.append(HCRIT)
        H2=H[:]
        X= logspace(-3, 3, 100)*1e5
        x, y, z, x2, y2=self.isolineas(H, X, steam_ph, start, rango)
        for i in range(len(H)-1, -1, -1):
            for j in range(len(X)-1, -1, -1):
                if bounds_ph(X[j], H[i], 0)!=0:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
                    del x2[i][j]
                    del y2[i][j]
            if len(x2[i])==0:
                del x2[i]
                del y2[i]
                H2.remove(H2[i])
        for i in range(len(H)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
            if len(x[i])==0:
                del x[i]
                del y[i]
                del z[i]
                H.remove(H[i])
        x_label, y_label, label, angle=self.labels([unidades.Enthalpy(i).config() for i in H2], x2, y2, self.Isoentalpica.posicion)
        self.isoentalpica=[x, y, z]
        if self.Isoentalpica.units:
            self.isoentalpica2=[x2, y2, x_label, y_label, ["H="+i+config.Configuracion("Enthalpy").text() for i in label], angle]
        else:
            self.isoentalpica2=[x2, y2, x_label, y_label, label, angle]

    def calcularIsobara(self, start=80, rango=5):
        """Método que actualiza los datos de isobaras"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isobaras...", None))
        if self.Isobara.Personalizar:
            P=self.Isobara.Lista
        else:
            P=arange(self.Isobara.inicio, self.Isobara.fin, self.Isobara.intervalo).tolist()
        if self.Isobara.Critica:
            P.append(PCRIT)
        X=linspace(0, 10000, 100)
        x, y, z, x2, y2=self.isolineas(X, P, steam_ps, start, rango)
        x=transpose(x)
        y=transpose(y)
        z=transpose(z)
        x2=transpose(x2)
        y2=transpose(y2)
        x=[list(i) for i in x]
        y=[list(i) for i in y]
        z=[list(i) for i in z]
        #Eliminamos puntos del grafico 3D fuera de los ejes
        for i in range(len(P)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
            if len(x[i])==0:
                del x[i]
                del y[i]
                del z[i]

        x_label, y_label, label, angle=self.labels([unidades.Pressure(i).config() for i in P], x2, y2, self.Isobara.posicion)
        self.isobara=[x, y, z]
        if self.Isobara.units:
            self.isobara2=[x2, y2, x_label, y_label, ["P="+i+config.Configuracion("Pressure").text() for i in label], angle]
        else:
            self.isobara2=[x2, y2, x_label, y_label, label, angle]

    def calcularIsoterma(self, start=75, rango=5):
        """Método que actualiza los datos de isotermas"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isotermas...", None))
        if self.Isoterma.Personalizar:
            T=self.Isoterma.Lista
        else:
            T=arange(self.Isoterma.inicio, self.Isoterma.fin, self.Isoterma.intervalo).tolist()
        if self.Isoterma.Critica:
            T.append(TCRIT)
        X= logspace(-3, 3, 200)*1e5
        x, y, z, x2, y2=self.isolineas(T, X, steam_pT, start, rango)
        
        #Añadimos puntos interiores de la campana de saturación
        X=linspace(1, 0, 50)
        for i , t in enumerate(T):
            if t<TCRIT:
                xi, yi, zi, xi2, yi2=self.isolineas(X, [t], steam_Tx, start+rango, 0)
                temp=x2[i]+xi2[0]
                temp.sort(reverse=True)
                indice=temp.index(xi2[0][0])
                for j in range(len(xi2)):
                    x2[i].insert(indice+j, xi2[j][0])
                    y2[i].insert(indice+j, yi2[j][0])
                    
        #Eliminamos puntos del grafico 3D fuera de los ejes
        for i in range(len(T)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
            if len(x[i])==0:
                del x[i]
                del y[i]
                del z[i]
        
        x_label, y_label, label, angle=self.labels([unidades.Temperature(i).config() for i in T], x2, y2, self.Isoterma.posicion)
        self.isoterma=[x, y, z]
        if self.Isoterma.units:
            self.isoterma2=[x2, y2, x_label, y_label, ["T="+i+config.Configuracion("Temperature").text() for i in label], angle]
        else:
            self.isoterma2=[x2, y2, x_label, y_label, label, angle]

    def calcularIsocora(self, start=95, rango=5):
        """Método que actualiza los datos de isocoras"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isocoras...", None))
        if self.Isocora.Personalizar:
            V=self.Isocora.Lista
        else:
            V=arange(self.Isocora.inicio, self.Isocora.fin, self.Isocora.intervalo).tolist()
        if self.Isocora.Critica:
            VCRIT=steam_pT(PCRIT, TCRIT).v
            V.append(VCRIT)
        V2=V[:]
        X= logspace(-3, 3, 300)*1e5
        x, y, z, x2, y2=self.isolineas(V, X, steam_pv, start, rango)
        for i in range(len(V)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if bounds_pv(X[j], V[i], 0)!=0:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
                    del x2[i][j]
                    del y2[i][j]
            if len(x2[i])==0:
                del x2[i]
                del y2[i]
                V2.remove(V2[i])
        for i in range(len(V)-1, -1, -1):
            for j in range(len(x[i])-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
            if len(x[i])==0:
                del x[i]
                del y[i]
                del z[i]
                V.remove(V[i])
        x_label, y_label, label, angle=self.labels([unidades.SpecificVolume(i).config() for i in V2], x2, y2, self.Isocora.posicion)
        self.isocora=[x, y, z]     
        if self.Isocora.units:
            self.isocora2=[x2, y2, x_label, y_label, ["v="+i+config.Configuracion("SpecificVolume").text() for i in label], angle]
        else:
            self.isocora2=[x2, y2, x_label, y_label, label, angle]
        
    def calcularIsoX(self, start=95, rango=5):
        """Método que actualiza los datos de isocalidad"""
        self.statusbar.showMessage(QtWidgets.QApplication.translate("SteamTables", "Calculando isocalidades...", None))
        if self.IsoX.Personalizar:
            X=self.IsoX.Lista
        else:
            X=arange(self.IsoX.inicio, self.IsoX.fin, self.IsoX.intervalo)
        T= linspace(273.15, TCRIT, 100)
        x, y, z, x2, y2=self.isolineas(X, T, steam_Tx, start, rango)
        for i in range(len(X)-1, -1, -1):
            for j in range(len(T)-1, -1, -1):
                if x[i][j]<self.xdata[0][0] or x[i][j]>self.xdata[-1][-1] or y[i][j]<self.ydata[0][0] or y[i][j]>self.ydata[-1][-1]:
                    del x[i][j]
                    del y[i][j]
                    del z[i][j]
        x_label, y_label, label, angle=self.labels(X, x2, y2, self.IsoX.posicion)
        self.isoX=[x, y, z]        
        if self.IsoX.units:
            self.isoX2=[x2, y2, x_label, y_label, ["x="+i for i in label], angle]
        else:
            self.isoX2=[x2, y2, x_label, y_label, label, angle]

    def calcularPropiedades(self):
        """Calcula y muestra las propiedades del punto especificado"""
        if self.variablesCalculo.currentIndex()==0:
            p=self.presion.value
            T=self.temperatura.value
            vapor=steam_pT(p, T)
        elif self.variablesCalculo.currentIndex()==1:
            p=self.presion.value
            h=self.entalpia.value
            vapor=steam_ph(p, h)
        elif self.variablesCalculo.currentIndex()==2:
            p=self.presion.value
            s=self.entropia.value
            vapor=steam_ps(p, s)            
        elif self.variablesCalculo.currentIndex()==3:
            p=self.presion.value
            v=self.volumen.value
            vapor=steam_pv(p, v)            
        elif self.variablesCalculo.currentIndex()==4:
            T=self.temperatura.value
            s=self.entropia.value
            vapor=steam_Ts(T, s)            
        elif self.variablesCalculo.currentIndex()==5:
            T=self.temperatura.value
            x=self.fraccionVapor.value
            vapor=steam_Tx(T, x)  
        self.punto=vapor
        self.mostrarPropiedades()
        self.botonAdd.setEnabled(True) 
                
    def mostrarPropiedades(self):
        """Muestra los valores de las propiedades de la sección de puntos especificados por el usuario"""
        if self.punto!=0:
            self.presion.setValue(self.punto.p)
            self.temperatura.setValue(self.punto.T)
            self.volumen.setValue(self.punto.v)
            self.entalpia.setValue(self.punto.h)
            self.entropia.setValue(self.punto.s)
            self.fraccionVapor.setValue(self.punto.x)
            self.energiaInterna.setValue(self.punto.u)
            self.energiaGibbs.setValue(self.punto.h-self.punto.T*self.punto.s)
            self.energiaHelmholtz.setValue(self.punto.u-self.punto.T*self.punto.s)
            self.densidad.setValue(self.punto.rho)
            self.cp.setValue(self.punto.cp)
            self.cv.setValue(self.punto.cv)
            self.conductividad.setValue(self.punto.k)
            self.viscosidad.setValue(self.punto.mu)
            if self.punto.region !='\x04' and  self.punto.region !='\x03':
                self.velocidad.setValue(self.punto.w)
            else:
                self.velocidad.clear()
            if self.punto.region =='\x01':
                self.region.setText("1")
            elif self.punto.region =='\x02':
                self.region.setText("2")
            elif self.punto.region =='\x03':
                self.region.setText("3")
            elif self.punto.region =='\x04':
                self.region.setText("4")
            elif self.punto.region =='\x05':
                self.region.setText("5")
                
    def actualizarConfiguracion(self):
        """Actualiza los diferentes parámetros que puedan cambiar al cerrar el dialogo de preferencias
            Factores de conversión si han cambiado las unidades
                Valores de la configuración de la tabla 3D
                Valores en la sección de puntos especificados por el usuario
            """
        #TODO: Añadir tareas al cambiar la configuración
        self.factores_conversion()
        if self.factorx2==0:  #El eje x es la temperatura
            xmax=unidades.Temperature(self.ejeX_max.value, config.Configuracion("Temperature").func())
            xmin=unidades.Temperature(self.ejeX_min.value, config.Configuracion("Temperature").func())
            self.ejeX_max.setValue(representacion(xmax.config()))
            self.ejeX_min.setValue(representacion(xmin.config()))
        else: #En cualquier otro caso basta con usar el factor de correción para ese eje
            xmax=float(self.ejeX_max.value)*self.factorx2
            xmin=float(self.ejeX_min.value)*self.factorx2
            self.ejeX_max.setValue(representacion(xmax/self.factorx2))
            self.ejeX_min.setValue(representacion(xmin/self.factorx2))
        if self.factory2==0:  
            ymax=unidades.Temperature(self.ejeY_max.value, config.Configuracion("Temperature").func())
            ymin=unidades.Temperature(self.ejeY_min.value, config.Configuracion("Temperature").func())
            self.ejeY_max.setValue(representacion(ymax.config()))
            self.ejeY_min.setValue(representacion(ymin.config()))
        else: 
            ymax=float(self.ejeY_max.value)*self.factory2
            ymin=float(self.ejeY_min.value)*self.factory2      
            self.ejeY_max.setValue(representacion(ymax/self.factory2))
            self.ejeY_min.setValue(representacion(ymin/self.factory2))
        if self.factory2==0:  
            self.ejeY_max.setValue(representacion(ymax.config))
            self.ejeY_min.setValue(representacion(ymin.config))
        else: 
            self.ejeY_max.setValue(representacion(ymax/self.factory2))
            self.ejeY_min.setValue(representacion(ymin/self.factory2))
            
        self.mostrarUnidades()

    def mostrarUnidades(self):
        """Muestra el texto de las unidades en función de la configuración"""
        self.presion.actualizar()
        self.temperatura.actualizar()
        self.volumen.actualizar()
        self.entropia.actualizar()
        self.entalpia.actualizar()
        self.energiaInterna.actualizar()
        self.energiaGibbs.actualizar()
        self.energiaHelmholtz.actualizar()
        self.densidad.actualizar()
        self.cp.actualizar()
        self.cv.actualizar()
        self.conductividad.actualizar()
        self.viscosidad.actualizar()
        self.velocidad.actualizar()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    locale = QtCore.QLocale.system().name()
    qtTranslator = QtCore.QTranslator()
    if qtTranslator.load("UI_steamTables_" + locale, "i18n"):
        app.installTranslator(qtTranslator)
    SteamTables= Ui_SteamTables()
    SteamTables.show()
    sys.exit(app.exec_())
