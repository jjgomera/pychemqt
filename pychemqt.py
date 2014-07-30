#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import urllib2
import shutil

if sys.version_info[0] == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser

from PyQt4 import QtCore, QtGui

path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(path)

conf_dir = os.path.expanduser('~') + os.sep+".pychemqt"+os.sep
os.environ["pychemqt"] = path + os.path.sep

app = QtGui.QApplication(sys.argv)
app.setOrganizationName("pychemqt")
app.setOrganizationDomain("pychemqt")
app.setApplicationName("pychemqt")

# Translation
locale = QtCore.QLocale.system().name()
myTranslator = QtCore.QTranslator()
if myTranslator.load("pychemqt_" + locale, os.environ["pychemqt"] + "i18n"):
    app.installTranslator(myTranslator)
qtTranslator = QtCore.QTranslator()
if qtTranslator.load("qt_" + locale,
   QtCore.QLibraryInfo.location(QtCore.QLibraryInfo.TranslationsPath)):
    app.installTranslator(qtTranslator)

# Check external modules
from tools.dependences import optional_modules
for module, use in optional_modules:
    try:
        __import__(module)
        os.environ[module] = "True"
    except ImportError:
        print(QtGui.QApplication.translate("pychemqt", "%s don't found, %s"
                                           % (module, use)).toUtf8())
        os.environ[module] = ""


class SplashScreen(QtGui.QSplashScreen):
    """Clase que crea una ventana de splash"""
    def __init__(self):
        QtGui.QSplashScreen.__init__(self, QtGui.QPixmap(os.environ["pychemqt"]
                                                         +"/images/splash.jpg"))
        self.show()
        QtGui.QApplication.flush()

    def showMessage(self, msg):
        """Método para mostrar mensajes en la parte inferior de la ventana de splash"""
        labelAlignment = QtCore.Qt.Alignment(QtCore.Qt.AlignBottom |
                                             QtCore.Qt.AlignRight | QtCore.Qt.AlignAbsolute)
        QtGui.QSplashScreen.showMessage(self, msg, labelAlignment, QtGui.QColor(QtCore.Qt.white))
        QtGui.QApplication.processEvents()

    def clearMessage(self):
        QtGui.QSplashScreen.clearMessage(self)
        QtGui.QApplication.processEvents()

splash=SplashScreen()

#Check config files
splash.showMessage(QtGui.QApplication.translate("pychemqt", "Checking config files..."))
from lib import firstrun
if not os.path.isdir(conf_dir):
    os.mkdir(conf_dir)

if not os.path.isfile(conf_dir + "pychemqtrc"):
    Preferences = firstrun.Preferences()
    Preferences.write(open(conf_dir + "pychemqtrc", "w"))

#FIXME: Hasta que no sepa como prescindir de este archivo sera necesario
if not os.path.isfile(conf_dir + "pychemqtrc_temporal"):
    Config=firstrun.config()
    Config.write(open(conf_dir + "pychemqtrc_temporal", "w"))

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Checking cost index..."))
if not os.path.isfile(conf_dir + "CostIndex.dat"):
        with open(os.environ["pychemqt"] + "dat/costindex.dat") as cost_index:
            lista = cost_index.readlines()[-1][:-1].split(" ")
            with open(conf_dir + "CostIndex.dat", "w") as archivo:
                for data in lista:
                    archivo.write(data + os.linesep)

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Checking currency data"))
if not os.path.isfile(conf_dir+"moneda.dat"):
    #Add default currency from system without internet connection
    #Check internet connection
    #Copy file to pychemqt config dir if no connection available
    from lib.firstrun import getdata
    try:
        getdata(conf_dir+"moneda.dat")
    except urllib2.URLError:
        shutil.copy(os.environ["pychemqt"]+"dat"+os.sep+"moneda.dat", conf_dir+"moneda.dat")
        print(QtGui.QApplication.translate("pychemqt", "Internet connection error, using archived currency rates"))

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Checking custom database..."))
from lib.sql import createDatabase

#Import internal libraries
splash.showMessage(QtGui.QApplication.translate("pychemqt", "Importing libraries..."))
from UI import texteditor, newComponent, flujo, charts, plots, viewComponents
from UI.widgets import createAction, ClickableLabel, TreeEquipment
from lib.config import conf_dir, getComponents
from lib.project import Project
from lib.EoS import K, H
from lib import unidades


splash.showMessage(QtGui.QApplication.translate("pychemqt", "Importing equipments..."))
from equipment import *

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Importing tools..."))
from tools import UI_confComponents, UI_confTransport, UI_confThermo, UI_confUnits, UI_confResolution, UI_databank, UI_unitConverter, UI_steamTables, UI_psychrometry
from UI.conversor_unidades import moneda

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Loading main window..."))
from UI.mainWindow import UI_pychemqt
pychemqt = UI_pychemqt()

splash.showMessage(QtGui.QApplication.translate("pychemqt", "Loading project files..."))
if pychemqt.Preferences.getboolean("General", 'Load_Last_Project'):
    filename = pychemqt.settings.value("LastFile").toStringList()
    for file in sys.argv[1:]:
        filename.append(file)
    for fname in filename:
        if fname and QtCore.QFile.exists(fname):
            splash.showMessage(QtGui.QApplication.translate("pychemqt", "Loading project files...")+"\n"+fname)
            pychemqt.fileOpen(fname)


splash.finish(pychemqt)
pychemqt.show()
sys.exit(app.exec_())
