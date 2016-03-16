#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


import os
import sys
import urllib.error
import shutil
import logging
from optparse import OptionParser

from PyQt5 import QtCore, QtGui, QtWidgets


path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(path)

conf_dir = os.path.expanduser('~') + os.sep+".pychemqt"+os.sep
os.environ["pychemqt"] = path + os.path.sep

app = QtWidgets.QApplication(sys.argv)
app.setOrganizationName("pychemqt")
app.setOrganizationDomain("pychemqt")
app.setApplicationName("pychemqt")


# Parse command line options
parser = OptionParser()
parser.add_option("--debug", action="store_true")
parser.add_option("-l", "--log", dest="loglevel", default="INFO")
(options, args) = parser.parse_args()

if options.debug:
    loglevel = "DEBUG"
else:
    loglevel = options.loglevel
loglevel = getattr(logging, loglevel.upper())

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
        print((QtWidgets.QApplication.translate("pychemqt", "%s don't found, %s"
                                           % (module, use))))
        os.environ[module] = ""

class SplashScreen(QtWidgets.QSplashScreen):
    """Clase que crea una ventana de splash"""
    def __init__(self):
        QtWidgets.QSplashScreen.__init__(self,
              QtGui.QPixmap(os.environ["pychemqt"] + "/images/splash.jpg"))
        self.show()
        QtWidgets.QApplication.flush()

    def showMessage(self, msg):
        """Método para mostrar mensajes en la parte inferior de la ventana de
        splash"""
        labelAlignment = QtCore.Qt.Alignment(QtCore.Qt.AlignBottom |
                                             QtCore.Qt.AlignRight |
                                             QtCore.Qt.AlignAbsolute)
        QtWidgets.QSplashScreen.showMessage(self, msg, labelAlignment,
                                        QtGui.QColor(QtCore.Qt.white))
        QtWidgets.QApplication.processEvents()

    def clearMessage(self):
        QtWidgets.QSplashScreen.clearMessage(self)
        QtWidgets.QApplication.processEvents()

splash = SplashScreen()

# Check config files
splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Checking config files..."))
from lib import firstrun
if not os.path.isdir(conf_dir):
    os.mkdir(conf_dir)

if not os.path.isfile(conf_dir + "pychemqtrc"):
    Preferences = firstrun.Preferences()
    Preferences.write(open(conf_dir + "pychemqtrc", "w"))

# FIXME: Hasta que no sepa como prescindir de este archivo sera necesario
if not os.path.isfile(conf_dir + "pychemqtrc_temporal"):
    Config = firstrun.config()
    Config.write(open(conf_dir + "pychemqtrc_temporal", "w"))

# Logging configuration
logging.basicConfig(filename=conf_dir+'pychemqt.log', filemode='w',
                    level=loglevel, datefmt='%d-%b-%Y %H:%M:%S',
                    format='[%(asctime)s.%(msecs)d] %(levelname)s: %(message)s')
logging.info(QtWidgets.QApplication.translate("pychemqt",
                                          "Starting pychemqt"))

splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Checking cost index..."))
if not os.path.isfile(conf_dir + "CostIndex.dat"):
        with open(os.environ["pychemqt"] + "dat/costindex.dat") as cost_index:
            lista = cost_index.readlines()[-1][:-1].split(" ")
            with open(conf_dir + "CostIndex.dat", "w") as archivo:
                for data in lista:
                    archivo.write(data + os.linesep)

splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Checking currency data"))
if not os.path.isfile(conf_dir+"moneda.dat"):
    from lib.firstrun import getrates
    try:
        getrates(conf_dir+"moneda.dat")
    except urllib.error.URLError:
        origen = os.environ["pychemqt"]+"dat"+os.sep+"moneda.dat"
        shutil.copy(origen, conf_dir+"moneda.dat")
        print((QtWidgets.QApplication.translate("pychemqt",
              "Internet connection error, using archived currency rates")))

splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Checking custom database..."))
from lib.sql import createDatabase

# Import internal libraries
splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Importing libraries..."))
from UI import texteditor, newComponent, flujo, plots, viewComponents
import plots as charts
from UI.widgets import createAction, ClickableLabel, TreeEquipment
from lib.config import conf_dir, getComponents
from lib.project import Project
from lib.EoS import K, H
from lib import unidades


splash.showMessage(QtWidgets.QApplication.translate("pychemqt",
                                                "Importing equipments..."))
from equipment import *

splash.showMessage(QtWidgets.QApplication.translate("pychemqt", "Importing tools..."))
from tools import UI_confComponents, UI_confTransport, UI_confThermo, UI_confUnits, UI_confResolution, UI_databank, UI_unitConverter, UI_steamTables, UI_psychrometry
from UI.conversor_unidades import moneda

splash.showMessage(QtWidgets.QApplication.translate("pychemqt", "Loading main window..."))
from UI.mainWindow import UI_pychemqt
pychemqt = UI_pychemqt()

splash.showMessage(QtWidgets.QApplication.translate("pychemqt", "Loading project files..."))
logging.info(QtWidgets.QApplication.translate("pychemqt", "Loading project files"))

pychemqt.show()

filename = []
if pychemqt.Preferences.getboolean("General", 'Load_Last_Project'):
    filename = pychemqt.settings.value("LastFile")
    if filename is None:
        filename = []
for file in args:
    filename.append(file)
for fname in filename:
    if fname and QtCore.QFile.exists(fname):
        splash.showMessage(QtWidgets.QApplication.translate("pychemqt", "Loading project files...")+"\n"+fname)
        logging.info(QtWidgets.QApplication.translate("pychemqt", "Loading project")+ ": %s" %fname)
        pychemqt.fileOpen(fname)
splash.finish(pychemqt)

sys.exit(app.exec_())
