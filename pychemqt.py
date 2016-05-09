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


import logging
from optparse import OptionParser
import os
import shutil
import sys
import urllib.error


# Add pychemqt folder to python path
path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(path)

# Define pychemqt environment
os.environ["pychemqt"] = path + os.path.sep
conf_dir = os.path.expanduser('~') + os.sep+".pychemqt"+os.sep

# Check mandatory external dependences
# PyQt5
try:
    from PyQt5 import QtCore, QtGui, QtWidgets
except ImportError as err:
    print("PyQt5 don't found, you need install it")
    raise err

# Qt application definition
app = QtWidgets.QApplication(sys.argv)
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


# scipy
try:
    import scipy
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Python scipy library don't found, you need install it")
    print(msg)
    raise err
else:
    mayor, minor, corr = map(int, scipy.version.version.split("."))
    if minor < 14:
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "scipy version too old, try to install a updated version")  # noqa
        raise ImportError(msg)

# numpy
try:
    import numpy
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Python numpy library don't found, you need install it")
    print(msg)
    raise err
else:
    mayor, minor, corr = map(int, numpy.version.version.split("."))
    if mayor < 1 or minor < 8:
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "numpy version too old, try to install a updated version")  # noqa
        raise ImportError(msg)

# matplotlib
try:
    import matplotlib
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Python matplotlib don't found, you need install it")
    print(msg)
    raise err
else:
    mayor, minor, corr = map(int, matplotlib.__version__.split("."))
    if mayor < 1 or minor < 4:
        msg = QtWidgets.QApplication.translate(
            "pychemqt", "matplotlib version too old, try to install a updated version")  # noqa
        raise ImportError(msg)

# python-graph
try:
    from pygraph.classes.graph import graph  # noqa
    from pygraph.algorithms.cycles import find_cycle  # noqa
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Python-graph don't found, you need install it")
    print(msg)
    raise err


# Check external optional modules
from tools.dependences import optional_modules  # noqa
for module, use in optional_modules:
    try:
        __import__(module)
        os.environ[module] = "True"
    except ImportError:
        print("%s don't found, %s" % (module, use))
        os.environ[module] = ""

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

# Logging configuration
if not os.path.isfile(conf_dir + "pychemqt.log"):
    os.mknod(conf_dir + "pychemqt.log")
fmt = '[%(asctime)s.%(msecs)d] %(levelname)s: %(message)s'
logging.basicConfig(filename=conf_dir+'pychemqt.log', filemode='w',
                    level=loglevel, datefmt='%d-%b-%Y %H:%M:%S', format=fmt)
logging.info(
    QtWidgets.QApplication.translate("pychemqt", "Starting pychemqt"))


class SplashScreen(QtWidgets.QSplashScreen):
    """Class to defne a splash screen to show progress of loading"""
    def __init__(self):
        QtWidgets.QSplashScreen.__init__(
            self,
            QtGui.QPixmap(os.environ["pychemqt"] + "/images/splash.jpg"))
        self.show()
        QtWidgets.QApplication.flush()

    def showMessage(self, msg):
        """Método para mostrar mensajes en la parte inferior de la ventana de
        splash"""
        align = QtCore.Qt.Alignment(QtCore.Qt.AlignBottom |
                                    QtCore.Qt.AlignRight |
                                    QtCore.Qt.AlignAbsolute)
        color = QtGui.QColor(QtCore.Qt.white)
        QtWidgets.QSplashScreen.showMessage(self, msg, align, color)
        QtWidgets.QApplication.processEvents()

    def clearMessage(self):
        QtWidgets.QSplashScreen.clearMessage(self)
        QtWidgets.QApplication.processEvents()

splash = SplashScreen()


# Checking config files
from lib import firstrun  # noqa
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking config files..."))

# Checking config folder
if not os.path.isdir(conf_dir):
    os.mkdir(conf_dir)

# Checking config file
if not os.path.isfile(conf_dir + "pychemqtrc"):
    Preferences = firstrun.Preferences()
    Preferences.write(open(conf_dir + "pychemqtrc", "w"))

# FIXME: Hasta que no sepa como prescindir de este archivo sera necesario
if not os.path.isfile(conf_dir + "pychemqtrc_temporal"):
    Config = firstrun.config()
    Config.write(open(conf_dir + "pychemqtrc_temporal", "w"))

# Checking costindex
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking cost index..."))
if not os.path.isfile(conf_dir + "CostIndex.dat"):
        with open(os.environ["pychemqt"] + "dat/costindex.dat") as cost_index:
            lista = cost_index.readlines()[-1][:-1].split(" ")
            with open(conf_dir + "CostIndex.dat", "w") as archivo:
                for data in lista:
                    archivo.write(data + os.linesep)

# Checking currency rates
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking currency data"))
if not os.path.isfile(conf_dir+"moneda.dat"):
    try:
        firstrun.getrates(conf_dir+"moneda.dat")
    except urllib.error.URLError:
        origen = os.environ["pychemqt"]+"dat"+os.sep+"moneda.dat"
        shutil.copy(origen, conf_dir+"moneda.dat")
        print(QtWidgets.QApplication.translate("pychemqt",
              "Internet connection error, using archived currency rates"))

# Checkin database with custom components
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking custom database..."))
from lib.sql import createDatabase  # noqa
if not os.path.isfile(conf_dir + "databank.db"):
    createDatabase(conf_dir + 'databank.db')


# Import internal libraries
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Importing libraries..."))
from lib import *  # noqa
from UI import *  # noqa
from equipment import UI_equipments, equipments  # noqa
from tools import *  # noqa
from plots import *  # noqa


splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Loading main window..."))
from UI.mainWindow import UI_pychemqt  # noqa
pychemqt = UI_pychemqt()

msg = QtWidgets.QApplication.translate("pychemqt", "Loading project files")
splash.showMessage(msg + "...")
logging.info(msg)

filename = []
if pychemqt.Preferences.getboolean("General", 'Load_Last_Project'):
    filename = pychemqt.lastFile
    if filename is None:
        filename = []
for file in args:
    filename.append(file)
for fname in filename:
    if fname and QtCore.QFile.exists(fname):
        msg = QtWidgets.QApplication.translate("pychemqt",
                                               "Loading project files...")
        splash.showMessage(msg + "\n" + fname)
        logging.info(msg + ": " + fname)
        pychemqt.fileOpen(fname)

pychemqt.show()
splash.finish(pychemqt)

sys.exit(app.exec_())
