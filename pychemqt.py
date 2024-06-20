#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=ungrouped-imports,unused-import,wrong-import-order
# pylint: disable=wrong-import-position, too-few-public-methods

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


import argparse
from configparser import ConfigParser
import logging
import os
import shutil
import sys
import urllib.error


# Parse command line options
desc = """pychemqt intended as a free software tool for calculation and \
design of chemical engineering unit operations."""
further = """For any suggestions, comments, bug ... you can contact me at \
https://github.com/jjgomera/pychemqt or by email jjgomera@gmail.com."""

parser = argparse.ArgumentParser(description=desc, epilog=further)
parser.add_argument("-l", "--log", dest="loglevel", default="INFO",
                    help="Set level of report in log file")
parser.add_argument("--debug", action="store_true",
                    help="Enable loglevel to debug, the more verbose option")
parser.add_argument("-n", "--nosplash", action="store_true",
                    help="Don't show the splash screen at start")
parser.add_argument("--style", help="Set qt style")
parser.add_argument("projectFile", nargs="*",
                    help="Optional pychemqt project files to load at startup")
args = parser.parse_args()


# Add pychemqt folder to python path
path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(path)

# Define pychemqt environment
os.environ["pychemqt"] = path + os.sep
conf_dir = os.path.expanduser("~") + os.sep + ".pychemqt" + os.sep

# Check mandatory external dependences
# qt
try:
    from tools.qt import QtCore, QtGui, QtWidgets
except ImportError as err:
    print("PyQt could not be found, you must install it" + os.linesep)
    print("PyQt5 and PyQt6 are supported")
    raise err

# Qt application definition
app = QtWidgets.QApplication(sys.argv)
app.setOrganizationName("pychemqt")
app.setOrganizationDomain("pychemqt")
app.setApplicationName("pychemqt")

# Qt style definition
if args.style is not None:
    style = QtWidgets.QStyleFactory.create(args.style)
    if style:
        app.setStyle(style)

# Add style options
app.setStyleSheet(
    "QDialogButtonBox {dialogbuttonbox-buttons-have-icons: true;}")


# Check qt configuration file
settings = QtCore.QSettings()
if not settings.contains("LastFile"):
    filename = QtCore.QVariant()
    settings.setValue("LastFile", filename)
    recentFiles = QtCore.QVariant()
    settings.setValue("RecentFiles", recentFiles)
    settings.setValue("Geometry", QtCore.QVariant())
    settings.setValue("MainWindow/State", QtCore.QVariant())


# Translation
locale = QtCore.QLocale.system().name()
myTranslator = QtCore.QTranslator()
if myTranslator.load("pychemqt_" + locale, os.environ["pychemqt"] + "i18n"):
    app.installTranslator(myTranslator)
qtTranslator = QtCore.QTranslator()
path = QtCore.QLibraryInfo.path(QtCore.QLibraryInfo.LibraryPath.TranslationsPath)
if qtTranslator.load("qt_" + locale, path):
    app.installTranslator(qtTranslator)


# scipy
try:
    import scipy
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "scipy could not be found, you must install it.")
    print(msg)
    raise err
mayor, minor = map(int, scipy.version.version.split(".")[:2])
if mayor == 0 and minor < 14:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Your version of scipy is too old, you must update it.")
    raise ImportError(msg)

# numpy
try:
    import numpy as np
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "numpy could not be found, you must install it.")
    print(msg)
    raise err
mayor, minor = map(int, np.version.version.split(".")[:2])
if mayor < 1 or minor < 8:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Your version of numpy is too old, you must update it.")
    raise ImportError(msg)

# matplotlib
try:
    import matplotlib
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "matplotlib could not be found, you must install it.")
    print(msg)
    raise err
mayor, minor = map(int, matplotlib.__version__.split(".")[:2])
if mayor < 1 or (mayor == 1 and minor < 4):
    msg = QtWidgets.QApplication.translate(
        "pychemqt",
        "Your version of matplotlib is too old, you must update it.")
    raise ImportError(msg)

# iapws
# Externalized version of iapws, to avoid duple maintenance
try:
    import iapws
except ImportError as err:
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "iapws could not be found, you must install it.")
    print(msg)
    raise err
if iapws.__version__ < "1.5.3":
    msg = QtWidgets.QApplication.translate(
        "pychemqt", "Your version of iapws is too old, you must update it.")
    raise ImportError(msg)


# TODO: Disable python-graph external dependence, functional mock up in
# project yet useless
# python-graph
# try:
    # from pygraph.classes.graph import graph
    # from pygraph.algorithms.cycles import find_cycle
# except ImportError as err:
    # msg = QtWidgets.QApplication.translate(
        # "pychemqt", "Python-graph don't found, you need install it")
    # print(msg)
    # raise err


# Check external optional modules
from tools.dependences import optional_modules  # noqa
for module, use in optional_modules:
    if module == "Qsci":
        # Special case for Qsci, a optional module from qt
        from tools.qt import Qsci
        if Qsci:
            os.environ[module] = "True"
        else:
            print(QtWidgets.QApplication.translate(
                "pychemqt", f"{module} could not be found, {use}"))
            os.environ[module] = ""
    else:
        try:
            __import__(module)
            os.environ[module] = "True"
        except ImportError:
            print(QtWidgets.QApplication.translate(
                "pychemqt", f"{module} could not be found, {use}"))
            os.environ[module] = ""
        else:
            # Check required version
            if module == "CoolProp":
                import CoolProp.CoolProp as CP
                version = CP.get_global_param_string("version")
                mayor, minor = map(int, version.split(".")[:2])
                if mayor < 6:
                    print(QtWidgets.QApplication.translate(
                        "pychemqt",
                        f"Find CoolProp {version} but CoolProp 6 required"))
                    os.environ[module] = ""


# Logging configuration
if args.debug:
    loglevel = "DEBUG"
else:
    loglevel = args.loglevel
loglevel = getattr(logging, loglevel.upper())

# Checking config folder
if not os.path.isdir(conf_dir):
    os.mkdir(conf_dir)

fmt = "[%(asctime)s.%(msecs)d] %(levelname)s: %(message)s"
logging.basicConfig(filename=conf_dir + "pychemqt.log", filemode="w",
                    level=loglevel, datefmt="%d-%b-%Y %H:%M:%S", format=fmt)
logging.info(QtWidgets.QApplication.translate("pychemqt", "Starting pychemqt"))


# Derive numpy error log to pychemqt log
class NumpyErrorLog():
    """Numpy error message catch and send to pychemqt log
    Use debug level for this messages"""
    @staticmethod
    def write(message):
        """Write error message to log file"""
        logging.debug(message)


np.seterrcall(NumpyErrorLog)
np.seterr(all='log')


class SplashScreen(QtWidgets.QSplashScreen):
    """Class to define a splash screen to show loading progress"""
    def __init__(self):
        QtWidgets.QSplashScreen.__init__(
            self,
            QtGui.QPixmap(os.path.join(
                os.environ["pychemqt"], "images", "splash.jpg")))
        QtWidgets.QApplication.processEvents()

    def showMessage(self, message):
        """Procedure to update message in splash"""
        align = (QtCore.Qt.AlignmentFlag.AlignBottom
                 | QtCore.Qt.AlignmentFlag.AlignRight
                 | QtCore.Qt.AlignmentFlag.AlignAbsolute)
        color = QtGui.QColor(QtCore.Qt.GlobalColor.white)
        QtWidgets.QSplashScreen.showMessage(self, message, align, color)
        QtWidgets.QApplication.processEvents()

    def clearMessage(self):
        """Clear message of splash screen"""
        QtWidgets.QSplashScreen.clearMessage(self)
        QtWidgets.QApplication.processEvents()


splash = SplashScreen()
if not args.nosplash:
    splash.show()


# Checking config files
from tools import firstrun  # noqa
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking config files..."))

# Checking config file
default_Preferences = firstrun.Preferences()
change = False
if not os.path.isfile(conf_dir + "pychemqtrc"):
    with open(conf_dir + "pychemqtrc", "w") as conf_file:
        default_Preferences.write(conf_file)
        Preferences = default_Preferences
        change = True
else:
    # Check Preferences options to find set new options
    Preferences = ConfigParser()
    Preferences.read(conf_dir + "pychemqtrc")
    for section in default_Preferences.sections():
        if not Preferences.has_section(section):
            Preferences.add_section(section)
            change = True
        for option in default_Preferences.options(section):
            if not Preferences.has_option(section, option):
                value = default_Preferences.get(section, option)
                Preferences.set(section, option, value)
                change = True
                logging.warning("Using default configuration option for "
                                "%s:%s, run preferences dialog for configure",
                                section, option)
    if change:
        with open(conf_dir + "pychemqtrc", "w") as conf_file:
            Preferences.write(conf_file)

# FIXME: This file might not to be useful but for now I use it to save project
# configuration data
if not os.path.isfile(conf_dir + "pychemqtrc_temporal"):
    Config = firstrun.config()
    with open(conf_dir + "pychemqtrc_temporal", "w") as conf_file:
        Config.write(conf_file)

# Checking costindex
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking cost index..."))
if not os.path.isfile(conf_dir + "CostIndex.dat"):
    orig = os.path.join(os.environ["pychemqt"], "dat", "costindex.dat")
    with open(orig) as cost_index:
        lista = cost_index.readlines()[-1].split(" ")
        with open(conf_dir + "CostIndex.dat", "w") as archivo:
            for data in lista:
                archivo.write(data.replace(os.linesep, "") + os.linesep)

# Checking currency rates
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking currency data"))
if not os.path.isfile(conf_dir + "moneda.dat"):
    # Exchange rates file don't available
    # Try to retrieve exchange rates from web service
    try:
        firstrun.getrates(conf_dir + "moneda.dat")
    except (urllib.error.URLError, urllib.error.HTTPError) as err:
        # Internet error, get hardcoded exchanges from pychemqt distribution
        # Possible outdated file, try to update each some commits
        logging.error(err)
        origen = os.path.join(os.environ["pychemqt"], "dat", "moneda.dat")
        shutil.copy(origen, conf_dir + "moneda.dat")
        print(QtWidgets.QApplication.translate("pychemqt",
                 "Internet connection error, using archived currency rates"))

# Checking database with custom components
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Checking custom database..."))
if not os.path.isfile(conf_dir + "databank.db"):
    firstrun.createDatabase(conf_dir + "databank.db")

# Import internal libraries
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Importing libraries..."))
import lib  # noqa
import UI  # noqa
import equipment  # noqa
import tools  # noqa
import plots  # noqa

# Load main program UI
splash.showMessage(QtWidgets.QApplication.translate(
    "pychemqt", "Loading main window..."))
from UI.mainWindow import UI_pychemqt  # noqa
pychemqt = UI_pychemqt()

# Load project files, opened in last pychemqt session and/or specified in
# command line
txt = QtWidgets.QApplication.translate("pychemqt", "Loading project files")
splash.showMessage(txt + "...")
logging.info(txt)

if change:
    lib.config.Preferences = Preferences

filename = []
if lib.config.Preferences.getboolean("General", "Load_Last_Project"):
    filename = pychemqt.lastFile
    if filename is None:
        filename = []
for file in args.projectFile:
    filename.append(file)
for fname in filename:
    if fname and QtCore.QFile.exists(fname):
        splash.showMessage(txt + "\n" + fname)
        logging.info(txt + ": " + fname)
        pychemqt.fileOpen(fname)


def exceptfunction(error, message, traceback):
    """Manage error message to avoid print to console"""
    sys.__excepthook__(error, message, traceback)


sys.excepthook = exceptfunction

# Finish splash and start qt main loop
pychemqt.show()
splash.finish(pychemqt)
sys.exit(app.exec())
