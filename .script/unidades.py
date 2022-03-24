#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Unidades first start to add to firstrun module

import sys
from pprint import pprint

from PyQt5 import QtCore, QtWidgets


# Qt application definition
app = QtWidgets.QApplication(sys.argv)
app.setOrganizationName("pychemqt")
app.setOrganizationDomain("pychemqt")
app.setApplicationName("pychemqt")
# Translations
locale = QtCore.QLocale.system().name()
myTranslator = QtCore.QTranslator()
if myTranslator.load("pychemqt_" + locale, "/home/jjgomera/pychemqt/i18n"):
    app.installTranslator(myTranslator)

from lib.unidades import _magnitudes  # noqa


# For get a fresh new list of magnitudes when we add some new, the list can be
# add at start of lib/firstrun.py file:
magnitudes = []
for magnitud, title, unit in _magnitudes:
    magnitudes.append(magnitud)
pprint(magnitudes, compact=True, indent=4, width=79)
