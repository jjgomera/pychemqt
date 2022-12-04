#!/usr/bin/python3
# -*- coding: utf-8 -*-

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

###############################################################################
# Module to support both PyQt5 and PyQt6 libraries in pychemqt project
# For future in mind let PyQt6 as default version and monkey patch the PyQt5
# library with the minor changes necessary to work
###############################################################################


try:
    from PyQt6 import QtWidgets, QtGui, QtCore, QtSvg, QtSvgWidgets

    try:
        from PyQt6 import Qsci
    except ImportError:
        Qsci = False

    # Define qt version, for check version if it's necessary different code
    __qt__ = 6

except ImportError:
    from PyQt5 import QtWidgets, QtGui, QtCore, QtSvg

    try:
        from PyQt5 import Qsci
    except ImportError:
        Qsci = False

    __qt__ = 5

    # Define changes
    QtCore.QLibraryInfo.path = QtCore.QLibraryInfo.location
    QtCore.QLibraryInfo.LibraryPath = QtCore.QLibraryInfo.LibraryLocation

    QtSvgWidgets = QtSvg
    QtGui.QAction = QtWidgets.QAction
    QtGui.QActionGroup = QtWidgets.QActionGroup
    QtGui.QShortcut = QtWidgets.QShortcut

__all__ = ["QtCore", "QtGui", "QtWidgets", "QtSvg", "QtSvgWidgets", "Qsci"]
