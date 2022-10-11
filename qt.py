#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


try:
    from PyQt6 import QtWidgets, QtGui, QtCore, QtSvg, QtSvgWidgets, Qsci

    # Define qt version if it's necessary
    __qt__ = 6

    # Define changes
    QtCore.QLibraryInfo.location = QtCore.QLibraryInfo.path
    QtSvg.QGraphicsSvgItem = QtSvgWidgets.QGraphicsSvgItem
    QtWidgets.QAction = QtGui.QAction

except ImportError:
    from PyQt5 import QtWidgets, QtGui, QtCore, QtSvg, Qsci

    __qt__ = 5
