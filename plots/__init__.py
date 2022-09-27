#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Module to define plots tools
###############################################################################

from PyQt5.QtWidgets import QApplication

from plots import drag, moody, standing
from plots.heatTransfer import chartHE


_all = {
    QApplication.translate("pychemqt", "Petro"): (standing.Standing_Katz, ),
    QApplication.translate("pychemqt", "Fluid Flow"): (moody.Moody, drag.Drag),
    QApplication.translate("pychemqt", "Heat Exchanger"): chartHE}

__all__ = ["moody", "standing", "chartHE"]

# List used in Preferences main window to add a subfolder of plot configuration
Pref = moody.Config, drag.Config, standing.Config
