#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# lib module
# module with general library functionality of pychemqt
###############################################################################


import glob
import os


files = glob.glob(os.path.join("lib", "*.py"))

__all__ = ["EoS", "mEoS"]
for file in files:
    fname, ext = os.path.splitext(os.path.basename(file))
    if fname != "__init__" and fname not in ["project"]:
        __all__.append(fname)

# Sort without upper-lower case classification
__all__.sort(key=lambda v: v.upper())
