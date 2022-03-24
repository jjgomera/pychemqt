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


# Main script for unit testing


import os

# Define pychemqt environment
os.environ["pychemqt"] = os.path.abspath('.')
os.environ["freesteam"] = "False"
os.environ["pybel"] = "False"
os.environ["CoolProp"] = "False"
os.environ["refprop"] = "False"
os.environ["ezodf"] = "False"
os.environ["openpyxl"] = "False"
os.environ["xlwt"] = "False"
os.environ["icu"] = "False"
os.environ["reportlab"] = "False"
os.environ["PyQt5.Qsci"] = "False"


# Don't print the numpy RuntimeWarning
from numpy import seterr
seterr("ignore")

import warnings
warnings.simplefilter("ignore")

from unittest import TextTestRunner, TestSuite
from test_lib import TestLib

suite = TestSuite()
suite.addTest(TestLib)

runner = TextTestRunner(failfast=True)
results = runner.run(suite)
