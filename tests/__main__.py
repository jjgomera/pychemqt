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


# Main script for unit testing


import os
import warnings
from unittest import TestSuite, TextTestRunner

import numpy as np

from test_lib import TestLib

# Define pychemqt environment
os.environ["pychemqt"] = os.path.abspath('.')
os.environ["freesteam"] = "False"
os.environ["openbabel"] = "False"
os.environ["CoolProp"] = "False"
os.environ["refprop"] = "False"
os.environ["ezodf"] = "False"
os.environ["openpyxl"] = "False"
os.environ["xlwt"] = "False"
os.environ["icu"] = "False"
os.environ["reportlab"] = "False"
os.environ["Qsci"] = "False"


# Don't print the numpy RuntimeWarning
np.seterr("ignore")

warnings.simplefilter("ignore")

suite = TestSuite()
suite.addTest(TestLib)

runner = TextTestRunner(failfast=False)
results = runner.run(suite)
