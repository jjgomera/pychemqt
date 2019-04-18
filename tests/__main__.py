#!/usr/bin/python3
# -*- coding: utf-8 -*-

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
