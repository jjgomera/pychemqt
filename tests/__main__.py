#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Main script for unit testing


import os
from unittest import TextTestRunner, TestSuite

from test_lib import TestLib


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


suite = TestSuite()
suite.addTest(TestLib)

runner = TextTestRunner()
results = runner.run(suite)
