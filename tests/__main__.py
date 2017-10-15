#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Initialization of environment
import os
os.environ["pychemqt"] = "/home/travis/build/jjgomera/pychemqt/"
os.environ["freesteam"] = "True"
os.environ["pybel"] = "True"
os.environ["CoolProp"] = "True"
os.environ["refprop"] = "True"
os.environ["ezodf"] = "True"
os.environ["openpyxl"] = "True"
os.environ["xlwt"] = "False"
os.environ["icu"] = "False"
os.environ["reportlab"] = "False"
os.environ["PyQt5.Qsci"] = "True"


# Main script for unit testing

from unittest import TextTestRunner, TestSuite

from test_lib import TestLib

suite = TestSuite()
suite.addTest(TestLib)

runner = TextTestRunner()
results = runner.run(suite)
