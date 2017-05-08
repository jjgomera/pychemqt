#!/usr/bin/python
# -*- coding: utf-8 -*-

from doctest import DocTestSuite
from unittest import TestSuite

from lib import adimensional, compuestos, petro, unidades


TestLib = TestSuite()
TestLib.addTest(DocTestSuite(adimensional))
TestLib.addTest(DocTestSuite(compuestos))
TestLib.addTest(DocTestSuite(petro))
TestLib.addTest(DocTestSuite(unidades))
