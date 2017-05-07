#!/usr/bin/python
# -*- coding: utf-8 -*-

from doctest import DocTestSuite
from unittest import TestSuite

from lib import adimensional, petro, unidades


TestLib = TestSuite()
TestLib.addTest(DocTestSuite(unidades))
TestLib.addTest(DocTestSuite(adimensional))
TestLib.addTest(DocTestSuite(petro))
