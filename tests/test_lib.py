#!/usr/bin/python
# -*- coding: utf-8 -*-

from doctest import DocTestSuite
from unittest import TestSuite

import lib

TestLib = TestSuite()
for mname in lib.__all__[2:]:
    module = lib.__getattribute__(mname)
    TestLib.addTest(DocTestSuite(module))
