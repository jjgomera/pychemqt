#!/usr/bin/python
# -*- coding: utf-8 -*-


import os
import sys
from doctest import DocTestSuite
from unittest import TestLoader, TestSuite

import lib
from lib import EoS

TestLib = TestSuite()
for mname in lib.__all__[2:]:
    module = lib.__getattribute__(mname)
    TestLib.addTest(DocTestSuite(module))

# Add lib.mEoS submodule test
loader = TestLoader()
path = os.path.dirname(os.path.realpath(sys.argv[0]))
tests = loader.discover(os.path.join(path, "lib", "mEoS"), pattern="*.py")
TestMEOS = TestSuite(tests)
TestLib.addTest(TestMEOS)

# Add lib.EoS submodule test
for module in EoS.__all__:
    TestLib.addTest(DocTestSuite(module))

# Add lib.EoS.Cubic submodule test
for module in EoS.Cubic._all:
    TestLib.addTest(DocTestSuite(module.__module__))
