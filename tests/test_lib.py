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


import os
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
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
tests = loader.discover(os.path.join(path, "lib", "mEoS"), pattern="*.py")
TestMEOS = TestSuite(tests)
TestLib.addTest(TestMEOS)

# Add lib.EoS submodule test
for module in EoS.__all__:
    TestLib.addTest(DocTestSuite(module))

# Add lib.EoS.Cubic submodule test
for module in EoS.Cubic._all:
    TestLib.addTest(DocTestSuite(module.__module__))
