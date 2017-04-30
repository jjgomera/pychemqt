#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
from unittest import TestLoader, TestSuite

loader = TestLoader()
path = os.path.dirname(os.path.realpath(sys.argv[0]))
tests = loader.discover(os.path.join(path, "lib", "mEoS"), pattern="*.py")
TestMEOS = TestSuite(tests)
