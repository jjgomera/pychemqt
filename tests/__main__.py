#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Main script for unit testing

from unittest import TextTestRunner, TestSuite

from tests.test_meos import TestMEOS
from tests.test_lib import TestLib

suite = TestSuite()
suite.addTest(TestMEOS)
suite.addTest(TestLib)

runner = TextTestRunner()
results = runner.run(suite)
