#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Main script for unit testing

from unittest import TextTestRunner, TestSuite

from test_lib import TestLib

suite = TestSuite()
suite.addTest(TestLib)

runner = TextTestRunner()
results = runner.run(suite)
